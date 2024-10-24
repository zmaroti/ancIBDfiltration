#!/usr/bin/python3.6
# import needed for the command part (when code is split to cmd/class)
import argparse
import re

# needed for class
import h5py
import numpy as np
import pandas as pd

class IBDScore:
    def readH5(self, h5Fn):
        f = h5py.File(h5Fn, 'r')

        # extract the map values only (we know that raw IBD marker indexes correspond to the indexes of this array)
        mapM = np.array(f['variants']['MAP'])

        # store the map data in a list for the given chromosome
        self.mapData.append(mapM)

        # this list is used to store an offset for the first indexes of each chromosome in a single
        # ndarray (markerIBDCounts) this way can have a global median/etc
        self.markerOffsets.append(self.totalMarkerCount) 

        # size of total marker counts
        self.totalMarkerCount += mapM.size 

    def calcStats(self):
        # calculate median, Q1, Q3, IQR of raw IBD count of all markers that is within the IBDs (means at least 1 IBD crosses it)
        withinIBDCounts = self.markerIBDCounts[self.markerIBDCounts > 0]
        self.median = np.median  (withinIBDCounts)
        self.Q1     = np.quantile(withinIBDCounts, 0.25)
        self.Q3     = np.quantile(withinIBDCounts, 0.75)
        self.IQR    = self.Q3 - self.Q1

    def addIBDCounts(self, start, end):
        self.markerIBDCounts[start:end] += 1

    def IBDScore(self, start, end):
        return np.sum(self.markerScores[start:end])

    def readRawIBDs(self, prefix_rawIBD):
        # IBD count on all chromosomes without length filtration
        self.AllIBDCounts = []

        # list to store raw IBDs per chromosome that has sufficient length
        self.RawIBDs = []

        # a numpy array for all IBD counts of each markers at different chromosomes
        self.markerIBDCounts = np.zeros(self.totalMarkerCount, dtype=np.int32)

        if self.verbose:
            print("Reading raw IBD data: ", end = " ", flush = True)

        # read all raw IBD files as a pandas dataframes
        for chrIdx, ch in enumerate(self.chs):
            ibdFn = f"{prefix_rawIBD}{ch}.tsv"

            if self.verbose:
                print(ibdFn, end = " ", flush = True)

            # read data, update marker IBD counts
            allIBDs = pd.read_csv(f"{ibdFn}", sep="\t",
                dtype={'Start': np.int32, 'End': np.int32, 'StartM': np.float32, 'EndM': np.float32, 'length': np.int32, 'lengthM': np.float32, 'ch': np.int32, 'iid1': str, 'iid2': str})

            # row count of all unfiltered IBDs
            self.AllIBDCounts.append(allIBDs.shape[0])

            # filter on length criteria
            rawIBDs = allIBDs.loc[allIBDs['lengthM'] * 100 > self.cM].copy().reset_index(drop=True)

            # add density column for debug
            if self.debug:
                rawIBDs['Density'] = rawIBDs['length'] / (100 * rawIBDs['lengthM'])

            # append raw IBDs by chromosome
            self.RawIBDs.append(rawIBDs)

            # offset in the marker count data for the given chromosome
            chrOffset = self.markerOffsets[chrIdx]
            starts = chrOffset + rawIBDs['Start']
            ends   = chrOffset + rawIBDs['End']

            # add IBD counts to all included markers (for all IBDs of the given chromosome)
            # Note it is granted that marker indices correspond in marker/raw IBD output (in case it is from same analysis)
            np.vectorize(self.addIBDCounts)(starts, ends)

        # calculate median, Q1, Q3, IQR of raw IBD count of all markers that is within the IBDs (means at least 1 IBD crosses it)
        self.calcStats()

        if self.verbose:
            print ("Done")

        # calculate median of raw IBD count of all markers that is within the IBDs (means at least 1 IBD crosses it)
        self.median = np.median(self.markerIBDCounts[self.markerIBDCounts > 0])

    def scoreRawIBDs(self):
        if self.verbose:
            print(f"Scoring IBDs by marker information (expected/observed IBD counts of the included markers)")

        # calculate the marker scores
        # IBD counts below Q1 are treated as Q1 (so scores are never treated too optimistic)
        minCounts = np.where(self.markerIBDCounts < self.Q1, self.Q1, self.markerIBDCounts)

        # it is granted that all markers included in any IBD we will have a positive IBDcount (see readRawIBDs)
        # some markers may have 0 IBD crossing, there score is meaningless.
        self.markerScores = self.median / np.where(self.markerIBDCounts > 0, minCounts, np.nan)

        # score IBDs, iterate through each chromosome

        for chrIdx, chrIBDs in enumerate(self.RawIBDs):
            # offset in the marker count data for the given chromosome
            chrOffset = self.markerOffsets[chrIdx]

            # scores for these IBD entries
            starts = chrOffset + chrIBDs['Start']
            ends   = chrOffset + chrIBDs['End']

            # get the sum of all marker score for each IBD fragments (between start/end markers)
            scores = np.vectorize(self.IBDScore)(starts, ends)

            # add the IBD data df a new column with the IBD scores
            self.RawIBDs[chrIdx]['score'] = scores


    # Although in our current settings we had seen that UGLY areas are continuous, it is not granted and we may have
    # ugly pattern of GOOD/BAD regions in an IBD
    # 1) IBD extends into bad are, but good part is >cM  -> trim BAD, recalc length
    #     >cM
    # [------------BAAAAD] 
    #              XXXXXX
    # remainder > cM => trim BAD, recalculate length
    # 2) same as above, but BAAAD at one side is dispersed -> trim all regions till we have >cM OK area, recalc length
    # [--BAAAD----BAAAAD-----------]
    #  XXXXXXXXXXXXXXXXX
    # 3) not all BAD has to be eliminated, here we have evidence that IBD is this large
    #      >cM                    >cM
    # [-----------------BAAAD-----------------------]  -> keep the IBD, no recalculatiuon
    # 4) a highly unlikely case
    #  <cM      >cM                  >cM              <cM
    # [--BAD---------------BAAAD-----------------BAD------]  -> trim the IBD (both sides), recalc the length
    #  XXXXX                                     XXXXXXXXX
    # 5) the total length is >cM, still each good segments individual lengths are less than that -> eliminiate the IBD
    # [------BAAAAAD---BAAAAAAAD---BAAAAAAAD-BAAAAAAD]
    #        XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    #
    # In case we trimmed an IBD, check if remainder < cM => remove IBD
    # #################################################################################################
    # the strategy should be to trim from both sides till we have a >cM stretch of OK IBD segment
    # throw out IBD where we don't have that, recaclulate the size of the OK part
    def trimLeft(self, markerM, starts, ends):
        # we have the start/end marker index of bad regions in starts, and ends
        # the aim is to find the first OK region from right to left that is longer > cM
        badIdx       = 0 # first idx in the bad regions
        mapMarkerIdx = 0 # first idx in the markerMap
        while badIdx < len(starts):
            if mapMarkerIdx < starts[badIdx]:
                okLenM = markerM[starts[badIdx] - 1] - markerM[mapMarkerIdx]

                if okLenM * 100 > self.cM:
                    return mapMarkerIdx, True, badIdx

            mapMarkerIdx = ends[badIdx] + 1
            badIdx = badIdx + 1

        # check if we have still good markers left
        if mapMarkerIdx < len(markerM):
            # we exited above loop since all badIdxs were consumed, so check if end of last bad region to end of IBD segment is >cM
            okLenM = markerM[len(markerM)-1] - markerM[ends[badIdx - 1] + 1]

            if okLenM * 100 > self.cM:
                return mapMarkerIdx, True, badIdx

        # either we consumed all markerIdx or last OK segment was <cM
        return mapMarkerIdx, False, badIdx

    def trimRight(self, markerM, starts, ends, lastBadIdx):
        # we have the start/end marker index of bad regions in starts, and ends
        # the aim is to find the first OK region from right to left that is longer > cM
        badIdx       = len(starts)  - 1
        mapMarkerIdx = len(markerM) - 1

        # lastBadIdx holds the last bad region checked, we don't have to check more regions as we did from left already
        while badIdx >= lastBadIdx:
            if mapMarkerIdx > ends[badIdx]:
                okLenM = markerM[mapMarkerIdx] - markerM[ends[badIdx] + 1]

                if okLenM * 100 > self.cM:
                    return mapMarkerIdx

            mapMarkerIdx = starts[badIdx] - 1
            badIdx = badIdx - 1

        # we stripped all bad regions, just return the last OK markerIdx
        return mapMarkerIdx

    # find which OK IBDs has potential length issues (containing markers in problematic regions)
    # calculate a new length by trimming problematic regions from both sides
    # drop out IBD where newLength is < cM
    # recalculate the start/end of IBDs that were trimmed, based on the new length
    def fixOkIBDsLength(self, X):
        # calc X*IQR (SD) threshold above median
        SDX = self.median + X * self.IQR

        if self.verbose:
            print(f"Median={self.median}, lower quartile={self.Q1}, upper quartile={self.Q3}, IQR={self.IQR}, {X}SD={SDX}")
            print(f"Fixing length of IBDs containing BAD markers (>{X} SD above median expected IBD counts)")

        self.lengthFixDroppedIBDs = []

        for chrIdx, okIBDs in enumerate(self.OkIBDs):
            allMarkerM = self.mapData[chrIdx]

            # IBD indexes to dropout due to length threshold
            dropOutIdxs = []

            # iterate this way, as we can keep original row index labels
            for IBDidx in okIBDs.index:
                # start, end of marker idx (allMarkerM - marker map data of current chromosome)
                start = okIBDs['Start'][IBDidx]
                end   = okIBDs['End'][IBDidx]

                # start, end count index of the given IBD's markers in the markerIBDCounts data (offset + index)
                startCount = self.markerOffsets[chrIdx] + start
                endCount   = self.markerOffsets[chrIdx] + end

                # the marker counts of all markers in the current IBD segment
                markerCounts = self.markerIBDCounts[startCount:endCount]

                # check which markers are above thresh; pad the array at the flanks
                badMarkers = np.concatenate(([False], markerCounts > SDX, [False]))

                # start positions of bad marker blocks (true in array) in the IBD
                starts = np.flatnonzero(~badMarkers[:-1] &  badMarkers[1:])

                # we have 5-6% bad regions so most IBDs are expected to be OK, where we can skip early
                # no bad regions -> nothing to do
                if len(starts) == 0:
                    continue

                # end positions of bad marker blocks in the IBD
                ends   = np.flatnonzero( badMarkers[:-1] & ~badMarkers[1:]) - 1

                # map data of each markers in the current IBD segment
                markerM = allMarkerM[start:end]

                # trim regions with bad markers at the left
                newStart, okIBD, lastBadIdx = self.trimLeft(markerM, starts, ends)

                # if the size of IBD is less than cM after trimming, drop the IBD as FP
                if not okIBD:
                    dropOutIdxs.append(IBDidx)
                    continue    # skip rest

                # left we have IBD >cM; but we have still bad regions from the right
                if lastBadIdx < len(starts):
                    newEnd = self.trimRight(markerM, starts, ends, lastBadIdx)

                # otherwise End and EndM stays the original
                else:
                    newEnd = len(markerM) - 1

                # NOTE, markerM is is a slice that contains only marker Map info of the markers of the current IBD
                newStartM, newEndM = markerM[newStart], markerM[newEnd]

                # update IBD data
                self.OkIBDs[chrIdx].at[IBDidx, 'Start' ]  = start + newStart # add the offset (chrosome index) of the IBD start
                self.OkIBDs[chrIdx].at[IBDidx, 'StartM']  = newStartM
                self.OkIBDs[chrIdx].at[IBDidx, 'End']     = start + newEnd   # add the offset (chrosome index) of the IBD start
                self.OkIBDs[chrIdx].at[IBDidx, 'EndM']    = newEndM
                self.OkIBDs[chrIdx].at[IBDidx, 'length']  = newEnd - newStart
                self.OkIBDs[chrIdx].at[IBDidx, 'lengthM'] = newEndM - newStartM

            # keep track on IBDs excluded due to length fix for debug/statistic
            self.lengthFixDroppedIBDs.append(self.OkIBDs[chrIdx].loc[dropOutIdxs])

            # we processed all IBDs for the current chromosome and we may have some IBD indexes to drop (due to length recalc)
            if (len(dropOutIdxs) > 0):
                self.OkIBDs[chrIdx] = okIBDs.drop(dropOutIdxs)

    def FilterRawIBDs(self, thresh):
        if self.verbose:
            print(f"Filtering IBDs with >{thresh} score")

        for rawIBDs in self.RawIBDs:
            # filter based on length and score
            okIBDs = rawIBDs.loc[rawIBDs['score'] > thresh]

            # append to list
            self.OkIBDs.append(okIBDs)

    def saveTsv(self, outFn, dfs):
        if self.verbose:
            if (len(dfs) > 0):
                print(f"Saving output to: {outFn}")
            else:
                print(f"Nothing to save to: {outFn}")

        if(len(dfs) > 0):
            df = pd.concat(dfs)
            df.to_csv(outFn, index=False, sep='\t', na_rep='NA')

    def meanCounts(self, start, end):
        return np.mean(self.chrMarkerIBDCounts[start:end])

    # FOR DEBUGING
    def SaveBADMarkerRegions(self, X, outPref):
        SDX = self.median + X * self.IQR
        badRegions = []

        for chrIdx, ch in enumerate(self.chs):
            # we have all offests in this list
            start = self.markerOffsets[chrIdx]
            # marker counts of the given chromosome are between two offsets
            if chrIdx + 1 < len (self.markerOffsets):
                end = self.markerOffsets[chrIdx + 1]
            # except the last one, here we use the total marker counts
            else:
                end = self.totalMarkerCount

            # the IBD counts of all markers of the given chromosome
            self.chrMarkerIBDCounts = self.markerIBDCounts[start:end]

            # vector of True/False padded with two falses at sides
            chrbadMarkers = np.concatenate(([False], self.chrMarkerIBDCounts > SDX, [False]))

            # places of switch from false->true (signing start of a BAD region)
            starts = np.flatnonzero(~chrbadMarkers[:-1] &  chrbadMarkers[1:])

            if len(starts) > 0:
                # if we have starts, we have ends as well (we padded with a false at the end)
                ends   = np.flatnonzero( chrbadMarkers[:-1] & ~chrbadMarkers[1:])

                # mean IBD count of markers included in the masked region
                means  = np.vectorize(self.meanCounts)(starts, ends)

                allMarkerM = self.mapData[chrIdx] # mapinfo (by markeridxs) for the current chromosome
                startMs    = allMarkerM[starts]
                endMs      = allMarkerM[ends-1]

                br = pd.DataFrame()
                br['Start']  = starts
                br['End']    = ends
                br['length'] = ends - starts
                br['StartM'] = startMs
                br['EndM']   = endMs
                br['lengthM']   = endMs - startMs
                br['mean (IBD count)'] = means
                br.insert(0, 'CH', ch)
                badRegions.append(br)

        # save the bad regions in a file
        self.saveTsv(f"{outPref}_badMarkerRegions.tsv", badRegions)

    def PrintStats(self, outPref, thresh):
        # we have to do this after filtering so all IBDs that fell out due to length recalc is also considered
        consensusIBDs   = []
        onlyDensityIBDs = []
        onlyScoreIBDs   = []
        markerInfo      = []
        excludedIBDs    = []

        if self.verbose:
            print("STATS in IBD counts by chromosome (score/density method)")
            print(f"CH\tall\t>{self.cM}cM\tscore\tdensity\tcons\tonly dens\tonly score\tsize excl\tscore excl\tlength fix excl")

        for chrIdx, ch in enumerate(self.chs):
            # we have all offests in this list
            start = self.markerOffsets[chrIdx]
            # marker counts of the given chromosome are between two offsets
            if chrIdx + 1 < len (self.markerOffsets):
                end = self.markerOffsets[chrIdx + 1]
            # except the last one, here we use the total marker counts
            else:
                end = self.totalMarkerCount

            # save the actual chr markerIndex mapM information as lifted over markers indexes will be different (hard for debug)
            mi = pd.DataFrame()
            mi['mapM'] = self.mapData[chrIdx]
            mi.insert(0, 'CH', ch)
            mi.insert(1, 'markerIdx', range(len(mi)))
            # the IBD counts of all markers of the given chromosome
            mi['IBDCount']    = self.markerIBDCounts[start:end]
            mi['markerScore'] = self.median
            mi['markerScore'] = mi['markerScore'].div(mi['IBDCount'])

            markerInfo.append(mi)

            rawIBDs        = self.RawIBDs[chrIdx]
            okIBDs         = self.OkIBDs[chrIdx]
            lengthDropIBDs = self.lengthFixDroppedIBDs[chrIdx]
            densityIBDs    = rawIBDs.loc[(rawIBDs['length'] / (100 * rawIBDs['lengthM']) > 220)]
            excluded       = rawIBDs.loc[rawIBDs['score'] <= thresh]

            oIdxs = set(okIBDs.index.values)
            dIdxs = set(densityIBDs.index.values)

            consIdxs  = list(oIdxs & dIdxs)
            onlyDens  = list(dIdxs - oIdxs)
            onlyScore = list(oIdxs - dIdxs)

            # print stats
            if self.verbose:
                sizeExcl  = self.AllIBDCounts[chrIdx] - rawIBDs.shape[0]
                scoreExcl = rawIBDs.shape[0] - lengthDropIBDs.shape[0] - okIBDs.shape[0]
                print(f"{ch}\t{self.AllIBDCounts[chrIdx]}\t{rawIBDs.shape[0]}\t{okIBDs.shape[0]}\t{densityIBDs.shape[0]}\t{len(consIdxs)}\t{len(onlyDens)}\t{len(onlyScore)}\t{sizeExcl}\t{scoreExcl}\t{lengthDropIBDs.shape[0]}")

            consensusIBDs.append(rawIBDs.iloc[consIdxs])
            onlyDensityIBDs.append(rawIBDs.iloc[onlyDens])
            onlyScoreIBDs.append(rawIBDs.iloc[onlyScore])
            excludedIBDs.append(excluded)

        # for DEBUGGING we save IBDs filtered by original density criteria, new method, consenus (both)
        # dropped IBDs due to length fix (trimming IBD part extending into mask area)
        # and the marker data the analysis is relied on including the experimental IBD counts at each markers
        if self.debug:
            self.saveTsv(f"{outPref}_consensus.tsv", consensusIBDs)
            self.saveTsv(f"{outPref}_onlyDensity.tsv", onlyDensityIBDs)
            self.saveTsv(f"{outPref}_onlyScore.tsv", onlyScoreIBDs)
            self.saveTsv(f"{outPref}_lengthFixDrop.tsv", self.lengthFixDroppedIBDs)
            self.saveTsv(f"{outPref}_markerInfo.tsv", markerInfo)
            self.saveTsv(f"{outPref}_excluded.tsv", excludedIBDs)
            self.SaveBADMarkerRegions(args.SD, args.outprefix)

    def Save(self, outPref, thresh):
        # save the result
        self.saveTsv(f"{outPref}.tsv", self.OkIBDs)

        # print stats, potentially save supporting data (in debug mode)
        self.PrintStats(outPref, thresh)

    def __init__(self, args):
        # chromosomes to load/analyse
        self.chs = args.ch
        # length threshold in cM
        self.cM = args.cM
        # number of markers
        self.totalMarkerCount = 0
        # list of map data (lengthM) for each markers (per chromosome)
        self.mapData = []
        # cumulative marker counts of chromosomes
        self.markerOffsets  = []
        # list to store ok IBDs (per chromosome)
        self.OkIBDs  = []

        self.verbose = args.verbose
        self.debug   = args.debug

        # read marker data
        if self.verbose:
            print("Reading marker data: ", end = " ", flush = True)

        for ch in self.chs:
            h5Fn = f"{args.h5prefix}{ch}.h5"

            if self.verbose:
                print(h5Fn, end = " ", flush = True)

            self.readH5(h5Fn)

        if self.verbose:
            print("Done")

        # read the raw IBDs and calculate IBD count per marker
        self.readRawIBDs(args.rawIBDprefix)

        # calculate stats and score IBDs by marker information
        self.scoreRawIBDs()

########################################################################################################################
# command part for parsing arguments and call the code
def getChList(str):
    chs = []

    for tag in  str.split(","):
        match = re.match(r"^(\d+)([-]{,1})(\d*)$", tag)

        if match:
            if match.group(2) == "-":
                chs.extend( range(int(match.group(1)), int(match.group(3)) + 1) )
            else:
                chs.append(int(match.group(1)))
        else:
            print (f"Invalid chromosome {tag}")
            quit()

    return (chs)

parser = argparse.ArgumentParser(description='TRUTH - True Ratio of Unbiased IBD Through Hypothesis, a tool to filter raw ancIBD segments based on the experimental spatial IBD segment density across the genome.')
parser.add_argument('h5prefix', help='prefix of the h5 file (h5dir/chr)')
parser.add_argument('rawIBDprefix', help='prefix of the raw IBD file (outdir/ch)')
parser.add_argument('outprefix', help='prefix of the output file. Default ALL_FILTERED')
parser.add_argument('--ch', action='store', dest='ch', type=str, required=False, default='1-22', help='Chromosome (INT) or range of chromosomes (1-22), default: 1-22')
parser.add_argument('--cm', action='store', dest = 'cM', type=int, required=False, default=8, help='Minimum length of IBD fragments in cM to be kept, default 8')
parser.add_argument('--sc', action='store', dest = 'sc', type=int, required=False, default=1760, help='Minimum IBD score (integer) to be kept, default 8 * 220 = 1760')
parser.add_argument('--sd', action='store', dest = 'SD', type=int, required=False, default=6, help='Minimum SD above the median (by IQR) to consider marker as BAD at IBD length recalculation, default 3')
parser.add_argument('--verbose', action='store_true', dest='verbose', default = False, help='Output some progress messages, and summary IBD stats per chromosome')
parser.add_argument('--debug', action='store_true', dest='debug', default = False, help='Save also the consensus/only density/only score filtered, excluded IBDs, the badregions and the detailed marker information scores in separate files.')
args = parser.parse_args()

if args.verbose == True:
    print(f"Chromosome(s):          {args.ch}")
    print(f"h5 file prefix:         {args.h5prefix}")
    print(f"Raw IBD data prefix:    {args.rawIBDprefix}")
    print(f"outprefix:              {args.outprefix}")
    print(f"min IBD length (cM):    {args.cM}")
    print(f"min IBD score:          {args.sc}")
    print(f"min SD (length recalc): {args.SD}")

# parse chromosome list
args.ch = getChList(args.ch)

# print the run parameters in case we are verbose
# read map data, and raw IBD tsv files
IBDs = IBDScore(args)

# filter the IBDs by length, and score (default is 8cM length and min 8 * 220 score)
IBDs.FilterRawIBDs(args.sc)

# fix length of IBDs extending into BAD marker areas (default use 6SD above the meadian)
IBDs.fixOkIBDsLength(args.SD)

# save the filtered IBDs for all chromosomes in a single file
# if debug is also provided then it saves difference compared to 220SNP/cM method + list of IBDs dropped due to eliminating" BAD markers from IBD"
IBDs.Save(args.outprefix, args.sc)
