TRUTH - True Ratio of Unbiased IBD Through Hypothesis, a tool to filter raw ancIBD segments based on experimental IBD density.
======================================================
This repository contains the python code to filter the raw ancIBD output for true positive IBD segments.

The approach is described in detail in our upcoming manuscript:
Oszkár Schütz*, Zoltán Maróti*, Balázs Tihanyi, Attila P. Kiss, Emil Nyerki, Alexandra Gînguță, Petra Kiss, Gergely I. B. Varga, Bence Kovács, Kitti Maár, Bernadett Ny. Kovacsóczy, Nikoletta Lukács, István Major, Antónia Marcsik, Eszter Patyi, Anna Szigeti, Zoltán Tóth, Dorottya Walter, Gábor Wilhelm, Réka Cs. Andrási, Zsolt Bernert, Luca Kis, Liana Loredana Oța, György Pálfi, Gábor Pintye, Dániel Pópity, Angela Simalcsik, Andrei Dorian Soficaru, Olga Spekker, Sándor Varga, Endre Neparáczki and Tibor Török: **Unveiling the Origins and Genetic Makeup of the 'Forgotten People': A Study of the Sarmatian-Period Population in the Carpathian Basin**

doi: (in progress)

### Dependencies
The code depends on numpy, pandas and the h5py python libraries which needs to be installed in your python environment.

### Usage
```sh
usage: truthFilter.py [-h] [--ch CH] [--cm CM] [--sc SC] [--sd SD] [--verbose]
                      [--debug]
                      h5prefix rawIBDprefix outprefix

TRUTH - True Ratio of Unbiased IBD Through Hypothesis, a tool to filter raw
ancIBD segments based on the experimental spatial IBD segment density across
the genome.

positional arguments:
  h5prefix      prefix of the h5 file (h5dir/chr)
  rawIBDprefix  prefix of the raw IBD file (outdir/ch)
  outprefix     prefix of the output file. Default ALL_FILTERED

optional arguments:
  -h, --help    show this help message and exit
  --ch CH       Chromosome (INT) or range of chromosomes (1-22), default: 1-22
  --cm CM       Minimum length of IBD fragments in cM to be kept, default 8
  --sc SC       Minimum IBD score (integer) to be kept, default 8 * 220 = 1760
  --sd SD       Minimum SD above the median (by IQR) to consider marker as BAD
                at IBD length recalculation, default 3
  --verbose     Output some progress messages, and summary IBD stats per
                chromosome
  --debug       Save also the consensus/only density/only score filtered,
                excluded IBDs, the badregions and the detailed marker
                information scores in separate files.
```

In case the ancIBD hdf5 files are under the **hdf5** directory as **ch(1-22).h5** files and the raw ancIBD output is in the **output** directory as **ch(1-22).tsv** files, then the tool invoked with the appropriate parameters:
```sh
truthFilter.py hdf5/chr output/ch FILTERED
```
will filter out >8cM, >1760 IBD informativity score IBD segments. IBDs extending into regions with 3SD higher IBD count compared to the median IBD count of the analysed markers are truncated for the low confidence mask areas, and kept only when the remaining high confidence region is above the length threshold. The results are saved in the FILTERED.tsv file. The optional parameters can be used to control the length, score, mask thresholds, and the list of chromosomes to be analyzed. The *--verbose* flag will print out progress and some summary stats on the filtered IBD segments. The *--debug* flag will result in the following additional files:
* FILTERED_consensus.tsv
* FILTERED_onlyDensity.tsv
* FILTERED_onlyScore.tsv
* FILTERED_lengthFixDrop.tsv
* FILTERED_markerInfo.tsv
* FILTERED_excluded.tsv
* FILTERED_badMarkerRegions.tsv

### Method description and performance analysis
A [PDF](Method_description_and_performance_analysis.pdf) is available for detailed description of the underlying biological hypothesis and the in-depth comparison with the density based approach implemented in the ancIBD package. The performance of both methods were evaluated in a large cohort of experimental, low coverage ancient shotgun WGS data (described in the accompanying manuscript).
