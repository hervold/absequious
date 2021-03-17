# absequious
NGS antibody sequence analysis tool

Absequious uses a Hidden Markov Model to identify antibody variable region domains and calculate high-level diversity statistics.  It relies on [HMMer](http://hmmer.org/) for alignments, and is therefore fast enough for analysis of NGS data.

## Features
- clone frequencies
- domain annotation (frameworks & CDRs)
- sequence quality annotations (frameshift, stop codons)

## Usage

```
python3 -m absequious aln sequences.fa
```

This will produce 2 files, `sequences.fa_reads.csv` and `sequences.fa_summ.csv`.  (A different base filename can be specified with the `--output_base` parameter.)  The "reads" file contains every input read, translated and with regions identified, while the "summary" file contains high-level diversity statistics and unique sequences.

## TODO
- clustering / binning
- liability annotations

## Technologies
- HMMER for domain annotations
