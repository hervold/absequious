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

## TODO
- clustering / binning
- liability annotations

## Technologies
- HMMER for domain annotations
