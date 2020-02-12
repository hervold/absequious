from Bio.Seq import Seq
from enum import Enum


class Dir(Enum):
    fwd = 0
    rev_comp = 1


def translate(dna):
    """
    translate DNA sequence to protein
    raises Bio.Data.CodonTable.TranslationError
    """
    return str(Seq(dna).translate())


def revcomp(dna):
    """reverse compliment"""
    return "".join(
        {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}[b] for b in dna.upper()[::-1]
    )


def translate_six(dna):
    """
    returns six translations, as tulples:
        (direction, offset, translation)
    """
    ret = []
    for rc in (False, True):
        for offset in range(3):
            # to avoid biopython warnings, we need to trim any partial codons, but we don't
            # want to @trim to be set to 0, or we'll get an empty string
            trim = -1 * ((len(dna) - offset) % 3) or None

            ret.append(
                (
                    Dir.rev_comp if rc else Dir.fwd,
                    offset,
                    translate(revcomp(dna)[offset:trim] if rc else dna[offset:trim]),
                )
            )
    return ret
