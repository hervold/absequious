import sys
from enum import Enum
import inspect
import os
from Bio.Seq import Seq
from pathlib import Path


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


def get_script_dir(follow_symlinks=True):
    """
    https://stackoverflow.com/questions/3718657/how-to-properly-determine-current-script-directory/22881871#22881871
    """
    if getattr(sys, "frozen", False):  # py2exe, PyInstaller, cx_Freeze
        path = Path(sys.executable).resolve()
    else:
        path = Path(inspect.getabsfile(get_script_dir))
    if follow_symlinks:
        path = path.resolve()
    return path.parent.parent
