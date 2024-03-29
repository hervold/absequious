import re
from enum import Enum

from . import AlnState, Unreachable
from .utils import revcomp


class ParseError(ValueError):
    pass


class NoAlignmentFound(ValueError):
    pass


class HMMAln:
    """
    attributes:
    - seq_id: best-matching target
    - 
    """

    _domain_re = re.compile(
        r"== domain .* score: ([-+]?\d*\.?\d+) bits;  conditional E-value: ([-+]?\d*\.?\d+)"
    )

    def __init__(self, input, dna_seq, translated):
        blocks = HMMAln.split_at_boundaries(
            input,
            [
                ("Domain annotation for each sequence", None),
                ("Alignments for each domain", "seq_table"),
                ("Internal pipeline statistics summary", "alignments"),
            ],
        )
        self._blocks = blocks
        self.seq_id, self.best_match = HMMAln.parse_seq_table(blocks["seq_table"])

        self.tgt_seq = translated[self.seq_id]
        comp, offs = self.seq_id.split(":")[-2:]
        if not offs.startswith("offset_"):
            raise ParseError(f'WTF offset "{offs}"')
        self.dna_seq = (dna_seq if comp == "fwd" else revcomp(dna_seq))[int(offs[7:]) :]
        self.tgt_len = len(self.tgt_seq)

        self.score_and_eval, self.annots = HMMAln.parse_aln(
            self.seq_id, blocks["alignments"]
        )

    @staticmethod
    def parse_seq_table(block):
        if (
            block
            and block[0] == "[No targets detected that satisfy reporting thresholds]"
        ):
            raise NoAlignmentFound()
        if (
            len(block) < 4
            or not block[0].startswith(">>")
            or not block[1].startswith("#")
            or not block[2].startswith("--")
        ):
            raise ParseError(
                "couldn't parse block starting with 'Domain annotation for each sequence': format",
                block,
            )
        seq_id = block[0][3:].strip()
        iden = lambda x: x
        best_match = dict(
            (label, conv(value))
            for (label, conv), value in zip(
                (
                    (None, iden),
                    (None, iden),
                    ("score", float),
                    ("bias", float),
                    ("c-Evalue", float),
                    ("i-Evalue", float),
                    ("hmm_from", int),
                    ("hmm_to", int),
                    (None, iden),
                    ("tgt_from", int),
                    ("tgt_to", int),
                    (None, iden),
                    ("env_from", int),
                    ("env_to", int),
                    (None, iden),
                    ("acc", float),
                ),
                block[3].split(),
            )
            if label
        )
        if len(best_match) != 11:
            raise ParseError(
                "couldn't parse block starting with 'Domain annotation for each sequence': table"
            )
        return seq_id, best_match

    @staticmethod
    def parse_aln(seq_id, block):
        score_and_eval = HMMAln._domain_re.match(block[0]).groups() if block else None
        if len(block) < 5 or not score_and_eval or not block[3].startswith(seq_id):
            raise ParseError(
                "couldn't parse block starting with 'Domain annotation for each sequence': format"
            )
        _, ref_st, ref_guide, ref_end = block[1].split()
        scores = block[2].strip()
        tgt_id, tgt_st, tgt_guide, tgt_end = block[3].split()
        if seq_id != tgt_id:
            raise ParseError(
                "inconsistent target names: '{}' != '{}'".format(tgt_id, seq_id)
            )
        annots = []
        for ref, score, tgt in zip(ref_guide, scores, tgt_guide):
            if ref == ".":
                annots.append((tgt, AlnState.insert))
            elif tgt.lower() == ref.lower():
                annots.append((tgt, AlnState.match_high))
            elif score == "+":
                annots.append((tgt, AlnState.match_low))
            elif tgt == "-":
                annots.append((tgt, AlnState.delete))
            else:
                annots.append((tgt, AlnState.mismatch))
        return score_and_eval, annots

    @staticmethod
    def split_at_boundaries(input, boundaries):
        """iterate over input, splitting at boundary strings.  return a dictionary mapping boundary name to subsequent lines of text"""

        # make a mutable copy
        boundaries = boundaries[:]

        ret = {}
        curr_bound, curr_nom, buff = "", None, []
        for line_ in input:
            line = line_.strip()
            if line.startswith(curr_bound):
                if curr_nom:
                    ret[curr_nom] = buff
                buff = []
                try:
                    curr_bound, curr_nom = boundaries.pop(0)
                except IndexError:
                    break
            else:
                if line:
                    buff.append(line)
        if buff and curr_nom:
            ret[curr_nom] = buff
        return ret


def annot_fmt(annots):
    return "".join(
        {
            AlnState.mismatch: "x",
            AlnState.insert: "i",
            AlnState.match_high: "M",
            AlnState.match_low: "M",
            AlnState.delete: "d",
        }[a]
        for _, a in annots
    )
