from enum import Enum
import re


class ParseError(ValueError):
    pass


class Unreachable(ValueError):
    pass


class AlnState(Enum):
    match_high = 0
    match_low = 1
    mismatch = 2
    insert = 3
    delete = 4


class HMMAln:
    """
    attributes:
    - seq_id: best-matching target
    - 
    """

    _domain_re = re.compile(
        r"== domain .* score: ([-+]?\d*\.?\d+) bits;  conditional E-value: ([-+]?\d*\.?\d+)"
    )

    def __init__(self, input):
        blocks = HMMAln.split_at_boundaries(
            input,
            [
                ("Domain annotation for each sequence", None),
                ("Alignments for each domain", "seq_table"),
                ("Internal pipeline statistics summary", "alignments"),
            ],
        )
        print(blocks)
        self.seq_id, self.best_match = HMMAln.parse_seq_table(blocks["seq_table"])
        self.score_and_eval, self.annots = HMMAln.parse_aln(
            self.seq_id, blocks["alignments"]
        )

    @staticmethod
    def parse_seq_table(block):
        if (
            len(block) < 4
            or not block[0].startswith(">>")
            or not block[1].startswith("#")
            or not block[2].startswith("--")
        ):
            raise ParseError(
                "couldn't parse block starting with 'Domain annotation for each sequence': format"
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
                    ("ali_from", int),
                    ("ali_to", int),
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
        score_and_eval = HMMAln._domain_re.match(block[0]) if block else None
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
                annots.append(AlnState.insert)
            elif tgt.lower() == ref.lower():
                annots.append(AlnState.match_high)
            elif score == "+":
                annots.append(AlnState.match_low)
            elif tgt == "-":
                annots.append(AlnState.delete)
            elif score == " ":
                annots.append(AlnState.mismatch)
            else:
                raise Unreachable(ref, score, tgt)
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
