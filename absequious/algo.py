from . import AlnState
from collections import Counter


def _update_ins_ct(ins_cts, aln_states, pos):
    last_pos = pos
    ct = 0
    for state in aln_states:
        if state == AlnState.insert:
            print("~~~ found ins")
            ct += 1
        else:
            if ct:
                ins_cts[last_pos] = max(ins_cts[last_pos], ct)
                ct = 0
            last_pos += 1
    if ct:
        ins_cts[last_pos] = max(ins_cts[last_pos], ct)


def insert_padding(alignments):
    """
    given a list of alignments, find the maximum number of insertions at each position
    """
    ctr = Counter()
    for aln in alignments:
        _update_ins_ct(ctr, aln.annots, aln.best_match["hmm_from"])
    return ctr
