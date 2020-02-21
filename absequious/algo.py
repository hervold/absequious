from . import AlnState, Unreachable
from collections import Counter


def _update_ins_ct(ins_cts, aln_states, pos):
    last_pos = pos
    ct = 0
    for state in aln_states:
        if state == AlnState.insert:
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


def multi_aln(padding_ctr, alignments):
    all_alns = []
    for aln in alignments:
        padded = []
        # pad the start
        for i in range(aln.best_match["hmm_from"]):
            padded.extend("-" for _ in range(max(1, padding_ctr.get(i, 1))))
        ref = 0
        for tgt, state in aln.annots:
            if state is AlnState.insert:
                padded.append(tgt)
            elif state is AlnState.delete:
                padded.append("-")
            elif state in (AlnState.match_high, AlnState.match_low):
                padded.append(tgt.upper())
            elif state is AlnState.mismatch:
                padded.append(tgt.lower())
            else:
                raise Unreachable()
        all_alns.append("".join(padded))
    return all_alns
