from . import AlnState, Unreachable
from collections import Counter

# length of each domain in reference space
DOMAIN_LENS = (
    ("H-FR1", 25),
    ("H-CDR1", 8),
    ("H-FR2", 17),
    ("H-CDR2", 8),
    ("H-FR3", 38),
)


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


def multi_aln(padding_ctr, alignments, DOMS=DOMAIN_LENS):
    def split(domain_lens, strip_start, annots, start_pos, acc):
        if domain_lens:
            dom_name, dom_len = domain_lens[0]
        else:
            dom_name, dom_len = "", float("inf")
        if strip_start:
            for i in range(min(strip_start, dom_len if dom_len > 1 else strip_start)):
                acc.extend("-" for _ in range(padding_ctr.get(i, 1)))

            if strip_start > dom_len:
                return [(dom_name, "".join(acc))] + split(
                    domain_lens[1:] if domain_lens else [],
                    strip_start - dom_len,
                    annots,
                    0,
                    [],
                )
            # else: dom_len > strip_start
            dom_len -= strip_start

        for pos in range(start_pos, len(annots)):
            tgt, state = annots[pos]
            if state is AlnState.insert:
                acc.append(tgt)
            elif state is AlnState.delete:
                acc.append("-")
            elif state in (AlnState.match_high, AlnState.match_low):
                acc.append(tgt.upper())
            elif state is AlnState.mismatch:
                acc.append(tgt.lower())
            else:
                raise Unreachable()

            if padding_ctr.get(pos, 1) > 1:
                acc.extend("-" for _ in range(padding_ctr[pos] - 1))

            # handle domain labels
            if state is not AlnState.insert:
                dom_len -= 1
                if dom_len < 1:
                    return [(dom_name, "".join(acc))] + split(
                        domain_lens[1:] if domain_lens else (), 0, annots, pos + 1, []
                    )
        return [(dom_name, "".join(acc))] if acc else []

    return [
        split(DOMS, aln.best_match["hmm_from"], aln.annots, 0, []) for aln in alignments
    ]
