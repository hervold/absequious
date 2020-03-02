from collections import Counter, OrderedDict
from itertools import zip_longest

import pandas as pd

from . import AlnState, Unreachable

# length of each domain in reference space
DOMAIN_LENS = (
    ("H-FR1", 25),
    ("H-CDR1", 8),
    ("H-FR2", 17),
    ("H-CDR2", 8),
    ("H-FR3", 39),
    ("H-CDR3", 15),
    ("H-FR4", 11),
)

# FIXME: grab this from hmmsearch output
HMM_LEN = 122


def _update_ins_ct(ins_cts, aln_states, pos):
    last_pos = pos
    ct = 0
    for _, state in aln_states:
        if state == AlnState.insert:
            ct += 1
        else:
            if ct:
                ins_cts[last_pos] = max(ins_cts[last_pos], ct)
                ct = 0
            last_pos += 1
    if ct:
        ins_cts[last_pos] = max(ins_cts[last_pos], ct)


def _insert_padding(ctr, alignments):
    """
    given a list of alignments, find the maximum number of insertions at each position
    """
    for aln in alignments:
        _update_ins_ct(ctr, aln.annots, aln.best_match["hmm_from"])
    return ctr


def insert_padding(alignments):
    return _insert_padding(Counter(), alignments)


def split_and_pad(padding_ctr, aln, DOMS=DOMAIN_LENS):
    """
    given @padding_ctr and HMMAln @aln, split the aligned sequence at domain boundaries and
    insert padding
    """

    def split(domain_lens, strip_start, annots, start_pos, acc):
        if domain_lens:
            dom_name, dom_len = domain_lens[0]
        else:
            dom_name, dom_len = "", float("inf")
        if strip_start:
            if padding_ctr is not None:
                for i in range(
                    min(strip_start, dom_len if dom_len > 1 else strip_start)
                ):
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

            if padding_ctr is not None and padding_ctr.get(pos, 1) > 1:
                acc.extend("-" for _ in range(padding_ctr[pos] - 1))

            # handle domain labels
            if state is not AlnState.insert:
                dom_len -= 1
                if dom_len < 1:
                    return [(dom_name, "".join(acc))] + split(
                        domain_lens[1:] if domain_lens else (), 0, annots, pos + 1, []
                    )
        return [(dom_name, "".join(acc))] if acc else []

    return split(DOMS, aln.best_match["hmm_from"], aln.annots, 0, [])


def multi_aln(padding_ctr, alignments, DOMS=DOMAIN_LENS):
    return [split_and_pad(padding_ctr, aln, DOMS) for aln in alignments]


def has_frameshift(aln):
    """
    heuristics for determining whether a frameshift is present.

    TODO: we should incorporate quality scores, if available
    """
    if aln.best_match["hmm_from"] > 2 and aln.best_match["tgt_from"] > 2:
        return True
    if (
        aln.best_match["hmm_to"] < HMM_LEN - 1
        and aln.best_match["tgt_to"] < aln.tgt_len - 1
    ):
        return True
    return False


def has_stop(aln):
    return "*" in aln.tgt_seq[aln.best_match["tgt_from"] : aln.best_match["tgt_to"]]


def report(alignments):
    padding_by_pos = None  # insert_padding([x for x in alignments if x is not None])
    padded_alns = [
        None if aln is None else split_and_pad(padding_by_pos, aln)
        for aln in alignments
    ]
    _df = OrderedDict(
        (k, [])
        for k in ["read"]
        + [nom for nom, _ in DOMAIN_LENS]
        + ["complete?", "frameshift?", "stop?",]
    )
    for aln, padded_aln in zip(alignments, padded_alns):
        if aln is not None:
            _df["read"].append(aln.seq_id)
            for (dom_name, _), (_, padded_seq) in zip_longest(
                DOMAIN_LENS, padded_aln, fillvalue=("", "")
            ):
                _df[dom_name].append(padded_seq)
            _df["complete?"].append(len(padded_aln) == 7)
            _df["frameshift?"].append(has_frameshift(aln))
            _df["stop?"].append(has_stop(aln))

    return pd.DataFrame(_df)


def summary(report, alns):
    failed = sum(1 for x in alns if x is None)
    tot = len(report) + failed
    return (
        ("failed", failed / (len(report) + failed), failed, tot),
        ("complete", report["complete?"].sum() / tot, report["complete?"].sum(), tot,),
        (
            "frameshift",
            report["frameshift?"].sum() / tot,
            report["frameshift?"].sum(),
            tot,
        ),
        ("stop_codon", report["stop?"].sum() / tot, report["stop?"].sum(), tot),
    )


def cdr3_freq(report, alns):
    """
    return a list of tuples of the form:
      (cdr3_sequene, fraction, count, total_reads)
    CDR3s with only 1 read are dropped
    """
    failed = sum(1 for x in alns if x is None)
    tot = len(report) + failed
    cdr3_cts = []
    report.groupby("H-CDR3").apply(
        lambda x: cdr3_cts.append(
            (x["H-CDR3"].values[0], x.count()[0] / tot, x.count()[0], tot)
        )
    )
    cdr3_cts.sort(key=(lambda row: row[1]), reverse=True)
    return [t for t in cdr3_cts if t[2] > 1]
