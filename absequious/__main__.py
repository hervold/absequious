import argparse
import csv
import multiprocessing
import subprocess
import sys
from contextlib import nullcontext
from io import StringIO
from multiprocessing import Pool
from pathlib import Path
from tempfile import NamedTemporaryFile, TemporaryDirectory

from Bio import SeqIO

from . import utils
from .algo import report
from .parse import HMMAln, NoAlignmentFound, annot_fmt

DEFAULT_HMM = Path(utils.get_script_dir()) / "data" / "ighv.hmm"


def trans6(rec, fout):
    t = {}
    for (comp, offset, seq) in utils.translate_six(str(rec.seq)):
        seq_id = f"{rec.id}:{comp.name}:offset_{offset}"
        t[seq_id] = seq
        fout.write(f">{seq_id}\n{seq}\n".encode("utf-8"))
    fout.flush()
    return t


def run_trans6(args):
    with open(args.filename) as fin, open(args.filename + ".trans.fa", "wb") as fout:
        for rec in SeqIO.parse(fin, "fasta"):
            _ = trans6(rec, fout)


def single_pipeline(t):
    rec, temp_dir = t
    with TemporaryDirectory() as temp_dir:
        with NamedTemporaryFile(dir=temp_dir, suffix=".trans.fa") as trans_f:
            translated = trans6(rec, trans_f)
            raw_aln = subprocess.run(
                ["hmmsearch", "--notextw", args.hmm, trans_f.name],
                stdout=subprocess.PIPE,
            )
            try:
                x = HMMAln(StringIO(raw_aln.stdout.decode("utf-8")), translated)
            except NoAlignmentFound:
                return None
            except Exception as e:
                print("~~~~ exception processing", rec.id)
                print(type(e))
                raise
            return x


def run_pipeline(args):
    with TemporaryDirectory() as temp_dir, open(args.input_filename) as fin:
        if not args.no_multiprocess:
            alns = []
            for rec in SeqIO.parse(fin, "fasta"):
                alns.append(single_pipeline((rec, temp_dir)))
        else:
            with Pool(multiprocessing.cpu_count()) as p:
                alns = p.map(
                    single_pipeline,
                    ((rec, temp_dir) for rec in SeqIO.parse(fin, "fasta")),
                )
    fail_ct = sum(1 for x in alns if x is None)
    df = report(alns)
    with (
        nullcontext(sys.stdout)
        if args.output_filename == "-"
        else open(args.output_filename, "w")
    ) as fout:
        df.to_csv(fout, index=False)
        fout.flush()
        fout.write("\nfail_ct: {:d}\n".format(fail_ct))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.set_defaults(func=lambda _: parser.print_help())
    parser.add_argument(
        "-T", "--no-multiprocess", action="store_false", help="disable multiprocessing"
    )

    subparsers = parser.add_subparsers()

    trans_args = subparsers.add_parser(
        "trans6",
        help="helper function: translate DNA Fasta file to 6 protein sequences",
    )
    trans_args.add_argument("filename")
    trans_args.set_defaults(func=run_trans6)

    aln_args = subparsers.add_parser(
        "aln", help="translate DNA and align resulting protein sequences to HMM"
    )
    aln_args.add_argument("input_filename")
    aln_args.add_argument("output_filename")
    aln_args.add_argument("--hmm", default=DEFAULT_HMM)
    aln_args.set_defaults(func=run_pipeline)

    args = parser.parse_args()
    args.func(args)
