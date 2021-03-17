import argparse
import csv
import multiprocessing
import subprocess
import sys
from io import StringIO
from multiprocessing import Pool
from pathlib import Path
from tempfile import NamedTemporaryFile, TemporaryDirectory

from Bio import SeqIO

from . import utils, algo
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
        reader = SeqIO.parse(fin, "fastq") if args.filename.lower().endswith(".fastq") else SeqIO.parse(fin, "fasta")
        for rec in reader:
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
                x = HMMAln(
                    StringIO(raw_aln.stdout.decode("utf-8")), str(rec.seq), translated
                )
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
            reader = SeqIO.parse(fin, "fastq") if args.input_filename.lower().endswith(".fastq") else SeqIO.parse(fin, "fasta")

            with Pool(multiprocessing.cpu_count()) as p:
                alns = p.map(
                    single_pipeline,
                    ((rec, temp_dir) for rec in reader)
                )

    df = algo.report(alns)
    summ = algo.summary(df, alns)
    freq = algo.full_seq_freq(df, alns)

    def dump_m(l, fout):
        for t in l:
            fout.write(",".join(map(str, t)))
            fout.write("\n")

    # default output_base to be the same as input_filename
    # eg: with input_filename "foo.fa" and no output_base specified, outputs "foo.fa_reads.csv"
    # and "foo.fa_summ.csv" are produced
    output_base = args.output_base
    if args.output_base is None:
        output_base = args.input_filename

    if output_base == "-":
        # write everything to stdout
        dump_m(summ, sys.stdout)
        sys.stdout.write("\n")
        dump_m(freq, sys.stdout)
        sys.stdout.write("\n")
        df.to_csv(sys.stdout, index=False)

    else: 
        with open(output_base + "_reads.csv", "w") as fout:
            df.to_csv(fout, index=False)
        with open(output_base + "_summ.csv", "w") as fout:
            dump_m(summ, fout)
            fout.write("\n")
            dump_m(freq, fout)


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
    aln_args.add_argument(
        "--output_base",
        help="root of output filename; given foo, we create foo_cdr3.csv and foo.csv",
    )
    aln_args.add_argument("--hmm", default=DEFAULT_HMM)
    aln_args.set_defaults(func=run_pipeline)

    args = parser.parse_args()
    args.func(args)
