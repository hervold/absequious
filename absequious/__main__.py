import argparse
import subprocess
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory, NamedTemporaryFile
import csv
from multiprocessing import Pool
import multiprocessing
from Bio import SeqIO

from . import utils
from .parse import HMMAln
from .algo import insert_padding, multi_aln

DEFAULT_HMM = Path(utils.get_script_dir()) / "data" / "ighv.hmm"


def trans6(rec, fout):
    for (comp, offset, seq) in utils.translate_six(str(rec.seq)):
        fout.write(
            ">{}:{}:offset_{}\n{}\n".format(rec.id, comp.name, offset, seq).encode(
                "utf-8"
            )
        )
    fout.flush()


def run_trans6(args):
    with open(args.filename) as fin, open(args.filename + ".trans.fa", "wb") as fout:
        for rec in SeqIO.parse(fin, "fasta"):
            trans6(rec, fout)


def single_pipeline(t):
    rec, temp_dir = t
    with TemporaryDirectory() as temp_dir, open(args.filename) as fin:
        with NamedTemporaryFile(dir=temp_dir, suffix=".trans.fa") as trans_f:
            trans6(rec, trans_f)
            raw_aln = subprocess.run(
                ["hmmsearch", args.hmm, trans_f.name], stdout=subprocess.PIPE
            )
            x = HMMAln(StringIO(raw_aln.stdout.decode("utf-8")))
            print(x.seq_id)
            return x


def run_pipeline(args):
    with TemporaryDirectory() as temp_dir, open(args.filename) as fin:
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
        print("///", alns[-1].best_match)
    padding_by_pos = insert_padding(alns)
    for doms, aln in zip(multi_aln(padding_by_pos, alns), alns):
        padded_seq = "|".join(s for _, s in doms)
        print(
            f'{aln.best_match["hmm_from"]:03d}:{aln.best_match["hmm_to"]:03d}  {padded_seq}'
        )


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
    aln_args.add_argument("filename")
    aln_args.add_argument("--hmm", default=DEFAULT_HMM)
    aln_args.set_defaults(func=run_pipeline)

    args = parser.parse_args()
    args.func(args)
