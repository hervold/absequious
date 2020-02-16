import argparse
import subprocess
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory

from Bio import SeqIO

from . import utils
from .parse import HMMAln

DEFAULT_HMM = Path(utils.get_script_dir()) / "data" / "ighv.hmm"


def trans6(fin, fout):
    for rec in SeqIO.parse(fin, "fasta"):
        for (comp, offset, seq) in utils.translate_six(str(rec.seq)):
            fout.write(">{}:{}:offset_{}\n{}\n".format(rec.id, comp.name, offset, seq))


def run_trans6(args):
    with open(args.filename) as fin, open(args.filename + ".trans.fa", "w") as fout:
        trans6(fin, fout)


def run_pipeline(args):
    with TemporaryDirectory() as temp_dir, open(args.filename) as fin:
        trans_fname = Path(temp_dir) / (Path(args.filename).name + ".trans.fa")
        with open(args.filename) as find, open(trans_fname, "w") as fout:
            trans6(fin, fout)
        raw_aln = subprocess.run(
            ["hmmsearch", args.hmm, trans_fname], stdout=subprocess.PIPE
        )
        aln = HMMAln(StringIO(raw_aln.stdout.decode("utf-8")))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.set_defaults(func=lambda _: parser.print_help())
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
