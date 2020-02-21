import argparse
import subprocess
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory, NamedTemporaryFile

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


def single_pipeline(rec, temp_dir):
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
    alns = []
    with TemporaryDirectory() as temp_dir, open(args.filename) as fin:
        for rec in SeqIO.parse(fin, "fasta"):
            alns.append(single_pipeline(rec, temp_dir))
            print("///", alns[-1].best_match)
    print("~~~ num:", len(alns))
    padding_by_pos = insert_padding(alns)
    print(alns[0]._blocks)
    print("~~~", padding_by_pos)
    for s in multi_aln(padding_by_pos, alns):
        print(s)


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
