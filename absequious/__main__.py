import argparse
from . import utils
from Bio import SeqIO


def trans6(args):
    with open(args.filename) as fin, open(args.filename + ".trans.fa", "w") as fout:
        for rec in SeqIO.parse(fin, "fasta"):
            for (comp, offset, seq) in utils.translate_six(str(rec.seq)):
                fout.write(
                    ">{}:{}:offset_{}\n{}\n".format(rec.id, comp.name, offset, seq)
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    trans_args = subparsers.add_parser(
        "trans6",
        help="helper function: translate DNA Fasta file to 6 protein sequences",
    )
    trans_args.add_argument("filename")
    trans_args.set_defaults(func=trans6)

    aln_args = subparsers.add_parser("aln")

    args = parser.parse_args()
    args.func(args)
