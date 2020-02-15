from absequious.parse import HMMAln


def test_parse(hmmsearch_output):
    x = HMMAln(hmmsearch_output)
