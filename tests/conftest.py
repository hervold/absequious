import inspect
import os
import sys
from pathlib import Path
import pytest

sys.path.append(".")
import absequious
import absequious.utils


@pytest.yield_fixture
def hmmsearch_output():
    fname = (
        absequious.utils.get_script_dir().parent
        / "tests"
        / "fixtures"
        / "KY199430_1_hmmsearch.txt"
    )
    with open(fname) as f:
        yield f


# "KY199430_1.fa"
