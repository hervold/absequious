import inspect
import os
import sys
from pathlib import Path
import pytest

sys.path.append(".")
import absequious


def get_script_dir(follow_symlinks=True):
    if getattr(sys, "frozen", False):  # py2exe, PyInstaller, cx_Freeze
        path = os.path.abspath(sys.executable)
    else:
        path = inspect.getabsfile(get_script_dir)
    if follow_symlinks:
        path = os.path.realpath(path)
    return os.path.dirname(path)


@pytest.yield_fixture
def hmmsearch_output():
    fname = Path(get_script_dir()) / "fixtures" / "KY199430_1_hmmsearch.txt"
    with open(fname) as f:
        yield f


# "KY199430_1.fa"
