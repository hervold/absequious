from enum import Enum


class AlnState(Enum):
    match_high = 0
    match_low = 1
    mismatch = 2
    insert = 3
    delete = 4


def get_script_dir(follow_symlinks=True):
    """
    https://stackoverflow.com/questions/3718657/how-to-properly-determine-current-script-directory/22881871#22881871
    """
    if getattr(sys, "frozen", False):  # py2exe, PyInstaller, cx_Freeze
        path = os.path.abspath(sys.executable)
    else:
        path = inspect.getabsfile(get_script_dir)
    if follow_symlinks:
        path = os.path.realpath(path)
    return os.path.dirname(path)
