import os

from pkg_resources import resource_filename

__all__ = ["get_data_file"]


def get_data_file(relative_path):
    """Get the full path to one of the reference files.

    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the repex folder).

    """
    fn = resource_filename('bgmol.systems', relative_path)

    if not os.path.exists(fn):
        raise ValueError("Sorry! %s does not exist. If you just added it, you'll have to re-install" % fn)

    return fn