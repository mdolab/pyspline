# Local modules
from .libspline import version


def fortranVersion():
    """
    This routine returns the git tag and commit hash that the current fortran library was compiled from.
    This may not match the current git version if you have not compiled since your last pull or change to the repo.

    Returns
    -------
    ver : string
        the value returned by ``git describe --dirty --always ---tags`` at compile time of the fortran library.
    """
    ver = version()
    return ver.decode("utf-8").strip()
