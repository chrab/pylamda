"""pylamda: read, write, verify, compare, combine and fix files in the LAMDA format.

The format is that of the Leiden Atomic and Molecular Database (LAMDA,
https://home.strw.leidenuniv.nl/~moldata/molformat.html). pylamda follows the standard
strictly: collision partners are identified by their ID only, and deviations are
reported as issues, never silently interpreted.

>>> import pylamda
>>> data = pylamda.read("co.dat")                    # doctest: +SKIP
>>> data.issues                                      # doctest: +SKIP
>>> pylamda.write(data, "co_normalised.dat")         # doctest: +SKIP
"""

from importlib.metadata import PackageNotFoundError, version

from .issues import CODES, Issue, LamdaReadError
from .model import (
    PARTNER_IDS,
    PARTNERS,
    STANDARD_PARTNERS,
    BlockSource,
    CollisionSet,
    LamdaData,
    SourceMap,
    partner_name,
    split_partner_line,
)
from .reader import read
from .references import find_dois
from .writer import write

try:
    __version__ = version("pylamda")
except PackageNotFoundError:  # not installed
    __version__ = "0+unknown"

__all__ = [
    "CODES",
    "PARTNERS",
    "PARTNER_IDS",
    "STANDARD_PARTNERS",
    "BlockSource",
    "CollisionSet",
    "Issue",
    "LamdaData",
    "LamdaReadError",
    "SourceMap",
    "__version__",
    "find_dois",
    "partner_name",
    "read",
    "split_partner_line",
    "write",
]
