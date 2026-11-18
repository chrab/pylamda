"""pylamda: read, write, verify, compare, combine and fix files in the LAMDA format.

The format is that of the Leiden Atomic and Molecular Database (LAMDA,
https://home.strw.leidenuniv.nl/~moldata/molformat.html). pylamda follows the standard
strictly: collision partners are identified by their ID only, and deviations are
reported as issues, never silently interpreted.

>>> import pylamda
>>> data = pylamda.read("co.dat")                    # doctest: +SKIP
>>> data.issues                                      # doctest: +SKIP
>>> pylamda.verify(data)                             # doctest: +SKIP
>>> pylamda.write(data, "co_normalised.dat")         # doctest: +SKIP
"""

from importlib.metadata import PackageNotFoundError, version

from .compare import LamdaDiff, LevelMatch, PartnerDiff, compare, critical_density, match_levels
from .issues import CODES, SEVERITIES, Issue, LamdaReadError
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
from .verify import verify
from .writer import write

try:
    __version__ = version("pylamda")
except PackageNotFoundError:  # not installed
    __version__ = "0+unknown"

__all__ = [
    "CODES",
    "PARTNERS",
    "PARTNER_IDS",
    "PartnerDiff",
    "SEVERITIES",
    "STANDARD_PARTNERS",
    "BlockSource",
    "CollisionSet",
    "Issue",
    "LevelMatch",
    "LamdaData",
    "LamdaDiff",
    "LamdaReadError",
    "SourceMap",
    "__version__",
    "compare",
    "critical_density",
    "find_dois",
    "match_levels",
    "partner_name",
    "read",
    "split_partner_line",
    "verify",
    "write",
]
