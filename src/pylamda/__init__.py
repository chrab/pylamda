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

from .bib import BibEntry, BibResult, bibtex, collect_references
from .combine import (
    AddedPartner,
    CombineError,
    DroppedPartner,
    ScaledPartner,
    add_partner,
    drop_partner,
    patch_scale_partner,
    scale_partner,
)
from .compare import LamdaDiff, LevelMatch, PartnerDiff, compare, critical_density, match_levels
from .fix import FixError, FixRecord, add_reference, fix, patch
from .issues import CODES, SEVERITIES, Issue, LamdaReadError, NotLamdaError
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
    "AddedPartner",
    "BibEntry",
    "BibResult",
    "BlockSource",
    "CollisionSet",
    "CombineError",
    "DroppedPartner",
    "FixError",
    "FixRecord",
    "Issue",
    "LevelMatch",
    "LamdaData",
    "LamdaDiff",
    "LamdaReadError",
    "NotLamdaError",
    "ScaledPartner",
    "SourceMap",
    "__version__",
    "add_partner",
    "add_reference",
    "bibtex",
    "collect_references",
    "compare",
    "critical_density",
    "drop_partner",
    "find_dois",
    "fix",
    "match_levels",
    "partner_name",
    "patch",
    "read",
    "patch_scale_partner",
    "scale_partner",
    "split_partner_line",
    "verify",
    "write",
]
