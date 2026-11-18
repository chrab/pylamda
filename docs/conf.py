"""Sphinx configuration of the pylamda documentation.

Build locally (from the repository root)::

    sphinx-build -W -b html docs docs/_build/html
"""

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as get_version

project = "pylamda"
author = "Christian Rab"
copyright = "2026, Christian Rab"
try:
    release = get_version("pylamda")
except PackageNotFoundError:  # built from a source tree without installing
    release = "0+unknown"
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx_rtd_theme",
    "sphinx_copybutton",
    "myst_parser",
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
root_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
language = "en"

# --- autodoc / napoleon (NumPy docstrings, types are written in the docstrings) ---
autodoc_typehints = "none"
autodoc_member_order = "bysource"
autodoc_default_options = {"members": True, "show-inheritance": True}
autodoc_class_signature = "separated"
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_ivar = True  # dataclass fields as :ivar:, no duplicate attribute entries
napoleon_use_param = True
napoleon_use_rtype = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
}

# --- Markdown (docs/design.md) ---
myst_heading_anchors = 3

# --- HTML ---
html_theme = "sphinx_rtd_theme"
html_theme_options = {"navigation_depth": 3}
html_show_sourcelink = False
html_title = f"pylamda {release}"
