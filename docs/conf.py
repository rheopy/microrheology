import os
import sys

sys.path.insert(0, os.path.abspath(".."))

project = "microrheology"
author = "Marco Caggioni"
release = "0.1.0"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
]

myst_enable_extensions = ["dollarmath", "colon_fence"]
templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_theme = "sphinx_rtd_theme"
autodoc_member_order = "bysource"
