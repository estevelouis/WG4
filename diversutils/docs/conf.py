from datetime import datetime

extensions = [
    "sphinx.ext.extlinks",
    "sphinx.ext.intersphinx"
]
templates_path = ["_templates"]
source_suffix = ".rst"
master_doc = "index"

project = "DiversUtils"
copyright = f"{datetime.now().year} Université Paris-Saclay"

exclude_patterns = ["_build"]

extlinks = {
    "git_tag": ("https://github.com/estevelouis/WG4/tree/%s", "%s"),
    "bug": ("https://github.com/estevelouis/WG4/issues/%s", "#%s"),
    "feature": ("https://github.com/estevelouis/WG4/issues/%s", "#%s"),
    "issue": ("https://github.com/estevelouis/WG4/issues/%s", "#%s"),
}

html_theme = "alabaster"
html_sidebars = {
    "**": [
        "about.html",
        "navigation.html",
        "relations.html",
        "searchbox.html",
        # "donate.html",
    ]
}
html_theme_options = {
    "description": "Measuring diversity in linguistic data",
    "github_user": "estevelouis",
    "github_repo": "WG4",
    "fixed_sidebar": True,
    # "tidelift_url": "https://tidelift.com/subscription/pkg/pypi-alabaster?utm_source=pypi-alabaster&utm_medium=referral&utm_campaign=docs",  # noqa
}
