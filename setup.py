"""Editable install so scripts/ can do `from cdk2moo.x import y`.

    pip install -e .

src-layout: the importable package lives under src/, which means an
un-installed checkout cannot be imported by accident from the repo root.
"""

from setuptools import find_packages, setup

setup(
    name="cdk2moo",
    version="0.1.0",
    description="CDK2 multi-objective molecular optimization proof-of-concept",
    python_requires=">=3.11",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    # Dependencies are managed by environment.yml, not here: most of them
    # (rdkit, vina) must come from conda-forge, not PyPI.
    install_requires=[],
)
