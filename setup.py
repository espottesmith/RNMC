#!/usr/bin/env python

from pathlib import Path

from setuptools import find_packages, setup

module_dir = Path(__file__).resolve().parent
base_dir = module_dir.parent.parent

with open(base_dir / "README.md") as f:
    long_desc = f.read()

setup(
    name="rnmc",
    setup_requires=["setuptools"],
    description="",
    long_description=long_desc,
    long_description_content_type="text/markdown",
    url="https://github.com/danielbarter/RNMC",
    author="Daniel Barter",
    author_email="danielbarter@gmail.com",
    license="modified BSD",
    packages=find_packages("src"),
    package_dir={"": "src"},
    zip_safe=False,
    install_requires=[
        "setuptools",
    ],
    tests_require=["pytest"],
    python_requires=">=3.6",
    version='0.1'
)