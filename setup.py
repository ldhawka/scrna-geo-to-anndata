from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="scrna-geo-to-anndata",
    version="0.1.0",
    author="Luvna Dhawka",
    author_email="ldhawka@unc.edu",  
    description="A Python tool for compiling and preprocessing single-cell RNA-seq data from GEO databases into AnnData objects",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ldhawka/scrna-geo-to-anndata",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    keywords="single-cell, RNA-seq, GEO, bioinformatics, scanpy, anndata",
)