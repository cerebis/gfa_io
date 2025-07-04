from setuptools import find_packages, setup

setup(
    name="gfa_io",
    version="0.3",
    packages=find_packages("src") ,
    package_dir={"": "src"},
    url="https://github.com/cerebis/gfa_io",
    license="GNU GPLv3",
    author="Matthew Z DeMaere",
    author_email="matt.demaere@gmail.com",
    description="A simple package for reading GFA files with a few utilities",
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU Affero General Public License v3",
        "Operating System :: POSIX :: Linux",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 4 - Beta",
    ],
    install_requires=["biopython", "networkx", "pandas", "tqdm"],
    entry_points={
        "console_scripts": ["gfa_utils=gfa_io.command_line:main"],
    },
)
