from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="cwi-ssn",
    version="0.1.0",
    author="CWI-SSN Contributors",
    description="Coverage-Weighted Identity-based Sequence Similarity Network",
    long_description=long_description,
    long_description_content_type="text/markdown",

    py_modules=[
        "cwi_ssn",
        "ssn_blastp",
        "ssn_build",
    ],

    install_requires=[
        "biopython<=1.85",
        "pandas>=2.0,<3.0",
        "networkx>=3.0,<4.0",
    ],

    entry_points={
        "console_scripts": [
            "cwi-ssn=cwi_ssn:main",
        ],
    },

    python_requires=">=3.8,<3.13",

    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.12",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)