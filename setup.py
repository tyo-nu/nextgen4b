import os
from setuptools import setup, find_packages

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "nextgen4b",
    version = "0.0.1.dev",
    author = "Ted Cybulski",
    author_email = "tcybulski@gmail.com",
    description = (""),
    license = "MIT",
    keywords = "bioinformatics",
    url = "http://github.com/tcyb/nextgen4b",
    packages=find_packages(),
    long_description=read('readme.md'),
    classifiers=[
        "Development Status :: 3 - Alpha"
    ],
)