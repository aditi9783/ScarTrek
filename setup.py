import os
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

version = os.getenv("TRAVIS_TAG", "0.0.1")  # default if TRAVIS_TAG doesn't exist ...
version = version if version else "0.0.1"  # ... or if TRAVIS_TAG is blank
version = "0.0.3"  # override for local branch testing
print("Using version %s" % version)

setuptools.setup(
    name="aditi9783",
    version=version,
    author="Aditi Gupta",
    author_email="aditi9783@gmail.com",
    description="Python package for identifying gene disrupting and restoring indels in whole-genome sequencing data of Mycobacterium tuberculosis samples",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/aditi9783/ScarTrek",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=2.7',
    install_requires=[
        # 'bowtie2',
        # 'pysam'
        # 'bcftools'
    ],
    tests_require=[
        'nose'
    ],
    entry_points={
        'console_scripts': ['find-scars=scartrek.find_scars:main'],
    },
)
