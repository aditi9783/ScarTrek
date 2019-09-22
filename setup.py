import os
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

print("My envs!")
print(os.environ["TRAVIS_TAG"])
print(os.environ["TRAVIS_SUDO"])

setuptools.setup(
    name="aditi9783",
    version="0.0.2",
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
)
