from distutils.core import setup
import os

os.umask(022)

setup(
    name='gamefft',
    version='0.1.0',
    author='Alioune Schurz',
    author_email='alioune.schurz@gmail.com',
    packages=['gamefft','gamefft.npsdb','gamefft.solvers','gamefft.solvers.cuda','gamefft.test'],
    scripts=['bin/csccp-solver-cli'],
    url='',
    license='LICENSE.txt',
    description='Programs for solving Chemical Substituent Core Combinatorial Problem',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy >= 1.7.1",
        "openbabel-python >= 1.3",
    ],
)

