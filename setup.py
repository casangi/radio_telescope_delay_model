import sys
import subprocess
import re
import os
import shutil
from typing import Any
from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.command.develop import develop
from distutils.command.clean import clean
from distutils.core import Command
from setuptools import setup, find_packages

_compile: Any = __import__('_compilers')

if _compile.is_windows:
    DLLNAME = 'libcalc11.dll'
    DLLNAME_WRAP = 'libcalc11_wrapper.dll'
else:
    DLLNAME = 'libcalc11.so'
    DLLNAME_WRAP = 'libcalc11_wrapper.so'

with open('README.md', "r") as fid:   #encoding='utf-8'
    long_description = fid.read()
    

def make_library():

    cmd = 'make -C ./calc11/src/ clean'
    if subprocess.call(cmd, shell=True, env=_compile.compiler_environ) != 0:
        raise IOError('Compilation failed')

    cmd = 'make -C ./calc11/src/ install'
    if subprocess.call(cmd, shell=True, env=_compile.compiler_environ) != 0:
            raise IOError('Compilation failed')


    
class MakeLibrary(Command):
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        make_library()

setup(
    name='radio_telescope_delay_model',
    version='0.0.1',
    description='',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='National Radio Astronomy Observatory',
    author_email='casa-feedback@nrao.edu',
    url='https://github.com/casangi/radio_telescope_delay_model',
    license='Apache-2.0',
    packages=find_packages(),
    install_requires=['numpy>=1.18.1'],
    cmdclass={'build': MakeLibrary}
)
