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
    #Check if windows work
    
    subprocess.call(('rm -rf calc11/lib/*.so'), shell=True)
    os.chdir('./calc11/src')
    subprocess.call(('rm -rf *.o','rm -rf *.so'), shell=True)
    
    COMPILER_F = 'gfortran'
    FLAGS_F = '-ffree-form -ffree-line-length-none -fPIC -I. -g -O2 -c'
    
    SOURCES_F = ['almacalc', 'cut1m', 'cwobm', 'ainit', 'csitm', 'adrvr', 'cctiu', 'catiu', 'cdtdb', 'cpepu', 'cnutm', 'cnutu', 'cnutu6', 'cdiuu', 'cm20u', 'crosu', 'cetdm', 'cptdm', 'cocem', 'hardisp', 'dsitu', 'astrm', 'catmm', 'caxom', 'cuvm', 'ctheu', 'cmatu', 'cmabd', 'cvecu', 'dkill']
    cmd  = COMPILER_F + ' ' + FLAGS_F + ' ' + ' '.join([ s + '.f ' for s in SOURCES_F])
    os.system('pwd')
    if subprocess.call(cmd, shell=True) != 0:
                    raise IOError('Compilation failed')
                    
    COMPILER_C = 'gcc'
    FLAGS_C = '-fPIC -c'
    SOURCES_C = ['almaout']
    cmd = COMPILER_C + ' ' + FLAGS_C + ' '  + ' '.join([s + '.c ' for s in SOURCES_C])
    if subprocess.call(cmd, shell=True) != 0:
                    raise IOError('Compilation failed')
    
    SOURCES = SOURCES_F + SOURCES_C
    FLAGS = '-shared -o'
    cmd = COMPILER_F + ' ' + FLAGS + ' ' + DLLNAME + ' ' + ' '.join([ s + '.o ' for s in SOURCES])
    if subprocess.call(cmd, shell=True) != 0:
                    raise IOError('Compilation failed')
                    
    cmd = COMPILER_F + ' -fPIC -shared -o calc11_wrapper.so calc11_wrapper.f90 -L. -lcalc11'
    if subprocess.call(cmd, shell=True) != 0:
                    raise IOError('Compilation failed')
                    
    subprocess.call(('mv calc11_wrapper.so ../lib'), shell=True)
    subprocess.call(('mv libcalc11.so ../lib'), shell=True)
    
    os.chdir('../../')
    
#    cmd = 'make -C ./calc11/src/ clean'
#    if subprocess.call(cmd, shell=True, env=_compile.compiler_environ) != 0:
#        raise IOError('Compilation failed')
#
#    cmd = 'make -C ./calc11/src/ install'
#    if subprocess.call(cmd, shell=True, env=_compile.compiler_environ) != 0:
#            raise IOError('Compilation failed')

    
class MakeLibrary(Command):
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        make_library()
        
class SharedLibrary(build_py):

    def run(self):
        make_library()
        build_py.run(self)

class DevelopLibrary(develop):

    def run(self):
        make_library()
        develop.run(self)

setup(
    name='radio_telescope_delay_model',
    version='0.0.1',
    description='Radio Telescope Delay Model',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='National Radio Astronomy Observatory',
    author_email='casa-feedback@nrao.edu',
    url='https://github.com/casangi/radio_telescope_delay_model',
    license='Apache-2.0',
    packages=find_packages(),
    install_requires=['numpy>=1.18.1'],
    cmdclass={'build_py': SharedLibrary,'make': MakeLibrary,'develop': DevelopLibrary}
)
