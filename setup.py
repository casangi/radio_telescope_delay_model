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
import platform

def call_command(cmd):
    try:
        return subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT).decode().strip()
    except:
        return None

def get_version(compiler):
    ver = call_command('gfortran -dumpversion')
    if ver and '.' not in ver:
        ver = call_command(compiler + ' -dumpfullversion')
    return ver
    
import re

def compare_version(version1, version2):
    def normalize(v):
        return [int(x) for x in re.sub(r'(\.0+)*$','', v).split(".")]
        
    return ((normalize(version1) > normalize(version2)) - (normalize(version1) < normalize(version2)))
    
if platform.system() == 'Window':
    DLLNAME = 'libcalc11.dll'
    DLLNAME_WRAP = 'libcalc11_wrapper.dll'
else:
    DLLNAME = 'libcalc11.so'
    DLLNAME_WRAP = 'libcalc11_wrapper.so'

with open('README.md', "r") as fid:   #encoding='utf-8'
    long_description = fid.read()
    

def make_library():

    #print(get_version('gcc'))
    #print(get_version('gfortran'))
    
    if get_version('gcc') is None:
        raise IOError('gcc version 4.8.5 or later is required.')
    
    if compare_version(get_version('gcc'),'4.8.5') < 0:
        raise IOError('gcc must be al least version 4.8.5.')
    
    if platform.system() == 'Windows':
        subprocess.call(('del ./calc11/lib/*.so'), shell=True)
        subprocess.call(('del ./calc11/src/*.o'), shell=True)
        subprocess.call(('del ./calc11/src/*.so'), shell=True)
    else:
        subprocess.call(('rm -rf ./calc11/lib/*.so'), shell=True)
        subprocess.call(('rm -rf ./calc11/src/*.o'), shell=True)
        subprocess.call(('rm -rf ./calc11/src/*.so'), shell=True)
        
    os.chdir('./calc11/src')
    
    cmd  = 'gfortran -ffree-form -ffree-line-length-none -fPIC -I. -g -O2 -c almacalc.f  cut1m.f  cwobm.f  ainit.f  csitm.f  adrvr.f  cctiu.f  catiu.f  cdtdb.f  cpepu.f  cnutm.f  cnutu.f  cnutu6.f  cdiuu.f  cm20u.f  crosu.f  cetdm.f  cptdm.f  cocem.f  hardisp.f  dsitu.f  astrm.f  catmm.f  caxom.f  cuvm.f  ctheu.f  cmatu.f  cmabd.f  cvecu.f  dkill.f calc11_wrapper.f90'
    if subprocess.call(cmd, shell=True) != 0:
        raise IOError('Compilation failed')
                    
    cmd = 'gcc -fPIC -c almaout.c '
    if subprocess.call(cmd, shell=True) != 0:
        raise IOError('Compilation failed')
    
    cmd = 'gfortran -shared -o libcalc11.so almacalc.o  cut1m.o  cwobm.o  ainit.o  csitm.o  adrvr.o  cctiu.o  catiu.o  cdtdb.o  cpepu.o  cnutm.o  cnutu.o  cnutu6.o  cdiuu.o  cm20u.o  crosu.o  cetdm.o  cptdm.o  cocem.o  hardisp.o  dsitu.o  astrm.o  catmm.o  caxom.o  cuvm.o  ctheu.o  cmatu.o  cmabd.o  cvecu.o  dkill.o  almaout.o calc11_wrapper.o'
    if subprocess.call(cmd, shell=True) != 0:
        raise IOError('Compilation failed')
    
    
    if platform.system() == 'Windows':
        subprocess.call(('move libcalc11.so ../lib'), shell=True)
    else:
        subprocess.call(('mv libcalc11.so ../lib'), shell=True)
    os.chdir('../../')

    
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
    version='0.0.4',
    description='Radio Telescope Delay Model',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='National Radio Astronomy Observatory',
    author_email='casa-feedback@nrao.edu',
    url='https://github.com/casangi/radio_telescope_delay_model',
    license='Apache-2.0',
    packages=find_packages(),
    install_requires=['numpy>=1.18.1'],
    include_package_data=True,
    zip_safe=False,
    cmdclass={'build_py': SharedLibrary,'make': MakeLibrary,'develop': DevelopLibrary},
    python_requires='>=3.6'
)
