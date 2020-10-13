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
    
    #print(platform.system())
    #print(os.environ.copy())
    
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
    #print(cmd)
    if subprocess.call(cmd, shell=True) != 0:
        raise IOError('Compilation failed')
                    
    cmd = 'gcc -fPIC -c almaout.c '
    #print(cmd)
    if subprocess.call(cmd, shell=True) != 0:
        raise IOError('Compilation failed')
    
    cmd = 'gfortran -shared -o libcalc11.so almacalc.o  cut1m.o  cwobm.o  ainit.o  csitm.o  adrvr.o  cctiu.o  catiu.o  cdtdb.o  cpepu.o  cnutm.o  cnutu.o  cnutu6.o  cdiuu.o  cm20u.o  crosu.o  cetdm.o  cptdm.o  cocem.o  hardisp.o  dsitu.o  astrm.o  catmm.o  caxom.o  cuvm.o  ctheu.o  cmatu.o  cmabd.o  cvecu.o  dkill.o  almaout.o calc11_wrapper.o'
    #print(cmd)
    if subprocess.call(cmd, shell=True) != 0:
        raise IOError('Compilation failed')
    
    
    if platform.system() == 'Windows':
        subprocess.call(('move libcalc11.so ../lib'), shell=True)
    else:
        subprocess.call(('mv libcalc11.so ../lib'), shell=True)
    os.chdir('../../')
    
    '''
    if get_version('gcc') is None:
        raise IOError('gcc version 4.8.5 or later is required.')
    
    if compare_version(get_version('gcc'),'4.8.5') < 0:
        raise IOError('gcc must be al least version 4.8.5.')
    
    if platform.system() == 'Windows':
        subprocess.call(('rm -rf calc11\lib\*.so'), shell=True)
        os.chdir('.\calc11\src')
        subprocess.call(('del -rf *.o','del -rf *.so'), shell=True)
    else:
        subprocess.call(('rm -rf calc11/lib/*.so'), shell=True)
        os.chdir('./calc11/src')
        subprocess.call(('rm -rf *.o','rm -rf *.so'), shell=True)
    
    COMPILER_F = 'gfortran'
    FLAGS_F = '-ffree-form -ffree-line-length-none -fPIC -I. -g -O2 -c'
    
    SOURCES_F = ['almacalc', 'cut1m', 'cwobm', 'ainit', 'csitm', 'adrvr', 'cctiu', 'catiu', 'cdtdb', 'cpepu', 'cnutm', 'cnutu', 'cnutu6', 'cdiuu', 'cm20u', 'crosu', 'cetdm', 'cptdm', 'cocem', 'hardisp', 'dsitu', 'astrm', 'catmm', 'caxom', 'cuvm', 'ctheu', 'cmatu', 'cmabd', 'cvecu', 'dkill']
    cmd  = COMPILER_F + ' ' + FLAGS_F + ' ' + ' '.join([ s + '.f ' for s in SOURCES_F])
    print(cmd)
    if subprocess.call(cmd, shell=True) != 0:
        raise IOError('Compilation failed')
                    
    COMPILER_C = 'gcc'
    FLAGS_C = '-fPIC -c'
    SOURCES_C = ['almaout']
    cmd = COMPILER_C + ' ' + FLAGS_C + ' '  + ' '.join([s + '.c ' for s in SOURCES_C])
    print(cmd)
    if subprocess.call(cmd, shell=True) != 0:
        raise IOError('Compilation failed')
    
    SOURCES = SOURCES_F + SOURCES_C
    FLAGS = '-shared -o'
    cmd = COMPILER_F + ' ' + FLAGS + ' ' + DLLNAME + ' ' + ' '.join([ s + '.o ' for s in SOURCES])
    print(cmd)
    if subprocess.call(cmd, shell=True) != 0:
        raise IOError('Compilation failed')
                    
    cmd = COMPILER_F + ' -fPIC -shared -o calc11_wrapper.so calc11_wrapper.f90 -L. -lcalc11'
    print(cmd)
    if subprocess.call(cmd, shell=True) != 0:
        raise IOError('Compilation failed')
    
    
    if platform.system() == 'Windows':
        subprocess.call(('move calc11_wrapper.so ../lib'), shell=True)
        subprocess.call(('move libcalc11.so ../lib'), shell=True)
    else:
        subprocess.call(('mv calc11_wrapper.so ../lib'), shell=True)
        subprocess.call(('mv libcalc11.so ../lib'), shell=True)
    
    os.chdir('../../')
    '''
    
    
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
