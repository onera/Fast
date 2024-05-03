#!/usr/bin/env python
from distutils.core import setup, Extension
import os

#=============================================================================
# FastS requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Fortran compiler: defined in config.py
# Numpy
# KCore library
#=============================================================================

# Write setup.cfg file
import KCore.Dist as Dist
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

# Test if fastc exists =======================================================
(fastcVersion, fastcIncDir, fastcLibDir) = Dist.checkFastC()

from KCore.config import *

# Compilation des fortrans ====================================================
if f77compiler == "None":
    print("Error: a fortran 77 compiler is required for compiling FastS.")
os.system('./MakeLinks')
args = Dist.getForArgs(); opt = ''
for c, v in enumerate(args): opt += 'FOPT'+str(c)+'='+v+' '
if f77compiler == 'gfortran': opt += 'FOPT3=-fdefault-real-8'

os.system("make -e FC="+f77compiler+" WDIR=FastS/Metric "+opt)
os.system("make -e FC="+f77compiler +" WDIR=FastS/Compute "+opt)
os.system('./RmLinks')
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs, include dirs and libraries =============================
libraryDirs = ["build/"+prod, kcoreLibDir, fastcLibDir]
includeDirs = [numpyIncDir, kcoreIncDir, fastcIncDir]
libraries = ["MetricF", "ComputeF", "kcore", "fastc"]
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs    

# Extensions ==================================================================
import srcs
listExtensions = []
listExtensions.append(
    Extension('FastS.fasts',
              sources=['FastS/fastS.cpp']+srcs.cpp_srcs,
              include_dirs=["FastS"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )
    
# setup ======================================================================
setup(
    name="FastS",
    version="4.0",
    description="Fast for structured grids.",
    author="Onera",
    package_dir={"":"."},
    packages=['FastS'],
    ext_modules=listExtensions
    )

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
