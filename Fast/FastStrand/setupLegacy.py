#!/usr/bin/env python
from distutils.core import setup, Extension
import os

#=============================================================================
# FastStrand requires:
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
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs, include dirs and libraries =============================
libraryDirs = ["build/"+prod, kcoreLibDir, fastcLibDir]
includeDirs = [numpyIncDir, kcoreIncDir, fastcIncDir]
libraries = ["kcore", "fastc"]
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs    

# Extensions ==================================================================
import srcs
listExtensions = []
listExtensions.append(
    Extension('FastStrand.faststrand',
              sources=['FastStrand/fastStrand.cpp']+srcs.cpp_srcs,
              include_dirs=["FastStrand"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )
    
# setup ======================================================================
setup(
    name="FastStrand",
    version="4.0",
    description="Fast for strand grids.",
    author="ONERA",
    package_dir={"":"."},
    packages=['FastStrand'],
    ext_modules=listExtensions
    )

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
