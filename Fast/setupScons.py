#!/usr/bin/env python
from distutils.core import setup, Extension
import os

#=============================================================================
# Fast requires:
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

# Test if connector exists =======================================================
(connectorVersion, connectorIncDir, connectorLibDir) = Dist.checkConnector()

from KCore.config import *

# Test if libmpi exists ======================================================
(mpi, mpiIncDir, mpiLibDir, mpiLibs) = Dist.checkMpi(additionalLibPaths, additionalIncludePaths)

# Compilation des fortrans ====================================================
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs, include dirs and libraries =============================
libraryDirs = ["build/"+prod, kcoreLibDir,connectorLibDir]
includeDirs = [numpyIncDir, kcoreIncDir,connectorIncDir]
libraries = ["fast", "kcore","connector"]
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

ADDITIONALCPPFLAGS=[]
if mpi:
    libraryDirs.append(mpiLibDir)
    includeDirs.append(mpiIncDir)
    ADDITIONALCPPFLAGS = ['-D_MPI']
    libraries += mpiLibs
    
# Extensions ==================================================================
listExtensions = []
listExtensions.append(
    Extension('Fast.fast',
              sources=['Fast/fast.cpp'],
              include_dirs=["Fast"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )
    
# setup ======================================================================
setup(
    name="Fast",
    version="3.0",
    description="Fast Navier-Stokes solver.",
    author="Onera",
    package_dir={"":"."},
    packages=['Fast'],
    ext_modules=listExtensions
    )

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
