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

# Test if xcore exists =======================================================
(xcoreVersion, xcoreIncDir, xcoreLibDir) = Dist.checkXCore()

# Test if connector exists ==================================================
(connectorVersion, connectorIncDir, connectorLibDir) = Dist.checkConnector()

# Test if fastc exists =====================================================
(fastcVersion, fastcIncDir, fastcLibDir) = Dist.checkFastC()

# Test if fasts exists =====================================================
(fastsVersion, fastsIncDir, fastsLibDir) = Dist.checkFastS()



from KCore.config import *

# Compilation des fortrans ====================================================
if f77compiler == "None":
    print("Error: a fortran 77 compiler is required for compiling Fast.")
args = Dist.getForArgs(); opt = ''
for c, v in enumerate(args): opt += 'FOPT'+str(c)+'='+v+' '
os.system("make -e FC="+f77compiler+" WDIR=Fast/Fortran "+opt)
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'

# Setting libraryDirs, include dirs and libraries =============================
libraryDirs = ["build/"+prod, kcoreLibDir, xcoreLibDir, connectorLibDir, fastcLibDir,fastsLibDir, fastlbmLibDir ]
includeDirs = [numpyIncDir, kcoreIncDir, xcoreIncDir, connectorIncDir, pythonIncDir,  fastcIncDir, fastsIncDir, fastlbmIncDir]
libraries = ["fastc", "fasts", "fastlbm", "kcore", "xcore", "connector"]
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# Extensions ==================================================================
import srcs
listExtensions = []
listExtensions.append(
    Extension('Fast.fast',
              sources=['Fast/fast.cpp']+srcs.cpp_srcs,
              include_dirs=["Fast"]+additionalIncludePaths+includeDirs,
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )
    
# setup ======================================================================
setup(
    name="Fast",
    version="3.6",
    description="Fast Navier-Stokes solver.",
    author="Onera",
    package_dir={"":"."},
    packages=['Fast'],
    ext_modules=listExtensions
    )

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
