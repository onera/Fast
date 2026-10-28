#=============================================================================
# FastS requires:
# ELSAPROD variable defined in environment
# C++ compiler
# Fortran compiler: defined in config.py
# Numpy
# KCore
#=============================================================================
import os
from setuptools import setup, Extension
import KCore.Dist as Dist

# Compiler settings must be set in installBase.json / installBaseUser.json
additionalLibPaths = Dist.getAdditionalLibPaths()
additionalIncludePaths = Dist.getAdditionalIncludePaths()
additionalLibs = Dist.getAdditionalLibs()

# Write setup.cfg file
Dist.writeSetupCfg()

# Write build info to file
Dist.writeBuildInfo("FastS/FastS.py")

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

# Test if xcore exists =======================================================
(xcoreVersion, xcoreIncDir, xcoreLibDir) = Dist.checkModuleCassiopee("XCore")

# Test if connector exists =====================================================
(connectorVersion, connectorIncDir, connectorLibDir) = Dist.checkModuleCassiopee("Connector")

# Test if fast exists =======================================================
(fastcVersion, fastcIncDir, fastcLibDir) = Dist.checkModuleFast("FastC")

# Test if libmpi exists ======================================================
(mpi, mpiIncDir, mpiLibDir, mpiLibs) = Dist.checkMpi()

prod = os.getenv("ELSAPROD")
if prod is None:
    raise ValueError("Environment variable ELSAPROD must be set.")

# Setting libraryDirs, include dirs and libraries =============================
libraryDirs = ["build/"+prod, kcoreLibDir, xcoreLibDir, connectorLibDir, fastcLibDir, '.']
includeDirs = [numpyIncDir, kcoreIncDir, xcoreIncDir, connectorIncDir, fastcIncDir]
libraries = ["fasts", "fastc", "connector", "xcore", "kcore"]

(ok, libs, paths) = Dist.checkFortranLibs()
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs()
libraryDirs += paths; libraries += libs
ADDITIONALCPPFLAGS = []
if mpi:
    libraryDirs.append(mpiLibDir)
    includeDirs.append(mpiIncDir)
    ADDITIONALCPPFLAGS = ['-D_MPI']
    libraries += mpiLibs

# Extensions ==================================================================
listExtensions = [
    Extension(
        'FastS.fasts',
        sources=['FastS/fastS.cpp'],
        include_dirs=[".","FastS"]+additionalIncludePaths+includeDirs,
        library_dirs=additionalLibPaths+libraryDirs,
        libraries=libraries+additionalLibs,
        extra_compile_args=Dist.getCppArgs(),
        extra_link_args=Dist.getLinkArgs()+['-p']
    )
]

# setup ======================================================================
setup(
    name="fast-fastS",
    version="4.0",
    description="Fast for structured grids.",
    author="ONERA",
    url="https://fast.onera.fr",
    packages=['FastS'],
    package_dir={"":"."},
    ext_modules=listExtensions
)
