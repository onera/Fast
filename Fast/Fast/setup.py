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
Dist.writeBuildInfo("Fast/Fast.py")

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkModuleCassiopee("KCore")

# Test if xcore exists =======================================================
(xcoreVersion, xcoreIncDir, xcoreLibDir) = Dist.checkModuleCassiopee("XCore")

# Test if connector exists =====================================================
(connectorVersion, connectorIncDir, connectorLibDir) = Dist.checkModuleCassiopee("Connector")

# Test if fast exists =======================================================
#(fastcVersion, fastcIncDir, fastcLibDir) = Dist.checkFastC()
(fastcVersion, fastcIncDir, fastcLibDir) = Dist.checkModuleFast("FastC")

# Test if fasts exists =====================================================
#(fastsVersion, fastsIncDir, fastsLibDir) = Dist.checkFastS()
(fastsVersion, fastsIncDir, fastsLibDir) = Dist.checkModuleFast("FastS")

# Test if fastp exists =====================================================
#(fastpVersion, fastpIncDir, fastpLibDir) = Dist.checkFastP()

# Test if fastlbm exists =====================================================
#(fastlbmVersion, fastlbmIncDir, fastlbmLibDir) = Dist.checkFastLBM()
#(fastaslbmVersion, fastaslbmIncDir, fastaslbmLibDir) = Dist.checkFastASLBM()


# Test if libmpi exists ======================================================
(mpi, mpiIncDir, mpiLibDir, mpiLibs) = Dist.checkMpi()

prod = os.getenv("ELSAPROD")
if prod is None:
    raise ValueError("Environment variable ELSAPROD must be set.")

# Setting libraryDirs, include dirs and libraries =============================
libraryDirs = ["build/"+prod, kcoreLibDir, xcoreLibDir, connectorLibDir, fastcLibDir, fastsLibDir]
includeDirs = [numpyIncDir, kcoreIncDir, xcoreIncDir, connectorIncDir, fastcIncDir, fastsIncDir]
libraries = ["fast", "fasts", "fastc", "connector", "xcore", "kcore"]
#libraryDirs += [fastaslbmLibDir]
#includeDirs += [fastaslbmIncDir]
#libraries += ["fastaslbm"]

#libraryDirs = ["build/"+prod, kcoreLibDir, xcoreLibDir, connectorLibDir, fastcLibDir,fastsLibDir, fastlbmLibDir]
#includeDirs = [numpyIncDir, kcoreIncDir, xcoreIncDir, connectorIncDir, fastcIncDir, fastsIncDir, fastlbmIncDir]
#libraries = ["fast", "fasts", "fastlbm", "fastc", "connector", "xcore", "kcore"]
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
        'Fast.fast',
        sources=['Fast/fast.cpp'],
        include_dirs=[".","Fast"]+additionalIncludePaths+includeDirs,
        library_dirs=additionalLibPaths+libraryDirs,
        libraries=libraries+additionalLibs,
        extra_compile_args=Dist.getCppArgs(),
        extra_link_args=Dist.getLinkArgs()
    )
]

# setup ======================================================================
setup(
    name="fast-fast",
    version="4.0",
    description="Fast NS and LBM solvers.",
    author="ONERA",
    url="https://fast.onera.fr",
    packages=['Fast'],
    package_dir={"":"."},
    ext_modules=listExtensions
)
