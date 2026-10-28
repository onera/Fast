# Installation de libfast pour etre accessible par les autres modules
# Si libfast.a existe, on la recopie
# Sinon, on cherche fast.so ou fast.pyd, on le recopie en libfast.so ou dll
import os
import sys
import shutil
import KCore.Dist as Dist

system = Dist.getSystem()[0]
if system == 'Windows':
    __EXTMODULE__ = '.pyd'
    __EXTSHARED__ = '.dll'
else:
    __EXTMODULE__ = '.so'
    __EXTSHARED__ = '.so'

# Check PYTHONPATH
prefix = Dist.parseInstallPathPrefix(*sys.argv)
installPath = os.path.join(prefix, Dist.getInstallPathSuffix())
libPath = os.path.dirname(os.path.dirname(installPath))
Dist.checkPythonPath(installPath)
Dist.checkLdLibraryPath(libPath)

prod = os.getenv("ELSAPROD")
if prod is None:
    raise ValueError("Environment variable ELSAPROD must be set.")

# La librarie statique existe?
installPathLocal = 'build/'+prod
a = os.access(installPathLocal+"/libfasts.a", os.F_OK)
if a:
    shutil.copyfile(installPathLocal+"/libfasts.a", libPath+"/libfasts.a")
else: # Essai en dynamique
    a = os.access(installPathLocal+"/fasts"+__EXTMODULE__, os.F_OK)
    if a:
        shutil.copyfile(installPathLocal+"/fasts"+__EXTMODULE__,
                        libPath+"/libfasts"+__EXTSHARED__)
    else:
        print("Error: fasts"+__EXTMODULE__+" can not be found.")

# Copie aussi les fichiers de config json
shutil.copyfile("FastS/buildInfo.json", os.path.join(installPath, "FastS/buildInfo.json"))