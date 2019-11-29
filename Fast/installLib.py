# Installation de libfast pour etre accessible par les autres modules
# Si libfast.a existe, on la recopie
# Sinon, on cherche fast.so ou fast.pyd, on le recopie en libfast.so ou dll
import os, shutil
import platform
system = platform.uname()[0]

if system == 'Windows':
    __EXTMODULE__ = '.pyd'
    __EXTSHARED__ = '.dll'
else:
    __EXTMODULE__ = '.so'
    __EXTSHARED__ = '.so'

import KCore.installPath as K
libPath = K.libPath
installPathLocal = K.installPath

# La librarie statique existe?
a = os.access(installPathLocal+"/Fast/libfast.a", os.F_OK)
if a:
    shutil.copyfile(installPathLocal+"/Fast/libfast.a", libPath+"/libfast.a")
else: # Essai en dynamique
    a = os.access(installPathLocal+"/Fast/fast"+__EXTMODULE__, os.F_OK)
    if a:
        shutil.copyfile(installPathLocal+"/Fast/fast"+__EXTMODULE__,
                        libPath+"/libfast"+__EXTSHARED__) 
    else:
        print("Error: fast"+__EXTMODULE__+" can not be found.")
