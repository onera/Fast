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
prod = os.getenv("ELSAPROD")
if prod is None: prod = 'xx'
installPathLocal = 'build/'+prod

# La librarie statique existe?
a = os.access(installPathLocal+"/libfast.a", os.F_OK)
if a:
    shutil.copyfile(installPathLocal+"/libfast.a", libPath+"/libfast.a")
else: # Essai en dynamique
    a = os.access(installPathLocal+"/fast"+__EXTMODULE__, os.F_OK)
    if a:
        shutil.copyfile(installPathLocal+"/fast"+__EXTMODULE__,
                        libPath+"/libfast"+__EXTSHARED__) 
    else:
        print("Error: fast"+__EXTMODULE__+" can not be found.")
