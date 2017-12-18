#==============================================================================
# Fichiers C++
#==============================================================================
cpp_srcs = ['Fast/checkVariables.cpp',
            'Fast/checkNumericsValue.cpp',
            'Fast/Eos/perfect.cpp']

#==============================================================================
# Fichiers fortran
#==============================================================================
for_srcs = ['Fast/Fortran/Cons2VelocityF.for',
            'Fast/Fortran/cpTimeStepConvF.for',
            'Fast/Fortran/cpTimeStepConvMotionF.for',
            'Fast/Fortran/cpTimeStepTurF.for',
            'Fast/Fortran/DenConvecF.for',
            'Fast/Fortran/HeatCoefF.for',
            'Fast/Fortran/HeatFluxF.for',
            'Fast/Fortran/MuF.for',         
            'Fast/Fortran/PressF.for',
            'Fast/Fortran/TempF.for',
            'Fast/Fortran/SetAllBValuesAtF.for',
            'Fast/Fortran/ExtrapAllBValuesF.for',
            'Fast/Fortran/ViscosityTensorF.for',
            'Fast/Fortran/ViscousFluxDensF.for']
