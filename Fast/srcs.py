#==============================================================================
# Fichiers C++
#==============================================================================
cpp_srcs = ["Fast/computePT.cpp",
            "Fast/gsdr3.cpp",
            "Fast/INTERP/interplbmns_.cpp",
            "Fast/INTERP/recuplbmns_.cpp",
            "Fast/INTERP/compute_sij.cpp",
            "Fast/FILTER/filtering_.cpp"]
#==============================================================================
# Fichiers fortran
#==============================================================================
for_srcs = ["Fast/FILTER/move_to_temp.for",
            "Fast/INTERP/copy_valuespara.for",
            "Fast/INTERP/interp_rk3para.for",
            'Fast/INTERP/shiftstk_para.for',
            "Fast/FILTER/filtrage_5.for"]
