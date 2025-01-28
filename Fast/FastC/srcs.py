#==============================================================================
# Fichiers C++
#==============================================================================
cpp_srcs = ['FastC/checkNumericsValue.cpp',
            'FastC/getRange.cpp',
            'FastC/Init/initNuma.cpp',
            'FastC/Metric/init_metric.cpp',
            'FastC/ALE/motionlaw.cpp',
            'FastC/Compute/souszones_list.cpp',
            'FastC/HPC_LAYER/distributeThreads.cpp',
            'FastC/Com/setInterpTransfersFast.cpp']

#==============================================================================
# Fichiers fortran
#==============================================================================
for_srcs = ['FastC/Metric/skmtr.for',
            'FastC/Metric/tijk_extrap.for',
            'FastC/Metric/dist_extrap.for',
            'FastC/Metric/error.for',
            'FastC/Metric/cp_tijk.for',
            'FastC/Metric/cp_vol.for',
            'FastC/Metric/nature_geom_dom.for',
            'FastC/Metric/move.for',
            'FastC/Metric/cptijk_ale.for',
            'FastC/Metric/cpventijk_ale.for',
            'FastC/Init/copynuma.for',
            'FastC/STAT/cpmys_rij.for',
            'FastC/BC/indice_cl_sdm.for',
            'FastC/Compute/init_ssiter_bloc.for',
            'FastC/Compute/skip_lu.for',
            'FastC/Compute/ssdom_lu_ijk.for',
            'FastC/Compute/shape_tab_mtr.for',
            'FastC/HPC_LAYER/indice_boucle_ssdom.for',
            'FastC/HPC_LAYER/indice_boucle_lu.for',
            'FastC/HPC_LAYER/indice_boucle_scater.for',
            'FastC/HPC_LAYER/crsdm.for',
            'FastC/HPC_LAYER/crsdm_scater.for',
            'FastC/HPC_LAYER/synchro_omp_scater.for',
            'FastC/HPC_LAYER/topo_scater.for',
            'FastC/HPC_LAYER/verrou.for',
            'FastC/HPC_LAYER/flush_integer.for',
            'FastC/HPC_LAYER/flush_real.for']
