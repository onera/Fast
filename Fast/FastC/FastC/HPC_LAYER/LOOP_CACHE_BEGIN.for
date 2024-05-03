#ifdef E_OMP_SOUS_DOMAIN
        !!!  a reordonner plus tard pour traiter d'abord les sousdomaine verrous
        !pour les thread voisins
        do kcache = 1, ijkv_sdm(3)
        do jcache = 1, ijkv_sdm(2)
        do icache = 1, ijkv_sdm(1)

#endif
