      enddo
      enddo
      enddo  !boucle patern bloc


 9999 continue

#if CHECK_SPLIT > 0
#include "FastC/HPC_LAYER/check_split1.for"
#endif

      enddo
      enddo
      enddo  !boucle bloc_thread

#if CHECK_SPLIT > 0
#include "FastC/HPC_LAYER/check_split2.for"
#endif
