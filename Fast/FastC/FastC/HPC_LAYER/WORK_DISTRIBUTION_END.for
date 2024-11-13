c      enddo
c      enddo
c      enddo  !boucle patern bloc


c 9999 continue

c#if CHECK_SPLIT > 0
c#include "../FastC/FastC/HPC_LAYER/check_split1.for"
c#endif

c      enddo
c      enddo
c      enddo  !boucle bloc_thread

c#if CHECK_SPLIT > 0
c#include "../FastC/FastC/HPC_LAYER/check_split2.for"
c#endif
