       lt2 = lt +  v2mtr
       lt3 = lt +  v3mtr
       ti(lt)  =rot(1,1)*ti0(lt) +rot(1,2)*ti0(lt2) +rot(1,3)*ti0(lt3)
       ti(lt2) =rot(2,1)*ti0(lt) +rot(2,2)*ti0(lt2) +rot(2,3)*ti0(lt3)
       ti(lt3) =rot(3,1)*ti0(lt) +rot(3,2)*ti0(lt2) +rot(3,3)*ti0(lt3)
