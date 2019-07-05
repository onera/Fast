            tcx  = tijk(lmtr,ic)*ci_mtr
            tcy  = tijk(lmtr,jc)*cj_mtr
            tcz  = tijk(lmtr,kc)*ck_mtr

            surf = sqrt(tcx*tcx + tcy*tcy +tcz*tcz)
            surf = max(surf,1e-30)
            s_1  = 1./surf

            tcx = tcx*s_1
            tcy = tcy*s_1
            tcz = tcz*s_1

            ventx      = ventijk(ldp ,1)
            venty      = ventijk(ldp ,2)
            ventz      = ventijk(ldp ,kc_vent)*ck_vent

           !on calcule l'inverse de la surface de la facette

            u  =  rop(il,2)
            v  =  rop(il,3)
            w  =  rop(il,4)

            qn = (u*tcx+v*tcy+w*tcz)

            ut = u-qn*tcx
            vt = v-qn*tcy
            wt = w-qn*tcz

            ua = 2.*ut -u
            va = 2.*vt -v
            wa = 2.*wt -w

            rop(ir,1) = rop(il,1)
            rop(ir,2) = ventx*c_ale + ua
            rop(ir,3) = venty*c_ale + va
            rop(ir,4) = ventz*c_ale + wa

            rop(ir,5) = rop(il,5)
