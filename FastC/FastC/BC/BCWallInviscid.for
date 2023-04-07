            tcx  = tijk(lmtr,ic)*ci_mtr
            tcy  = tijk(lmtr,jc)*cj_mtr
            tcz  = tijk(lmtr,kc)*ck_mtr

            surf = sqrt(tcx*tcx + tcy*tcy +tcz*tcz)
            surf = max(surf,1e-30)
            s_1  = 1./surf

            tcx = tcx*s_1
            tcy = tcy*s_1
            tcz = tcz*s_1

#if __DEBUG__
            ! encadre le depassement intentionnel
            if (c_ale > 0) then
#endif
            ventx      = ventijk(ldp ,1)
            venty      = ventijk(ldp ,2)
            ventz      = ventijk(ldp ,kc_vent)*ck_vent
#if __DEBUG__
            else
            ventx = 0.
            venty = 0.
            ventz = 0.
            endif
#endif

            ! on calcule l inverse de la surface de la facette
            
            u  =  rop(il,2) - ventx*c_ale
            v  =  rop(il,3) - venty*c_ale
            w  =  rop(il,4) - ventz*c_ale

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
