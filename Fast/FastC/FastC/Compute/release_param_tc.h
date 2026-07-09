    for (E_Int nopass = 1; nopass <= npass_transfer; nopass++)
    {
      if      (nopass == 1)
       { 
         RELEASESHAREDN( pyParam_int_tc1, param_int_tc1); 
         if ( pyParam_real_tc1 != Py_None) { RELEASESHAREDN( pyParam_real_tc1, param_real_tc1);}
       }
      else if (nopass == 2)
       { 
         RELEASESHAREDN( pyParam_int_tc2, param_int_tc2); 
         if ( pyParam_real_tc2 != Py_None) { RELEASESHAREDN( pyParam_real_tc2, param_real_tc2);}
       }
      else if (nopass == 3)
       { 
         RELEASESHAREDN( pyParam_int_tc3, param_int_tc3); 
         if ( pyParam_real_tc3 != Py_None) { RELEASESHAREDN( pyParam_real_tc3, param_real_tc3);}
       }
      else if (nopass == 4) 
       { 
         RELEASESHAREDN( pyParam_int_tc4, param_int_tc4); 
         if ( pyParam_real_tc4 != Py_None) { RELEASESHAREDN( pyParam_real_tc4, param_real_tc4);}
       }
      else if (nopass == 5) 
       { 
         RELEASESHAREDN( pyParam_int_tc5, param_int_tc5); 
         if ( pyParam_real_tc5 != Py_None) { RELEASESHAREDN( pyParam_real_tc5, param_real_tc5);}
       }
      else if (nopass == 6) 
       { 
         RELEASESHAREDN( pyParam_int_tc6, param_int_tc6); 
         if ( pyParam_real_tc6 != Py_None) { RELEASESHAREDN( pyParam_real_tc6, param_real_tc6);}
       }
      else if (nopass == 7) 
       { 
         RELEASESHAREDN( pyParam_int_tc7, param_int_tc7); 
         if ( pyParam_real_tc7 != Py_None) { RELEASESHAREDN( pyParam_real_tc7, param_real_tc7);}
       }
      else if (nopass == 8) 
       { 
         RELEASESHAREDN( pyParam_int_tc8, param_int_tc8); 
         if ( pyParam_real_tc8 != Py_None) { RELEASESHAREDN( pyParam_real_tc8, param_real_tc8);}
       }
      else if (nopass == 9) 
       { 
         RELEASESHAREDN( pyParam_int_tc9, param_int_tc9); 
         if ( pyParam_real_tc9 != Py_None) { RELEASESHAREDN( pyParam_real_tc9, param_real_tc9);}
       }
      else if (nopass == 10) 
       { 
         RELEASESHAREDN( pyParam_int_tc10, param_int_tc10); 
         if ( pyParam_real_tc10 != Py_None) { RELEASESHAREDN( pyParam_real_tc10, param_real_tc10);}
       }
      else if (nopass == 11) 
       { 
         RELEASESHAREDN( pyParam_int_tc11, param_int_tc11); 
         if ( pyParam_real_tc11 != Py_None) { RELEASESHAREDN( pyParam_real_tc11, param_real_tc11);}
       }
      else if (nopass == 12) 
       { 
         RELEASESHAREDN( pyParam_int_tc12, param_int_tc12); 
         if ( pyParam_real_tc12 != Py_None) { RELEASESHAREDN( pyParam_real_tc12, param_real_tc12);}
       }
      else if (nopass == 13) 
       { 
         RELEASESHAREDN( pyParam_int_tc13, param_int_tc13); 
         if ( pyParam_real_tc13 != Py_None) { RELEASESHAREDN( pyParam_real_tc13, param_real_tc13);}
       }
      else if (nopass == 14) 
       { 
         RELEASESHAREDN( pyParam_int_tc14, param_int_tc14); 
         if ( pyParam_real_tc14 != Py_None) { RELEASESHAREDN( pyParam_real_tc14, param_real_tc14);}
       }
      else if (nopass == 15) 
       { 
         RELEASESHAREDN( pyParam_int_tc15, param_int_tc15); 
         if ( pyParam_real_tc15 != Py_None) { RELEASESHAREDN( pyParam_real_tc15, param_real_tc15);}
       }
      else {printf("transfert : npass > 15 pas codee."); return NULL;}
    }
