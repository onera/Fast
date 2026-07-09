    for (E_Int nopass = 1; nopass <= npass_transfer; nopass++)
    {
      //printf("transfert : npass %d %d \n",nopass, npass_transfer ); fflush(0);
      char str_real[16];char str_int[15];
      if      (nopass == 1)
       { strcpy(str_int, "param_int_tc1"); strcpy(str_real, "param_real_tc1");
         pyParam_int_tc1  = PyDict_GetItemString(work, str_int );
         pyParam_real_tc1 = PyDict_GetItemString(work, str_real); 
         K_NUMPY::getFromNumpyArray(pyParam_int_tc1 , param_int_tc1 );  int_tc[ nopass-1 ]= param_int_tc1 -> begin();
         if ( pyParam_real_tc1 != Py_None) { K_NUMPY::getFromNumpyArray(pyParam_real_tc1, param_real_tc1); real_tc[ nopass-1 ]= param_real_tc1-> begin(); }
         else{ real_tc[nopass-1] = NULL;}
       }
      else if (nopass == 2)
       { strcpy(str_int, "param_int_tc2"); strcpy(str_real, "param_real_tc2");
         pyParam_int_tc2  = PyDict_GetItemString(work, str_int );
         pyParam_real_tc2 = PyDict_GetItemString(work, str_real); 
         K_NUMPY::getFromNumpyArray(pyParam_int_tc2 , param_int_tc2 );  int_tc[ nopass-1 ]= param_int_tc2 -> begin();
         if ( pyParam_real_tc2 != Py_None) { K_NUMPY::getFromNumpyArray(pyParam_real_tc2, param_real_tc2); real_tc[ nopass-1 ]= param_real_tc2-> begin(); }
         else{ real_tc[nopass-1] = NULL;}
       }
      else if (nopass == 3)
       { strcpy(str_int, "param_int_tc3"); strcpy(str_real, "param_real_tc3");
         pyParam_int_tc3  = PyDict_GetItemString(work, str_int );
         pyParam_real_tc3 = PyDict_GetItemString(work, str_real); 
         K_NUMPY::getFromNumpyArray(pyParam_int_tc3 , param_int_tc3 );  int_tc[ nopass-1 ]= param_int_tc3 -> begin();
         if ( pyParam_real_tc3 != Py_None) { K_NUMPY::getFromNumpyArray(pyParam_real_tc3, param_real_tc3); real_tc[ nopass-1 ]= param_real_tc3-> begin(); }
         else{ real_tc[nopass-1] = NULL;}
       }
      else if (nopass == 4) 
       { strcpy(str_int, "param_int_tc4"); strcpy(str_real, "param_real_tc4");
         pyParam_int_tc4  = PyDict_GetItemString(work, str_int );
         pyParam_real_tc4 = PyDict_GetItemString(work, str_real); 
         K_NUMPY::getFromNumpyArray(pyParam_int_tc4 , param_int_tc4 );  int_tc[ nopass-1 ]= param_int_tc4 -> begin();
         if ( pyParam_real_tc4 != Py_None) { K_NUMPY::getFromNumpyArray(pyParam_real_tc4, param_real_tc4); real_tc[ nopass-1 ]= param_real_tc4-> begin(); }
         else{ real_tc[nopass-1] = NULL;}
       }
      else if (nopass == 5) 
       { strcpy(str_int, "param_int_tc5"); strcpy(str_real, "param_real_tc5");
         pyParam_int_tc5  = PyDict_GetItemString(work, str_int );
         pyParam_real_tc5 = PyDict_GetItemString(work, str_real); 
         K_NUMPY::getFromNumpyArray(pyParam_int_tc5 , param_int_tc5 );  int_tc[ nopass-1 ]= param_int_tc5 -> begin();
         if ( pyParam_real_tc5 != Py_None) { K_NUMPY::getFromNumpyArray(pyParam_real_tc5, param_real_tc5); real_tc[ nopass-1 ]= param_real_tc5-> begin(); }
         else{ real_tc[nopass-1] = NULL;}
       }
      else if (nopass == 6) 
       { strcpy(str_int, "param_int_tc6"); strcpy(str_real, "param_real_tc6");
         pyParam_int_tc6  = PyDict_GetItemString(work, str_int );
         pyParam_real_tc6 = PyDict_GetItemString(work, str_real); 
         K_NUMPY::getFromNumpyArray(pyParam_int_tc6 , param_int_tc6 );  int_tc[ nopass-1 ]= param_int_tc6 -> begin();
         if ( pyParam_real_tc6 != Py_None) { K_NUMPY::getFromNumpyArray(pyParam_real_tc6, param_real_tc6); real_tc[ nopass-1 ]= param_real_tc6-> begin(); }
         else{ real_tc[nopass-1] = NULL;}
       }
      else if (nopass == 7) 
       { strcpy(str_int, "param_int_tc7"); strcpy(str_real, "param_real_tc7");
         pyParam_int_tc7  = PyDict_GetItemString(work, str_int );
         pyParam_real_tc7 = PyDict_GetItemString(work, str_real); 
         K_NUMPY::getFromNumpyArray(pyParam_int_tc7 , param_int_tc7 );  int_tc[ nopass-1 ]= param_int_tc7 -> begin();
         if ( pyParam_real_tc7 != Py_None) { K_NUMPY::getFromNumpyArray(pyParam_real_tc7, param_real_tc7); real_tc[ nopass-1 ]= param_real_tc7-> begin(); }
         else{ real_tc[nopass-1] = NULL;}
       }
      else if (nopass == 8) 
       { strcpy(str_int, "param_int_tc8"); strcpy(str_real, "param_real_tc8");
         pyParam_int_tc8  = PyDict_GetItemString(work, str_int );
         pyParam_real_tc8 = PyDict_GetItemString(work, str_real); 
         K_NUMPY::getFromNumpyArray(pyParam_int_tc8 , param_int_tc8 );  int_tc[ nopass-1 ]= param_int_tc8 -> begin();
         if ( pyParam_real_tc8 != Py_None) { K_NUMPY::getFromNumpyArray(pyParam_real_tc8, param_real_tc8); real_tc[ nopass-1 ]= param_real_tc8-> begin(); }
         else{ real_tc[nopass-1] = NULL;}
       }
      else if (nopass == 9) 
       { strcpy(str_int, "param_int_tc9"); strcpy(str_real, "param_real_tc9");
         pyParam_int_tc9  = PyDict_GetItemString(work, str_int );
         pyParam_real_tc9 = PyDict_GetItemString(work, str_real); 
         K_NUMPY::getFromNumpyArray(pyParam_int_tc9 , param_int_tc9 );  int_tc[ nopass-1 ]= param_int_tc9 -> begin();
         if ( pyParam_real_tc9 != Py_None) { K_NUMPY::getFromNumpyArray(pyParam_real_tc9, param_real_tc9); real_tc[ nopass-1 ]= param_real_tc9-> begin(); }
         else{ real_tc[nopass-1] = NULL;}
       }
      else if (nopass == 10) 
       { strcpy(str_int, "param_int_tc10"); strcpy(str_real, "param_real_tc10");
         pyParam_int_tc10  = PyDict_GetItemString(work, str_int );
         pyParam_real_tc10 = PyDict_GetItemString(work, str_real); 
         K_NUMPY::getFromNumpyArray(pyParam_int_tc10 , param_int_tc10 );  int_tc[ nopass-1 ]= param_int_tc10 -> begin();
         if ( pyParam_real_tc10 != Py_None) { K_NUMPY::getFromNumpyArray(pyParam_real_tc10, param_real_tc10); real_tc[ nopass-1 ]= param_real_tc10-> begin(); }
         else{ real_tc[nopass-1] = NULL;}
       }
      else if (nopass == 11) 
       { strcpy(str_int, "param_int_tc11"); strcpy(str_real, "param_real_tc11");
         pyParam_int_tc11  = PyDict_GetItemString(work, str_int );
         pyParam_real_tc11 = PyDict_GetItemString(work, str_real); 
         K_NUMPY::getFromNumpyArray(pyParam_int_tc11 , param_int_tc11 );  int_tc[ nopass-1 ]= param_int_tc11 -> begin();
         if ( pyParam_real_tc11 != Py_None) { K_NUMPY::getFromNumpyArray(pyParam_real_tc11, param_real_tc11); real_tc[ nopass-1 ]= param_real_tc11-> begin(); }
         else{ real_tc[nopass-1] = NULL;}
       }
      else if (nopass == 12) 
       { strcpy(str_int, "param_int_tc12"); strcpy(str_real, "param_real_tc12");
         pyParam_int_tc12  = PyDict_GetItemString(work, str_int );
         pyParam_real_tc12 = PyDict_GetItemString(work, str_real); 
         K_NUMPY::getFromNumpyArray(pyParam_int_tc12 , param_int_tc12 );  int_tc[ nopass-1 ]= param_int_tc12 -> begin();
         if ( pyParam_real_tc12 != Py_None) { K_NUMPY::getFromNumpyArray(pyParam_real_tc12, param_real_tc12); real_tc[ nopass-1 ]= param_real_tc12-> begin(); }
         else{ real_tc[nopass-1] = NULL;}
       }
      else if (nopass == 13) 
       { strcpy(str_int, "param_int_tc13"); strcpy(str_real, "param_real_tc13");
         pyParam_int_tc13  = PyDict_GetItemString(work, str_int );
         pyParam_real_tc13 = PyDict_GetItemString(work, str_real); 
         K_NUMPY::getFromNumpyArray(pyParam_int_tc13 , param_int_tc13 );  int_tc[ nopass-1 ]= param_int_tc13 -> begin();
         if ( pyParam_real_tc13 != Py_None) { K_NUMPY::getFromNumpyArray(pyParam_real_tc13, param_real_tc13); real_tc[ nopass-1 ]= param_real_tc13-> begin(); }
         else{ real_tc[nopass-1] = NULL;}
       }
      else if (nopass == 14) 
       { strcpy(str_int, "param_int_tc14"); strcpy(str_real, "param_real_tc14");
         pyParam_int_tc14  = PyDict_GetItemString(work, str_int );
         pyParam_real_tc14 = PyDict_GetItemString(work, str_real); 
         K_NUMPY::getFromNumpyArray(pyParam_int_tc14 , param_int_tc14 );  int_tc[ nopass-1 ]= param_int_tc14 -> begin();
         if ( pyParam_real_tc14 != Py_None) { K_NUMPY::getFromNumpyArray(pyParam_real_tc14, param_real_tc14); real_tc[ nopass-1 ]= param_real_tc14-> begin(); }
         else{ real_tc[nopass-1] = NULL;}
       }
      else if (nopass == 15) 
       { strcpy(str_int, "param_int_tc15"); strcpy(str_real, "param_real_tc15");
         pyParam_int_tc15  = PyDict_GetItemString(work, str_int );
         pyParam_real_tc15 = PyDict_GetItemString(work, str_real); 
         K_NUMPY::getFromNumpyArray(pyParam_int_tc15 , param_int_tc15 );  int_tc[ nopass-1 ]= param_int_tc15 -> begin();
         if ( pyParam_real_tc15 != Py_None) { K_NUMPY::getFromNumpyArray(pyParam_real_tc15, param_real_tc15); real_tc[ nopass-1 ]= param_real_tc15-> begin(); }
         else{ real_tc[nopass-1] = NULL;}
       }
      else {printf("transfert : npass > 15 pas codee."); return NULL;}
    }
    if (npass_transfer == 0){int_tc[0] = NULL; real_tc[0] = NULL; }
