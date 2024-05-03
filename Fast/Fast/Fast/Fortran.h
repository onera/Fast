// Les protos fortran
extern "C"
{
  void copy_valuespara_(E_Int* ipt_para_int, E_Int* ind_loop, E_Int* ind_loop_, E_Float* iptro, E_Float* stock, E_Int& ind, E_Int& sol, E_Int& taille);

  void interp_rk3para_( E_Int* ipt_param_int, E_Int* ind_loop, E_Int* ind_loop_, E_Float* iptro, E_Float* stock, E_Int& nstep, E_Int& nitrun, E_Int& taille);

  void shiftstk_para_( E_Int* ipt_param_int, E_Int* ind_loop, E_Int* ind_loop_, E_Float* stock, E_Int& sol_i, E_Int& sol_d, E_Int& taille);

  void fill_ghostcellspara_( E_Int* param_int, E_Int* ind_loop, E_Int* ind_loop_, E_Float* ropDnr, E_Float* ropRcv);

  void move_to_temp_( E_Int nd, E_Int nidom, E_Int* param_intt, E_Float* iptro, E_Float* iptroflt); //E_Int ng, E_Int nx, E_Int ny, E_Int nz );

  void filtrage_5_( E_Int nd, E_Int nidom, E_Int* param_intt, E_Float* iptro, E_Float* iptroflt);

}
