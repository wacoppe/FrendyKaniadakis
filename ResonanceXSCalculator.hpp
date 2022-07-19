#ifndef RESONANCE_XS_CALCULATOR_H
#define RESONANCE_XS_CALCULATOR_H

#include "Config/CommonData.hpp"
#include "CommonUtils/VectorCopier.hpp"
#include "CommonUtils/TabInterpolator.hpp"
#include "CommonUtils/TabAdjuster.hpp"
#include "NuclearDataUtils/NuclearDataObject.hpp"
#include "MathUtils/MathUtils.hpp"

#ifndef NO_AMUR_MODE
#include "ReconResonance/RmatrixLimited/run/RmatrixLimitedCalculator.hpp"
#endif //NO_AMUR_MODE

namespace frendy
{
  class ResonanceXSCalculator
  {
    protected:
      frendy::ErrorManager    err_obj;
      frendy::VectorClearer   clr_obj;
      frendy::VectorCopier    cp_obj;
      frendy::TabInterpolator ti_obj;
      frendy::TabAdjuster     ta_obj;
      frendy::MathUtils       math_obj;

#ifdef DEBUG_MODE
  public:
#endif
      Real8 k_part, third;
      static const int add_ele_no = 100000;

      frendy::NuclearDataObject        nucl_data_obj;
      frendy::ResonanceDataContainer   reso_data_obj;
#ifndef NO_AMUR_MODE
      frendy::RmatrixLimitedCalculator r_matrix_calc_obj;
#endif //NO_AMUR_MODE

      Integer calc_reso_xs_flg;
      Integer nucl_data_set_flg;

      Real8  temp, temp_nucl;

      ///////////// Declaracao de Variaveis para serem armazenadas ////////////////

      // PosINAC_FRENDY - V3

      Real8 psi_eff;

      // PosINAC_FRENDY - V5

      Real8 w_x_eff;
      Real8 w_y_eff;
      Real8 w_eff;
      Real8 x_eff;

      // PosINAC_FRENDY - V6

      Real8 l_max_eff;
      Real8 m_max_eff;
      Real8 m_eff;
      Real8 E0_eff;

      // PosINAC_FRENDY - V8

      Real8 qsi_eff;
      Real8 fxqsi_eff;

      // PosINAC_FRENDY - V10

      Real8 gtt_eff;

      Real8 smr_gg_eff;
   

      /////////////////////////////////////////////////////////////////////////////

      Real8  error_value, error_max_value, error_int_value;
      vector<vector<vector<Real8> > >            ene_array;
      vector<vector<vector<vector<Real8> > > >   sig_array;

      //For R-matrix limited
      vector<vector<vector<Integer> > >          react_type_list;
      vector<vector<vector<Real8> > >            q_array;

      int                  nis;
      vector<int>          ner;
      vector<vector<int> > nls;

      Integer reso_flg;
      vector<Real8> abn;
      vector<Integer> fis_width_flg;
      vector<vector<Integer> > reso_region_flg, xs_formula_flg, nro, radius_calc_flg, self_shielding_flg;
      vector<vector<Real8> >   spi, ap;
      vector<vector<vector<Real8> > > awri;

      vector<vector<vector<Integer> > > nbt_nro, int_nro;
      vector<vector<vector<Real> > >    e_int_nro, ap_nro;

      Real8 rpi;
      vector<vector<vector<Real8> > >          wave_no_part, chan_rad;
      vector<vector<vector<vector<Real8> > > > gx_er, sl_er, pl_er;
      vector<vector<Real8> >                   q_wei, q_abs;

      Real8         xs_pot;
      vector<Real8> ene_potential, xs_potential;

      void calc_resonance_xs();

      void check_ene_data(int i, int j, vector<Real8>& ene_data);

      void calc_reso_xs_scat_radius_only(int i, int j, Real8& ene_val, vector<Real8>& sig_val);
      void calc_reso_xs_bw_single(int i, int j, Real8& ene_val, vector<Real8>& sig_val);
      void calc_reso_xs_bw_multi( int i, int j, Real8& ene_val, vector<Real8>& sig_val);
      void calc_reso_xs_rm(       int i, int j, Real8& ene_val, vector<Real8>& sig_val);
      void calc_reso_xs_adler(    int i, int j, Real8& ene_val, vector<Real8>& sig_val);
      void calc_reso_xs_r_matrix( int i, int j, Real8& ene_val, vector<Real8>& sig_val);

      virtual void calc_reso_xs_unreso_a(int i, int j, Real8& ene_val, vector<Real8>& sig_val);
      virtual void calc_reso_xs_unreso_c(int i, int j, Real8& ene_val, vector<Real8>& sig_val);

      //For calc_reso_xs_rm
      void calc_complex_matrix(Real8 mat_r[3][3], Real8 mat_i[3][3],
                               Real8& ene_val, Real8& er, Real8& gn, Real8& gg, 
                               Real8& gfa_root, Real8& gfb_root, Real8& pl, Real8& pl_er_inv, Integer& gf);
      Integer get_ex_chan_spin_chk(int kchannl, Integer& kngtv, Integer& kpstv, 
                                   int jjl, int jj, int numj);

      void calc_rho(int k, Real8& ene, Real8& wave_no, Real8& rho, Real8& rho_h);
      void calc_rho(int i, int j, int k, Real8& ene, Real8& wave_no, Real8& rho, Real8& rho_h);
      void calc_channel_radius(int i, int j, Real8& ene, Real8& radii);

      void add_middle_energy_grid(int i, int j, vector<Real8>& ene_data, vector<vector<Real8> >& sig_data);
      Integer check_energy_grid_distance(int i, int j, int ele_no,
                                         vector<Real8>& ene_data, vector<vector<Real8> >& sig_data,
                                         Real8& mid_ene, vector<Real8>& mid_sig);
      void insert_middle_energy_grid(int ele_no, vector<Real8>& ene_data, vector<vector<Real8> >& sig_data,
                                      Real8& mid_ene, vector<Real8>& mid_sig);
      void add_xs_at_each_grid(vector<Real8>& new_ene,      vector<vector<Real8> >& new_sig, 
                               vector<Real8>& new_ene_part, vector<vector<Real8> >& new_sig_part);

      void delete_overlap_grid(vector<Real8>& ene_data, vector<vector<Real8> >& sig_data);
      
      void unify_energy_grid();
      void calc_unreso_xs();

      void calc_const_value_bw(int i, int j);
      void calc_const_value_rm(int i, int j);
      void calc_const_value_adler(int i, int j);
      void calc_const_value_unreso(int i, int j);

      void set_point_quadrature_weight_and_abscissa();

      void check_react_type_list(int i, int j);

      void input_data_check();

      void clear_all();
      
      //For set_multi_array_data
      int                             i_pre;
      int                             j_pre;
      int                             nls_multi;
      vector<int>                     nrs_multi;
      vector<int>                     njs_multi;
      vector<vector<int> >            nlj_multi;
      Integer                         reso_region_flg_multi;
      Integer                         xs_formula_flg_multi;
      Integer                         fis_width_flg_multi;
      Integer                         radius_calc_flg_multi;
      Integer                         li_multi;
      Integer                         ne_multi;
      Integer                         nro_multi;
      vector<Integer>                 nbt_nro_multi;
      vector<Integer>                 int_nro_multi;
      vector<Integer>                 l_multi;
      vector<Integer>                 lrx_multi;
      vector<vector<Integer> >        muf_multi;
      vector<Real>                    e_int_nro_multi;
      vector<Real>                    ap_nro_multi;
      vector<vector<vector<Real> > >  gf_b_multi;
      Real8                           abn_multi;
      Real8                           spi_multi;
      Real8                           den_multi;
      Real8                           ap_multi;
      vector<Real8>                   wave_no_part_multi;
      vector<Real8>                   qx_multi;
      vector<Real8>                   awri_multi;
      vector<Real8>                   apl_multi;
      vector<Real8>                   chan_rad_multi;
      vector<Real8>                   at_multi;
      vector<Real8>                   af_multi;
      vector<Real8>                   ac_multi;
      vector<vector<Real8> >          er_multi;
      vector<vector<Real8> >          gn_multi;
      vector<vector<Real8> >          gg_multi;
      vector<vector<Real8> >          gf_multi;
      vector<vector<Real8> >          gfa_multi;
      vector<vector<Real8> >          gfb_multi;
      vector<vector<Real8> >          gfa_root_multi;
      vector<vector<Real8> >          gfb_root_multi;
      vector<vector<Real8> >          aj_multi;
      vector<vector<Real8> >          sl_er_multi;
      vector<vector<Real8> >          pl_er_inv_multi;
      vector<vector<Real8> >          gx_er_multi;
      vector<vector<Real8> >          dx_multi;
      vector<vector<Real8> >          amun_multi;
      vector<vector<Real8> >          amux_multi;
      vector<vector<Real8> >          amuf_multi;
      vector<vector<Real8> >          gnox_multi;
      vector<vector<Real8> >          ggx_multi;
      vector<vector<vector<Real8> > > det_multi;
      vector<vector<vector<Real8> > > dwt_multi;
      vector<vector<vector<Real8> > > grt_multi;
      vector<vector<vector<Real8> > > git_multi;
      vector<vector<vector<Real8> > > def_multi;
      vector<vector<vector<Real8> > > dwf_multi;
      vector<vector<vector<Real8> > > grf_multi;
      vector<vector<vector<Real8> > > gif_multi;
      vector<vector<vector<Real8> > > dec_multi;
      vector<vector<vector<Real8> > > dwc_multi;
      vector<vector<vector<Real8> > > grc_multi;
      vector<vector<vector<Real8> > > gic_multi;
      vector<vector<vector<Real8> > > es_multi;
      vector<vector<vector<Real8> > > d_c_multi;
      vector<vector<vector<Real8> > > gx_c_multi;
      vector<vector<vector<Real8> > > gno_c_multi;
      vector<vector<vector<Real8> > > gg_c_multi;
      vector<vector<vector<Real8> > > gf_c_multi;
      
      void set_multi_array_data(int i, int j);
      void set_multi_array_data_bw(int i, int j);
      void set_multi_array_data_rm(int i, int j);
      void set_multi_array_data_adler(int i, int j);
      void set_multi_array_data_r_matrix(int i, int j);

      virtual void set_multi_array_data_unreso_a(int i, int j);
      virtual void set_multi_array_data_unreso_c(int i, int j);
      
      virtual void clear_multi_array_data();

    public:
      //constructor
      ResonanceXSCalculator(void);

      //destructor
      virtual ~ResonanceXSCalculator(void);

      virtual void clear();
      void         clear_reso_data();

      void calc_reso_xs_each_case(int i, int j, Real8& ene_val, vector<Real8>& sig_val);

      void set_nucl_data_obj(frendy::NuclearDataObject& data_obj);
      void set_resonance_grid(vector<vector<vector<Real8> > >& real_vec);

      void set_temp(Real8 real_val);


      ////////////////////////////////////////////////////////////////////////////////////////
      /////////////// Declaracao de objetos .set que armazenarao as variaveis ////////////////
      ////////////////////////////////////////////////////////////////////////////////////////
    
      // PosINAC_FRENDY - V3

      void set_psi(Real8 real_val);

      // PosINAC_FRENDY - V5

      void set_w_x(Real8 real_val);   
      void set_w_y(Real8 real_val);
      void set_w(Real8 real_val);
      void set_x(Real8 real_val);

      // PosINAC_FRENDY - V6

      void set_l_max(Real8 real_val);
      void set_m_max(Real8 real_val);
      void set_m(Real8 real_val);
      void set_E0(Real8 real_val);

      // PosINAC_FRENDY - V8

      void set_qsi(Real8 real_val);
      void set_fxqsi(Real8 real_val);

      // PosINAC_FRENDY - V10

      void set_gtt(Real8 real_val);

      void set_smr_gg(Real8 real_val);


      ////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////


      void set_err(Real8 real_val);
      void set_err_max(Real8 real_val);
      void set_err_int(Real8 real_val);

      void calc_shift_penetration_factor(Integer& l_val, Real8& rho, Real8& sl, Real8& pl);
      void calc_phase_shift(Integer& l_val, Real8& rho, Real8& phi);
      void calc_penetrability_factor(Integer& l, Real8& rho, Real8& amun, Real8& vl);
      void calc_width_fluctuation_factor( Real8& gn, Real8& gf, Real8& gg, int mu, int nu, int lamda,
                                          vector<Real8>& sig_val, Real8 df );

      vector<vector<vector<Integer> > >        get_resonance_react_type_list();
      vector<vector<vector<Real8> > >          get_resonance_q_val();
      vector<vector<vector<Real8> > >          get_resonance_grid();
      vector<vector<vector<vector<Real8> > > > get_resonance_xs();

      Real8 get_temp();

      /////////////// Declaração de objetos .get que "chamarao" as variaveis armazenadas ////////////////

      // PosINAC_FRENDY - V3

      Real8 get_psi();

      // PosINAC_FRENDY - V5

      Real8 get_w_x();
      Real8 get_w_y();
      Real8 get_w();
      Real8 get_x();

      // PosINAC_FRENDY - V6
      
      Real8 get_l_max();
      Real8 get_m_max();
      Real8 get_m();
      Real8 get_E0();

      // PosINAC_FRENDY - V8

      Real8 get_qsi();
      Real8 get_fxqsi();


      // PosINAC_FRENDY - V10

      Real8 get_gtt();

      Real8 get_smr_gg();




      ///////////////////////////////////////////////////////////////////////////////////////



      Real8 get_err();
      Real8 get_err_max();
      Real8 get_err_int();

      frendy::NuclearDataObject      get_nucl_data_obj();
      frendy::ResonanceDataContainer get_reso_data_obj();

      Real8 get_xs_potential();
  };
}

#endif //RESONANCE_XS_CALCULATOR_H
