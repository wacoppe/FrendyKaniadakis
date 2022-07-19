#include "ReconResonance/ResonanceXSCalculator.hpp"
#include <math.h>                                      // PosINAC_FRENDY - V3
#include "Faddeeva.cc" //adicionei para tentar usar a funcao erro complexa do pacote do MIT 20220622

using namespace frendy;

//constructor
ResonanceXSCalculator::ResonanceXSCalculator(void)
{
  clear_all();
  rpi    = sqrt(M_PI);
  k_part = (sqrt(2.0*amu_n*amu*ev)*1.0e-12)/h_bar;
  third  = 1.0/3.0;

  set_point_quadrature_weight_and_abscissa();
}

//destructor
ResonanceXSCalculator::~ResonanceXSCalculator(void)
{
  clear_all();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_resonance_xs()
{
  if( calc_reso_xs_flg != 0 )
  {
    return;
  }

  calc_reso_xs_flg     = 1;
  if( reso_flg != 1 )
  {
    clr_obj.clear_vec_array3_real8(ene_array);
    clr_obj.clear_vec_array4_real8(sig_array);
    return;
  }

  input_data_check();

  int i_max = nis;
  sig_array.resize(i_max);
  if( i_max != static_cast<int>(ene_array.size()) )
  {
    string class_name = "ResonanceXSCalculator";
    string func_name  = "calc_resonance_xs()";

    ostringstream oss01, oss02;
    oss01 << nis;
    oss02 << static_cast<int>(ene_array.size());
    string str_data01 = "NIS           : " + oss01.str();
    string str_data02 = "Size of array : " + oss02.str();

    vector<string> err_com;
    err_com.push_back(str_data01);
    err_com.push_back(str_data02);
    err_com.push_back("NIS and energy vector size is different.");
    err_com.push_back("Please check the energy vector.");

    err_obj.output_runtime_error(class_name, func_name, err_com);
  }

  for(int i=0; i<i_max; i++)
  {
    int j_max = ner[i];
    sig_array[i].resize(j_max);
    if( j_max != static_cast<int>(ene_array[i].size()) )
    {
      string class_name = "ResonanceXSCalculator";
      string func_name  = "calc_resonance_xs()";

      ostringstream oss01, oss02, oss03, oss04;
      oss01 << i;
      oss02 << i_max;
      oss03 << ner[i];
      oss04 << static_cast<int>(ene_array[i].size());
      string str_data01 = "NIS           : " + oss01.str() + " / " + oss02.str();
      string str_data02 = "NER           : " + oss03.str();
      string str_data03 = "Size of array : " + oss04.str();

      vector<string> err_com;
      err_com.push_back(str_data01);
      err_com.push_back(str_data02);
      err_com.push_back(str_data03);
      err_com.push_back("NER and energy vector size is different.");
      err_com.push_back("Please check the energy vector.");

      err_obj.output_runtime_error(class_name, func_name, err_com);
    }

    for(int j=0; j<j_max; j++)
    {
      vector<Real8>          ene_data = ene_array[i][j];
      vector<vector<Real8> > sig_data;

      check_ene_data(i, j, ene_data);

      int k_max = static_cast<int>(ene_data.size());
      sig_data.resize(k_max);
      for(int k=0; k<k_max; k++)
      {
        calc_reso_xs_each_case(i, j, ene_data[k], sig_data[k]);
      }
      
      int add_grid_chk = 0;
      if( reso_region_flg[i][j] == 2 && self_shielding_flg[i][j] == 1 )
      {
        add_grid_chk = -1;
      }
      else if( reso_region_flg[i][j] == 0 )
      {
        add_grid_chk = -1;
        ene_array[i][j].resize(0);
        clr_obj.clear_vec_array2_real8(sig_array[i][j]);
      }

      if( add_grid_chk == 0 )
      {
        add_middle_energy_grid(i, j, ene_data, sig_data);
        delete_overlap_grid(ene_data, sig_data);
      }

      ene_array[i][j] = ene_data;
      sig_array[i][j] = sig_data;
      ene_data.clear();
      clr_obj.clear_vec_array2_real8(sig_data);

      //This option is removed since the addition of the resonance cross section
      //to MT=19 (first chance fission) is changed from this class to the
      //ReconstructXSMerger class (update_xs_tot_sc_fis_rad function).
      //check_react_type_list(i, j);
    }
  }
  
  unify_energy_grid();
  calc_unreso_xs();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::check_ene_data(int i, int j, vector<Real8>& ene_data)
{
  vector<Real8> ene_data_new;
  Real el = reso_data_obj.get_lower_ene_limit()[i][j] * (1.0 - min_ene_dif);
  Real eh = reso_data_obj.get_upper_ene_limit()[i][j] * (1.0 + min_ene_dif);

  int k_max    = static_cast<int>(ene_data.size());
  int warn_flg = 0;
  for(int k=0; k<k_max; k++)
  {
    if( ene_data[k] > el && ene_data[k] < eh )
    {
      ene_data_new.push_back(ene_data[k]);
    }
    else if( warn_flg == 0 )
    {
      warn_flg = 1;

      string class_name = "ResonanceXSCalculator";
      string func_name  = "check_ene_data(int i, int j, vector<Real8>& ene_data)";

      ostringstream oss01, oss02, oss03, oss04, oss05, oss06;
      oss01 << i;
      oss02 << nis;
      oss03 << j;
      oss04 << ner[i];
      oss05 << el;
      oss06 << eh;
      string str_data01 = "Number of isotopes (NIS)               : " + oss01.str() + " / " + oss02.str();
      string str_data02 = "Number of resonance energy range (NER) : " + oss03.str() + " / " + oss04.str();
      string str_data03 = "Lower limit for an energy range (EL)   : " + oss05.str();
      string str_data04 = "Upper limit for an energy range (EH)   : " + oss06.str();

      vector<string> err_com;
      err_com.push_back(str_data01);
      err_com.push_back(str_data02);
      err_com.push_back(str_data03);
      err_com.push_back(str_data04);
      err_com.push_back("The set energy grid data is not within the energy range.");

      err_obj.output_caution(class_name, func_name, err_com);
    }
  }

  ene_data = ene_data_new;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_reso_xs_each_case(int i, int j, Real8& ene_val, vector<Real8>& sig_val)
{
  set_multi_array_data(i, j);

  sig_val.clear();
  sig_val.resize(xs_type_no);
  for(int l=0; l<xs_type_no; l++)
  {
    sig_val[l] = 0.0;
  }

  if( reso_flg != 1 )
  {
    return;
  }

  xs_pot = 0.0;
  if( reso_region_flg_multi == 1 )
  {
    switch(xs_formula_flg_multi)
    {
      case 1: //Single Level Breit-Wigner
        calc_reso_xs_bw_single(i, j, ene_val, sig_val);
        break;
      case 2: //Multi Level Breit-Wigner
        calc_reso_xs_bw_multi(i, j, ene_val, sig_val);
        break;
      case 3: //Reich-Moore
        calc_reso_xs_rm(i, j, ene_val, sig_val);
        break;
      case 4: //Adler-Adler
        calc_reso_xs_adler(i, j, ene_val, sig_val);
        break;
      case 7: //R-Matrix Limited Format
        calc_reso_xs_r_matrix(i, j, ene_val, sig_val);
        break;
      default:
        string class_name = "ResonanceXSCalculator";
        string func_name  = "calc_reso_xs_each_case(int i, int j, Real8& ene_val, vector<Real8>& sig_val)";

        ostringstream oss01, oss02, oss03, oss04, oss05, oss06;
        oss01 << i;
        oss02 << nis;
        oss03 << j;
        oss04 << ner[i];
        oss05 << reso_region_flg_multi;
        oss06 << xs_formula_flg_multi;
        string str_data01 = "Number of isotopes (NIS)               : " + oss01.str() + " / " + oss02.str();
        string str_data02 = "Number of resonance energy range (NER) : " + oss03.str() + " / " + oss04.str();
        string str_data03 = "Resonance region flg (LRU)             : " + oss05.str();
        string str_data04 = "XS formula flg (LRF)                   : " + oss06.str();

        vector<string> err_com;
        err_com.push_back(str_data01);
        err_com.push_back(str_data02);
        err_com.push_back(str_data03);
        err_com.push_back(str_data04);
        err_com.push_back("This xs_formula_flg (LRF) value is not applicable in this program.");
        err_com.push_back("Supported xs_formula_flg (LRF) value is 1 - 4 and 7.");

        err_obj.output_runtime_error(class_name, func_name, err_com);
        break;
    }

    if( xs_formula_flg_multi != 7 )
    {
      Real8 radii = 0.0;
      calc_channel_radius(i, j, ene_val, radii);
      xs_pot = 4.0 * M_PI * radii * radii;
    }
    else
    {
      vector<vector<Real> > scat_radius_true;
      scat_radius_true = reso_data_obj.get_r_matrix_data_obj(i, j).get_scat_radius_true();

      Real xs_pot_tot = 0.0;
      Real data_no    = 0.0;
      for(int l=0; l<static_cast<int>(scat_radius_true.size()); l++)
      {
        for(int m=0; m<static_cast<int>(scat_radius_true[l].size()); m++)
        {
          if( scat_radius_true[l][m] > min_value )
          {
            xs_pot_tot += 4.0 * M_PI * scat_radius_true[l][m] * scat_radius_true[l][m];
            data_no    += 1.0;
          }
        }
      }
      clr_obj.clear_vec_array2_real(scat_radius_true);

      xs_pot = 0.0;
      if( data_no > min_value )
      {
        xs_pot = xs_pot_tot / data_no;
      }
    }
  }
  else if( reso_region_flg_multi == 2 ) //(Unresolved resonance)
  {
    if( xs_formula_flg_multi == 1 )
    {
      //Case A(fis_width_flg(lfw)=0), Case B(fis_width_flg(lfw)=1)
      calc_reso_xs_unreso_a(i, j, ene_val, sig_val);
    }
    else //xs_formula_flg_multi==2
    {
      //Case C
      calc_reso_xs_unreso_c(i, j, ene_val, sig_val);
    }

    if( self_shielding_flg[i][j] != 0 )
    {
      for(int l=0; l<xs_type_no; l++)
      {
        sig_val[l] = 0.0;
      }
    }
  }
  else //if( reso_region_flg_multi == 0 )
  {
    calc_reso_xs_scat_radius_only(i, j, ene_val, sig_val);
    xs_pot = sig_val[scatter_xs];
  }

  ene_potential.push_back(ene_val);
  xs_potential.push_back(xs_pot);

  int l_max = static_cast<int>(sig_val.size());
  for(int l=0; l<l_max; l++)
  {
    if( sig_val[l] < 0.0 )
    {
      sig_val[l] = 0.0;
    }
  }

  //Multiply the abundance
  if( abn_multi > 0.0 )
  {
    for(int l=0; l<l_max; l++)
    {
      sig_val[l] = abn_multi * sig_val[l];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_reso_xs_scat_radius_only(int i, int j, 
                                                          Real8& ene_val, vector<Real8>& sig_val)
{
  Real8 s_p = 4.0 * M_PI * ap_multi * ap_multi;
  sig_val[total_xs]   = s_p;
  sig_val[scatter_xs] = s_p;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_reso_xs_bw_single(int i, int j, Real8& ene_val, vector<Real8>& sig_val)
{
  //Single Level Breit Wigner

  Real8   x, x2, x3, x4, x5, qsi, qsi2, qsi3, qsi4, qsi5, fxqsi;   //Declaração de variáveis para cálculo de Tempo_eff

  
  Real8   tbk_coef, spifac, pifac;
  Real8   wave_no,    rho,    rho_h,    sl,    pl, phi;
  Real8   wave_no_qx, rho_qx, rho_h_qx, sl_qx, pl_qx;
  Real8   qx, e_qx, cos2p, sin2p, sinsq, delta;
  Real8   er, gn, gg, gf, gx, gj, edelt, gne, gtt;
  Real8   pl_er_inv;
  Real8   ex, ax, ay, psi, chi, smr;
  Real8   comfac, add;
  Integer l_val;

  if( temp > 0.0 )
  {
    tbk_coef = 1.0 / (temp*boltzmann_const * 4.0*ene_val);
  }
  else
  {
    tbk_coef = 0.0;
    delta    = 0.0;
  }

  spifac  = den_multi;
  wave_no = wave_no_part_multi[0] * sqrt(fabs(ene_val));
  pifac   = 4.0*M_PI / (wave_no * wave_no);
  
  int l_max = nls_multi;
  for(int l=0; l<l_max; l++)
  {
    qx   = qx_multi[l];
    e_qx = ene_val + qx;

    l_val = l_multi[l];
    calc_rho(l, ene_val, wave_no, rho, rho_h);
    calc_shift_penetration_factor(l_val, rho, sl, pl);
    calc_phase_shift(l_val, rho_h, phi);
    
    if( lrx_multi[l] != 0 )
    {
      wave_no_qx = wave_no_part_multi[l] * sqrt(fabs(e_qx));
      calc_rho(l, e_qx, wave_no_qx, rho_qx, rho_h_qx);
      calc_shift_penetration_factor(l_val, rho_qx, sl_qx, pl_qx);
    }
    else
    {
      pl_qx = 0.0;
    }
    
    cos2p = cos(2.0*phi);
    sin2p = sin(2.0*phi);
    sinsq = sin(phi);
    sinsq = sinsq * sinsq;
    if( temp > 0.0 )
    {
      delta = sqrt(awri_multi[l]*tbk_coef); 
      //delta = sqrt(reso_data_obj.get_mass_isotope()[i][j][l]*tbk_coef); // PosINAC_FRENDY - V4
                                                                           //teste para verificar a origem de awri_multi[l] 
    }
    
    int m_max = nrs_multi[l];

    for(int m=0; m<m_max; m++)
    {
      er        = er_multi[l][m];
      gn        = gn_multi[l][m];
      gg        = gg_multi[l][m];
      gf        = gf_multi[l][m];
      gx        = gg + gf;
      gj        = (aj_multi[l][m]+0.5)*spifac;
      pl_er_inv = pl_er_inv_multi[l][m];
      edelt     = 2.0*(ene_val - (er + 0.5*gn*(sl_er_multi[l][m]-sl)*pl_er_inv));
      gne       = gn * pl * pl_er_inv;
      gtt       = gne + gx + gx_er_multi[l][m]*pl_qx;

      if( temp > 0.0 )
      {
        ex    = edelt/gtt;          // definição de x eq (5.2.13)
        ay    = 0.5*gtt*delta;      // equivale a qsi * 0.5, ou seja h  (5.2.25)   

        
        ///////////////////////////////////////////////////////////////////////////
        ////////////////// MODIFICAÇÕES SUGERIDAS PELO AQUILINO ///////////////////
        ///////////////////////////////////////////////////////////////////////////
        //// p/ gerar dados s/ Kaniadakis, deixe essa parte do codigo comentada ///
        ///////////////////////////////////////////////////////////////////////////
        /*
        qsi = ay / 0.5 ;
        x = fabs(ex) ;

        qsi2 = pow( qsi, 2.0 );
        qsi3 = pow( qsi, 3.0 );
        qsi4 = pow( qsi, 4.0 );
        qsi5 = pow( qsi, 5.0 );

        x2 = pow( x, 2.0 );
        x3 = pow( x, 3.0 );
        x4 = pow( x, 4.0 );
        x5 = pow( x, 5.0 );
        
        fxqsi = (-0.0095668) + 0.0042817*x + 1.2301*qsi + 0.00025628*x2 + (-0.1449)*x*qsi 
        + (-1.2704)*qsi2 + (-0.000061751)*x3 + 0.0062072*x2*qsi + 0.92537*x*qsi2 + 0.76292*qsi3
        + 0.0000029873*x4 + (-0.00010748)*x3*qsi + (-0.023709)*x2*qsi2 + (-2.1902)*x*qsi3 + 6.4376*qsi4
        + (-0.000000041619)*x5 + 0.00000080913*x4*qsi + 0.00015675*x3*qsi2 + 0.023925*x2*qsi3 + 1.7734*x*qsi4
        + (-9.4704)*qsi5;

        ay = fxqsi * 0.5;
        */
        /////////////////////////////////////////////////////////////////////////
        ////////// FIM DAS MODIFICAÇÕES SUGERIDAS PELO AQUILINO /////////////////
        /////////////////////////////////////////////////////////////////////////
        

        ax = ay*ex;                                 //ax é igual a U da eq (5.2.24)

        math_obj.calc_cerfc(ax, ay, psi, chi);
        Real8 cerfc_coef = rpi*ay;

        if (m==641)                 // quando m=1 entao trata-se da segunda regiao de ressonancia
        {
          set_w_x(psi);           // PosINAC_FRENDY - V5          
                                  // psi e chi antes de multiplicar por cerfc_coef 
                                  // equivalem as partes reais e imaginarias da funcao erro
          set_w_y(chi);           // PosINAC_FRENDY - V5    
          set_w(psi+chi);         // Soma da parte imaginaria e real da funcao w (omega)
                                  // PosINAC_FRENDY - V5

          //set_x(ex);              // PosINAC_FRENDY - V6
          
        }

        

        //////////////////////////////////////////////////////
        ///////////IMPLEMENTAÇÃO DE KANI ANALÍTICO////////////
        //////////////////////////////////////////////////////
        /*
          
        
          Real8  Kappa, Bk, f1, f2, x1, chi1;

          chi1 = ay / 0.5 ;
          x1 = ex * 1.0 ;
          Kappa = 0.1  ;
          complex <double> i {0.,1.};
          complex <double> aa, d1, d2, d3, d4, d5, p1, p2, p3, theta, D, om1, om2, f3, omg, Sg, z, p11;
          
                   
          Bk=(1)*(pow(2*Kappa,1.5))*(1+3*Kappa/2)*tgamma(1/(2*Kappa)+0.75)/tgamma(1/(2*Kappa)-0.75);  //Bk que vou tentar imprimir
          
          f1 = (chi1*(sqrt(M_PI))*Bk)/4; //parte mais à esquerda da função
          
          f2 = exp((pow(chi1,2)-pow(chi1,2)*pow(x1,2))/4); //#Exponencial mais externa da solução
              
          aa=(-1.*i*pow(chi1,2)*x2);  //#Termos presentes em p1, p2 e p3
          
          d1= pow(chi1,4); //#Termos presentes em p1, p2 e p3
                  
          d2 = (2*pow(chi1,2)*pow(Kappa,2));
                
          d3 = sqrt(d1 - d2);

          d4 = aa + d3;

          d5 = aa - d3;
                  
          p1 = (d4)/(2.0*chi1); //#Eq. 3.92 da tese
                
          p2 = (d5)/(2.0*chi1);  //#Eq. 3.93 da tese

          p3 = (d3)/(2.0*chi1);  //#Eq. 3.93 da tese
          
          theta = (x1/2)*d3; //#definição de theta da tese (Eq. 3.95)
          
          D = ((2-(2*erf(chi1/2)))/(1-pow(Kappa,2)))*cos(theta); //# D(Eq. 3.99 da tese)
          
          om1 = (sin(theta))*((Faddeeva::erf(p1)*(pow(Kappa,2)))-Faddeeva::erf(p1)+(Faddeeva::erf(p2)*(pow(Kappa,2)))-Faddeeva::erf(p2)); //# Omega 1(Eq. 3.90 da tese)
          
          om2 = (cos(theta))*((2.0*Faddeeva::erf(p3)*(pow(Kappa,2)))-2.0*Faddeeva::erf(p3)-(Faddeeva::erf(p1)*(pow(Kappa,2)))+Faddeeva::erf(p1)+(Faddeeva::erf(p2)*(pow(Kappa,2)))-Faddeeva::erf(p2)); //# Omega 1(Eq. 3.90 da tese)
          
          f3 = d3/((-pow(chi1,2))+(2*pow(Kappa,2))); //#primeiro trecho dentro de omega g

          omg =  f3*(exp(pow(-Kappa,2)/2))*((i*om1)+(om2));   //# Omega geral(Eq. 3.100 da tese)
          
          Sg = (f1*f2)*(D+omg); //# Solução geral (Eq. 3.98)
          
          //cerfc_coef = real(Sg);

          psi = real(Sg);

        */
        //////////////////////////////////////////////////////
        ////////FIM DA IMPLEMENTAÇÃO DE KANI ANALÍTICO////////
        //////////////////////////////////////////////////////

        psi *= cerfc_coef;       // deve ser comentado para implementação de KANI ANALITICO
        chi *= cerfc_coef;
        smr  = pifac*gj*gne/(gtt*gtt);
        
        sig_val[scatter_xs]   +=  smr * ((cos2p*gtt - gx)*psi + sin2p*gtt*chi);
        sig_val[fission_xs]   +=  smr * gf * psi;
        sig_val[radiation_xs] +=  smr * gg * psi;

        comfac = pifac*gj*gne/(edelt*edelt + gtt*gtt);
        add    = comfac * (gne*cos2p - 2.0*gx*sinsq + edelt*sin2p);

        
        /////////////////////// Armazenamento de variaveis /////////////////////

        if (m==2)
        {
          set_psi(psi);                                           // PosINAC_FRENDY - V3
          set_qsi(qsi);                                           // PosINAC_FRENDY - V8
          set_fxqsi(fxqsi);                                       // PosINAC_FRENDY - V8
          set_E0((er + 0.5*gn*(sl_er_multi[l][m]-sl)*pl_er_inv));  // PosINAC_FRENDY - V8
          set_x(ex);              // PosINAC_FRENDY - V
          set_gtt(gtt);
          set_smr_gg(smr * gg);
        }

        set_l_max(l_max*1.0);     // PosINAC_FRENDY - V6 
        set_m_max(m_max*1.0);     // PosINAC_FRENDY - V6 

        /*
        // Para plotar variáveis com diferentes valores de m (vairias ressonancias)        
        int m_ref = int( ( ene_val - 55.5    )/         0.000275        );    // PosINAC_FRENDY - V7
        /////                     ene_data[0])/(ene_data[1]-ene_data[0]);

        if (m<=m_ref)         //Define em qual iteração as variáveis serão salvas para plot   //  
        {
          set_m(m*1.0);
          set_E0((er + 0.5*gn*(sl_er_multi[l][m]-sl)*pl_er_inv));
        }
        */
    
        ////////////////// Fim do Armazenamento de Variaveis //////////////////
      
      }
      else
      {
        comfac = pifac*gj*gne/(edelt*edelt + gtt*gtt);
        add    = comfac * (gne*cos2p - 2.0*gx*sinsq + edelt*sin2p);
        
        sig_val[scatter_xs]   +=  add;
        sig_val[fission_xs]   +=  comfac * gf;
        sig_val[radiation_xs] +=  comfac * gg;
      }
    }
    sig_val[scatter_xs] += static_cast<Real8>(2*l_val+1)*pifac*sinsq;
  }
  sig_val[total_xs] = sig_val[scatter_xs] + sig_val[fission_xs] + sig_val[radiation_xs];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_reso_xs_bw_multi(int i, int j, Real8& ene_val, vector<Real8>& sig_val)
{
  //Multi-Level Breit Wigner
  if( ( temp > 0.0 ) || ( temp_nucl > min_ene_val ) )
  {
    string class_name = "ResonanceXSCalculator";
    string func_name  = "calc_reso_xs_bw_multi(int i, int j, Real8& ene_val, vector<Real8>& sig_val)";

    ostringstream oss01, oss02, oss03, oss04, oss05, oss06;
    oss01 << i;
    oss02 << nis;
    oss03 << j;
    oss04 << ner[i];
    oss05 << temp;
    oss06 << temp_nucl;
    string str_data01 = "Number of isotpes (NIS)                : " + oss01.str() + " / " + oss02.str();
    string str_data02 = "Number of resonance energy range (NER) : " + oss03.str() + " / " + oss04.str();
    string str_data03 = "Temperature                            : " + oss05.str();
    string str_data04 = "Temperature in the nuclear data        : " + oss06.str();

    vector<string> err_com;
    err_com.push_back(str_data01);
    err_com.push_back(str_data02);
    err_com.push_back(str_data03);
    err_com.push_back(str_data04);

    if( temp > 0.0 )
    {
      err_com.push_back("Temperature is larger than 0[K].");
    }
    if( temp_nucl > min_ene_val )
    {
      err_com.push_back("Temperature in the nuclear data file is larger than 0[K].");
    }
    err_com.push_back("Multi-Level Brite Wigner is only applicable T=0[K].");
    err_obj.output_runtime_error(class_name, func_name, err_com);
  }
  
  Real8   pifac;
  Real8   wave_no,    rho,    rho_h,    sl,    pl, phi;
  Real8   wave_no_qx, rho_qx, rho_h_qx, sl_qx, pl_qx;
  Real8   qx, e_qx, cos2p, sin2p, fl, ajmin, ajmax, aj, sum, diff;
  Real8   gn, gg, gf, gne, gtt, ax, comfac, add;
  Real8   pl_er_inv;
  Integer l_val;
  int     nj, aj_no;
  vector<Real8> gj;
  vector<vector<Real8> > sigj;
  gj.clear();
  clr_obj.clear_vec_array2_real8(sigj);

  wave_no = wave_no_part_multi[0] * sqrt(fabs(ene_val));
  pifac  = M_PI / (wave_no * wave_no);
  
  int l_max = nls_multi;
  for(int l=0; l<l_max; l++)
  {
    qx   = qx_multi[l];
    e_qx = ene_val + qx;

    l_val = l_multi[l];
    calc_rho(l, ene_val, wave_no, rho, rho_h);
    calc_shift_penetration_factor(l_val, rho, sl, pl);
    calc_phase_shift(l_val, rho_h, phi);
    
    if( lrx_multi[l] != 0 )
    {
      wave_no_qx = wave_no_part_multi[l] * sqrt(fabs(e_qx));
      calc_rho(l, e_qx, wave_no_qx, rho_qx, rho_h_qx);
      calc_shift_penetration_factor(l_val, rho_qx, sl_qx, pl_qx);
    }
    else
    {
      pl_qx = 0.0;
    }
    
    cos2p = 1.0 - cos(2.0*phi);
    sin2p = sin(2.0*phi);
    
    fl    = static_cast<Real8>(l_val);
    ajmin = fabs(fabs(spi_multi-fl) - 0.5);
    ajmax = spi_multi+fl+0.5;
    nj    = static_cast<int>(round(ajmax - ajmin + 1.0));
    
    aj    = ajmin;
    sum   = 0.0;
    gj.resize(nj);
    for(int m=0; m<nj; m++)
    {
      gj[m] = (2.0*aj+1.0)*den_multi;
      aj  += 1.0;
      sum += gj[m];
    }
    diff = 2.0*fl+1.0-sum;
    
    sigj.resize(nj);
    for(int m=0; m<nj; m++)
    {
      sigj[m].resize(2);
      for(int n=0; n<2; n++)
      {
        sigj[m][n] = 0.0;
      }
    }
    
    int m_max = nrs_multi[l];
    for(int m=0; m<m_max; m++)
    {
      aj_no           = static_cast<int>(round(aj_multi[l][m] - ajmin));
      gn              = gn_multi[l][m];
      gg              = gg_multi[l][m];
      gf              = gf_multi[l][m];
      pl_er_inv       = pl_er_inv_multi[l][m];
      gne             = gn * pl * pl_er_inv;
      gtt             = gne + (gg + gf) + gx_er_multi[l][m]*pl_qx;
      ax              = 2.0*(ene_val - (er_multi[l][m] + 0.5*gn*(sl_er_multi[l][m]-sl)*pl_er_inv))/gtt;
      comfac          = 2.0*gne/(gtt*(1.0+ax*ax));
      sigj[aj_no][0] += comfac;
      sigj[aj_no][1] += comfac * ax;
      
      comfac *= gj[aj_no]/gtt;
      sig_val[fission_xs]   +=  comfac * gf;
      sig_val[radiation_xs] +=  comfac * gg;
    }
    
    add = 0.0;
    for(int m=0; m<nj; m++)
    {
      add += gj[m]*((cos2p - sigj[m][0])*(cos2p - sigj[m][0]) + (sin2p + sigj[m][1])*(sin2p + sigj[m][1]));
    }
    sig_val[scatter_xs] += add + 2.0*diff*cos2p;
    
    gj.clear();
    clr_obj.clear_vec_array2_real8(sigj);
  }
  sig_val[scatter_xs]   *= pifac;
  sig_val[fission_xs]   *= 2.0*pifac;
  sig_val[radiation_xs] *= 2.0*pifac;
  sig_val[total_xs]     = sig_val[scatter_xs] + sig_val[fission_xs] + sig_val[radiation_xs];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_reso_xs_rm(int i, int j, Real8& ene_val, vector<Real8>& sig_val)
{
  //Reich-Moore
  if( ( temp > 0.0 ) || ( temp_nucl > min_ene_val ) )
  {
    string class_name = "ResonanceXSCalculator";
    string func_name  = "calc_reso_xs_rm(int i, int j, Real8& ene_val, vector<Real8>& sig_val)";

    ostringstream oss01, oss02, oss03, oss04, oss05, oss06;
    oss01 << i;
    oss02 << nis;
    oss03 << j;
    oss04 << ner[i];
    oss05 << temp;
    oss06 << temp_nucl;
    string str_data01 = "Number of isotopes (NIS)               : " + oss01.str() + " / " + oss02.str();
    string str_data02 = "Number of resonance energy range (NER) : " + oss03.str() + " / " + oss04.str();
    string str_data03 = "Temperature                            : " + oss05.str();
    string str_data04 = "Temperature in the nuclear data        : " + oss06.str();

    vector<string> err_com;
    err_com.push_back(str_data01);
    err_com.push_back(str_data02);
    err_com.push_back(str_data03);
    err_com.push_back(str_data04);
    if( temp > 0.0 )
    {
      err_com.push_back("Temperature is larger than 0[K].");
    }
    if( temp_nucl > min_ene_val )
    {
      err_com.push_back("Temperature in the nuclear data file is larger than 0[K].");
    }
    err_com.push_back("Reich-Moore can only applicable T=0[K].");
    err_obj.output_runtime_error(class_name, func_name, err_com);
  }
  
  Real8   wave_no, pifac;
  Real8   apl, rho, rho_h, sl, pl, phi, p1, p2, phi_sq;
  Real8   fl, ajmin, ajmax, ajc, gj, aj_ori, aj;
  Real8   t1, t2, t3, t4, uur, uui, termf, termt, termn;
  Real8   dd, dd_sq, rr, ss, ss_sq, amag, rri, ssi, xx;
  Integer gf, l_val, kngtv, kpstv, ex_chan_spin_chk;
  int     nj, jjl;
  
  vector<Real8>  aj_ori_vec, pl_er_inv_vec, er, gn, gg, gfa_root, gfb_root;
  Real8 mat_r[3][3], mat_i[3][3];
  
  wave_no = wave_no_part_multi[0] * sqrt(fabs(ene_val));
  pifac   = M_PI / (wave_no * wave_no);
  
  gf = 0;
  int l_max = nls_multi;
  for(int l=0; l<l_max; l++)
  {
    l_val = l_multi[l];
    apl   = apl_multi[l];
    calc_rho(l, ene_val, wave_no, rho, rho_h);
    if( apl != 0.0 )
    {
      rho_h = wave_no*apl;
      
      if( radius_calc_flg_multi == 1 )
      {
        rho = wave_no*apl;
      }
    }
    calc_shift_penetration_factor(l_val, rho, sl, pl);
    calc_phase_shift(l_val, rho_h, phi);
    
    p1     = cos(2.0*phi);
    p2     = sin(2.0*phi); 
    phi_sq = phi*phi;
    
    fl    = static_cast<Real8>(l_val);
    ajmin = fabs(fabs(spi_multi - fl) - 0.5);
    ajmax = spi_multi + fl + 0.5;
    nj    = static_cast<int>(round(ajmax - ajmin + 1.0));
    ajc   = ajmin - 1.0;
    
    jjl = 1;
    if( l_val!=0 )
    {
      if( fl > (spi_multi - 0.5) && fl <= spi_multi )
      {
        jjl = 0;
      }
    }
    
    int n_max = nrs_multi[l];
    aj_ori_vec    = aj_multi[l];
    pl_er_inv_vec = pl_er_inv_multi[l];
    er            = er_multi[l];
    gn            = gn_multi[l];
    gg            = gg_multi[l];
    gfa_root      = gfa_root_multi[l];
    gfb_root      = gfb_root_multi[l];
    for(int m=0; m<nj; m++)
    {
      ajc += 1.0;
      gj   = (2.0*ajc+1.0)*den_multi;
      
      for(int kchanl=0; kchanl<2; kchanl++)
      {
        kngtv = 0;
        kpstv = 0;
        for(int ii=0; ii<3; ii++)
        {
          for(int jj=0; jj<3; jj++)
          {
            mat_r[ii][jj] = 0.0;
            mat_i[ii][jj] = 0.0;
          }
        }
        
        for(int n=0; n<n_max; n++)
        {
          aj_ori = aj_ori_vec[n];
          aj     = fabs(aj_ori);
          if( fabs(aj - ajc) <= 0.25 )
          {
            if( aj_ori < 0.0 )
            {
              kngtv++;
            }
            else if( aj_ori > 0.0 )
            {
              kpstv++;
            }
            
            if( (kchanl == 0 && aj_ori < 0.0) || (kchanl == 1 && aj_ori > 0.0) )
            {
              //skip
            }
            else
            {
              calc_complex_matrix(mat_r, mat_i, ene_val, er[n], gn[n], gg[n], gfa_root[n], gfb_root[n],
                                  pl, pl_er_inv_vec[n], gf);
            }
          }
        }
        
        ex_chan_spin_chk = get_ex_chan_spin_chk(kchanl, kngtv, kpstv, jjl, m, nj);
        if( ex_chan_spin_chk != 0 )
        {
          if( gf != 0 )
          {
            mat_r[0][0] += 1.0;
            mat_r[1][1] += 1.0;
            mat_r[2][2] += 1.0;
            
            mat_r[1][0] = mat_r[0][1];
            mat_r[2][0] = mat_r[0][2];
            mat_r[2][1] = mat_r[1][2];
            
            mat_i[1][0] = mat_i[0][1];
            mat_i[2][0] = mat_i[0][2];
            mat_i[2][1] = mat_i[1][2];
            math_obj.calc_inv_comp_matrix_3x3(mat_r, mat_i);
            
            dd    = 2.0*mat_r[0][0] - 1.0;
            ss    = 2.0*mat_i[0][0];
            t1    = mat_r[0][1];
            t2    = mat_i[0][1];
            t3    = mat_r[0][2];
            t4    = mat_i[0][2];
            termf = 4.0*gj*(t1*t1 + t2*t2 + t3*t3 + t4*t4);
            uur  =  1.0 - p1*dd - p2*ss;
            uui  = -1.0 * p2*dd + p1*ss;
            termt = 2.0 * gj *  uur; 
            termn = gj * (uur*uur + uui*uui);
          }
          else
          {
            dd    = mat_r[0][0];
            rr    = 1.0 + dd;
            ss    = mat_i[0][0];
            ss_sq = ss*ss;
            amag  = 1.0 / (rr*rr + ss_sq);
            rri   =  2.0 * rr*amag - 1.0;
            ssi   = -2.0 * ss*amag;
            uur   =  1.0 - p1*rri - p2*ssi;
            uui   = -1.0 * p2*rri + p1*ssi;
            
            if( fabs(dd) < 3.0e-4 && phi_sq < 9.0e-8 )
            {
              dd_sq = dd*dd;
              xx = (2.0*dd + 2.0*(dd_sq + ss_sq + phi_sq + p2*ss) - 2.0*phi_sq*(dd_sq + ss_sq)) * amag;
              termt = 2.0 * gj * xx;
              termn = gj * (xx*xx + uui*uui);
            }
            else
            {
              termt = 2.0*gj*uur;
              termn = gj * (uur*uur + uui*uui);
            }
            termf = 0.0;
          }
        
          if( ex_chan_spin_chk == 2 )
          {
            termn += 2.0 * gj * (1.0 - p1);
            termt += 2.0 * gj * (1.0 - p1);
          }
          
          sig_val[total_xs]     += termt;
          sig_val[scatter_xs]   += termn;
          sig_val[fission_xs]   += termf;
          sig_val[radiation_xs] += termt - termf - termn;
        } //if( ex_chan_spin_chk > 0 )
      } //for(int kchanl=0; kchanl<2; kchanl++)
    } //for(int m=0; m<nj; m++)
  } //for(int l=0; l<l_max; l++)
  
  sig_val[total_xs]     *= pifac;
  sig_val[scatter_xs]   *= pifac;
  sig_val[fission_xs]   *= pifac;
  sig_val[radiation_xs] *= pifac;

  aj_ori_vec.clear();
  pl_er_inv_vec.clear();
  er.clear();
  gn.clear();
  gg.clear();
  gfa_root.clear();
  gfb_root.clear();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_complex_matrix(Real8 mat_r[3][3], Real8 mat_i[3][3],
                                                Real8& ene_val, Real8& er, Real8& gn, Real8& gg, 
                                                Real8& gfa_root, Real8& gfb_root, Real8& pl, Real8& pl_er_inv, 
                                                Integer& gf)
{
  Real8 a_val[3];
  
  a_val[0] = sqrt(gn*pl*pl_er_inv);
  a_val[1] = gfa_root;
  a_val[2] = gfb_root;
  
  Real8 diff = er - ene_val;
  Real8 den = 1.0 / (diff*diff + 0.25*gg*gg);
  Real8 de2 = 0.5 * diff*den;
  Real8 gg4 = 0.25* gg*den;
  if( gfa_root != 0.0 || gfb_root != 0.0 )
  {
    for(int ii=0; ii<3; ii++)
    {
      for(int jj=ii; jj<3; jj++)
      {
        mat_r[ii][jj] += gg4*a_val[ii]*a_val[jj];
        mat_i[ii][jj] -= de2*a_val[ii]*a_val[jj];
      }
    }
    gf = 1;
  }
  else
  {
    mat_r[0][0] += gg4*a_val[0]*a_val[0];
    mat_i[0][0] -= de2*a_val[0]*a_val[0];
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

Integer ResonanceXSCalculator::get_ex_chan_spin_chk(int kchanl, Integer& kngtv, Integer& kpstv, 
                                                    int jjl, int jj, int numj)
{
  Integer ex_chan_spin_chk = 0;
  jj++; //a[i] in C++ is a[i+1] in Fortran
  
  if( kchanl == 0 )
  {
    if( kpstv > 0 )
    {
      if( kngtv == 0 )
      {
        if( jj > jjl && jj < numj )
        {
          ex_chan_spin_chk = 2;
        }
        else
        {
          ex_chan_spin_chk = 1;
        }
      }
      else if( kngtv > 0 )
      {
        ex_chan_spin_chk = 1;
      }
    }
    else if( kpstv == 0 )
    {
      if( kngtv == 0 )
      {
        if( jj > jjl && jj < numj )
        {
          ex_chan_spin_chk = 2;
        }
        else
        {
          ex_chan_spin_chk = 1;
        }
      }
    }
  }
  else if( kchanl == 1 )
  {
    if( kpstv > 0 )
    {
      if( kngtv > 0 )
      {
        ex_chan_spin_chk = 1;
      }
    }
    else if( kpstv == 0 )
    {
      if( kngtv > 0 )
      {
        if( jj > jjl && jj < numj )
        {
          ex_chan_spin_chk = 2;
        }
        else
        {
          ex_chan_spin_chk = 1;
        }
      }
    }
  }
  return ex_chan_spin_chk;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_reso_xs_adler(int i, int j, Real8& ene_val, vector<Real8>& sig_val)
{
  //Adler-Adler
  Real8 delta, e_root, e_sq, e_inv, e_sq_inv, e_cub_inv, wave_no, c_val, omg, snf, csf;
  Real8 bakt, bakf, bakc;
  Integer li;
  
  if( temp > 0.0 )
  {
    delta = sqrt( awri_multi[0] / (4.0*boltzmann_const*temp*ene_val) );
  }
  else
  {
    delta = 0.0;
  }
  
  e_root     = sqrt(fabs(ene_val));
  e_sq       = ene_val * ene_val;
  e_inv      = 1.0 / ene_val;
  e_sq_inv   = e_inv*e_inv;
  e_cub_inv  = e_inv*e_sq_inv;
  wave_no    = wave_no_part_multi[0] * e_root;
  c_val      = chan_rad_multi[0] / e_root;
  omg        = 2.0*wave_no*ap_multi;
  snf        = sin(omg);
  csf        = cos(omg);
  
  bakt = 0.0;
  bakf = 0.0;
  bakc = 0.0;
  li = li_multi;
  if( li < 5 )
  {
    string class_name = "ResonanceXSCalculator";
    string func_name  = "calc_reso_xs_adler(int i, int j, Real8& ene_val, vector<Real8>& sig_val)";

    ostringstream oss01, oss02, oss03, oss04, oss05;
    oss01 << i;
    oss02 << nis;
    oss03 << j;
    oss04 << ner[i];
    oss05 << li;
    string str_data01 = "Number of isotopes (NIS)               : " + oss01.str() + " / " + oss02.str();
    string str_data02 = "Number of resonance energy range (NER) : " + oss03.str() + " / " + oss04.str();
    string str_data03 = "adler_calc_flg (LI)                    : " + oss05.str();

    vector<string> err_com;
    err_com.push_back(str_data01);
    err_com.push_back(str_data02);
    err_com.push_back(str_data03);
    err_com.push_back("adler_calc_flg(LI) is less than 5.");
    err_com.push_back("This program can not treat LI < 5  at MF02MT151.");
    err_com.push_back("(reso_region_flg_multi(LRU)=1, xs_formula_flg(LRF)=4)");
    err_obj.output_runtime_error(class_name, func_name, err_com);
  }
  
  //Background total cross section
  if ( li == 5 || li == 7 )
  {
    bakt = (at_multi[0] + at_multi[1]*e_inv + at_multi[2]*e_sq_inv + at_multi[3]*e_cub_inv
          + at_multi[4]*ene_val + at_multi[5]*e_sq) * c_val;
  }
  
  //Background fission cross section
  if( li == 6 || li == 7 )
  {
    bakf = (af_multi[0] + af_multi[1]*e_inv + af_multi[2]*e_sq_inv + af_multi[3]*e_cub_inv
         +  af_multi[4]*ene_val + af_multi[5]*e_sq) * c_val;
  }
  
  //Background capture cross section
  bakc = (ac_multi[0] + ac_multi[1]*e_inv + ac_multi[2]*e_sq_inv + ac_multi[3]*e_cub_inv
       +  ac_multi[4]*ene_val + ac_multi[5]*e_sq) * c_val;
  
  //Calculate resonance contribution
  Real8 det, dwt, grt, git, def, dwf, grf, gif, dec, dwc, grc, gic;
  Real8 x_t, psi_t, chi_t, x_f, psi_f, chi_f, x_c, psi_c, chi_c;
  Real8 ax, ay;
  int l_max = nls_multi;
  for(int l=0; l<l_max; l++)
  {
    int m_max = njs_multi[l];
    for(int m=0; m<m_max; m++)
    {
      int n_max = nlj_multi[l][m];
      for(int n=0; n<n_max; n++)
      {
        det = det_multi[l][m][n];
        dwt = dwt_multi[l][m][n];
        grt = grt_multi[l][m][n];
        git = git_multi[l][m][n];
        x_t = det - ene_val / dwt;
        
        def = def_multi[l][m][n];
        dwf = dwf_multi[l][m][n];
        grf = grf_multi[l][m][n];
        gif = gif_multi[l][m][n];
        x_f = def - ene_val / dwf;
        
        dec = dec_multi[l][m][n];
        dwc = dwc_multi[l][m][n];
        grc = grc_multi[l][m][n];
        gic = gic_multi[l][m][n];
        x_c = dec - ene_val / dwc;
        
        if( temp <= 0.0 )
        {
          psi_t = 1.0 / (1.0 + x_t*x_t);
          chi_t = x_t * psi_t;
          
          psi_f = 1.0 / (1.0 + x_f*x_f);
          chi_f = x_f * psi_f;
          
          psi_c = 1.0 / (1.0 + x_c*x_c);
          chi_c = x_c * psi_c;
        }
        else
        {
          ay     = dwt * delta;
          ax     = ay * x_t;
          math_obj.calc_cerfc(ax, ay, psi_t, chi_t);
          psi_t  = rpi*ay*psi_t;
          chi_t  = rpi*ay*chi_t;
          
          ay     = dwf * delta;
          ax     = ay * x_f;
          math_obj.calc_cerfc(ax, ay, psi_f, chi_f);
          psi_f  = rpi*ay*psi_f;
          chi_f  = rpi*ay*chi_f;
          
          ay     = dwc * delta;
          ax     = ay * x_c;
          math_obj.calc_cerfc(ax, ay, psi_c, chi_c);
          psi_c  = rpi*ay*psi_c;
          chi_c  = rpi*ay*chi_c;
        }
        
        sig_val[total_xs]     += ((grt*csf + git*snf) * psi_t + (git*csf - grt*snf) * chi_t) / dwt;
        sig_val[fission_xs]   += (grf*psi_f + gif*chi_f) / dwf;
        sig_val[radiation_xs] += (grc*psi_c + gic*chi_c) / dwc;
      }
    }
  }
  
  //Add background cross section
  sig_val[total_xs]     = (sig_val[total_xs]     + bakt) * c_val + 2.0*c_val*(1.0-csf)/e_root;
  sig_val[fission_xs]   = (sig_val[fission_xs]   + bakf) * c_val;
  sig_val[radiation_xs] = (sig_val[radiation_xs] + bakc) * c_val;
  
  if( li != 6 )
  {
    sig_val[scatter_xs] = sig_val[total_xs] - (sig_val[fission_xs] + sig_val[radiation_xs]);
  }
  else
  {
    sig_val[total_xs] = sig_val[fission_xs] + sig_val[radiation_xs];
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_reso_xs_r_matrix(int i, int j, Real8& ene_val, vector<Real8>& sig_val)
{
#ifndef NO_AMUR_MODE
  r_matrix_calc_obj.calc_resonance_xs( ene_val, sig_val );
#else
  string class_name = "ResonanceXSCalculator";
  string func_name  = "calc_reso_xs_r_matrix(int i, int j, Real8& ene_val, vector<Real8>& sig_val)";

  ostringstream oss01, oss02, oss03, oss04, oss05, oss06;
  oss01 << i+1;
  oss02 << j+1;
  string str_data01 = "Number of isotopes (NIS)               : " + oss01.str();
  string str_data02 = "Number of resonance energy range (NER) : " + oss02.str();

  vector<string> err_com;
  err_com.push_back(str_data01);
  err_com.push_back(str_data02);
  err_com.push_back("NO_AMUR_MODE flag is used in this calculation.");
  err_com.push_back("In this case, R-matrix limited formula cannot be treated.");
  err_com.push_back("If user wants to process R-matrix limited, please remove this flag.");

  err_obj.output_runtime_error(class_name, func_name, err_com);
#endif //NO_AMUR_MODE
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_reso_xs_unreso_a(int i, int j, Real8& ene_val, vector<Real8>& sig_val)
{
  Real8   e_root, wave_no, const_val, rho, rho_h, vl, phi, sph, spot;
  Real8   dx, aj, amun, gnox, ggx, gfx, gj, gnx;
  Integer fis_width_flg_val, l_val;
  Integer muf;
  vector<Real8>   sig_tmp;
  vector<Real>    es, gf;
  vector<Integer> nbt_gf, int_gf;
  
  xs_pot = 0.0;

  if( nro_multi == 0 && ap_multi < 0.0)
  {
    string class_name = "ResonanceXSCalculator";
    string func_name  = "calc_reso_xs_unreso_a(int i, int j, Real8& ene_val, vector<Real8>& sig_val)";

    ostringstream oss01, oss02, oss03, oss04, oss05, oss06;
    oss01 << i;
    oss02 << nis;
    oss03 << j;
    oss04 << ner[i];
    oss05 << nro_multi;
    oss06 << ap_multi;
    string str_data01 = "Number of isotopes (NIS)               : " + oss01.str() + " / " + oss02.str();
    string str_data02 = "Number of resonance energy range (NER) : " + oss03.str() + " / " + oss04.str();
    string str_data03 = "scat_radius_ene_dependence_flg (NRO)   : " + oss05.str();
    string str_data04 = "scat_radius (AP)                       : " + oss06.str();

    vector<string> err_com;
    err_com.push_back(str_data01);
    err_com.push_back(str_data02);
    err_com.push_back(str_data03);
    err_com.push_back(str_data04);
    err_com.push_back("scat_radius(AP) at reso_region_flg(LRU)=2 is less than 0.0.");
    err_com.push_back("scat_radius(AP) or scat_radius_ene_dependence_flg(NRO) value is not appropriate.");

    err_obj.output_runtime_error(class_name, func_name, err_com);
  }
  
  sig_tmp.resize(xs_type_no); //total, scatter, fission, capture
  
  int k_max = static_cast<int>(es_multi.size());
  es.resize(k_max);
  for(int k=0; k<k_max; k++)
  {
    es[k] = es_multi[k][0][0];
  }
  nbt_gf.resize(1);
  int_gf.resize(1);
  nbt_gf[0] = k_max;
  int_gf[0] = int_lin_lin; //linear-linear
  
  fis_width_flg_val = fis_width_flg_multi;
  
  e_root  = sqrt(fabs(ene_val));
  wave_no = wave_no_part_multi[0] * e_root;
  const_val = 2.0*M_PI*M_PI / (wave_no*wave_no);
  
  int l_max = nls_multi;
  for(int l=0; l<l_max; l++)
  {
    l_val = l_multi[l];
    calc_rho(l, ene_val, wave_no, rho, rho_h);
    calc_phase_shift(l_val, rho_h, phi);
    sph = sin(phi);
    
    int m_max = njs_multi[l];
    for(int m=0; m<m_max; m++)
    {
      dx   = dx_multi[l][m];
      aj   = aj_multi[l][m];
      amun = amun_multi[l][m];
      gnox = gnox_multi[l][m];
      ggx  = ggx_multi[l][m];
      
      if( fis_width_flg_val == 0 ) //Case A
      {
        muf  = 1;
        gfx  = 0.0;
      }
      else //Case B
      {
        muf = muf_multi[l][m];
        gf  = gf_b_multi[l][m];
        ti_obj.interpolation_tab1(ene_val, gfx, nbt_gf, int_gf, es, gf);
      }
      calc_penetrability_factor(l_val, rho, amun, vl);
      vl = vl*e_root;
      
      gnx      = gnox * vl;
      gj       = const_val*gnx*(2.0*aj + 1.0) * den_multi / dx;
      
      calc_width_fluctuation_factor( gnx, gfx, ggx, static_cast<int>(round(amun)),
                                     static_cast<int>(muf), 1, sig_tmp, 0.0 );
      
      sig_val[scatter_xs]   += sig_tmp[scatter_xs]  *gj*gnx - 2.0*gj*sph*sph;
      sig_val[fission_xs]   += sig_tmp[fission_xs]  *gj*gfx;
      sig_val[radiation_xs] += sig_tmp[radiation_xs]*gj*ggx;
    }
    spot = sph / wave_no;
    spot = 4.0*M_PI*(2.0*static_cast<Real8>(l_val) + 1.0) * spot*spot;
    sig_val[scatter_xs] += spot; //add potential scatter xs
    xs_pot              += spot;
  }
  sig_val[total_xs] = sig_val[scatter_xs] + sig_val[fission_xs] + sig_val[radiation_xs];
  
  sig_tmp.clear();
  es.clear();
  gf.clear();
  nbt_gf.clear();
  int_gf.clear();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_reso_xs_unreso_c(int i, int j, Real8& ene_val, vector<Real8>& sig_val)
{
  Real8   e_root, wave_no, const_val, rho, rho_h, vl, phi, sph, spot;
  Real8   aj, amux, amun, amuf, gj, gnx;
  Real8   es_m, es_p, d_m, d_val, d_p, gx_m, gx, gx_p, gn_m, gn, gn_p, gg_m, gg, gg_p, gf_m, gf, gf_p;
  Integer l_val;
  int mu, nu, lamda, es_pos;

  Integer int_val = int_lin_lin;
  
  vector<Real8> sig_tmp;
  vector<Real8> es_vec;

  xs_pot = 0.0;

  if( nro_multi == 0 && ap_multi < 0.0)
  {
    string class_name = "ResonanceXSCalculator";
    string func_name  = "calc_reso_xs_unreso_c(int i, int j, Real8& ene_val, vector<Real8>& sig_val)";

    ostringstream oss01, oss02, oss03, oss04, oss05, oss06;
    oss01 << i;
    oss02 << nis;
    oss03 << j;
    oss04 << ner[i];
    oss05 << nro_multi;
    oss06 << ap_multi;
    string str_data01 = "Number of istopes (NIS)                : " + oss01.str() + " / " + oss02.str();
    string str_data02 = "Number of resonance energy range (NER) : " + oss03.str() + " / " + oss04.str();
    string str_data03 = "scat_radius_ene_dependence_flg (NRO)   : " + oss05.str();
    string str_data04 = "scat_radius (AP)                       : " + oss06.str();

    vector<string> err_com;
    err_com.push_back(str_data01);
    err_com.push_back(str_data02);
    err_com.push_back(str_data03);
    err_com.push_back(str_data04);
    err_com.push_back("scat_radius(AP) at reso_region_flg(LRU)=2 is less than 0.0.");
    err_com.push_back("scat_radius (AP) or scat_radius_ene_dependence_flg(NRO) value is not appropriate.");

    err_obj.output_runtime_error(class_name, func_name, err_com);
  }
  
  sig_tmp.resize(xs_type_no); //total, scatter, fission, capture
  
  e_root  = sqrt(fabs(ene_val));
  wave_no = wave_no_part_multi[0] * e_root;
  const_val = 2.0*M_PI*M_PI / (wave_no*wave_no);
  
  int l_max = nls_multi;
  for(int l=0; l<l_max; l++)
  {
    l_val = l_multi[l];
    calc_rho(l, ene_val, wave_no, rho, rho_h);
    calc_phase_shift(l_val, rho_h, phi);
    sph = sin(phi);
    
    int m_max = njs_multi[l];
    for(int m=0; m<m_max; m++)
    {
      aj   = aj_multi[l][m];
      amux = amux_multi[l][m];
      amun = amun_multi[l][m];
      amuf = amuf_multi[l][m];
      
      es_vec = es_multi[l][m];
      es_pos = -1;
      for(int n=1; n<static_cast<int>(es_vec.size()); n++)
      {
        if( es_vec[n] > ene_val )
        {
          es_pos = n;
          break;
        }
      }
      if( es_pos < 0 )
      {
        es_pos = static_cast<int>(es_vec.size()) - 1;
      }
      
      es_m = es_vec[es_pos-1];
      es_p = es_vec[es_pos];
      d_m  = d_c_multi[l][m][es_pos-1];
      d_p  = d_c_multi[l][m][es_pos];
      gx_m = gx_c_multi[l][m][es_pos-1];
      gx_p = gx_c_multi[l][m][es_pos];
      gn_m = gno_c_multi[l][m][es_pos-1];
      gn_p = gno_c_multi[l][m][es_pos];
      gg_m = gg_c_multi[l][m][es_pos-1];
      gg_p = gg_c_multi[l][m][es_pos];
      gf_m = gf_c_multi[l][m][es_pos-1];
      gf_p = gf_c_multi[l][m][es_pos];

      //Modify es, d, gx, gn, gg and gf
      if( es_m < min_value )
      {
        es_m = 1.0E-12;
      }
      if( es_p < min_value )
      {
        es_p = 1.0E-12;
      }
      if( d_m < min_value )
      {
        d_m = 1.0E-12;
      }
      if( d_p < min_value )
      {
        d_p = 1.0E-12;
      }
      if( gx_m < min_value )
      {
        gx_m = 1.0E-12;
      }
      if( gx_p < min_value )
      {
        gx_p = 1.0E-12;
      }
      if( gn_m < min_value )
      {
        gn_m = 1.0E-12;
      }
      if( gn_p < min_value )
      {
        gn_p = 1.0E-12;
      }
      if( gg_m < min_value )
      {
        gg_m = 1.0E-12;
      }
      if( gg_p < min_value )
      {
        gg_p = 1.0E-12;
      }
      if( gf_m < min_value )
      {
        gf_m = 1.0E-12;
      }
      if( gf_p < min_value )
      {
        gf_p = 1.0E-12;
      }
      
      ti_obj.interpolation_1d(int_val, ene_val, d_val, es_m, d_m,  es_p, d_p);
      ti_obj.interpolation_1d(int_val, ene_val, gx,    es_m, gx_m, es_p, gx_p);
      ti_obj.interpolation_1d(int_val, ene_val, gn,    es_m, gn_m, es_p, gn_p);
      ti_obj.interpolation_1d(int_val, ene_val, gg,    es_m, gg_m, es_p, gg_p);
      ti_obj.interpolation_1d(int_val, ene_val, gf,    es_m, gf_m, es_p, gf_p);
      
      if( gx < 1.0e-8 )
      {
        gx = 0.0;
      }
      if( gf < 1.0e-8 )
      {
        gf = 0.0;
      }
      
      calc_penetrability_factor(l_val, rho, amun, vl);
      vl    = vl*e_root;
      gnx   = gn * vl;
      gj    = const_val*gnx*(2.0*aj + 1.0) * den_multi / d_val;
      
      mu    = static_cast<int>(round(amun));
      nu    = static_cast<int>(round(amuf));
      lamda = static_cast<int>(round(amux));
      calc_width_fluctuation_factor( gnx, gf, gg, mu, nu, lamda, sig_tmp, gx );
      
      sig_val[scatter_xs]   += sig_tmp[scatter_xs]  *gj*gnx - 2.0*gj*sph*sph;
      sig_val[fission_xs]   += sig_tmp[fission_xs]  *gj*gf;
      sig_val[radiation_xs] += sig_tmp[radiation_xs]*gj*gg;
    }
    spot = sph / wave_no;
    spot = 4.0*M_PI*(2.0*static_cast<Real8>(l_val) + 1.0) * spot*spot;
    sig_val[scatter_xs] += spot; //add potential scatter xs
    xs_pot              += spot;
  }
  sig_val[total_xs] = sig_val[scatter_xs] + sig_val[fission_xs] + sig_val[radiation_xs];
  
  sig_tmp.clear();
  es_vec.clear();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_rho(int k, Real8& ene, Real8& wave_no, Real8& rho, Real8& rho_h)
{
  if( nro_multi == 0 )
  {
    rho_h = wave_no * ap_multi;
    
    if( radius_calc_flg_multi == 0 )
    {
      rho = wave_no * chan_rad_multi[k];
    }
    else
    {
      rho = rho_h;
    }
  }
  else
  {
    ti_obj.interpolation_tab1(ene, rho_h, nbt_nro_multi, int_nro_multi, e_int_nro_multi, ap_nro_multi);
    rho_h *= wave_no;
    
    if( radius_calc_flg_multi == 0 )
    {
      rho = wave_no * chan_rad_multi[k];
    }
    else if( radius_calc_flg_multi == 1 )
    {
      if( reso_region_flg_multi != 2 )
      {
        rho = rho_h;
      }
      else
      {
        if( ap_multi > 0.0 )
        {
          rho = ap_multi;
        }
        else
        {
          rho = rho_h;
        }
      }
    }
    else
    {
      rho = wave_no * ap_multi;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_rho(int i, int j, int k, 
                                     Real8& ene, Real8& wave_no, Real8& rho, Real8& rho_h)
{
  if( nro[i][j] == 0 )
  {
    rho_h = wave_no * ap[i][j];
    
    if( radius_calc_flg[i][j] == 0 )
    {
      rho = wave_no * chan_rad[i][j][k];
    }
    else
    {
      rho = rho_h;
    }
  }
  else
  {
    ti_obj.interpolation_tab1(ene, rho_h, nbt_nro[i][j], int_nro[i][j], e_int_nro[i][j], ap_nro[i][j]);
    rho_h *= wave_no;
    
    if( radius_calc_flg[i][j] == 0 )
    {
      rho = wave_no * chan_rad[i][j][k];
    }
    else if( radius_calc_flg[i][j] == 1 )
    {
      if( reso_region_flg[i][j] != 2 )
      {
        rho = rho_h;
      }
      else
      {
        if( ap[i][j] > 0.0 )
        {
          rho = ap[i][j];
        }
        else
        {
          rho = rho_h;
        }
      }
    }
    else
    {
      rho = wave_no * ap[i][j];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_shift_penetration_factor(Integer& l, Real8& rho, Real8& sl, Real8& pl)
{
  if( l==0 )
  {
    sl = 0.0;
    pl = rho;
  }
  else if( l==1 )
  {
    Real8 rho_2 = rho*rho;
    Real8 den   = 1.0 / (1.0+rho_2);
    sl = -1.0*den;
    pl = rho*rho_2*den;
  }
  else if( l==2 )
  {
    Real8 rho_2 = rho  *rho;
    Real8 rho_4 = rho_2*rho_2;
    Real8 den   = 1.0 / (9.0+3.0*rho_2+rho_4);
    sl = -1.0*(18.0+3.0*rho_2)*den;
    pl = (rho*rho_4)*den;
  }
  else if( l==3 )
  {
    Real8 rho_2 = rho  *rho;
    Real8 rho_4 = rho_2*rho_2;
    Real8 rho_6 = rho_2*rho_4;
    Real8 den   = 1.0 / (225.0+45.0*rho_2+6.0*rho_4+rho_6);
    sl = -1.0*(675.0+90.0*rho_2+6.0*rho_4)*den;
    pl = (rho*rho_6)*den;
  }
  else if( l==4 )
  {
    Real8 rho_2 = rho  *rho;
    Real8 rho_4 = rho_2*rho_2;
    Real8 rho_6 = rho_2*rho_4;
    Real8 rho_8 = rho_4*rho_4;
    Real8 den   = 1.0 / (11025.0+1575.0*rho_2+135.0*rho_4+10.0*rho_6+rho_8);
    sl = -1.0*(44100.0+4725.0*rho_2+270.0*rho_4+10.0*rho_6)*den;
    pl = (rho*rho_8)*den;
  }
  else
  {
    string class_name = "ResonanceXSCalculator";
    string func_name  = "calc_shift_penetration_factor(Integer& l, Real8& rho, Real8& sl, Real8& pl)";

    ostringstream oss01, oss02, oss03, oss04;
    oss01 << l;
    oss02 << rho;
    oss03 << sl;
    oss04 << pl;
    string str_data01 = "l value   : " + oss01.str();
    string str_data02 = "rho value : " + oss02.str();
    string str_data03 = "sl value  : " + oss03.str();
    string str_data04 = "pl value  : " + oss04.str();

    vector<string> err_com;
    err_com.push_back(str_data01);
    err_com.push_back(str_data02);
    err_com.push_back(str_data03);
    err_com.push_back(str_data04);
    err_com.push_back("This l value can not treat in this program.");
    err_com.push_back("In this program, the maximum l value is 4.");
    err_com.push_back("Please check or modify this function.");

    err_obj.output_runtime_error(class_name, func_name, err_com);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_phase_shift(Integer& l, Real8& rho, Real8& phi)
{
  if( l==0 )
  {
    phi = rho;
  }
  else if( l==1 )
  {
    phi = rho - atan(rho);
  }
  else if( l==2 )
  {
    phi = rho - atan(3.0*rho/(3.0-rho*rho));
  }
  else if( l==3 )
  {
    Real8 rho_2 = rho*rho;
    phi = rho - atan(rho*(15.0-rho_2)/(15.0-6.0*rho_2));
  }
  else if( l==4 )
  {
    Real8 rho_2 = rho  *rho;
    Real8 rho_4 = rho_2*rho_2;
    phi = rho - atan(rho*(105.0-10.0*rho_2)/(105.0-45.0*rho_2+rho_4));
  }
  else
  {
    string class_name = "ResonanceXSCalculator";
    string func_name  = "calc_phase_shift(Integer& l, Real8& rho, Real8& phi)";

    ostringstream oss01, oss02, oss03;
    oss01 << l;
    oss02 << rho;
    oss03 << phi;
    string str_data01 = "l value   : " + oss01.str();
    string str_data02 = "rho value : " + oss02.str();
    string str_data03 = "phi value : " + oss03.str();

    vector<string> err_com;
    err_com.push_back(str_data01);
    err_com.push_back(str_data02);
    err_com.push_back(str_data03);
    err_com.push_back("This l value can not treat in this program.");
    err_com.push_back("In this program, the maximum l value is 4.");
    err_com.push_back("Please check or modify this function.");

    err_obj.output_runtime_error(class_name, func_name, err_com);
  }

  if( phi/rho < min_value )
  {
    phi = 0.0;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_channel_radius(int i, int j, Real8& ene, Real8& radii)
{
  if( nro[i][j] == 0 )
  {
    if( radius_calc_flg[i][j] == 0 )
    {
      radii = chan_rad[i][j][0];
    }
    else
    {
      radii = ap[i][j];
    }
  }
  else
  {
    Real8 ap_e;
    ti_obj.interpolation_tab1(ene, ap_e, nbt_nro[i][j], int_nro[i][j], e_int_nro[i][j], ap_nro[i][j]);

    if( radius_calc_flg[i][j] == 0 )
    {
      radii = chan_rad[i][j][0];
    }
    else if( radius_calc_flg[i][j] == 1 )
    {
      if( ap[i][j] > 0.0 )
      {
        radii = ap[i][j];
      }
      else
      {
        ti_obj.interpolation_tab1(ene, radii, nbt_nro[i][j], int_nro[i][j], e_int_nro[i][j], ap_nro[i][j]);
      }
    }
    else
    {
      radii = ap[i][j];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_penetrability_factor(Integer& l, Real8& rho, Real8& amun, Real8& vl)
{
  if( l==0 )
  {
    vl = 1.0;
  }
  else if( l==1 )
  {
    Real8 rho_2 = rho*rho;
    vl = rho_2 / (1.0 + rho_2);
  }
  else if( l==2 )
  {
    Real8 rho_2 = rho   * rho;
    Real8 rho_4 = rho_2 * rho_2;
    vl = rho_4 / (9.0 + 3.0*rho_2 + rho_4);
  }
  else if( l==3 )
  {
    Real8 rho_2 = rho   * rho;
    Real8 rho_4 = rho_2 * rho_2;
    Real8 rho_6 = rho_4 * rho_2;
    vl = rho_6 / (225.0 + 45.0*rho_2 + 6.0*rho_4 + rho_6);
  }
  else if( l==4 )
  {
    Real8 rho_2 = rho   * rho;
    Real8 rho_4 = rho_2 * rho_2;
    Real8 rho_6 = rho_4 * rho_2;
    Real8 rho_8 = rho_4 * rho_4;
    vl = rho_8 / (11025.0 + 1575.0*rho_2 + 135.0*rho_4 + 10.0*rho_6 + rho_8);
  }
  else
  {
    string class_name = "ResonanceXSCalculator";
    string func_name  = "calc_penetrability_factor(Integer& l, Real8& rho, Real8& amun, Real8& vl)";

    ostringstream oss01, oss02, oss03;
    oss01 << l;
    oss02 << rho;
    oss03 << amun;
    string str_data01 = "l value                                     : " + oss01.str();
    string str_data02 = "rho value                                   : " + oss02.str();
    string str_data03 = "No. of degrees of freedom in neutron width  : " + oss03.str();

    vector<string> err_com;
    err_com.push_back(str_data01);
    err_com.push_back(str_data02);
    err_com.push_back(str_data03);
    err_com.push_back("This l value can not treat in this program.");
    err_com.push_back("In this program, the maximum l value is 4.");
    err_com.push_back("Please check or modify this function.");

    err_obj.output_runtime_error(class_name, func_name, err_com);
  }

  vl *= amun;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_width_fluctuation_factor(Real8& gn, Real8& gf, Real8& gg, 
                                                          int mu, int nu, int lamda, 
                                                          vector<Real8>& sig_val, Real8 df )
{
  for(int i=0; i<xs_type_no; i++)
  {
    sig_val[i] = 0.0;
  }
  
  if( gn <= 0.0 || gg <= 0.0 || gf < 0.0 )
  {
    return;
  }
  else if( gf > 0.0 && df < 0.0 )
  {
    return;
  }
  
  Real8 real_val, coef_i_sc, coef_i_cap, coef_j_sc, coef_j_cap, coef_j_fis, den_i, den_j;
  mu--; //Fortran:[1-10], C++[0-9]
  nu--;
  lamda--;
  if( gf == 0.0 )
  {
    if( df == 0.0 )
    {
      for(int i=0; i<10; i++)
      {
        real_val = q_wei[mu][i]*q_abs[mu][i] / (gn*q_abs[mu][i] + gg);
        sig_val[scatter_xs]   += real_val * q_abs[mu][i];
        sig_val[radiation_xs] += real_val;
      }
    }
    else if( df > 0.0 )
    {
      for(int i=0; i<10; i++)
      {
        den_i      = gn*q_abs[mu][i] + gg;
        coef_i_cap = q_wei[mu][i]*q_abs[mu][i];
        coef_i_sc  = coef_i_cap *q_abs[mu][i];
        for(int k=0; k<10; k++)
        {
          real_val = q_wei[lamda][k] / (den_i + df*q_abs[lamda][k]);
          sig_val[scatter_xs]   += real_val * coef_i_sc;
          sig_val[radiation_xs] += real_val * coef_i_cap;
        }
      }
    }
  }
  else if( gf > 0.0 )
  {
    if( df == 0.0 )
    {
      for(int i=0; i<10; i++)
      {
        den_i      = gn*q_abs[mu][i] + gg;
        coef_i_cap = q_wei[mu][i]*q_abs[mu][i];
        coef_i_sc  = coef_i_cap *q_abs[mu][i];
        for(int j=0; j<10; j++)
        {
          real_val = q_wei[nu][j] / (den_i + gf*q_abs[nu][j]);
          sig_val[scatter_xs]   += real_val * coef_i_sc;
          sig_val[fission_xs]   += real_val * coef_i_cap * q_abs[nu][j];
          sig_val[radiation_xs] += real_val * coef_i_cap;
        }
      }
    }
    else if( df > 0.0 )
    {
      for(int i=0; i<10; i++)
      {
        den_i      = gn*q_abs[mu][i] + gg;
        coef_i_cap = q_wei[mu][i]*q_abs[mu][i];
        coef_i_sc  = coef_i_cap *q_abs[mu][i];
        for(int j=0; j<10; j++)
        {
          den_j      = den_i + gf*q_abs[nu][j];
          coef_j_sc  = coef_i_sc  * q_wei[nu][j];
          coef_j_cap = coef_i_cap * q_wei[nu][j];
          coef_j_fis = coef_i_cap * q_wei[nu][j] * q_abs[nu][j];
          for(int k=0; k<10; k++)
          {
            real_val = q_wei[lamda][k] / (den_j + df*q_abs[lamda][k]);
            sig_val[scatter_xs]   += real_val * coef_j_sc;
            sig_val[fission_xs]   += real_val * coef_j_fis;
            sig_val[radiation_xs] += real_val * coef_j_cap;
          }
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::set_multi_array_data(int i, int j)
{
  if( i_pre != i || j_pre != j )
  {
    react_type_list[i][j].resize(xs_type_no);
    react_type_list[i][j][0] =   1; //Total
    react_type_list[i][j][1] =   2; //Elastic scatter
    react_type_list[i][j][2] =  18; //Fission
    react_type_list[i][j][3] = 102; //Radiation

    q_array[i][j].resize(xs_type_no);
    q_array[i][j][0] = 0.0; //Total
    q_array[i][j][1] = 0.0; //Elastic scatter
    q_array[i][j][2] = 0.0; //Fission
    q_array[i][j][3] = 0.0; //Radiation

    clear_multi_array_data();
    i_pre = i;
    j_pre = j;
    
    if( fabs(abn[i] - 1.0) < min_ene_dif )
    {
      abn_multi  = -1.0;
    }
    else
    {
      abn_multi  = abn[i];
    }
    reso_region_flg_multi = reso_region_flg[i][j];
    xs_formula_flg_multi  = xs_formula_flg[i][j];
    radius_calc_flg_multi = radius_calc_flg[i][j];
    ap_multi              = ap[i][j];
    nro_multi             = nro[i][j];
    if( nro_multi != 0 )
    {
      nbt_nro_multi   = nbt_nro[i][j];
      int_nro_multi   = int_nro[i][j];
      e_int_nro_multi = e_int_nro[i][j];
      ap_nro_multi    = ap_nro[i][j];
    }
    
    if( reso_region_flg_multi == 1 )
    {
      switch(xs_formula_flg_multi)
      {
        case 1: //Single Level Breit-Wigner
          set_multi_array_data_bw(i, j);
          break;
        case 2: //Multi Level Breit-Wigner
          set_multi_array_data_bw(i, j);
          break;
        case 3: //Reich-Moore
          set_multi_array_data_rm(i, j);
          break;
        case 4: //Adler-Adler
          set_multi_array_data_adler(i, j);
          break;
        case 7: //R-Matrix Limited Format
          set_multi_array_data_r_matrix(i, j);
          break;
        default:
          string class_name = "ResonanceXSCalculator";
          string func_name  = "set_multi_array_data(int i, int j)";

          ostringstream oss01, oss02, oss03, oss04, oss05, oss06;
          oss01 << i;
          oss02 << nis;
          oss03 << j;
          oss04 << ner[i];
          oss05 << reso_region_flg_multi;
          oss06 << xs_formula_flg_multi;
          string str_data01 = "Number of isotopes (NIS)                : " + oss01.str() + " / " + oss02.str();
          string str_data02 = "Number of resonance energy region (NER) : " + oss03.str() + " / " + oss04.str();
          string str_data03 = "Resonance region flg (LRU)              : " + oss05.str();
          string str_data04 = "XS formula flg (LRF)                    : " + oss06.str();

          vector<string> err_com;
          err_com.push_back(str_data01);
          err_com.push_back(str_data02);
          err_com.push_back(str_data03);
          err_com.push_back(str_data04);
          err_com.push_back("This xs_formula_flg (LRF) value is not applicable in this program.");
          err_com.push_back("Supported xs_formula_flg (LRF) value is 1 - 4 and 7.");

          err_obj.output_runtime_error(class_name, func_name, err_com);
          break;
      }
    }
    else if( reso_region_flg_multi==2 ) //(Unresolved resonance)
    {
      if( xs_formula_flg_multi == 1 )
      {
        //Case A(fis_width_flg(lfw)=0), Case B(fis_width_flg(lfw)=1)
        set_multi_array_data_unreso_a(i, j);
      }
      else //xs_formula_flg_multi==2
      {
        //Case C
        set_multi_array_data_unreso_c(i, j);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::set_multi_array_data_bw(int i, int j)
{
  BreitWignerDataContainer bw_data_obj = reso_data_obj.get_bw_data_obj(i,j);

  spi_multi          = spi[i][j];
  wave_no_part_multi = wave_no_part[i][j];
  chan_rad_multi     = chan_rad[i][j];
  nls_multi          = nls[i][j];
  cp_obj.copy_vec_array1_real8(qx_multi, bw_data_obj.get_q_value());
  l_multi            = reso_data_obj.get_l_value()[i][j];
  lrx_multi          = bw_data_obj.get_competitive_width_flg();
  cp_obj.copy_vec_array2_real8(er_multi, bw_data_obj.get_ene_reso());
  cp_obj.copy_vec_array2_real8(gn_multi, bw_data_obj.get_gam_width_n());
  cp_obj.copy_vec_array2_real8(gg_multi, bw_data_obj.get_gam_width_rad());
  cp_obj.copy_vec_array2_real8(gf_multi, bw_data_obj.get_gam_width_fis());
  cp_obj.copy_vec_array2_real8(aj_multi, bw_data_obj.get_j_value_abs());
  sl_er_multi        = sl_er[i][j];
  gx_er_multi        = gx_er[i][j];
  bw_data_obj.clear();

  if( xs_formula_flg_multi == 1 )
  {
    den_multi  = 1.0/(2.0*spi_multi+1.0);
    awri_multi = awri[i][j];
  }
  else
  {
    den_multi  = 1.0/(4.0*spi_multi+2.0);
  }

  int k_max = static_cast<int>(er_multi.size());
  nrs_multi.resize(k_max);
  for(int k=0; k<k_max; k++)
  {
    nrs_multi[k] = static_cast<int>(er_multi[k].size());
  }

  k_max = static_cast<int>(pl_er[i][j].size());
  pl_er_inv_multi.resize(k_max);
  for(int k=0; k<k_max; k++)
  {
    int l_max = static_cast<int>(pl_er[i][j][k].size());
    pl_er_inv_multi[k].resize(l_max);
    for(int l=0; l<l_max; l++)
    {
      pl_er_inv_multi[k][l] = 1.0/pl_er[i][j][k][l];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::set_multi_array_data_rm(int i, int j)
{
  ReichMooreDataContainer rm_data_obj = reso_data_obj.get_rm_data_obj(i,j);

  spi_multi          = spi[i][j];
  den_multi          = 1.0/(4.0*spi_multi+2.0);
  wave_no_part_multi = wave_no_part[i][j];
  chan_rad_multi     = chan_rad[i][j];
  nls_multi          = nls[i][j];
  l_multi            = reso_data_obj.get_l_value()[i][j];
  cp_obj.copy_vec_array1_real8(apl_multi, rm_data_obj.get_scat_radius_l());
  cp_obj.copy_vec_array2_real8(aj_multi,  rm_data_obj.get_j_value_abs());
  cp_obj.copy_vec_array2_real8(er_multi,  rm_data_obj.get_ene_reso());
  cp_obj.copy_vec_array2_real8(gn_multi,  rm_data_obj.get_gam_width_n());
  cp_obj.copy_vec_array2_real8(gg_multi,  rm_data_obj.get_gam_width_rad());
  cp_obj.copy_vec_array2_real8(gfa_multi, rm_data_obj.get_gam_width_fis_a());
  cp_obj.copy_vec_array2_real8(gfb_multi, rm_data_obj.get_gam_width_fis_b());
  rm_data_obj.clear();

  int k_max = static_cast<int>(er_multi.size());
  nrs_multi.resize(k_max);
  for(int k=0; k<k_max; k++)
  {
    nrs_multi[k] = static_cast<int>(er_multi[k].size());
  }

  k_max = static_cast<int>(gfa_multi.size());
  gfa_root_multi.resize(k_max);
  gfb_root_multi.resize(k_max);
  for(int k=0; k<k_max; k++)
  {
    int l_max = static_cast<int>(gfa_multi[k].size());
    gfa_root_multi[k].resize(l_max);
    gfb_root_multi[k].resize(l_max);
    for(int l=0; l<l_max; l++)
    {
      if( gfa_multi[k][l] > 0.0 )
      {
        gfa_root_multi[k][l] = sqrt(gfa_multi[k][l]);
      }
      else if( gfa_multi[k][l] < 0.0 )
      {
        gfa_root_multi[k][l] = -1.0 * sqrt(-1.0*gfa_multi[k][l]);
      }
      else
      {
        gfa_root_multi[k][l] = 0.0;
      }

      if( gfb_multi[k][l] > 0.0 )
      {
        gfb_root_multi[k][l] = sqrt(gfb_multi[k][l]);
      }
      else if( gfb_multi[k][l] < 0.0 )
      {
        gfb_root_multi[k][l] = -1.0 * sqrt(-1.0*gfb_multi[k][l]);
      }
      else
      {
        gfb_root_multi[k][l] = 0.0;
      }
    }
  }

  k_max = static_cast<int>(pl_er[i][j].size());
  pl_er_inv_multi.resize(k_max);
  for(int k=0; k<k_max; k++)
  {
    int l_max = static_cast<int>(pl_er[i][j][k].size());
    pl_er_inv_multi[k].resize(l_max);
    for(int l=0; l<l_max; l++)
    {
      pl_er_inv_multi[k][l] = 1.0/pl_er[i][j][k][l];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::set_multi_array_data_adler(int i, int j)
{
  AdlerAdlerDataContainer adler_data_obj = reso_data_obj.get_adler_data_obj(i,j);

  awri_multi         = awri[i][j];
  wave_no_part_multi = wave_no_part[i][j];
  chan_rad_multi     = chan_rad[i][j];
  li_multi           = adler_data_obj.get_adler_calc_flg();
  cp_obj.copy_vec_array1_real8(at_multi,  adler_data_obj.get_back_ground_tot());
  cp_obj.copy_vec_array1_real8(af_multi,  adler_data_obj.get_back_ground_fis());
  cp_obj.copy_vec_array1_real8(ac_multi,  adler_data_obj.get_back_ground_rad());
  cp_obj.copy_vec_array3_real8(det_multi, adler_data_obj.get_ene_reso_tot());
  cp_obj.copy_vec_array3_real8(dwt_multi, adler_data_obj.get_gam_width_half_tot());
  cp_obj.copy_vec_array3_real8(grt_multi, adler_data_obj.get_symmetrical_data_tot());
  cp_obj.copy_vec_array3_real8(git_multi, adler_data_obj.get_non_symmetrical_data_tot());
  cp_obj.copy_vec_array3_real8(def_multi, adler_data_obj.get_ene_reso_fis());
  cp_obj.copy_vec_array3_real8(dwf_multi, adler_data_obj.get_gam_width_half_fis());
  cp_obj.copy_vec_array3_real8(grf_multi, adler_data_obj.get_symmetrical_data_fis());
  cp_obj.copy_vec_array3_real8(gif_multi, adler_data_obj.get_non_symmetrical_data_fis());
  cp_obj.copy_vec_array3_real8(dec_multi, adler_data_obj.get_ene_reso_rad());
  cp_obj.copy_vec_array3_real8(dwc_multi, adler_data_obj.get_gam_width_half_rad());
  cp_obj.copy_vec_array3_real8(grc_multi, adler_data_obj.get_symmetrical_data_rad());
  cp_obj.copy_vec_array3_real8(gic_multi, adler_data_obj.get_non_symmetrical_data_rad());
  adler_data_obj.clear();

  nls_multi = nls[i][j];

  int k_max = static_cast<int>(det_multi.size());
  njs_multi.resize(k_max);
  nlj_multi.resize(k_max);
  for(int k=0; k<k_max; k++)
  {
    int l_max = static_cast<int>(det_multi[k].size());
    njs_multi[k] = l_max;
    nlj_multi[k].resize(l_max);
    for(int l=0; l<l_max; l++)
    {
      nlj_multi[k][l] = static_cast<Integer>(det_multi[k][l].size());
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::set_multi_array_data_r_matrix(int i, int j)
{
#ifndef NO_AMUR_MODE
  r_matrix_calc_obj.set_reso_data_obj( nucl_data_obj.get_reso_data_obj() );

  r_matrix_calc_obj.set_resonance_data(i, j);

  react_type_list[i][j].clear();
  react_type_list[i][j] = r_matrix_calc_obj.get_react_type_list();

  q_array[i][j].clear();
  q_array[i][j] = r_matrix_calc_obj.get_q_val();
#endif //NO_AMUR_MODE
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::set_multi_array_data_unreso_a(int i, int j)
{
  UnresolvedResonanceDataContainer unreso_data_obj = reso_data_obj.get_unreso_data_obj(i,j);

  fis_width_flg_multi = fis_width_flg[i];
  spi_multi           = spi[i][j];
  den_multi           = 1.0/(4.0*spi_multi+2.0);
  es_multi            = unreso_data_obj.get_ene_unreso();
  wave_no_part_multi  = wave_no_part[i][j];
  chan_rad_multi      = chan_rad[i][j];
  ne_multi            = static_cast<Integer>(es_multi.size());
  nls_multi           = nls[i][j];
  l_multi             = reso_data_obj.get_l_value()[i][j];

  cp_obj.copy_vec_array2_real8(aj_multi,   unreso_data_obj.get_j_value_abs());
  cp_obj.copy_vec_array2_real8(amun_multi, unreso_data_obj.get_freedom_n());
  
  if( fis_width_flg_multi != 0 ) //Case B
  {
    muf_multi  = unreso_data_obj.get_freedom_fis_int();
    gf_b_multi = unreso_data_obj.get_ave_gam_width_fis();
  }

  int k_max = static_cast<int>(aj_multi.size());
  njs_multi.resize(k_max);
  dx_multi.resize(k_max);
  gnox_multi.resize(k_max);
  ggx_multi.resize(k_max);
  for(int k=0; k<k_max; k++)
  {
    njs_multi[k] = static_cast<int>(aj_multi[k].size());
    
    int l_max = njs_multi[k];
    dx_multi[k].resize(l_max);
    gnox_multi[k].resize(l_max);
    ggx_multi[k].resize(l_max);
    for(int l=0; l<l_max; l++)
    {
      dx_multi[k][l]   = unreso_data_obj.get_level_spacing()[k][l][0];
      gnox_multi[k][l] = unreso_data_obj.get_ave_gam_width_n()[k][l][0];
      ggx_multi[k][l]  = unreso_data_obj.get_ave_gam_width_rad()[k][l][0];
    }
  }
  unreso_data_obj.clear();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::set_multi_array_data_unreso_c(int i, int j)
{
  UnresolvedResonanceDataContainer unreso_data_obj = reso_data_obj.get_unreso_data_obj(i,j);

  spi_multi          = spi[i][j];
  den_multi          = 1.0/(4.0*spi_multi+2.0);
  wave_no_part_multi = wave_no_part[i][j];
  chan_rad_multi     = chan_rad[i][j];
  nls_multi          = nls[i][j];
  l_multi            = reso_data_obj.get_l_value()[i][j];
  cp_obj.copy_vec_array2_real8(aj_multi,    unreso_data_obj.get_j_value_abs());
  cp_obj.copy_vec_array2_real8(amux_multi,  unreso_data_obj.get_freedom_comp());
  cp_obj.copy_vec_array2_real8(amun_multi,  unreso_data_obj.get_freedom_n());
  cp_obj.copy_vec_array2_real8(amuf_multi,  unreso_data_obj.get_freedom_fis());
  cp_obj.copy_vec_array3_real8(es_multi,    unreso_data_obj.get_ene_unreso());
  cp_obj.copy_vec_array3_real8(d_c_multi,   unreso_data_obj.get_level_spacing());
  cp_obj.copy_vec_array3_real8(gx_c_multi,  unreso_data_obj.get_ave_gam_width_comp());
  cp_obj.copy_vec_array3_real8(gno_c_multi, unreso_data_obj.get_ave_gam_width_n());
  cp_obj.copy_vec_array3_real8(gg_c_multi,  unreso_data_obj.get_ave_gam_width_rad());
  cp_obj.copy_vec_array3_real8(gf_c_multi,  unreso_data_obj.get_ave_gam_width_fis());
  unreso_data_obj.clear();

  int k_max = static_cast<int>(d_c_multi.size());
  njs_multi.resize(k_max);
  for(int k=0; k<k_max; k++)
  {
    njs_multi[k] = static_cast<int>(d_c_multi[k].size());
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::add_middle_energy_grid( int i, int j,
                                         vector<Real8>& ene_data, vector<vector<Real8> >& sig_data )
{
  int k_max = static_cast<int>(ene_data.size());
  if( k_max == 0 )
  {
    return;
  }

  int                    l;
  Integer                int_chk;
  Real8                  mid_ene;
  vector<Real8>          mid_sig;
  vector<Real8>          new_ene, new_ene_part;
  vector<vector<Real8> > new_sig, new_sig_part;
  new_ene.clear();
  new_ene_part.clear();
  clr_obj.clear_vec_array2_real8(new_sig);
  clr_obj.clear_vec_array2_real8(new_sig_part);

  new_ene.push_back(ene_data[0]);
  new_sig.push_back(sig_data[0]);
  for(int k=1; k<k_max; k++)
  {
    new_ene_part.push_back(ene_data[k-1]);
    new_sig_part.push_back(sig_data[k-1]);
    new_ene_part.push_back(ene_data[k]);
    new_sig_part.push_back(sig_data[k]);

    l = 1;
    int_chk = check_energy_grid_distance(i, j, l, new_ene_part, new_sig_part, mid_ene, mid_sig);
    while( int_chk < 0 || l < static_cast<int>(new_ene_part.size()-1))
    {
      if( int_chk >= 0 )
      {
        l++;
      }
      else
      {
        insert_middle_energy_grid(l, new_ene_part, new_sig_part, mid_ene, mid_sig);
      }
      int_chk = check_energy_grid_distance(i, j, l, new_ene_part, new_sig_part, mid_ene, mid_sig);
    }
    add_xs_at_each_grid(new_ene, new_sig, new_ene_part, new_sig_part);
  }
  ene_data = new_ene;
  sig_data = new_sig;
  new_ene.clear();
  clr_obj.clear_vec_array2_real8(new_sig);
  mid_sig.clear();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

Integer ResonanceXSCalculator::check_energy_grid_distance( int i, int j, int ele_no, 
                                                           vector<Real8>& ene_data,
                                                           vector<vector<Real8> >& sig_data,
                                                           Real8& mid_ene, vector<Real8>& mid_sig )
{
  Integer chk_flg = 0;
  Real8 err, err_max, err_int, delta_ene;
  vector<Real8> delta_sig;

  int l_max  = static_cast<int>(sig_data[ele_no].size());
  delta_sig.resize(l_max);
  
  if( ene_data[ele_no-1] < 0.5 )
  {
    err     = 0.2 * error_value;
    err_max = 0.2 * error_max_value;
    //err_int = 0.2 * error_int_value;
    err_int = error_int_value;
  }
  else
  {
    err     = error_value;
    err_max = error_max_value;
    err_int = error_int_value;
  }
  
  mid_ene = 0.5*(ene_data[ele_no] + ene_data[ele_no-1]);

  //If energy distance is so small, add middle energy grid is skipped.
  if(ene_data[ele_no] - mid_ene <= min_ene_dif*mid_ene || mid_ene - ene_data[ele_no-1] <= min_ene_dif*mid_ene)
  {
    return chk_flg;
  }

  calc_reso_xs_each_case(i, j, mid_ene, mid_sig);
  
  chk_flg = 0;
  delta_ene = fabs(ene_data[ele_no] - ene_data[ele_no-1]);
  for(int l=1; l<l_max; l++)
  {
    delta_sig[l] = fabs(mid_sig[l] - 0.5*(sig_data[ele_no][l] + sig_data[ele_no-1][l]));
    if( delta_sig[l] < err * mid_sig[l] )
    {
      chk_flg++;
    }
    else if( delta_sig[l] < min_value && mid_sig[l] < min_value )
    {
      chk_flg++;
    }
  }

  //int dif_pos = 0;
  //Real r_val1 = 0.0;
  //Real r_val2 = 0.0;
  if( chk_flg == 3 )
  {
    if( ele_no > 2 )
    {
      if( delta_ene > 4.1*(ene_data[ele_no-1] - ene_data[ele_no-2]) )
      {
        chk_flg = -3;
      }
    }
  }
  else
  {
    for(int l=1; l<l_max; l++)
    {
      if( delta_sig[l] > err_max * fabs(mid_sig[l]) )
      {
        chk_flg = -1;
        break;
      }
      else if( delta_sig[l]*delta_ene >= 2.0*err_int*mid_ene )
      {
        //dif_pos = l;
        //r_val1  = delta_sig[l];
        //r_val2  = 2.0*err_int*mid_ene/delta_ene;
        chk_flg = -2;
        break;
      }
    }
  }
  delta_sig.clear();
  
  return chk_flg;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::insert_middle_energy_grid(int ele_no, 
                                                      vector<Real8>& ene_data, vector<vector<Real8> >& sig_data,
                                                      Real8& mid_ene, vector<Real8>& mid_sig)
{
  ta_obj.add_table_data_at_given_position(ene_data, sig_data, mid_ene, mid_sig, ele_no);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::add_xs_at_each_grid
                              ( vector<Real8>& new_ene,      vector<vector<Real8> >& new_sig, 
                                vector<Real8>& new_ene_part, vector<vector<Real8> >& new_sig_part )
{
  int k_max = static_cast<int>(new_ene_part.size());
  for(int k=1; k<k_max; k++)
  {
    new_ene.push_back(new_ene_part[k]);
    new_sig.push_back(new_sig_part[k]);
  }
  new_ene_part.clear();
  clr_obj.clear_vec_array2_real8(new_sig_part);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::delete_overlap_grid(vector<Real8>& ene_data, vector<vector<Real8> >& sig_data)
{
  vector<Real8>          sorted_ene;
  vector<vector<Real8> > sorted_sig;
  sorted_ene.clear();
  clr_obj.clear_vec_array2_real8(sorted_sig);
  
  int i_max = static_cast<int>(ene_data.size());
  if( i_max > 0 )
  {
    Real8 ene_pre, ene_cur;
    ene_pre = ene_data[0];
    sorted_ene.push_back(ene_pre);
    sorted_sig.push_back(sig_data[0]);

    for(int i=1; i<static_cast<int>(ene_data.size()); i++)
    {
      ene_cur = ene_data[i];

      if( fabs(ene_cur - ene_pre) > min_ene_dif*ene_pre )
      {
        sorted_ene.push_back(ene_cur);
        sorted_sig.push_back(sig_data[i]);
        
        ene_pre = ene_cur;
      }
    }
    ene_data = sorted_ene;
    sig_data = sorted_sig;
    
    sorted_ene.clear();
    clr_obj.clear_vec_array2_real8(sorted_sig);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::unify_energy_grid()
{
  vector<Real8> unify_grid;
  unify_grid.clear();
  
  //Copy all energy grid
  int i_max = static_cast<int>(ene_array.size());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(ene_array[i].size());
    for(int j=0; j<j_max; j++)
    {
      int k_max = static_cast<int>(ene_array[i][j].size());
      for(int k=0; k<k_max; k++)
      {
        unify_grid.push_back(ene_array[i][j][k]);
      }
    }
  }
  
  //Sort energy grid
  i_max = static_cast<int>(unify_grid.size());
  if( i_max > 0 )
  {
    vector<Real8> new_grid;
    new_grid.clear();
  
    sort(unify_grid.begin(), unify_grid.end());

    Real8 ene_pre, ene_cur;
    ene_pre = unify_grid[0];
    new_grid.push_back(ene_pre);

    for(int i=1; i<i_max; i++)
    {
      ene_cur = unify_grid[i];

      if( fabs(ene_cur - ene_pre) > ene_pre*min_ene_dif )
      {
        new_grid.push_back(ene_cur);
        ene_pre = ene_cur;
      }
    }
    unify_grid = new_grid;
    new_grid.clear();
  }
  
  //Calc reso_grid at each energy grid
  int ene_max = unify_grid.size();

  i_max = static_cast<int>(ene_array.size());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(ene_array[i].size());
    for(int j=0; j<j_max; j++)
    {
      int     k_min   = 0;
      int     ele_no  = 0;
      int     ele_max = static_cast<int>(ene_array[i][j].size());
      Real8   ene_min = ene_array[i][j][0];
      Real8   ene_val = ene_min;
      for(int k=0; k<ene_max; k++)
      {
        if( ene_min - unify_grid[k] <= ene_min*min_ene_dif )
        {
          k_min = k;
          break;
        }
      }

      vector<Real8>          new_ene;
      vector<vector<Real8> > new_sig;
      for(int k=k_min; k<ene_max; k++)
      {
        new_ene.push_back(unify_grid[k]);
        if( fabs(unify_grid[k] - ene_val) <= ene_val*min_ene_dif )
        {
          new_sig.push_back(sig_array[i][j][ele_no]);
          ele_no++;
          if( ele_no == ele_max )
          {
            break;
          }
          ene_val = ene_array[i][j][ele_no];
        }
        else
        {
          vector<Real8> sig_val;
          calc_reso_xs_each_case(i, j, unify_grid[k], sig_val);
          new_sig.push_back(sig_val);
        }
      }

      ene_array[i][j] = new_ene;
      sig_array[i][j] = new_sig;
      new_ene.clear();
      clr_obj.clear_vec_array2_real8(new_sig);
    }
  }
  unify_grid.clear();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_unreso_xs()
{
  int i_max = static_cast<int>(ene_array.size());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(ene_array[i].size());
    for(int j=0; j<j_max; j++)
    {
      if( reso_region_flg[i][j] == 2 )
      {
        if( self_shielding_flg[i][j] != 0 )
        {
          clear_multi_array_data();
          set_multi_array_data(i, j);

          int k_max = static_cast<int>(ene_array[i][j].size());
          vector<Real8> sig_val;
          sig_val.resize(xs_type_no);

          if( xs_formula_flg[i][j] == 1 ) //Case A or B
          {
            for(int k=0; k<k_max; k++)
            {
              for(int l=0; l<xs_type_no; l++)
              {
                sig_val[l] = 0.0;
              }
              calc_reso_xs_unreso_a(i, j, ene_array[i][j][k], sig_val);
              sig_array[i][j][k] = sig_val;
            }
          }
          else //Case C
          {
            for(int k=0; k<k_max; k++)
            {
              for(int l=0; l<xs_type_no; l++)
              {
                sig_val[l] = 0.0;
              }
              calc_reso_xs_unreso_c(i, j, ene_array[i][j][k], sig_val);
              sig_array[i][j][k] = sig_val;
            }
          }
          sig_val.clear();
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::input_data_check()
{
  if( error_int_value < 0.0 )
  {
    error_int_value = error_value / 20000.0;
  }

  if( error_max_value < 0.0 )
  {
    error_max_value = error_value * 10.0;
  }
  else if( error_max_value < error_value )
  {
    string class_name = "ResonanceXSCalculator";
    string func_name  = "input_data_check()";

    ostringstream oss01, oss02;
    oss01 << error_value;
    oss02 << error_max_value;
    string str_data01 = "Error value         : " + oss01.str();
    string str_data02 = "Maximum error value : " + oss02.str();

    vector<string> err_com;
    err_com.push_back(str_data01);
    err_com.push_back(str_data02);
    err_com.push_back("Maximum error is less than error value.");

    err_obj.output_caution(class_name, func_name, err_com);
  }

  if( ( nucl_data_set_flg == 0 ) || ( error_value < min_ene_dif ) || ( error_max_value < min_ene_dif ) ||
      ( error_int_value < min_ene_dif ) || ( temp < 0.0 ) || ( temp < temp_nucl ) )
  {
    string class_name = "ResonanceXSCalculator";
    string func_name  = "input_data_check()";

    ostringstream oss01, oss02, oss03, oss04, oss05, oss06;
    oss01 << error_value;
    oss02 << error_max_value;
    oss03 << error_int_value;
    oss04 << min_ene_dif;
    oss05 << temp;
    oss06 << temp_nucl;
    string str_data01 = "Error value                     : " + oss01.str();
    string str_data02 = "Maximum error value             : " + oss02.str();
    string str_data03 = "Integral error value            : " + oss03.str();
    string str_data04 = "Minimum difference              : " + oss04.str();
    string str_data05 = "Temperature                     : " + oss05.str();
    string str_data06 = "Temperature in the nuclear data : " + oss06.str();

    vector<string> err_com;
    err_com.push_back(str_data01);
    err_com.push_back(str_data02);
    err_com.push_back(str_data03);
    err_com.push_back(str_data04);
    err_com.push_back(str_data05);
    err_com.push_back(str_data06);

    if( nucl_data_set_flg == 0 )
    {
      err_com.push_back(str_data06);
    }

    if( error_value < 0.0 )
    {
      err_com.push_back("Please set the positive error value.");
    }
    else if( error_value < min_ene_dif )
    {
      err_com.push_back("Error value is too small (less than minimum difference).");
      err_com.push_back("Please check the error value.");
    }

    if( error_max_value < min_ene_dif )
    {
      err_com.push_back("Maximum rror value is too small (less than minimum difference).");
      err_com.push_back("Please check the error value.");
    }

    if( error_int_value < min_ene_dif )
    {
      err_com.push_back("Integral rror value is too small (less than minimum difference).");
      err_com.push_back("Please check the error value.");
    }

    if( temp < 0.0 )
    {
      err_com.push_back("The temperature[K] is less than 0.0.");
      err_com.push_back("Please check the temperature.");
    }
    else if( temp < temp_nucl )
    {
      err_com.push_back("The input temperature is larger than the temperature in the nuclear data file.");
    }
    err_obj.output_runtime_error(class_name, func_name, err_com);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::clear_all()
{

  rpi = 0.0;
  clr_obj.clear_vec_array2_real8(q_wei);
  clr_obj.clear_vec_array2_real8(q_abs);

  clear();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::clear()
{
  temp            =    0.0;
  temp_nucl       = -100.0;
  error_value     = -1.0;
  error_max_value = -1.0;
  error_int_value = -1.0;

  clear_reso_data();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::clear_reso_data()
{
#ifndef NO_AMUR_MODE
  r_matrix_calc_obj.clear();
#endif //NO_AMUR_MODE

  clear_multi_array_data();
  
  calc_reso_xs_flg  = 0;
  nucl_data_set_flg = 0;

  clr_obj.clear_vec_array3_real8(ene_array);
  clr_obj.clear_vec_array4_real8(sig_array);

  reso_flg = 0;

  nis = 0;
  abn.clear();
  ner.clear();
  fis_width_flg.clear();
  clr_obj.clear_vec_array2_int(reso_region_flg);
  clr_obj.clear_vec_array2_int(xs_formula_flg);
  clr_obj.clear_vec_array2_int(nro);
  clr_obj.clear_vec_array2_int(radius_calc_flg);
  clr_obj.clear_vec_array2_int(nls);
  clr_obj.clear_vec_array2_int(self_shielding_flg);
  clr_obj.clear_vec_array2_real8(spi);
  clr_obj.clear_vec_array2_real8(ap);
  clr_obj.clear_vec_array3_real8(awri);
  
  clr_obj.clear_vec_array3_int(nbt_nro);
  clr_obj.clear_vec_array3_int(int_nro);
  clr_obj.clear_vec_array3_real(e_int_nro);
  clr_obj.clear_vec_array3_real(ap_nro);

  clr_obj.clear_vec_array3_real8(wave_no_part);
  clr_obj.clear_vec_array3_real8(chan_rad);
  clr_obj.clear_vec_array4_real8(gx_er);
  clr_obj.clear_vec_array4_real8(sl_er);
  clr_obj.clear_vec_array4_real8(pl_er);
  
  nucl_data_obj.clear();
  reso_data_obj.clear();

  ene_potential.clear();
  xs_potential.clear();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::clear_multi_array_data()
{
  i_pre     = -1;
  j_pre     = -1;
  nls_multi =  0;
  nrs_multi.clear();
  njs_multi.clear();
  int i_max = static_cast<int>(nlj_multi.size());
  for(int i=0; i<i_max; i++)
  {
    nlj_multi[i].clear();
  }
  nlj_multi.clear();
  
  reso_region_flg_multi = 0;
  xs_formula_flg_multi  = 0;
  fis_width_flg_multi   = 0;
  radius_calc_flg_multi = 0;
  li_multi              = 0;
  ne_multi              = 0;
  nro_multi             = 0;
  nbt_nro_multi.clear();
  int_nro_multi.clear();
  l_multi.clear();
  lrx_multi.clear();
  clr_obj.clear_vec_array2_int(muf_multi);
  
  e_int_nro_multi.clear();
  ap_nro_multi.clear();
  clr_obj.clear_vec_array3_real(gf_b_multi);
  
  abn_multi = 0.0;
  spi_multi = 0.0;
  den_multi = 0.0;
  ap_multi  = 0.0;
  wave_no_part_multi.clear();
  qx_multi.clear();
  awri_multi.clear();
  apl_multi.clear();
  chan_rad_multi.clear();
  at_multi.clear();
  af_multi.clear();
  ac_multi.clear();
  clr_obj.clear_vec_array2_real8(er_multi);
  clr_obj.clear_vec_array2_real8(gn_multi);
  clr_obj.clear_vec_array2_real8(gg_multi);
  clr_obj.clear_vec_array2_real8(gf_multi);
  clr_obj.clear_vec_array2_real8(gfa_multi);
  clr_obj.clear_vec_array2_real8(gfb_multi);
  clr_obj.clear_vec_array2_real8(gfa_root_multi);
  clr_obj.clear_vec_array2_real8(gfb_root_multi);
  clr_obj.clear_vec_array2_real8(aj_multi);
  clr_obj.clear_vec_array2_real8(sl_er_multi);
  clr_obj.clear_vec_array2_real8(pl_er_inv_multi);
  clr_obj.clear_vec_array2_real8(gx_er_multi);
  clr_obj.clear_vec_array2_real8(dx_multi);
  clr_obj.clear_vec_array2_real8(amun_multi);
  clr_obj.clear_vec_array2_real8(amux_multi);
  clr_obj.clear_vec_array2_real8(amuf_multi);
  clr_obj.clear_vec_array2_real8(gnox_multi);
  clr_obj.clear_vec_array2_real8(ggx_multi);
  clr_obj.clear_vec_array3_real8(det_multi);
  clr_obj.clear_vec_array3_real8(dwt_multi);
  clr_obj.clear_vec_array3_real8(grt_multi);
  clr_obj.clear_vec_array3_real8(git_multi);
  clr_obj.clear_vec_array3_real8(def_multi);
  clr_obj.clear_vec_array3_real8(dwf_multi);
  clr_obj.clear_vec_array3_real8(grf_multi);
  clr_obj.clear_vec_array3_real8(gif_multi);
  clr_obj.clear_vec_array3_real8(dec_multi);
  clr_obj.clear_vec_array3_real8(dwc_multi);
  clr_obj.clear_vec_array3_real8(grc_multi);
  clr_obj.clear_vec_array3_real8(gic_multi);
  clr_obj.clear_vec_array3_real8(es_multi);
  clr_obj.clear_vec_array3_real8(d_c_multi);
  clr_obj.clear_vec_array3_real8(gx_c_multi);
  clr_obj.clear_vec_array3_real8(gno_c_multi);
  clr_obj.clear_vec_array3_real8(gg_c_multi);
  clr_obj.clear_vec_array3_real8(gf_c_multi);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::set_nucl_data_obj(frendy::NuclearDataObject& data_obj)
{
  clear_reso_data();
  nucl_data_set_flg = 1;

  nucl_data_obj = data_obj;

  reso_flg = nucl_data_obj.get_general_data_obj().get_reso_flg();
  temp_nucl = static_cast<Real8>(nucl_data_obj.get_general_data_obj().get_temp());
  if( temp_nucl < -1.0*min_ene_val )
  {
    string class_name = "ResonanceXSCalculator";
    string func_name  = "set_endf6_parser_no_cov(Endf6ParserNoCov& parser_obj)";

    ostringstream oss01, oss02;
    oss01 << temp_nucl;
    string str_data01 = "Temperature in the nuclear data : " + oss01.str();

    vector<string> err_com;
    err_com.push_back(str_data01);
    err_com.push_back("The temperature in the nuclear data file is less than 0[K].");
    err_com.push_back("Please check the nuclear data file at MF01MT451.");

    err_obj.output_runtime_error(class_name, func_name, err_com);
  }

  reso_data_obj    = nucl_data_obj.get_reso_data_obj();
  fis_width_flg    = reso_data_obj.get_fis_width_flg();
  reso_region_flg  = reso_data_obj.get_reso_region_flg();
  xs_formula_flg   = reso_data_obj.get_xs_formula_flg();
  nro              = reso_data_obj.get_scat_radius_ene_dependence_flg();
  radius_calc_flg  = reso_data_obj.get_radius_calc_flg();
  nbt_nro          = reso_data_obj.get_scat_radius_tab_range_data();
  int_nro          = reso_data_obj.get_scat_radius_tab_int_data();
  e_int_nro        = reso_data_obj.get_scat_radius_tab_ene_data();
  ap_nro           = reso_data_obj.get_scat_radius_tab_data();

  int i_max = static_cast<int>(reso_region_flg.size());
  nis = i_max;
  ner.resize(i_max);
  nls.resize(i_max);
  abn.resize(i_max);
  spi.resize(i_max);
  ap.resize(i_max);
  awri.resize(i_max);
  wave_no_part.resize(i_max);
  chan_rad.resize(i_max);
  gx_er.resize(i_max);
  sl_er.resize(i_max);
  pl_er.resize(i_max);
  self_shielding_flg.resize(i_max);
  sig_array.resize(i_max);
  q_array.resize(i_max);
  react_type_list.resize(i_max);
  for(int i=0; i<i_max; i++)
  {
    abn[i] = static_cast<Real8>(reso_data_obj.get_abundance_isotope()[i]);
    
    int j_max = static_cast<int>(reso_region_flg[i].size());
    ner[i] = j_max;
    nls[i].resize(j_max);
    spi[i].resize(j_max);
    ap[i].resize(j_max);
    awri[i].resize(j_max);
    wave_no_part[i].resize(j_max);
    chan_rad[i].resize(j_max);
    gx_er[i].resize(j_max);
    sl_er[i].resize(j_max);
    pl_er[i].resize(j_max);
    self_shielding_flg[i].resize(j_max);
    sig_array[i].resize(j_max);
    q_array[i].resize(j_max);
    react_type_list[i].resize(j_max);
    for(int j=0; j<j_max; j++)
    {
      spi[i][j] = static_cast<Real8>(reso_data_obj.get_spin_data()[i][j]);
      ap[i][j]  = static_cast<Real8>(reso_data_obj.get_scat_radius()[i][j]);

      cp_obj.copy_vec_array1_real8(awri[i][j], reso_data_obj.get_mass_isotope()[i][j]);
      nls[i][j] = static_cast<int>(awri[i][j].size());

      self_shielding_flg[i][j] = 0;
      if( reso_region_flg[i][j] == 1 )
      {
        if( xs_formula_flg[i][j]==1 || xs_formula_flg[i][j]==2 )
        {
          calc_const_value_bw(i, j);
        }
        else if( xs_formula_flg[i][j]==3 )
        {
          calc_const_value_rm(i, j);
        }
        else if( xs_formula_flg[i][j]==4 )
        {
          nls[i][j] = static_cast<int>(reso_data_obj.get_l_value()[i][j].size());
          calc_const_value_adler(i, j);
        }
        else if( xs_formula_flg[i][j]==7 )
        {
          nls[i][j] = 0;
        }
      }
      else if( reso_region_flg[i][j] == 2 )
      {
        calc_const_value_unreso(i, j);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_const_value_bw(int i, int j)
{
  BreitWignerDataContainer bw_data_obj = reso_data_obj.get_bw_data_obj(i,j);

  int k_max = nls[i][j];
  Real8 awri_val, wave_no, rho, rho_h, er;
  Real8 gt, gn, gg, gf;
  wave_no_part[i][j].resize(k_max);
  chan_rad[i][j].resize(k_max);
  gx_er[i][j].resize(k_max);
  sl_er[i][j].resize(k_max);
  pl_er[i][j].resize(k_max);
  for(int k=0; k<k_max; k++)
  {
    awri_val              = awri[i][j][k];
    wave_no_part[i][j][k] = k_part*awri_val/(awri_val+1.0);
    chan_rad[i][j][k]     = 0.123 * pow(awri_val*amu_n,third)+0.08;
  }
  
  for(int k=0; k<k_max; k++)
  {
    Integer l_val   = reso_data_obj.get_l_value()[i][j][k];
    Integer lrx_val = bw_data_obj.get_competitive_width_flg()[k];
    int     l_max   = static_cast<int>(bw_data_obj.get_ene_reso()[k].size());

    gx_er[i][j][k].resize(l_max);
    sl_er[i][j][k].resize(l_max);
    pl_er[i][j][k].resize(l_max);
    for(int l=0; l<l_max; l++)
    {
      er = static_cast<Real8>(bw_data_obj.get_ene_reso()[k][l]);
      wave_no = wave_no_part[i][j][k] * sqrt(fabs(er));
      calc_rho(i, j, k, er, wave_no, rho, rho_h);
      calc_shift_penetration_factor(l_val, rho, sl_er[i][j][k][l], pl_er[i][j][k][l]);
      if( lrx_val != 0 )
      {
        gt = bw_data_obj.get_gam_width_tot()[k][l];
        gn = bw_data_obj.get_gam_width_n()[k][l];
        gg = bw_data_obj.get_gam_width_rad()[k][l];
        gf = bw_data_obj.get_gam_width_fis()[k][l];
        gx_er[i][j][k][l] = (gt - (gn + gg + gf)) / pl_er[i][j][k][l];
      }
      else
      {
        gx_er[i][j][k][l] = 0.0;
      }
    } 
  }
  bw_data_obj.clear();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_const_value_rm(int i, int j)
{
  ReichMooreDataContainer rm_data_obj = reso_data_obj.get_rm_data_obj(i,j);

  int k_max = nls[i][j];
  Real8 awri_val, apl, wave_no, rho, rho_h, er;
  wave_no_part[i][j].resize(k_max);
  chan_rad[i][j].resize(k_max);
  gx_er[i][j].resize(k_max);
  sl_er[i][j].resize(k_max);
  pl_er[i][j].resize(k_max);
  for(int k=0; k<k_max; k++)
  {
    awri_val              = awri[i][j][k];
    wave_no_part[i][j][k] = k_part*awri_val/(awri_val+1.0);
    chan_rad[i][j][k]     = 0.123 * pow(awri_val*amu_n,third)+0.08;
  }
  
  for(int k=0; k<k_max; k++)
  {
    Integer l_val   = reso_data_obj.get_l_value()[i][j][k];
    int     l_max   = static_cast<int>(rm_data_obj.get_ene_reso()[k].size());
    gx_er[i][j][k].resize(l_max);
    sl_er[i][j][k].resize(l_max);
    pl_er[i][j][k].resize(l_max);

    apl = static_cast<Real8>(rm_data_obj.get_scat_radius_l()[k]);
    for(int l=0; l<l_max; l++)
    {
      er      = static_cast<Real8>(rm_data_obj.get_ene_reso()[k][l]);
      wave_no = wave_no_part[i][j][k] * sqrt(fabs(er));
      calc_rho(i, j, k, er, wave_no, rho, rho_h);
      if( radius_calc_flg[i][j] == 1 && fabs(apl) > min_ene_dif )
      {
        rho = wave_no * apl;
      }
      calc_shift_penetration_factor(l_val, rho, sl_er[i][j][k][l], pl_er[i][j][k][l]);
      gx_er[i][j][k][l] = 0.0;
    } 
  }
  rm_data_obj.clear();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_const_value_adler(int i, int j)
{
  int k_max = nls[i][j];
  Real8 awri_val, arat;
  wave_no_part[i][j].resize(k_max);
  chan_rad[i][j].resize(k_max);
  for(int k=0; k<k_max; k++)
  {
    awri_val = awri[i][j][k];
    arat     = awri_val/(awri_val+1.0);
    wave_no_part[i][j][k] = k_part*arat;
    chan_rad[i][j][k]     = 6.5099897e+5 / (arat*arat);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::calc_const_value_unreso(int i, int j)
{
  int k_max = nls[i][j];
  Real8 awri_val;
  wave_no_part[i][j].resize(k_max);
  chan_rad[i][j].resize(k_max);
  for(int k=0; k<k_max; k++)
  {
    awri_val = awri[i][j][k];
    wave_no_part[i][j][k] = k_part*awri_val/(awri_val+1.0);
    chan_rad[i][j][k]     = 0.123 * pow(awri_val*amu_n,third)+0.08;
  }

  self_shielding_flg[i][j] = reso_data_obj.get_unreso_data_obj(i,j).get_self_shielding_flg();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::set_point_quadrature_weight_and_abscissa()
{
  clr_obj.clear_vec_array2_real8(q_wei);
  clr_obj.clear_vec_array2_real8(q_abs);

  q_wei.resize(xs_type_no);
  q_abs.resize(xs_type_no);
  for(int i=0; i<4; i++)
  {
    q_wei[i].resize(10);
    q_abs[i].resize(10);
  }
  
  q_wei[0][0] = 1.11204130E-01;
  q_wei[0][1] = 2.35467980E-01;
  q_wei[0][2] = 2.84409870E-01;
  q_wei[0][3] = 2.24191270E-01;
  q_wei[0][4] = 1.09676680E-01;
  q_wei[0][5] = 3.04937890E-02;
  q_wei[0][6] = 4.29308740E-03;
  q_wei[0][7] = 2.58270470E-04;
  q_wei[0][8] = 4.90319650E-06;
  q_wei[0][9] = 1.40792060E-08;
  
  q_abs[0][0] = 3.00134650E-03;
  q_abs[0][1] = 7.85928860E-02;
  q_abs[0][2] = 4.32824150E-01;
  q_abs[0][3] = 1.33452670E+00;
  q_abs[0][4] = 3.04818460E+00;
  q_abs[0][5] = 5.82631980E+00;
  q_abs[0][6] = 9.94526560E+00;
  q_abs[0][7] = 1.57821280E+01;
  q_abs[0][8] = 2.39968240E+01;
  q_abs[0][9] = 3.62162080E+01;


  q_wei[1][0] = 3.37734180E-02;
  q_wei[1][1] = 7.99321710E-02;
  q_wei[1][2] = 1.28359370E-01;
  q_wei[1][3] = 1.76526160E-01;
  q_wei[1][4] = 2.13470430E-01;
  q_wei[1][5] = 2.11549650E-01;
  q_wei[1][6] = 1.33651860E-01;
  q_wei[1][7] = 2.26306590E-02;
  q_wei[1][8] = 1.63136380E-05;
  q_wei[1][9] = 2.74538300E-31;

  q_abs[1][0] = 1.32192030E-02;
  q_abs[1][1] = 7.23496240E-02;
  q_abs[1][2] = 1.90894730E-01;
  q_abs[1][3] = 3.95288420E-01;
  q_abs[1][4] = 7.40834430E-01;
  q_abs[1][5] = 1.34982930E+00;
  q_abs[1][6] = 2.52979830E+00;
  q_abs[1][7] = 5.23848940E+00;
  q_abs[1][8] = 1.38217720E+01;
  q_abs[1][9] = 7.56475250E+01;


  q_wei[2][0] = 3.33762140E-04;
  q_wei[2][1] = 1.85061080E-02;
  q_wei[2][2] = 1.23099460E-01;
  q_wei[2][3] = 2.99189230E-01;
  q_wei[2][4] = 3.34314750E-01;
  q_wei[2][5] = 1.77666570E-01;
  q_wei[2][6] = 4.26958940E-02;
  q_wei[2][7] = 4.07605750E-03;
  q_wei[2][8] = 1.17661150E-04;
  q_wei[2][9] = 5.09895460E-07;

  q_abs[2][0] = 1.00044880E-03;
  q_abs[2][1] = 2.61976290E-02;
  q_abs[2][2] = 1.44274720E-01;
  q_abs[2][3] = 4.44842230E-01;
  q_abs[2][4] = 1.01606150E+00;
  q_abs[2][5] = 1.94210660E+00;
  q_abs[2][6] = 3.31508850E+00;
  q_abs[2][7] = 5.26070920E+00;
  q_abs[2][8] = 7.99894140E+00;
  q_abs[2][9] = 1.20720690E+01;


  q_wei[3][0] = 1.76237880E-03;
  q_wei[3][1] = 2.15177490E-02;
  q_wei[3][2] = 8.09798490E-02;
  q_wei[3][3] = 1.87979980E-01;
  q_wei[3][4] = 3.01563350E-01;
  q_wei[3][5] = 2.96160910E-01;
  q_wei[3][6] = 1.07756490E-01;
  q_wei[3][7] = 2.51719140E-03;
  q_wei[3][8] = 8.96303880E-10;
  q_wei[3][9] = 0.00000000E+00;

  q_abs[3][0] = 1.32192030E-02;
  q_abs[3][1] = 7.23496240E-02;
  q_abs[3][2] = 1.90894730E-01;
  q_abs[3][3] = 3.95288420E-01;
  q_abs[3][4] = 7.40834430E-01;
  q_abs[3][5] = 1.34982930E+00;
  q_abs[3][6] = 2.52979830E+00;
  q_abs[3][7] = 5.23848940E+00;
  q_abs[3][8] = 1.38217720E+01;
  q_abs[3][9] = 7.56475250E+01;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::check_react_type_list(int i, int j)
{
  int fis_total_flg = -1;
  int fis_first_flg = -1;

  int k_max = static_cast<int>(react_type_list[i][j].size());
  for(int k=0; k<k_max; k++)
  {
    if( react_type_list[i][j][k] == 18 )
    {
      fis_total_flg = k;
    }
    else if( react_type_list[i][j][k] == 19 )
    {
      fis_first_flg = k;
    }
  }

  if( fis_total_flg >= 0 || fis_first_flg >= 0 )
  {
    if( fis_total_flg >= 0 && fis_first_flg >= 0 )
    {
      //Do not add the fission data
    }
    else //if( fis_total_flg < 0 || fis_first_flg < 0 )
    {
      if( fis_total_flg >= 0 )
      {
        react_type_list[i][j].push_back(19);
        sig_array[i][j].push_back(sig_array[i][j][fis_total_flg]);
        q_array[i][j].push_back(q_array[i][j][fis_total_flg]);
      }
      else if( fis_first_flg >= 0 )
      {
        react_type_list[i][j].push_back(18);
        sig_array[i][j].push_back(sig_array[i][j][fis_first_flg]);
        q_array[i][j].push_back(q_array[i][j][fis_first_flg]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::set_resonance_grid(vector<vector<vector<Real8> > >& real_vec)
{
  calc_reso_xs_flg = 0;
  ene_array        = real_vec;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::set_temp(Real8 real_val)
{
  calc_reso_xs_flg = 0;
  temp             = real_val;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Definicao dos objetos set_ //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

// PosINAC_FRENDY - V3

void ResonanceXSCalculator::set_psi(Real8 real_val)    // define o objeto que seta o valor de psi_eff    
{                                                                               
  psi_eff           = real_val;                              
}

// PosINAC_FRENDY - V5

void ResonanceXSCalculator::set_w_x(Real8 real_val)    // define o objeto que seta o valor de w_x_eff    
{                                                                               
  w_x_eff           = real_val;                              
}

void ResonanceXSCalculator::set_w_y(Real8 real_val)    // define o objeto que seta o valor de w_y_eff    
{                                                                               
  w_y_eff           = real_val;                              
}

void ResonanceXSCalculator::set_w(Real8 real_val)    // define o objeto que seta o valor de w_eff    
{                                                                               
  w_eff           = real_val;                              
}

void ResonanceXSCalculator::set_x(Real8 real_val)    // define o objeto que seta o valor de x_eff    
{                                                                               
  x_eff           = real_val;                              
}

// PosINAC_FRENDY - V6

void ResonanceXSCalculator::set_l_max(Real8 real_val)    // define o objeto que seta o valor de l_max_eff    
{                                                                               
  l_max_eff           = real_val;                              
}

void ResonanceXSCalculator::set_m_max(Real8 real_val)    // define o objeto que seta o valor de m_max_eff    
{                                                                               
  m_max_eff           = real_val;                              
}

void ResonanceXSCalculator::set_m(Real8 real_val)    // define o objeto que seta o valor de m_max_eff    
{                                                                               
  m_eff           = real_val;                              
}

void ResonanceXSCalculator::set_E0(Real8 real_val)    // define o objeto que seta o valor de m_max_eff    
{                                                                               
  E0_eff           = real_val;                              
}

// PosINAC_FRENDY - V8

void ResonanceXSCalculator::set_qsi(Real8 real_val)    // define o objeto que seta o valor de qsi_eff    
{                                                                               
  qsi_eff           = real_val;                              
}

void ResonanceXSCalculator::set_fxqsi(Real8 real_val)    // define o objeto que seta o valor de fxqsi_eff    
{                                                                               
  fxqsi_eff           = real_val;                              
}

// PosINAC_FRENDY - V10

void ResonanceXSCalculator::set_gtt(Real8 real_val)    // define o objeto que seta o valor de qsi_eff    
{                                                                               
  gtt_eff           = real_val;                              
}

void ResonanceXSCalculator::set_smr_gg(Real8 real_val)    // define o objeto que seta o valor de smr*gg   
{                                                                               
  smr_gg_eff           = real_val;                              
}








//////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Fim da definicao de objetos set_ ///////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::set_err(Real8 real_val)
{
  calc_reso_xs_flg = 0;
  error_value      = real_val;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::set_err_max(Real8 real_val)
{
  calc_reso_xs_flg = 0;
  error_max_value  = real_val;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

void ResonanceXSCalculator::set_err_int(Real8 real_val)
{
  calc_reso_xs_flg = 0;
  error_int_value = real_val;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//Getter

vector<vector<vector<Integer> > > ResonanceXSCalculator::get_resonance_react_type_list()
{
  calc_resonance_xs();
  return react_type_list;
}

vector<vector<vector<Real8> > >   ResonanceXSCalculator::get_resonance_q_val()
{
  calc_resonance_xs();
  return q_array;
}

vector<vector<vector<Real8> > >          ResonanceXSCalculator::get_resonance_grid()
{
  calc_resonance_xs();
  return ene_array;
}

vector<vector<vector<vector<Real8> > > > ResonanceXSCalculator::get_resonance_xs()
{
  calc_resonance_xs();
  return sig_array;
}

Real8 ResonanceXSCalculator::get_temp()
{
  return temp;
}



//////////////////////////////////////////////////////////////////////////
/////////////////////// Definicao dos objtos get_ ////////////////////////
//////////////////////////////////////////////////////////////////////////

// PosINAC_FRENDY - V3

Real8 ResonanceXSCalculator::get_psi()       // definição do objeto que retorna psi_eff          
{                                                          
  return psi_eff;                                        
}

// PosINAC_FRENDY - V5

Real8 ResonanceXSCalculator::get_w_x()       // definição do objeto que retorna w_x_eff          
{                                                          
  return w_x_eff;                                        
}

Real8 ResonanceXSCalculator::get_w_y()       // definição do objeto que retorna w_y_eff          
{                                                          
  return w_y_eff;                                        
}

Real8 ResonanceXSCalculator::get_w()       // definição do objeto que retorna w_eff          
{                                                          
  return w_eff;                                        
}

// PosINAC_FRENDY - V5

Real8 ResonanceXSCalculator::get_x()       // definição do objeto que retorna x_eff          
{                                                          
  return x_eff;                                        
}

// PosINAC_FRENDY - V6

Real8 ResonanceXSCalculator::get_l_max()       // definição do objeto que retorna l_max_eff          
{                                                          
  return l_max_eff;                                        
}

Real8 ResonanceXSCalculator::get_m_max()       // definição do objeto que retorna m_max_eff          
{                                                          
  return m_max_eff;                                        
}

Real8 ResonanceXSCalculator::get_m()       // definição do objeto que retorna m_eff          
{                                                          
  return m_eff;                                        
}

Real8 ResonanceXSCalculator::get_E0()       // definição do objeto que retorna E0_eff          
{                                                          
  return E0_eff;                                        
}

// PosINAC_FRENDY - V8

Real8 ResonanceXSCalculator::get_qsi()       // definição do objeto que retorna qsi_eff          
{                                                          
  return qsi_eff;                                        
}

Real8 ResonanceXSCalculator::get_fxqsi()       // definição do objeto que retorna fxqsi_eff          
{                                                          
  return fxqsi_eff;                                        
}

// PosINAC_FRENDY - V10

Real8 ResonanceXSCalculator::get_gtt()       // definição do objeto que retorna gtt_eff          
{                                                          
  return gtt_eff;                                        
}


Real8 ResonanceXSCalculator::get_smr_gg()       // definição do objeto que retorna smr_gg_eff          
{                                                          
  return smr_gg_eff;                                        
}



//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////



Real8 ResonanceXSCalculator::get_err()
{
  return error_value;
}

Real8 ResonanceXSCalculator::get_err_max()
{
  return error_max_value;
}

Real8 ResonanceXSCalculator::get_err_int()
{
  return error_int_value;
}


NuclearDataObject ResonanceXSCalculator::get_nucl_data_obj()
{
  //Add error and temperature data in GeneralDataContainer
  GeneralDataContainer general_data_obj = nucl_data_obj.get_general_data_obj();
  general_data_obj.set_error_value(static_cast<Real>(error_value));
  general_data_obj.set_temp(static_cast<Real>(temp));
  general_data_obj.set_special_derived_mat_flg(1);
  nucl_data_obj.set_general_data_obj(general_data_obj);
  general_data_obj.clear();

  //Add temperature data in UnresolvedCrossSectionDataContainer
  ResonanceDataContainer              reso_data_obj      = nucl_data_obj.get_reso_data_obj();
  UnresolvedCrossSectionDataContainer unreso_xs_data_obj = reso_data_obj.get_unreso_xs_data_obj();
  unreso_xs_data_obj.set_temp(static_cast<Real>(temp));
  reso_data_obj.set_unreso_xs_data_obj(unreso_xs_data_obj);
  nucl_data_obj.set_reso_data_obj(reso_data_obj);
  unreso_xs_data_obj.clear();
  reso_data_obj.clear();

  return nucl_data_obj;
}

ResonanceDataContainer ResonanceXSCalculator::get_reso_data_obj()
{
  return reso_data_obj;
}

Real8 ResonanceXSCalculator::get_xs_potential()
{
  calc_resonance_xs();

  int i_max = static_cast<int>(ene_potential.size());
  if( i_max == 0 )
  {
    return 0.0;
  }
  else if( i_max == 1 )
  {
    return xs_potential[0];
  }

  vector<Real8> ene_vec, xs_vec;
  ene_vec = ene_potential;
  sort(ene_vec.begin(), ene_vec.end());

  xs_vec.resize(i_max);
  for(int i=0; i<i_max; i++)
  {
    Real8 ene_val = ene_vec[i];
    Real8 ene_chk = fabs(min_ene_dif * ene_val);
    for(int j=0; j<i_max; j++)
    {
      if( fabs(ene_potential[j] - ene_val) < ene_chk )
      {
        xs_vec[i] = xs_potential[j];
        break;
      }
    }
  }

  Real ene_integ, xs_integ, xs_mid;
  ene_integ = 0.0;
  xs_integ  = 0.0;
  if( ene_vec[0] < min_value )
  {
    ene_vec[0] = min_value;
  }
  for(int i=1; i<i_max; i++)
  {
    if( ene_vec[i] < min_value )
    {
      ene_vec[i] = min_value;
    }
    ene_integ +=  log(ene_vec[i]) - log(ene_vec[i-1]);

    xs_mid     = 0.5*(xs_vec[i-1] + xs_vec[i]);
    xs_integ  += (log(ene_vec[i]) - log(ene_vec[i-1])) * xs_mid;
  }

  Real xs_pot_data = 0.0;
  if( fabs(ene_integ) > min_ene_dif )
  {
    xs_pot_data = xs_integ / ene_integ;
  }

  return xs_pot_data;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////

