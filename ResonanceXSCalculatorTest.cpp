#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#define DEBUG_MODE

#include <boost/test/unit_test.hpp>

#include "EndfUtils/Endf6Converter/Endf6Converter.hpp"
#include "ReconResonance/ResonanceXSCalculator.hpp"
#include "ReconResonance/ResonanceEnergyGridLinearizer.hpp"

using namespace frendy;

////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(BreitWigner_check01)
{
  ErrorManager err_obj;
  err_obj.set_err_mes_opt(err_obj.err_mes_debug);
  //err_obj.set_err_mes_opt(err_obj.err_mes_default);

  string file_name = "../EndfUtils/MFxxMTyyyParser/for_test/n_4331_43-Tc-099.dat";
  
  Endf6ParserNoCov cp_obj;
  cp_obj.set_file_name(file_name);

  Endf6Converter    conv_obj;
  NuclearDataObject nucl_data_obj;
  conv_obj.convert_endf_format_to_frendy(cp_obj, nucl_data_obj);

  ResonanceEnergyGridLinearizer sre_obj;
  sre_obj.set_nucl_data_obj(nucl_data_obj);

  vector<vector<vector<Real8> > > sre_grid  = sre_obj.get_resonance_grid();

  ResonanceXSCalculator rxs_obj;
  rxs_obj.set_temp(100.0);
  rxs_obj.set_err(1.0e-2);
  rxs_obj.set_err_max(2.0e-1);
  rxs_obj.set_err_int(5.0e-3);
  
  vector<vector<vector<Real8> > >          reso_grid;
  vector<vector<vector<vector<Real8> > > > reso_xs;

  BOOST_CHECK_LE(fabs((rxs_obj.get_temp()    - 100.0 )/100.0 ), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err()     - 1.0e-2)/1.0e-2), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err_max() - 2.0e-1)/2.0e-1), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err_int() - 5.0e-3)/5.0e-3), 1.0e-5);

  rxs_obj.set_temp(0.0);
  rxs_obj.set_nucl_data_obj(nucl_data_obj);
  rxs_obj.set_resonance_grid(sre_grid);
  
//*
  ofstream fout;
  fout.open("./comp_njoy/calc_data_bw.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);

  vector<Real8>          ene_data;
  vector<vector<Real8> > sig_data;
  ene_data.resize(5);
  sig_data.resize(5);
  ene_data[0] = 1.0e+0;
  ene_data[1] = 5.0e+0;
  ene_data[2] = 1.0e+1;
  ene_data[3] = 5.0e+1;
  ene_data[4] = 1.0e+2;

  ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  MF02MT151Converter mf02mt151conv;
  MF02MT151Parser par_obj;
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  int i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      fout << static_cast<int>(ene_data.size()) << endl;
      for(int k=0; k<static_cast<int>(ene_data.size()); k++)
      {
        fout << ene_data[k] << " ";
      }
      fout << endl;

      fout << rxs_obj.get_temp() << endl;
      fout << "+0 " << par_obj.get_NAPS()[i][j] << " " << par_obj.get_NLS()[i][j] << endl;
      fout << par_obj.get_AWRI()[i][j][0] << " " << par_obj.get_AP()[i][j] << " " 
           << par_obj.get_SPI()[i][j] << endl;

      if( i==0 && j==0 )
      {
        for(int k=0; k<static_cast<int>(par_obj.get_NLS()[i][j]); k++)
        {
          fout << par_obj.get_NRS_LRU01_LRF02()[i][j][k] <<  " " 
               << par_obj.get_L_LRU01_LRF02()[i][j][k]   <<  " "
               << par_obj.get_QX_LRU01_LRF02()[i][j][k]  <<  " "  
               << par_obj.get_LRX_LRU01_LRF02()[i][j][k] << endl;

          for(int l=0; l<static_cast<int>(par_obj.get_NRS_LRU01_LRF02()[i][j][k]); l++)
          {
            fout << par_obj.get_ER_LRU01_LRF02()[i][j][k][l] << " "  
                 << par_obj.get_AJ_LRU01_LRF02()[i][j][k][l] << " "  << endl;
            fout << par_obj.get_GN_LRU01_LRF02()[i][j][k][l] << " "
                 << par_obj.get_GG_LRU01_LRF02()[i][j][k][l] << " "
                 << par_obj.get_GF_LRU01_LRF02()[i][j][k][l] << " "
                 << rxs_obj.gx_er[i][j][k][l] << " " << endl;
            fout << rxs_obj.sl_er[i][j][k][l] << " " << rxs_obj.pl_er[i][j][k][l] << endl;
          }
        }
        fout << endl;
        fout << "======================================================================" << endl;
        fout << endl;
      }
    }
  }
  fout.close();
  
//*
  rxs_obj.calc_resonance_xs();
  reso_grid = rxs_obj.get_resonance_grid();
  reso_xs   = rxs_obj.get_resonance_xs();

  fout.open("./comp_njoy/calc_result_TC99_frendy.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);
  i_max = static_cast<int>(reso_grid.size());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(reso_grid[i].size());
    for(int j=0; j<j_max; j++)
    {
      fout << "i :" << i+1 << "/" << i_max << ", j :" << j+1 << "/" << j_max << endl;
      int k_max = static_cast<int>(reso_grid[i][j].size());
      for(int k=0; k<k_max; k++)
      {
        fout << k << "\t" << reso_grid[i][j][k] << "\t" 
             << reso_xs[i][j][k][0] << "\t" << reso_xs[i][j][k][1] << "\t" 
             << reso_xs[i][j][k][2] << "\t" << reso_xs[i][j][k][3] << "\t" << endl;
      }
      fout << endl;
      fout << "====================================================================================================" << endl;
      fout << endl;
    }
  }
  fout.close();
// */ 
}

////////////////////////////////////////////////////////////////////////////////////////////////////

//*
BOOST_AUTO_TEST_CASE(BreitWigner_check02)
{
  string file_name = "../EndfUtils/MFxxMTyyyParser/for_test/n_9434_94-Pu-238.dat";
  
  Endf6ParserNoCov cp_obj;
  cp_obj.set_file_name(file_name);

  Endf6Converter    conv_obj;
  NuclearDataObject nucl_data_obj;
  conv_obj.convert_endf_format_to_frendy(cp_obj, nucl_data_obj);

  ResonanceEnergyGridLinearizer sre_obj;
  sre_obj.set_nucl_data_obj(nucl_data_obj);

  vector<vector<vector<Real8> > > sre_grid  = sre_obj.get_resonance_grid();

  ResonanceXSCalculator rxs_obj;
  rxs_obj.set_temp(100.0);
  rxs_obj.set_err(1.0e-2);
  rxs_obj.set_err_max(2.0e-1);
  rxs_obj.set_err_int(5.0e-3);
  
  vector<vector<vector<Real8> > >          reso_grid;
  vector<vector<vector<vector<Real8> > > > reso_xs;

  BOOST_CHECK_LE(fabs((rxs_obj.get_temp()    - 100.0 )/100.0 ), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err()     - 1.0e-2)/1.0e-2), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err_max() - 2.0e-1)/2.0e-1), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err_int() - 5.0e-3)/5.0e-3), 1.0e-5);

  rxs_obj.set_temp(0.0);
  rxs_obj.set_nucl_data_obj(nucl_data_obj);
  rxs_obj.set_resonance_grid(sre_grid);

  ofstream fout;
  fout.open("./comp_njoy/calc_data_bw_Pu238.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);

  vector<Real8>          ene_data;
  vector<vector<Real8> > sig_data;
  ene_data.resize(5);
  sig_data.resize(5);
  ene_data[0] = 1.0e+0;
  ene_data[1] = 5.0e+0;
  ene_data[2] = 1.0e+1;
  ene_data[3] = 5.0e+1;
  ene_data[4] = 1.0e+2;

  ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  MF02MT151Converter mf02mt151conv;
  MF02MT151Parser par_obj;
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  int i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      if( i==0 && j==0 )
      {
        fout << static_cast<int>(ene_data.size()) << endl;
        for(int k=0; k<static_cast<int>(ene_data.size()); k++)
        {
          fout << ene_data[k] << " ";
        }
        fout << endl;

        fout << rxs_obj.get_temp() << endl;
        fout << "+0 " << par_obj.get_NAPS()[i][j] << " " << par_obj.get_NLS()[i][j] << endl;
        fout << par_obj.get_AWRI()[i][j][0] << " " << par_obj.get_AP()[i][j] << " " 
             << par_obj.get_SPI()[i][j] << endl;
        for(int k=0; k<static_cast<int>(par_obj.get_NLS()[i][j]); k++)
        {
          fout << par_obj.get_NRS_LRU01_LRF02()[i][j][k] <<  " " 
               << par_obj.get_L_LRU01_LRF02()[i][j][k]   <<  " "
               << par_obj.get_QX_LRU01_LRF02()[i][j][k]  <<  " "  
               << par_obj.get_LRX_LRU01_LRF02()[i][j][k] << endl;

          for(int l=0; l<static_cast<int>(par_obj.get_NRS_LRU01_LRF02()[i][j][k]); l++)
          {
            fout << par_obj.get_ER_LRU01_LRF02()[i][j][k][l] << " "  
                 << par_obj.get_AJ_LRU01_LRF02()[i][j][k][l] << " "  << endl;
            fout << par_obj.get_GN_LRU01_LRF02()[i][j][k][l] << " "
                 << par_obj.get_GG_LRU01_LRF02()[i][j][k][l] << " "
                 << par_obj.get_GF_LRU01_LRF02()[i][j][k][l] << " "
                 << rxs_obj.gx_er[i][j][k][l] << " " << endl;
            fout << rxs_obj.sl_er[i][j][k][l] << " " << rxs_obj.pl_er[i][j][k][l] << endl;
          }
        }
        fout << endl;
        fout << "======================================================================" << endl;
        fout << endl;
      }
    }
  }
  fout.close();
}
// */ 

////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(BreitWigner_check03)
{
  vector<string> file_name;  // CHAMA DADOS NUCLEARES PARA SEREM UTILIZADOS E ARMAZENA NO VETOR file_name
  
  file_name.push_back("./for_test/test_mf02mt151_lru01lrf01.dat");                          // n=0
  file_name.push_back("../EndfUtils/MFxxMTyyyParser/for_test/n_4331_43-Tc-099.dat");        // n=1
  file_name.push_back("./for_test/test_mf02mt151_lru01lrf01_Pu238.dat");                    // n=2
  file_name.push_back("../EndfUtils/MFxxMTyyyParser/for_test/n_9434_94-Pu-238.dat");        // n=3
  file_name.push_back("./for_test/test_mf02mt151_lru01lrf01_U235.dat");                     // n=4  // PosInac_FRENDY-V6
  file_name.push_back("../EndfUtils/MFxxMTyyyParser/for_test/n_9228_92-U-235.dat");         // n=5  // PosInac_FRENDY-V6
  file_name.push_back("./for_test/test_mf02mt151_lru01lrf01_Pu239.dat");                    // n=6  // PosInac_FRENDY-V9
  file_name.push_back("./for_test/test_mf02mt151_lru01lrf01_Xe135.dat");                    // n=7  // PosInac_FRENDY-V9
  file_name.push_back("./for_test/test_mf02mt151_lru01lrf01_Fe056.dat");                    // n=8  // PosInac_FRENDY-V9
  file_name.push_back("./for_test/test_mf02mt151_lru01lrf01_Gd155.dat");                    // n=9  // PosInac_FRENDY-V9
  file_name.push_back("./for_test/test_mf02mt151_lru01lrf01_Gd157.dat");                    // n=10  // PosInac_FRENDY-V9
  
  
  //file_name.push_back("./for_test/test_mf02mt151_lru01lrf01_H001.dat");                     // n=11  // PosInac_FRENDY-V9
  //file_name.push_back("./for_test/test_mf02mt151_lru01lrf01_U238.dat");                     // n=7  // PosInac_FRENDY-V7
  //file_name.push_back("../EndfUtils/MFxxMTyyyParser/for_test/n_9237_92-U-238.dat.dat");     // n=8  // PosInac_FRENDY-V7
  
  //file_name.push_back("./for_test/test_mf02mt151_lru01lrf01_O016.dat");                     // n=11  // PosInac_FRENDY-V9

  
  
  vector<string> ref_file;  // APENAS PARA REFERÊNCIA
  ref_file.push_back("./comp_njoy/calc_result_bw.dat");
  ref_file.push_back("./comp_njoy/calc_result_bw.dat");
  ref_file.push_back("./comp_njoy/calc_result_bw_Pu238.dat");
  //ref_file.push_back("./comp_njoy/calc_result_bw_Pu238.dat");
  
  ofstream fout;
  ifstream fin;
  for(int n=0; n<static_cast<int>(file_name.size()); n++)
  {
    Endf6ParserNoCov cp_obj;
    cp_obj.set_file_name(file_name[n]);

    Endf6Converter    conv_obj;
    NuclearDataObject nucl_data_obj;
    conv_obj.convert_endf_format_to_frendy(cp_obj, nucl_data_obj);
 
    ResonanceXSCalculator rxs_obj;
    rxs_obj.set_nucl_data_obj(nucl_data_obj); 

    
    vector<Real8>          ene_data;
    vector<vector<Real8> > sig_data;



    int h_max = 1000;       // h_max e o numero de pontos              // PosInac_FRENDY-V1
    //int h_max = 300000;                                             //     PosInac_FRENDY-V7
    ene_data.resize(h_max); // define o numero de vetores de ene_data  // PosInac_FRENDY-V1
    sig_data.resize(h_max); // define o numero de vetores de sig_data  // PosInac_FRENDY-V1
    

    //ene_data[0] = 2.05;     // primeiro valor de ene_data // PosInac_FRENDY-V1
    //ene_data[h_max] = 3.71; // ultimo valor de ene_data   // PosInac_FRENDY-V1
    //ene_data[0] = 0.0005;
    //ene_data[h_max] = 400; 
    //ene_data[0] = 331.5902-0.8;            //       PosInac_FRENDY-V7 
    //ene_data[h_max] = 331.5902+0.8;        //     PosInac_FRENDY-V7
    //ene_data[0] = pow(10,0);            //       PosInac_FRENDY-V9
    //ene_data[h_max] = 5*pow(10,3);        //     PosInac_FRENDY-V9
    //ene_data[0] = 54.85246;
    //ene_data[h_max] = 61.66754;
    ene_data[0] = 9.0;
    ene_data[h_max] = 11.0; 

  
    
    for ( int h=1 ; h < h_max ; h++ ) // define os outros valores de ene_val // PosInac_FRENDY-V1
    {
      ene_data[h] = ene_data[h-1] + ((ene_data[h_max] - ene_data[0]) / h_max);
    }
    
    //Definição para range de pontos de energias máximo e mínimo muito grande //     PosInac_FRENDY-V10
    /*
    ene_data[0] = pow(10,-2);                                                //       PosInac_FRENDY-V10
    ene_data[1000] = pow(10,-1);                                             //     PosInac_FRENDY-V10
    ene_data[2000] = pow(10,0);                                              //     PosInac_FRENDY-V10
    ene_data[3000] = pow(10,1);                                              //     PosInac_FRENDY-V10
    ene_data[4000] = pow(10,2);                                              //     PosInac_FRENDY-V10    
    ene_data[5000] = pow(10,3);                                              //     PosInac_FRENDY-V10    
    ene_data[6000] = pow(10,4);                                              //     PosInac_FRENDY-V10    
    ene_data[7000] = pow(10,5);                                              //     PosInac_FRENDY-V10
    ene_data[8000] = pow(10,6);                                              //     PosInac_FRENDY-V10
    ene_data[9000] = pow(10,7);                                              //     PosInac_FRENDY-V10    
    ene_data[10000] = pow(10,8);                                             //     PosInac_FRENDY-V10    
    for ( int p=0 ; p<10 ; p++ ) // qual a casa decimal                      //     PosInac_FRENDY-V10
    {
      int k = h_max/10;
      for (int g= 1; g < k; g++) // qual a posição em cada casa decimal
      {
        int h = k*p+g;
        ene_data[h] = ene_data[h-1] + ((ene_data[(p+1)*k] - ene_data[p*k]) / k);
      }
    }
    */


  /* não serão necessários
  
    ene_data.resize(5);
    sig_data.resize(5);   
    ene_data[0] = 1.0e+0;
    ene_data[1] = 5.0e+0;
    ene_data[2] = 1.0e+1;
    ene_data[3] = 5.0e+1;
    ene_data[4] = 1.0e+2;
  */  
    
    //*   

    // OUTPUT

    if( n==0 )
    {
      fout.open("./comp_njoy/calc_result_bw_Tc99_frendy.dat");
    }
    else if( n==2 )
    {
      fout.open("./comp_njoy/calc_result_bw_Pu238_frendy.dat");
    }
    else if( n==4 )        // PosInac_FRENDY-V7
    {
      fout.open("./comp_njoy/calc_result_bw_u235_frendy.dat");
    }
    else if( n==6 )        // PosInac_FRENDY-V9
    {
      fout.open("./comp_njoy/calc_result_bw_Pu239_frendy.dat");
    }
    else if( n==7 )        // PosInac_FRENDY-V9
    {
      fout.open("./comp_njoy/calc_result_bw_Xe135_frendy.dat");
    }
    else if( n==8 )        // PosInac_FRENDY-V9
    {
      fout.open("./comp_njoy/calc_result_bw_Fe056_frendy.dat");
    }
    else if( n==9 )        // PosInac_FRENDY-V9
    {
      fout.open("./comp_njoy/calc_result_bw_Gd155_frendy.dat");
    }
    else if( n==10 )        // PosInac_FRENDY-V9
    {
      fout.open("./comp_njoy/calc_result_bw_Gd157_frendy.dat");
    }

       
    
    /*

    else if( n==11 )        // PosInac_FRENDY-V9
    {
      fout.open("./comp_njoy/calc_result_bw_O016_frendy.dat");
    }

    else if( n==6 )        // PosInac_FRENDY-V7
    {
      fout.open("./comp_njoy/calc_result_bw_u238_frendy.dat");
    }

    else if( n==11 )        // PosInac_FRENDY-V9
    {
      fout.open("./comp_njoy/calc_result_bw_H001_frendy.dat");
    } 
    */

    fout.precision(15);
    fout.setf(ios::scientific);
    fout.setf(ios::showpos);
    fout.setf(ios::showpoint);

    ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
    MF02MT151Converter mf02mt151conv;
    MF02MT151Parser par_obj;
    mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

    int i_max = static_cast<int>(par_obj.get_NIS());
    for(int i=0; i<i_max; i++)
    {
      int j_max = static_cast<int>(par_obj.get_NER()[i]);
      for(int j=0; j<j_max; j++)
      {
        if( i==0 && j==0 )
        {
          int k_max = static_cast<int>(ene_data.size());     // k_max e o tamanho do ene_data
          //
          if( n==2 || n==6 || n==9 || n==10)                   //PosINAC_FRENDY - V9
          //if( n==0 || n==2 || n==4 || n==6 || n==7 || n==8 || n==9 || n==10 || n==11 || n==12 )                   //PosINAC_FRENDY - V9
          //if( n==2 )        // Fará apenas os calculos para U238    //PosINAC_FRENDY - V8
          {
            fout << "Single-Level_Breite_Wigner" << endl;
            //for(int t=0; t<6; t++)                           // t define qnts temepraturas serao
            //for(int t=1; t<11; t++)                             // PosINAC_FRENDY - V7
            for(int t=3; t<6; t++)
            {
              /*
              if(n==40785623)           // definindo valores de T partindo de qsi //PosINAC_FRENDY - V10
              {                   // n==10 --> Gd157
                Real8 Gamma = 0.1712328222612969;
                Real8 Kb    = 8.6173324*(pow(10,-5));
                Real8 A     = 1.55576*(pow(10,2));
                Real8 E0    = 58.26; 

                Real8 Numerador = pow(Gamma,2)*A/(4*Kb*E0);

                temp_val = Numerador/pow((0.05*static_cast<Real8>(t)),2);
              }
              else */
              
              Real8 temp_val = static_cast<Real8>(t)*500.0;  // define o range entre as temperaturas
              //temp_val = static_cast<Real8>(t)*100.0;
              
              
              fout << "Temp do meio:       " << temp_val << "[K]"   << endl;
              fout << "Massa do isotopo:   " <<  reso_data_obj.get_mass_isotope()[i][j][0] << endl;      // PosINAC_FRENDY - V4
              //fout << "Valores de i e j    : " <<  i  << " , " <<  j  << endl;                         // PosINAC_FRENDY - V4
              //fout << "l_max: " << rxs_obj.get_l_max() << endl; // PosINAC_FRENDY - V6
              //fout << "m_max: " << rxs_obj.get_m_max() << endl; // PosINAC_FRENDY - V6


                          
              rxs_obj.set_temp(temp_val);
              for(int k=0; k<k_max; k++)
              {
                rxs_obj.calc_reso_xs_each_case(i, j, ene_data[k], sig_data[k]);  // sig_data entra como sig_val

                if (k==0)                                   // PosINAC_FRENDY - V3
                {
                  //fout << "Gamma: " << rxs_obj.get_gtt() << endl; // PosINAC_FRENDY - V6
                  //fout << "E0: " << rxs_obj.get_E0() << endl; // PosINAC_FRENDY - V6
                  //fout << "qsi: " << rxs_obj.get_qsi() << endl; // PosINAC_FRENDY - V10
                  fout << endl;
                  fout << "ene_val                ";         // PosINAC_FRENDY - V3
                  //fout << "   sig_val[total_xs]   ";         // PosINAC_FRENDY - V3
                  //fout << "  sig_val[scatter_xs]  ";         // PosINAC_FRENDY - V3
                  //fout << "  sig_val[fission_xs]  ";         // PosINAC_FRENDY - V3
                  fout << " sig_val[radiation_xs] ";         // PosINAC_FRENDY - V3
                  fout << "          psi          ";         // PosINAC_FRENDY - V3
                  //fout << "        w_x (im)       ";         // PosINAC_FRENDY - V5
                  //fout << "       w_y (real)      ";         // PosINAC_FRENDY - V5
                  //fout << "      w (real+im)      ";         // PosINAC_FRENDY - V5
                  
                  //fout << "           m           ";         // PosINAC_FRENDY - V6
                  //fout << "          E0           ";         // PosINAC_FRENDY - V6
                  //fout << "         qsi           ";         // PosINAC_FRENDY - V8
                  //fout << "        fxqsi          ";         // PosINAC_FRENDY - V8
                  //fout << "    (qsi-fxqsi)/qsi    ";         // PosINAC_FRENDY - V8
                  //fout << "           x           ";         // PosINAC_FRENDY - V6
                  
                  fout << "        smr*gg         ";
                  

                  

                  fout << endl << endl;
                }

                fout << ene_data[k] << " ";
                //fout << sig_data[k][0] << " ";           //plota os valores de sig_val[total_xs] // PosINAC_FRENDY - V7
                //fout << sig_data[k][1] << " ";           //plota os valores de sig_val[scatter_xs] // PosINAC_FRENDY - V7
                //fout << sig_data[k][2] << " ";           //plota os valores de sig_val[fission_xs] // PosINAC_FRENDY - V7
                fout << sig_data[k][3] << " ";           //plota os valores de sig_val[radiation_xs] // PosINAC_FRENDY - V7
                
                /*
                for(int l=0; l<4; l++)                   // l define qual dos xs ele vai pegar
                {
                  fout << sig_data[k][l] << " ";         //plota os valores de XS
                }
                */

                fout << rxs_obj.get_psi() << " ";        // PosINAC_FRENDY - V3
                //fout << rxs_obj.get_w_x() << " ";        // PosINAC_FRENDY - V5
                //fout << rxs_obj.get_w_y() << " ";        // PosINAC_FRENDY - V5
                //fout << rxs_obj.get_w()   << " ";        // PosINAC_FRENDY - V5
                
                //fout << rxs_obj.get_m()   << " ";        // PosINAC_FRENDY - V6
                //fout << rxs_obj.get_E0()  << " ";          // PosINAC_FRENDY - V7
                //fout << rxs_obj.get_qsi() << " ";          // PosINAC_FRENDY - V8
                //fout << rxs_obj.get_fxqsi() << " ";        // PosINAC_FRENDY - V8
                //fout << (rxs_obj.get_qsi() - rxs_obj.get_fxqsi())/rxs_obj.get_qsi() << " "; // PosINAC_FRENDY - V8
                //fout << rxs_obj.get_x()   << " ";          // PosINAC_FRENDY - V5
                fout << rxs_obj.get_smr_gg() << " ";
                
               

                if (k == (k_max-1))                      // PosINAC_FRENDY - V3
                {
                  fout << endl << endl;                  // PosINAC_FRENDY - V3
                  
                  fout << "______________________ ";     // PosINAC_FRENDY - V3
                  fout << "______________________ ";     // PosINAC_FRENDY - V3
                  fout << "______________________ ";     // PosINAC_FRENDY - V3
                  fout << "______________________ ";     // PosINAC_FRENDY - V3
                  //fout << "______________________ ";     // PosINAC_FRENDY - V3
                  //fout << "______________________ ";     // PosINAC_FRENDY - V3
                  //fout << "_______________________";     // PosINAC_FRENDY - V5
                  //fout << "_______________________";     // PosINAC_FRENDY - V5
                  //fout << "_______________________";     // PosINAC_FRENDY - V5
                  //fout << "_______________________";     // PosINAC_FRENDY - V6

                  fout << endl << endl;                  // PosINAC_FRENDY - V3
                }

                fout << endl;
              }
              fout << endl;
            }
          } /*  não trabalharemos com multi-level
          else if( n==1 || n==3 || n==5 || n==7 )
          {
            fout << "Multi-Level Breite Wigner" << endl;
            rxs_obj.set_temp(0.0);
            for(int k=0; k<k_max; k++)
            {
              rxs_obj.calc_reso_xs_each_case(i, j,ene_data[k], sig_data[k]);
              fout << ene_data[k] << " ";
              for(int l=0; l<4; l++)
              {
                fout << sig_data[k][l] << " ";
              }
              fout << endl;
            } 
            fout << endl; 
          } 
          */
        }
      }
    }


     fout.close(); // PosINAC_FRENDY - V9   o fout.close precisa ser efetuado para todos os n's

    /* // PosINAC_FRENDY - V9 

    if( n==1 || n==3 || n==5 || n==7 )
    {
      fout.close();
    }
    */ 
    
    Real8         tmp_real;
    vector<Real8> sig_ref;
    sig_ref.resize(4);
    
    fin.open(ref_file[n].c_str());

    reso_data_obj = rxs_obj.get_reso_data_obj();
    mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

    i_max = static_cast<int>(par_obj.get_NIS());
    for(int i=0; i<i_max; i++)
    {
      int j_max = static_cast<int>(par_obj.get_NER()[i]);
      for(int j=0; j<j_max; j++)
      {
        if( i==0 && j==0 )
        {
          int k_max = static_cast<int>(ene_data.size());
          if( n==0 || n==2 || n==4 || n==6)
          {
            string line_data;
            while(!fin.eof())
            {
              getline(fin, line_data);
              if( line_data.find("Single-Level Breite Wigner") != string::npos )
              {
                break;
              }
            }                        /* ****** tirando o checador *******  
            for(int t=0; t<5; t++)
            {
              getline(fin, line_data); //TEMP: X.X
              Real8 temp_val = static_cast<Real8>(t)*100.0;
              rxs_obj.set_temp(temp_val);
              for(int k=0; k<k_max; k++)
              {
                rxs_obj.calc_reso_xs_each_case(i, j, ene_data[k], sig_data[k]);
                fin >> tmp_real >> sig_ref[0] >> sig_ref[1] >> sig_ref[2] >> sig_ref[3];
                for(int l=0; l<4; l++)
                {
                  //**********  !!! caution !!!    ********** 
                  //BOOST_CHECK_LE(fabs(sig_data[k][l] - sig_ref[l]), fabs(1.0e-6*sig_ref[l]+1.0e-10));
                  BOOST_CHECK_LE(fabs(sig_data[k][l] - sig_ref[l]), fabs(5.0e-2*sig_ref[l]+1.0e-10));
                }
              }
              getline(fin, line_data);
              getline(fin, line_data);
            }                              */
          }
          else if( n==1 || n==3 || n==5 || n==7 )
          {
            string line_data;
            while(!fin.eof())
            {
              getline(fin, line_data);
              if( line_data.find("Multi-Level Breite Wigner") != string::npos )
              {
                break;
              }
            }
            rxs_obj.set_temp(0.0);                   /* ****** tirando o checador ******* 
            for(int k=0; k<k_max; k++)
            {
              rxs_obj.calc_reso_xs_each_case(i, j,ene_data[k], sig_data[k]);
              fin >> tmp_real >> sig_ref[0] >> sig_ref[1] >> sig_ref[2] >> sig_ref[3];
              
              for(int l=0; l<4; l++)
              {
                BOOST_CHECK_LE(fabs(sig_data[k][l] - sig_ref[l]), fabs(5.0e-5*sig_ref[l]+1.0e-10));
              }
            }                            */
          }
        }
      }
    }
    fin.close();
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(Coef_check01)
{
  vector<Real8> real_vec;
  real_vec.resize(33);
  real_vec[ 0] = 1.0e-5;
  real_vec[ 1] = 2.0e-5;
  real_vec[ 2] = 5.0e-5;
  real_vec[ 3] = 1.0e-4;
  real_vec[ 4] = 2.0e-4;
  real_vec[ 5] = 5.0e-4;
  real_vec[ 6] = 1.0e-3;
  real_vec[ 7] = 2.0e-3;
  real_vec[ 8] = 5.0e-3;
  real_vec[ 9] = 1.0e-2;
  real_vec[10] = 2.0e-2;
  real_vec[11] = 5.0e-2;
  real_vec[12] = 1.0e-1;
  real_vec[13] = 2.0e-1;
  real_vec[14] = 5.0e-1;
  real_vec[15] = 1.0e+0;
  real_vec[16] = 2.0e+0;
  real_vec[17] = 5.0e+0;
  real_vec[18] = 1.0e+1;
  real_vec[19] = 2.0e+1;
  real_vec[20] = 5.0e+1;
  real_vec[21] = 1.0e+2;
  real_vec[22] = 2.0e+2;
  real_vec[23] = 5.0e+2;
  real_vec[24] = 1.0e+3;
  real_vec[25] = 2.0e+3;
  real_vec[26] = 5.0e+3;
  real_vec[27] = 1.0e+4;
  real_vec[28] = 2.0e+4;
  real_vec[29] = 5.0e+4;
  real_vec[30] = 1.0e+5;
  real_vec[31] = 2.0e+5;
  real_vec[32] = 5.0e+5;

  ResonanceXSCalculator rxs_obj;
  
  vector<vector<Real8> > sl,     pl,     phi,     w_r,     w_i;
  vector<vector<Real8> > sl_ref, pl_ref, phi_ref, w_r_ref, w_i_ref;
  sl.resize(5);
  pl.resize(5);
  phi.resize(5); 
  sl_ref.resize(5);
  pl_ref.resize(5);
  phi_ref.resize(5); 
  int i_max = static_cast<int>(real_vec.size());
  int j_max = static_cast<int>(real_vec.size());
  for(Integer l=0; l<5; l++)
  {
    sl[l].resize(i_max);
    pl[l].resize(i_max);
    phi[l].resize(i_max);
    sl_ref[l].resize(i_max);
    pl_ref[l].resize(i_max);
    phi_ref[l].resize(i_max);
    for(int i=0; i<i_max; i++)
    {
      rxs_obj.calc_shift_penetration_factor(l, real_vec[i], sl[l][i], pl[l][i]);
      rxs_obj.calc_phase_shift(l, real_vec[i], phi[l][i]);
    }
  }
  
  MathUtils math_obj;
  w_r.resize(i_max);
  w_i.resize(i_max);
  w_r_ref.resize(i_max);
  w_i_ref.resize(i_max);
  for(int i=0; i<i_max; i++)
  {
    w_r[i].resize(j_max);
    w_i[i].resize(j_max);
    w_r_ref[i].resize(j_max);
    w_i_ref[i].resize(j_max);
    for(int j=0; j<j_max; j++)
    {
      math_obj.calc_cerfc(real_vec[i], real_vec[j], w_r[i][j], w_i[i][j]);
    }
  }
  
  //*
  ofstream fout;
  fout.open("./comp_njoy/coef_check_frendy.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);
  
  fout << "sl, pl, phi" << endl;
  for(int i=0; i<i_max; i++)
  {
    fout << "  " << real_vec[i];
    for(int l=0; l<5; l++)
    {
      fout << " " << sl[l][i];
    }
    for(int l=0; l<5; l++)
    {
      fout << " " << pl[l][i];
    }
    for(int l=0; l<5; l++)
    {
      fout << " " << phi[l][i];
    }
    fout << endl;
  }
  fout << endl;
  
  fout << "w_r" << endl;
  for(int i=0; i<i_max; i++)
  {
    fout << "  " << real_vec[i];
    for(int j=0; j<j_max; j++)
    {
      fout << "  " << w_r[i][j];
    }
    fout << endl;
  }
  fout << endl;
  
  fout << "w_i" << endl;
  for(int i=0; i<i_max; i++)
  {
    fout << "  " << real_vec[i];
    for(int j=0; j<j_max; j++)
    {
      fout << "  " << w_i[i][j];
    }
    fout << endl;
  }
  fout.close();
  // */ 
  
  ifstream fin;
  string   line_data;
  fin.open("./comp_njoy/coef_check.dat");
  getline(fin, line_data); // sl, pl, phi
  for(int i=0; i<i_max; i++)
  {
    fin >> real_vec[i];
    fin >> sl_ref[0][i]  >> sl_ref[1][i]  >> sl_ref[2][i]  >> sl_ref[3][i]  >> sl_ref[4][i];
    fin >> pl_ref[0][i]  >> pl_ref[1][i]  >> pl_ref[2][i]  >> pl_ref[3][i]  >> pl_ref[4][i];
    fin >> phi_ref[0][i] >> phi_ref[1][i] >> phi_ref[2][i] >> phi_ref[3][i] >> phi_ref[4][i];
  }
  getline(fin, line_data);
  
  getline(fin, line_data);
  getline(fin, line_data); //w_r
  for(int i=0; i<i_max; i++)
  {
    fin >> real_vec[i];
    for(int j=0; j<j_max; j++)
    {
      fin >> w_r_ref[i][j];
    }
  }
  getline(fin, line_data);
  
  getline(fin, line_data);
  getline(fin, line_data); //w_i
  for(int i=0; i<i_max; i++)
  {
    fin >> real_vec[i];
    for(int j=0; j<j_max; j++)
    {
      fin >> w_i_ref[i][j];
    }
  }
  getline(fin, line_data);
  fin.close();
  
  for(int i=0; i<i_max; i++)
  {
    for(int l=0; l<5; l++)
    {
      BOOST_CHECK_LE(fabs(sl[l][i]  - sl_ref[l][i]), fabs(5.0e-5*sl_ref[l][i] +1.0e-10));
      BOOST_CHECK_LE(fabs(pl[l][i]  - pl_ref[l][i]), fabs(5.0e-5*pl_ref[l][i] +1.0e-10));
      if( phi_ref[l][i] > min_value )
      {
        BOOST_CHECK_LE(fabs(phi[l][i] - phi_ref[l][i]),fabs(5.0e-5*phi_ref[l][i]+5.0e-10));
      }
      else
      {
        BOOST_CHECK_LE(fabs(phi[l][i] - phi_ref[l][i]),1.0e-6);
      }
    }
    
    for(int j=0; j<j_max; j++)
    {
      //**********  !!! caution !!!    ********** 
      //BOOST_CHECK_LE(fabs(w_r[i][j] - w_r_ref[i][j]),fabs(5.0e-1*w_r_ref[i][j]+1.0e-10));
      BOOST_CHECK_LE(fabs(w_i[i][j] - w_i_ref[i][j]),fabs(1.0e-2*w_i_ref[i][j]+1.0e-10));
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(Reich_Moore_check01)
{
  string file_name = "../EndfUtils/MFxxMTyyyParser/for_test/n_9228_92-U-235.dat";
  
  Endf6ParserNoCov cp_obj;
  cp_obj.set_file_name(file_name);

  Endf6Converter    conv_obj;
  NuclearDataObject nucl_data_obj;
  conv_obj.convert_endf_format_to_frendy(cp_obj, nucl_data_obj);

  ResonanceEnergyGridLinearizer sre_obj;
  sre_obj.set_nucl_data_obj(nucl_data_obj);

  vector<vector<vector<Real8> > > sre_grid  = sre_obj.get_resonance_grid();

  ResonanceXSCalculator rxs_obj;
  rxs_obj.set_temp(100.0);
  rxs_obj.set_err(1.0e-2);
  rxs_obj.set_err_max(2.0e-1);
  rxs_obj.set_err_int(5.0e-3);
  
  vector<vector<vector<Real8> > >          reso_grid;
  vector<vector<vector<vector<Real8> > > > reso_xs;

  BOOST_CHECK_LE(fabs((rxs_obj.get_temp()    - 100.0 )/100.0 ), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err()     - 1.0e-2)/1.0e-2), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err_max() - 2.0e-1)/2.0e-1), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err_int() - 5.0e-3)/5.0e-3), 1.0e-5);

  rxs_obj.set_temp(0.0);
  rxs_obj.set_nucl_data_obj(nucl_data_obj);
  rxs_obj.set_resonance_grid(sre_grid);

  vector<Real8>          ene_data;
  vector<vector<Real8> > sig_data;
  ene_data.resize(5);
  sig_data.resize(5);
  ene_data[0] = 1.0e+0;
  ene_data[1] = 5.0e+0;
  ene_data[2] = 1.0e+1;
  ene_data[3] = 5.0e+1;
  ene_data[4] = 1.0e+2;


  //*
  ofstream fout;
  fout.open("./comp_njoy/calc_data_rm.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);

  ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  MF02MT151Converter mf02mt151conv;
  MF02MT151Parser par_obj;
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  int i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      fout << static_cast<int>(ene_data.size()) << endl;
      for(int k=0; k<static_cast<int>(ene_data.size()); k++)
      {
        fout << ene_data[k] << " ";
      }
      fout << endl;

      fout << par_obj.get_NAPS()[i][j]    << " " << par_obj.get_NLS()[i][j] << endl;
      fout << par_obj.get_AWRI()[i][j][0] << " " << par_obj.get_AP()[i][j] << " " 
           << par_obj.get_SPI()[i][j] << endl;

      if( i==0 && j==0 )
      {
        for(int k=0; k<static_cast<int>(par_obj.get_NLS()[i][j]); k++)
        {
          fout << par_obj.get_NRS_LRU01_LRF03()[i][j][k]  <<  " " 
               << par_obj.get_NLSC_LRU01_LRF03()[i][j]    <<  " " 
               << par_obj.get_L_LRU01_LRF03()[i][j][k]    <<  " "
               << par_obj.get_APL_LRU01_LRF03()[i][j][k]  << endl;

          for(int l=0; l<static_cast<int>(par_obj.get_NRS_LRU01_LRF03()[i][j][k]); l++)
          {
            fout << par_obj.get_ER_LRU01_LRF03()[i][j][k][l]  << " "  
                 << par_obj.get_AJ_LRU01_LRF03()[i][j][k][l]  << " "  << endl;
            fout << par_obj.get_GN_LRU01_LRF03()[i][j][k][l]  << " "
                 << par_obj.get_GG_LRU01_LRF03()[i][j][k][l]  << " "
                 << par_obj.get_GFA_LRU01_LRF03()[i][j][k][l] << " "
                 << par_obj.get_GFB_LRU01_LRF03()[i][j][k][l] << " "
                 << rxs_obj.gx_er[i][j][k][l] << " " << endl;
            fout << rxs_obj.sl_er[i][j][k][l] << " " << rxs_obj.pl_er[i][j][k][l] << endl;
          }
        }
        fout << endl;
        fout << "======================================================================" << endl;
        fout << endl;
      }
    }
  }
  fout.close();
  
  //*
  rxs_obj.calc_resonance_xs();
  reso_grid = rxs_obj.get_resonance_grid();
  reso_xs   = rxs_obj.get_resonance_xs();

  fout.open("./comp_njoy/calc_result_U235_frendy.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);
  i_max = static_cast<int>(reso_grid.size());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(reso_grid[i].size());
    for(int j=0; j<j_max; j++)
    {
      fout << "i :" << i+1 << "/" << i_max << ", j :" << j+1 << "/" << j_max << endl;
      int k_max = static_cast<int>(reso_grid[i][j].size());
      for(int k=0; k<k_max; k++)
      {
        fout << k << "\t" << reso_grid[i][j][k] << "\t"
             << reso_xs[i][j][k][0] << "\t" << reso_xs[i][j][k][1] << "\t"
             << reso_xs[i][j][k][2] << "\t" << reso_xs[i][j][k][3] << "\t" << endl;
      }
      fout << endl;
      fout << "====================================================================================================" << endl;
      fout << endl;
    }
  }
  fout.close();
//  */

  //*
  fout.open("./comp_njoy/calc_result_rm_frendy.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);
  //ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  //MF02MT151Converter mf02mt151conv;
  //MF02MT151Parser par_obj;
  //mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);
  i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      if( i==0 && j==0 )
      {
        fout << "Reich-Moore" << endl;
        rxs_obj.set_temp(0.0);
        int k_max = static_cast<int>(ene_data.size());
        for(int k=0; k<k_max; k++)
        {
          rxs_obj.calc_reso_xs_each_case(i, j,ene_data[k], sig_data[k]);
          fout << ene_data[k] << " ";
          for(int l=0; l<4; l++)
          {
            fout << sig_data[k][l] << " ";
          }
          fout << endl;
        }
      }
    }
  }
  fout << endl;
  fout.close();
  // */ 
  
  ifstream fin;
  fin.open("./comp_njoy/calc_result_rm.dat");
  
  Real8         tmp_real;
  vector<Real8> sig_ref;
  sig_ref.resize(4);
  
  //ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  //MF02MT151Converter mf02mt151conv;
  //MF02MT151Parser par_obj;
  //mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);
  //int i_max = static_cast<int>(par_obj.get_NIS());

  reso_data_obj = rxs_obj.get_reso_data_obj();
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      if( i==0 && j==0 )
      {
        string   line_data;
        getline(fin, line_data); // Reich-Moore
        rxs_obj.set_temp(0.0);
        int k_max = static_cast<int>(ene_data.size());
        for(int k=0; k<k_max; k++)
        {
          rxs_obj.calc_reso_xs_each_case(i, j,ene_data[k], sig_data[k]);
          fin >> tmp_real >> sig_ref[0] >> sig_ref[1] >> sig_ref[2] >> sig_ref[3];
          for(int l=0; l<4; l++)
          {
            BOOST_CHECK_LE(fabs(sig_data[k][l] - sig_ref[l]),fabs(5.0e-5*sig_ref[l]+1.0e-10));
          }
        }
        getline(fin, line_data);
      }
    }
  }
  fin.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(Reich_Moore_check03)
{
  string file_name = "../EndfUtils/MFxxMTyyyParser/for_test/n_2631_26-Fe-056.dat";
  
  Endf6ParserNoCov cp_obj;
  cp_obj.set_file_name(file_name);

  Endf6Converter    conv_obj;
  NuclearDataObject nucl_data_obj;
  conv_obj.convert_endf_format_to_frendy(cp_obj, nucl_data_obj);

  ResonanceEnergyGridLinearizer sre_obj;
  sre_obj.set_nucl_data_obj(nucl_data_obj);

  vector<vector<vector<Real8> > > sre_grid  = sre_obj.get_resonance_grid();

  ResonanceXSCalculator rxs_obj;
  rxs_obj.set_temp(100.0);
  rxs_obj.set_err(1.0e-2);
  rxs_obj.set_err_max(2.0e-1);
  rxs_obj.set_err_int(5.0e-3);
  
  vector<vector<vector<Real8> > >          reso_grid;
  vector<vector<vector<vector<Real8> > > > reso_xs;

  BOOST_CHECK_LE(fabs((rxs_obj.get_temp()    - 100.0 )/100.0 ), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err()     - 1.0e-2)/1.0e-2), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err_max() - 2.0e-1)/2.0e-1), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err_int() - 5.0e-3)/5.0e-3), 1.0e-5);

  rxs_obj.set_temp(0.0);
  rxs_obj.set_nucl_data_obj(nucl_data_obj);
  rxs_obj.set_resonance_grid(sre_grid);

  vector<Real8>          ene_data;
  vector<vector<Real8> > sig_data;
  ene_data.resize(3);
  sig_data.resize(3);
  ene_data[0] = 1150.90;
  ene_data[1] = 1151.00;
  ene_data[2] = 1151.10;


  //*
  ofstream fout;
  fout.open("./comp_njoy/calc_data_rm_Fe56.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);

  ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  MF02MT151Converter mf02mt151conv;
  MF02MT151Parser par_obj;
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  int i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      fout << static_cast<int>(ene_data.size()) << endl;
      for(int k=0; k<static_cast<int>(ene_data.size()); k++)
      {
        fout << ene_data[k] << " ";
      }
      fout << endl;

      fout << par_obj.get_NAPS()[i][j]    << " " << par_obj.get_NLS()[i][j] << endl;
      fout << par_obj.get_AWRI()[i][j][0] << " " << par_obj.get_AP()[i][j] << " " 
           << par_obj.get_SPI()[i][j] << endl;

      if( i==0 && j==0 )
      {
        for(int k=0; k<static_cast<int>(par_obj.get_NLS()[i][j]); k++)
        {
          fout << par_obj.get_NRS_LRU01_LRF03()[i][j][k]  <<  " " 
               << par_obj.get_NLSC_LRU01_LRF03()[i][j]    <<  " " 
               << par_obj.get_L_LRU01_LRF03()[i][j][k]    <<  " "
               << par_obj.get_APL_LRU01_LRF03()[i][j][k]  << endl;

          for(int l=0; l<static_cast<int>(par_obj.get_NRS_LRU01_LRF03()[i][j][k]); l++)
          {
            fout << par_obj.get_ER_LRU01_LRF03()[i][j][k][l]  << " "  
                 << par_obj.get_AJ_LRU01_LRF03()[i][j][k][l]  << " "  << endl;
            fout << par_obj.get_GN_LRU01_LRF03()[i][j][k][l]  << " "
                 << par_obj.get_GG_LRU01_LRF03()[i][j][k][l]  << " "
                 << par_obj.get_GFA_LRU01_LRF03()[i][j][k][l] << " "
                 << par_obj.get_GFB_LRU01_LRF03()[i][j][k][l] << " "
                 << rxs_obj.gx_er[i][j][k][l] << " " << endl;
            fout << rxs_obj.sl_er[i][j][k][l] << " " << rxs_obj.pl_er[i][j][k][l] << endl;
          }
        }
        fout << endl;
        fout << "======================================================================" << endl;
        fout << endl;
      }
    }
  }
  fout.close();
  
  //*
  rxs_obj.calc_resonance_xs();
  reso_grid = rxs_obj.get_resonance_grid();
  reso_xs   = rxs_obj.get_resonance_xs();

  fout.open("./comp_njoy/calc_result_Fe56_frendy.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);
  i_max = static_cast<int>(reso_grid.size());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(reso_grid[i].size());
    for(int j=0; j<j_max; j++)
    {
      fout << "i :" << i+1 << "/" << i_max << ", j :" << j+1 << "/" << j_max << endl;
      int k_max = static_cast<int>(reso_grid[i][j].size());
      for(int k=0; k<k_max; k++)
      {
        fout << k << "\t" << reso_grid[i][j][k] << "\t"
             << reso_xs[i][j][k][0] << "\t" << reso_xs[i][j][k][1] << "\t"
             << reso_xs[i][j][k][2] << "\t" << reso_xs[i][j][k][3] << "\t" << endl;
      }
      fout << endl;
      fout << "====================================================================================================" << endl;
      fout << endl;
    }
  }
  fout.close();
//  */

  ifstream fin;
  fin.open("./comp_njoy/calc_result_rm_Fe56.dat");
  
  Real8         tmp_real;
  vector<Real8> sig_ref;
  sig_ref.resize(4);
  
  //ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  //MF02MT151Converter mf02mt151conv;
  //MF02MT151Parser par_obj;
  //mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);
  //int i_max = static_cast<int>(par_obj.get_NIS());

  reso_data_obj = rxs_obj.get_reso_data_obj();
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      if( i==0 && j==0 )
      {
        string   line_data;
        getline(fin, line_data); // Reich-Moore
        rxs_obj.set_temp(0.0);
        int k_max = static_cast<int>(ene_data.size());
        for(int k=0; k<k_max; k++)
        {
          rxs_obj.calc_reso_xs_each_case(i, j,ene_data[k], sig_data[k]);
          fin >> tmp_real >> sig_ref[0] >> sig_ref[1] >> sig_ref[2] >> sig_ref[3];
          for(int l=0; l<4; l++)
          {
            BOOST_CHECK_LE(fabs(sig_data[k][l] - sig_ref[l]),fabs(5.0e-5*sig_ref[l]+1.0e-10));
          }
        }
        getline(fin, line_data);
      }
    }
  }
  fin.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(Reich_Moore_check04)
{
  string file_name = "../EndfUtils/MFxxMTyyyParser/for_test/n_9437_94-Pu-239.dat";
  
  Endf6ParserNoCov cp_obj;
  cp_obj.set_file_name(file_name);

  Endf6Converter    conv_obj;
  NuclearDataObject nucl_data_obj;
  conv_obj.convert_endf_format_to_frendy(cp_obj, nucl_data_obj);

  ResonanceEnergyGridLinearizer sre_obj;
  sre_obj.set_nucl_data_obj(nucl_data_obj);

  vector<vector<vector<Real8> > > sre_grid  = sre_obj.get_resonance_grid();

  ResonanceXSCalculator rxs_obj;
  rxs_obj.set_err(1.0e-2);
  rxs_obj.set_err_max(2.0e-1);
  rxs_obj.set_err_int(5.0e-3);
  
  vector<vector<vector<Real8> > >          reso_grid;
  vector<vector<vector<vector<Real8> > > > reso_xs;

  rxs_obj.set_temp(0.0);
  rxs_obj.set_nucl_data_obj(nucl_data_obj);
  rxs_obj.set_resonance_grid(sre_grid);

  vector<Real8>          ene_data;
  vector<vector<Real8> > sig_data;
  ene_data.resize(3);
  sig_data.resize(3);
  ene_data[0] = 1.2e-5;
  ene_data[1] = 1.3e-5;
  ene_data[2] = 1.4e-4;


  //*
  ofstream fout;
  fout.open("./comp_njoy/calc_data_rm_Pu239.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);

  ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  MF02MT151Converter mf02mt151conv;
  MF02MT151Parser par_obj;
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  int i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      fout << static_cast<int>(ene_data.size()) << endl;
      for(int k=0; k<static_cast<int>(ene_data.size()); k++)
      {
        fout << ene_data[k] << " ";
      }
      fout << endl;

      fout << par_obj.get_NAPS()[i][j]    << " " << par_obj.get_NLS()[i][j] << endl;
      fout << par_obj.get_AWRI()[i][j][0] << " " << par_obj.get_AP()[i][j] << " " 
           << par_obj.get_SPI()[i][j] << endl;

      if( i==0 && j==0 )
      {
        for(int k=0; k<static_cast<int>(par_obj.get_NLS()[i][j]); k++)
        {
          fout << par_obj.get_NRS_LRU01_LRF03()[i][j][k]  <<  " " 
               << par_obj.get_NLSC_LRU01_LRF03()[i][j]    <<  " " 
               << par_obj.get_L_LRU01_LRF03()[i][j][k]    <<  " "
               << par_obj.get_APL_LRU01_LRF03()[i][j][k]  << endl;

          for(int l=0; l<static_cast<int>(par_obj.get_NRS_LRU01_LRF03()[i][j][k]); l++)
          {
            fout << par_obj.get_ER_LRU01_LRF03()[i][j][k][l]  << " "  
                 << par_obj.get_AJ_LRU01_LRF03()[i][j][k][l]  << " "  << endl;
            fout << par_obj.get_GN_LRU01_LRF03()[i][j][k][l]  << " "
                 << par_obj.get_GG_LRU01_LRF03()[i][j][k][l]  << " "
                 << par_obj.get_GFA_LRU01_LRF03()[i][j][k][l] << " "
                 << par_obj.get_GFB_LRU01_LRF03()[i][j][k][l] << " "
                 << rxs_obj.gx_er[i][j][k][l] << " " << endl;
            fout << rxs_obj.sl_er[i][j][k][l] << " " << rxs_obj.pl_er[i][j][k][l] << endl;
          }
        }
        fout << endl;
        fout << "======================================================================" << endl;
        fout << endl;
      }
    }
  }
  fout.close();
  
  //*
  rxs_obj.calc_resonance_xs();
  reso_grid = rxs_obj.get_resonance_grid();
  reso_xs   = rxs_obj.get_resonance_xs();

  fout.open("./comp_njoy/calc_result_Pu239_frendy.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);
  i_max = static_cast<int>(reso_grid.size());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(reso_grid[i].size());
    for(int j=0; j<j_max; j++)
    {
      fout << "i :" << i+1 << "/" << i_max << ", j :" << j+1 << "/" << j_max << endl;
      int k_max = static_cast<int>(reso_grid[i][j].size());
      for(int k=0; k<k_max; k++)
      {
        fout << k << "\t" << reso_grid[i][j][k] << "\t"
             << reso_xs[i][j][k][0] << "\t" << reso_xs[i][j][k][1] << "\t"
             << reso_xs[i][j][k][2] << "\t" << reso_xs[i][j][k][3] << "\t" << endl;
      }
      fout << endl;
      fout << "====================================================================================================" << endl;
      fout << endl;
    }
  }
  fout.close();
//  */

  ifstream fin;
  fin.open("./comp_njoy/calc_result_rm_Pu239.dat");
  
  Real8         tmp_real;
  vector<Real8> sig_ref;
  sig_ref.resize(4);
  
  //ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  //MF02MT151Converter mf02mt151conv;
  //MF02MT151Parser par_obj;
  //mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);
  //int i_max = static_cast<int>(par_obj.get_NIS());

  reso_data_obj = rxs_obj.get_reso_data_obj();
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      if( i==0 && j==0 )
      {
        string   line_data;
        getline(fin, line_data); // Reich-Moore
        rxs_obj.set_temp(0.0);
        int k_max = static_cast<int>(ene_data.size());
        for(int k=0; k<k_max; k++)
        {
          rxs_obj.calc_reso_xs_each_case(i, j,ene_data[k], sig_data[k]);
          fin >> tmp_real >> sig_ref[0] >> sig_ref[1] >> sig_ref[2] >> sig_ref[3];
          for(int l=0; l<4; l++)
          {
            BOOST_CHECK_LE(fabs(sig_data[k][l] - sig_ref[l]),fabs(5.0e-5*sig_ref[l]+1.0e-10));
          }
        }
        getline(fin, line_data);
      }
    }
  }
  fin.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(Reich_Moore_check05)
{
  string file_name = "../EndfUtils/MFxxMTyyyParser/for_test/n_1837_18-Ar-040.dat";
  
  Endf6ParserNoCov cp_obj;
  cp_obj.set_file_name(file_name);

  Endf6Converter    conv_obj;
  NuclearDataObject nucl_data_obj;
  conv_obj.convert_endf_format_to_frendy(cp_obj, nucl_data_obj);

  ResonanceEnergyGridLinearizer sre_obj;
  sre_obj.set_nucl_data_obj(nucl_data_obj);

  vector<vector<vector<Real8> > > sre_grid  = sre_obj.get_resonance_grid();

  ResonanceXSCalculator rxs_obj;
  rxs_obj.set_err(1.0e-2);
  rxs_obj.set_err_max(2.0e-1);
  rxs_obj.set_err_int(5.0e-3);
  
  vector<vector<vector<Real8> > >          reso_grid;
  vector<vector<vector<vector<Real8> > > > reso_xs;

  rxs_obj.set_temp(0.0);
  rxs_obj.set_nucl_data_obj(nucl_data_obj);
  rxs_obj.set_resonance_grid(sre_grid);

  vector<Real8>          ene_data;
  vector<vector<Real8> > sig_data;
  ene_data.resize(4);
  sig_data.resize(4);
  ene_data[0] = 1.121538139648438e+06;
  ene_data[1] = 1.121539194824219e+06;
  ene_data[2] = 1.121540250000000e+06;
  ene_data[3] = 1.121541484375000e+06;

  //*
  ofstream fout;
  fout.open("./comp_njoy/calc_data_rm_Ar040.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);

  ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  MF02MT151Converter mf02mt151conv;
  MF02MT151Parser par_obj;
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  int i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      fout << static_cast<int>(ene_data.size()) << endl;
      for(int k=0; k<static_cast<int>(ene_data.size()); k++)
      {
        fout << ene_data[k] << " ";
      }
      fout << endl;

      fout << par_obj.get_NAPS()[i][j]    << " " << par_obj.get_NLS()[i][j] << endl;
      fout << par_obj.get_AWRI()[i][j][0] << " " << par_obj.get_AP()[i][j] << " " 
           << par_obj.get_SPI()[i][j] << endl;

      if( i==0 && j==0 )
      {
        for(int k=0; k<static_cast<int>(par_obj.get_NLS()[i][j]); k++)
        {
          fout << par_obj.get_NRS_LRU01_LRF03()[i][j][k]  <<  " " 
               << par_obj.get_NLSC_LRU01_LRF03()[i][j]    <<  " " 
               << par_obj.get_L_LRU01_LRF03()[i][j][k]    <<  " "
               << par_obj.get_APL_LRU01_LRF03()[i][j][k]  << endl;

          for(int l=0; l<static_cast<int>(par_obj.get_NRS_LRU01_LRF03()[i][j][k]); l++)
          {
            fout << par_obj.get_ER_LRU01_LRF03()[i][j][k][l]  << " "  
                 << par_obj.get_AJ_LRU01_LRF03()[i][j][k][l]  << " "  << endl;
            fout << par_obj.get_GN_LRU01_LRF03()[i][j][k][l]  << " "
                 << par_obj.get_GG_LRU01_LRF03()[i][j][k][l]  << " "
                 << par_obj.get_GFA_LRU01_LRF03()[i][j][k][l] << " "
                 << par_obj.get_GFB_LRU01_LRF03()[i][j][k][l] << " "
                 << rxs_obj.gx_er[i][j][k][l] << " " << endl;
            fout << rxs_obj.sl_er[i][j][k][l] << " " << rxs_obj.pl_er[i][j][k][l] << endl;
          }
        }
        fout << endl;
        fout << "======================================================================" << endl;
        fout << endl;
      }
    }
  }
  fout.close();
  
  //*
  rxs_obj.calc_resonance_xs();
  reso_grid = rxs_obj.get_resonance_grid();
  reso_xs   = rxs_obj.get_resonance_xs();

  fout.open("./comp_njoy/calc_result_Ar040_frendy.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);
  i_max = static_cast<int>(reso_grid.size());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(reso_grid[i].size());
    for(int j=0; j<j_max; j++)
    {
      fout << "i :" << i+1 << "/" << i_max << ", j :" << j+1 << "/" << j_max << endl;
      int k_max = static_cast<int>(reso_grid[i][j].size());
      for(int k=0; k<k_max; k++)
      {
        fout << k << "\t" << reso_grid[i][j][k] << "\t"
             << reso_xs[i][j][k][0] << "\t" << reso_xs[i][j][k][1] << "\t"
             << reso_xs[i][j][k][2] << "\t" << reso_xs[i][j][k][3] << "\t" << endl;
      }
      fout << endl;
      fout << "====================================================================================================" << endl;
      fout << endl;
    }
  }
  fout.close();
//  */

  ifstream fin;
  fin.open("./comp_njoy/calc_result_rm_Ar040.dat");
  
  Real8         tmp_real;
  vector<Real8> sig_ref;
  sig_ref.resize(4);
  
  //ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  //MF02MT151Converter mf02mt151conv;
  //MF02MT151Parser par_obj;
  //mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);
  //int i_max = static_cast<int>(par_obj.get_NIS());

  reso_data_obj = rxs_obj.get_reso_data_obj();
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      if( i==0 && j==0 )
      {
        string   line_data;
        getline(fin, line_data); // Reich-Moore
        rxs_obj.set_temp(0.0);
        int k_max = static_cast<int>(ene_data.size());
        for(int k=0; k<k_max; k++)
        {
          rxs_obj.calc_reso_xs_each_case(i, j,ene_data[k], sig_data[k]);
          fin >> tmp_real >> sig_ref[0] >> sig_ref[1] >> sig_ref[2] >> sig_ref[3];
          for(int l=0; l<4; l++)
          {
            BOOST_CHECK_LE(fabs(sig_data[k][l] - sig_ref[l]),fabs(5.0e-5*sig_ref[l]+1.0e-10));
          }
        }
        getline(fin, line_data);
      }
    }
  }
  fin.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

//*
BOOST_AUTO_TEST_CASE(AdlerAdler_check01)
{
  string file_name = "./for_test/test_mf02mt151_lru01lrf04.dat";
  
  Endf6ParserNoCov cp_obj;
  cp_obj.set_file_name(file_name);

  Endf6Converter    conv_obj;
  NuclearDataObject nucl_data_obj;
  conv_obj.convert_endf_format_to_frendy(cp_obj, nucl_data_obj);

  ResonanceEnergyGridLinearizer sre_obj;
  sre_obj.set_nucl_data_obj(nucl_data_obj);

  vector<vector<vector<Real8> > > sre_grid  = sre_obj.get_resonance_grid();

  ResonanceXSCalculator rxs_obj;
  rxs_obj.set_temp(100.0);
  rxs_obj.set_err(1.0e-2);
  rxs_obj.set_err_max(2.0e-1);
  rxs_obj.set_err_int(5.0e-3);
  
  vector<vector<vector<Real8> > >          reso_grid;
  vector<vector<vector<vector<Real8> > > > reso_xs;

  BOOST_CHECK_LE(fabs((rxs_obj.get_temp()    - 100.0 )/100.0 ), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err()     - 1.0e-2)/1.0e-2), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err_max() - 2.0e-1)/2.0e-1), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err_int() - 5.0e-3)/5.0e-3), 1.0e-5);

  rxs_obj.set_temp(0.0);
  rxs_obj.set_nucl_data_obj(nucl_data_obj);
  rxs_obj.set_resonance_grid(sre_grid);

  ofstream fout;
  fout.open("./comp_njoy/calc_data_aa.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);

  vector<Real8>          ene_data;
  vector<vector<Real8> > sig_data;
  ene_data.resize(5);
  sig_data.resize(5);
  ene_data[0] = 1.0e+0;
  ene_data[1] = 5.0e+0;
  ene_data[2] = 1.0e+1;
  ene_data[3] = 5.0e+1;
  ene_data[4] = 1.0e+2;

  ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  MF02MT151Converter mf02mt151conv;
  MF02MT151Parser par_obj;
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  int i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      if( i==0 && j==0 )
      {
        fout << static_cast<int>(ene_data.size()) << endl;
        for(int k=0; k<static_cast<int>(ene_data.size()); k++)
        {
          fout << ene_data[k] << " ";
        }
        fout << endl;

        fout << rxs_obj.get_temp() << endl;
        fout << par_obj.get_NLS()[i][j] << " " << par_obj.get_LI_LRU01_LRF04()[i][j] << endl;
        fout << par_obj.get_AWRI()[i][j][0] << " " << par_obj.get_AP()[i][j] << endl;
        
        if( par_obj.get_LI_LRU01_LRF04()[i][j] != 6 )
        {
          for(int k=0; k<6; k++)
          {
            fout << par_obj.get_AT_LRU01_LRF04()[i][j][k] << " ";
          }
          fout << endl;
        }
        if( par_obj.get_LI_LRU01_LRF04()[i][j] != 5 )
        {
          for(int k=0; k<6; k++)
          {
            fout << par_obj.get_AF_LRU01_LRF04()[i][j][k] << " ";
          }
          fout << endl;
        }
        for(int k=0; k<6; k++)
        {
          fout << par_obj.get_AC_LRU01_LRF04()[i][j][k] << " ";
        }
        fout << endl;
        
        for(int k=0; k<static_cast<int>(par_obj.get_NLS()[i][j]); k++)
        {
          fout << par_obj.get_NJS_LRU01_LRF04()[i][j][k] << endl;

          for(int l=0; l<static_cast<int>(par_obj.get_NJS_LRU01_LRF04()[i][j][k]); l++)
          {
            fout << par_obj.get_NLJ_LRU01_LRF04()[i][j][k][l] << endl;
            
            for(int m=0; m<static_cast<int>(par_obj.get_NLJ_LRU01_LRF04()[i][j][k][l]); m++)
            {
              fout << par_obj.get_DET_LRU01_LRF04()[i][j][k][l][m] << " "
                   << par_obj.get_DWT_LRU01_LRF04()[i][j][k][l][m] << " "
                   << par_obj.get_GRT_LRU01_LRF04()[i][j][k][l][m] << " "
                   << par_obj.get_GIT_LRU01_LRF04()[i][j][k][l][m] << " "
                   << par_obj.get_DEF_LRU01_LRF04()[i][j][k][l][m] << " "
                   << par_obj.get_DWF_LRU01_LRF04()[i][j][k][l][m] << " "
                   << par_obj.get_GRF_LRU01_LRF04()[i][j][k][l][m] << " "
                   << par_obj.get_GIF_LRU01_LRF04()[i][j][k][l][m] << " "
                   << par_obj.get_DEC_LRU01_LRF04()[i][j][k][l][m] << " "
                   << par_obj.get_DWC_LRU01_LRF04()[i][j][k][l][m] << " "
                   << par_obj.get_GRC_LRU01_LRF04()[i][j][k][l][m] << " "
                   << par_obj.get_GIC_LRU01_LRF04()[i][j][k][l][m] << endl;
            }
          }
        }
        fout << endl;
        fout << "======================================================================" << endl;
        fout << endl;
      }
    }
  }
  fout.close();
}
// */ 

////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(AdlerAdler_check02)
{
  string file_name = "./for_test/test_mf02mt151_lru01lrf04.dat";
  
  Endf6ParserNoCov cp_obj;
  cp_obj.set_file_name(file_name);

  Endf6Converter    conv_obj;
  NuclearDataObject nucl_data_obj;
  conv_obj.convert_endf_format_to_frendy(cp_obj, nucl_data_obj);
  
  ResonanceXSCalculator rxs_obj;
  rxs_obj.set_nucl_data_obj(nucl_data_obj);
  
  vector<Real8>          ene_data;
  vector<vector<Real8> > sig_data;
  ene_data.resize(5);
  sig_data.resize(5);
  ene_data[0] = 1.0e+0;
  ene_data[1] = 5.0e+0;
  ene_data[2] = 1.0e+1;
  ene_data[3] = 5.0e+1;
  ene_data[4] = 1.0e+2;

  //*
  ofstream fout;
  fout.open("./comp_njoy/calc_result_aa_frendy.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);

  ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  MF02MT151Converter mf02mt151conv;
  MF02MT151Parser par_obj;
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  int i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      if( i==0 && j==0 )
      {
        int k_max = static_cast<int>(ene_data.size());
        fout << "Adler-Adler" << endl;
        for(int t=0; t<5; t++)
        {
          Real8 temp_val = static_cast<Real8>(t)*100.0;
          fout << "Temp: " << temp_val << "[K]" << endl;
          rxs_obj.set_temp(temp_val);
          for(int k=0; k<k_max; k++)
          {
            rxs_obj.calc_reso_xs_each_case(i, j, ene_data[k], sig_data[k]);
            fout << ene_data[k] << " ";
            for(int l=0; l<4; l++)
            {
              fout << sig_data[k][l] << " ";
            }
            fout << endl;
          }
          fout << endl;
        }
      }
    }
  }
  fout.close();
  // */ 
  
  
  ifstream fin;
  fin.open("./comp_njoy/calc_result_aa.dat");
  
  string        line_data;
  Real8         tmp_real;
  vector<Real8> sig_ref;
  sig_ref.resize(4);
  
  //ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  //MF02MT151Converter mf02mt151conv;
  //MF02MT151Parser par_obj;
  //mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);
  //int i_max = static_cast<int>(par_obj.get_NIS());

  reso_data_obj = rxs_obj.get_reso_data_obj();
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      if( i==0 && j==0 )
      {
        int k_max = static_cast<int>(ene_data.size());
        getline(fin, line_data);  //Adler-Adler
        for(int t=0; t<5; t++)
        {
          Real8 temp_val = static_cast<Real8>(t)*100.0;
          rxs_obj.set_temp(temp_val);
          getline(fin, line_data); // TEMP: X.X
          for(int k=0; k<k_max; k++)
          {
            rxs_obj.calc_reso_xs_each_case(i, j, ene_data[k], sig_data[k]);
            fin >> tmp_real >> sig_ref[0] >> sig_ref[1] >> sig_ref[2] >> sig_ref[3];
            for(int l=0; l<4; l++)
            {
              BOOST_CHECK_LE(fabs(sig_data[k][l] - sig_ref[l]),fabs(5.0e-5*sig_ref[l]+1.0e-10));
            }
          }
          getline(fin, line_data);
          getline(fin, line_data);
        }
      }
    }
  }
  fin.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

//*
BOOST_AUTO_TEST_CASE(Unresolved_CaseC_check01)
{
  string file_name = "../EndfUtils/MFxxMTyyyParser/for_test/n_9228_92-U-235.dat";
  
  Endf6ParserNoCov cp_obj;
  cp_obj.set_file_name(file_name);

  Endf6Converter    conv_obj;
  NuclearDataObject nucl_data_obj;
  conv_obj.convert_endf_format_to_frendy(cp_obj, nucl_data_obj);

  ResonanceEnergyGridLinearizer sre_obj;
  sre_obj.set_nucl_data_obj(nucl_data_obj);

  vector<vector<vector<Real8> > > sre_grid = sre_obj.get_resonance_grid();

  ResonanceXSCalculator rxs_obj;
  rxs_obj.set_temp(100.0);
  rxs_obj.set_err(1.0e-2);
  rxs_obj.set_err_max(2.0e-1);
  rxs_obj.set_err_int(5.0e-3);
  
  vector<vector<vector<Real8> > >          reso_grid;
  vector<vector<vector<vector<Real8> > > > reso_xs;

  BOOST_CHECK_LE(fabs((rxs_obj.get_temp()    - 100.0 )/100.0 ), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err()     - 1.0e-2)/1.0e-2), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err_max() - 2.0e-1)/2.0e-1), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err_int() - 5.0e-3)/5.0e-3), 1.0e-5);

  rxs_obj.set_temp(0.0);
  rxs_obj.set_nucl_data_obj(nucl_data_obj);
  rxs_obj.set_resonance_grid(sre_grid);

  ofstream fout;
  fout.open("./comp_njoy/calc_data_case_c.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);

  vector<Real8>          ene_data;
  vector<vector<Real8> > sig_data;
  ene_data.resize(5);
  sig_data.resize(5);
  ene_data[0] = 1.0e+0;
  ene_data[1] = 5.0e+0;
  ene_data[2] = 1.0e+1;
  ene_data[3] = 5.0e+1;
  ene_data[4] = 1.0e+2;

  ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  MF02MT151Converter mf02mt151conv;
  MF02MT151Parser par_obj;
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  int i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      if( i==0 && j==1 )
      {
        fout << static_cast<int>(ene_data.size()) << endl;
        for(int k=0; k<static_cast<int>(ene_data.size()); k++)
        {
          fout << ene_data[k] << " ";
        }
        fout << endl;

        fout << rxs_obj.get_temp() << endl;
        fout << par_obj.get_LFW()[i]     << " " << par_obj.get_NRO()[i][j] << " "
             << par_obj.get_NAPS()[i][j] << " " << par_obj.get_NLS()[i][j] << endl;
        fout << par_obj.get_AWRI()[i][j][0] << " " << par_obj.get_SPI()[i][j] << " " 
             << par_obj.get_AP()[i][j] << endl;
        
        for(int k=0; k<static_cast<int>(par_obj.get_NLS()[i][j]); k++)
        {
          fout << par_obj.get_NJS_LRU02()[i][j][k] << " " << par_obj.get_L_LRU02()[i][j][k] << endl;

          for(int l=0; l<static_cast<int>(par_obj.get_NJS_LRU02()[i][j][k]); l++)
          {
            fout << par_obj.get_AJ_LRU02_C()[i][j][k][l]   << " " << par_obj.get_INT_LRU02_C()[i][j][k][l] << " " 
                 << par_obj.get_NE_LRU02_C()[i][j][k][l]   << endl;
            fout << par_obj.get_AMUX_LRU02_C()[i][j][k][l] << " " << par_obj.get_AMUN_LRU02_C()[i][j][k][l] << " " 
                 << par_obj.get_AMUG_LRU02_C()[i][j][k][l] << " " << par_obj.get_AMUF_LRU02_C()[i][j][k][l] << endl;
            
            for(int m=0; m<static_cast<int>(par_obj.get_NE_LRU02_C()[i][j][k][l]); m++)
            {
              fout << par_obj.get_ES_LRU02_C()[i][j][k][l][m] << " ";
            }
            fout << endl;
            
            for(int m=0; m<static_cast<int>(par_obj.get_NE_LRU02_C()[i][j][k][l]); m++)
            {
              fout << par_obj.get_D_LRU02_C()[i][j][k][l][m] << " ";
            }
            fout << endl;
            
            for(int m=0; m<static_cast<int>(par_obj.get_NE_LRU02_C()[i][j][k][l]); m++)
            {
              fout << par_obj.get_GX_LRU02_C()[i][j][k][l][m] << " ";
            }
            fout << endl;
            
            for(int m=0; m<static_cast<int>(par_obj.get_NE_LRU02_C()[i][j][k][l]); m++)
            {
              fout << par_obj.get_GNO_LRU02_C()[i][j][k][l][m] << " ";
            }
            fout << endl;
            
            for(int m=0; m<static_cast<int>(par_obj.get_NE_LRU02_C()[i][j][k][l]); m++)
            {
              fout << par_obj.get_GG_LRU02_C()[i][j][k][l][m] << " ";
            }
            fout << endl;
            
            for(int m=0; m<static_cast<int>(par_obj.get_NE_LRU02_C()[i][j][k][l]); m++)
            {
              fout << par_obj.get_GF_LRU02_C()[i][j][k][l][m] << " ";
            }
            fout << endl;
          }
        }
        fout << endl;
        fout << "======================================================================" << endl;
        fout << endl;
      }
    }
  }
  fout.close();
}
// */ 

////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(Unresolved_CaseC_check02)
{
  string file_name = "../EndfUtils/MFxxMTyyyParser/for_test/n_9228_92-U-235.dat";
  
  Endf6ParserNoCov cp_obj;
  cp_obj.set_file_name(file_name);

  Endf6Converter    conv_obj;
  NuclearDataObject nucl_data_obj;
  conv_obj.convert_endf_format_to_frendy(cp_obj, nucl_data_obj);
  
  ResonanceXSCalculator rxs_obj;
  rxs_obj.set_nucl_data_obj(nucl_data_obj);
  
  vector<Real8>          ene_data;
  vector<vector<Real8> > sig_data;
  ene_data.resize(5);
  sig_data.resize(5);
  ene_data[0] = 1.0e+0;
  ene_data[1] = 5.0e+0;
  ene_data[2] = 1.0e+1;
  ene_data[3] = 5.0e+1;
  ene_data[4] = 1.0e+2;

  //*
  ofstream fout;
  fout.open("./comp_njoy/calc_result_case_c_frendy.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);

  ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  MF02MT151Converter mf02mt151conv;
  MF02MT151Parser par_obj;
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  int i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      if( i==0 && j==1 )
      {
        int k_max = static_cast<int>(ene_data.size());
        fout << "Unresolved Resonance (Case C)" << endl;
        for(int t=0; t<5; t++)
        {
          Real8 temp_val = static_cast<Real8>(t)*100.0;
          fout << "Temp: " << temp_val << "[K]" << endl;
          rxs_obj.set_temp(temp_val);
          for(int k=0; k<k_max; k++)
          {
            sig_data[k].resize(4);
            for(int l=0; l<4; l++)
            {
              sig_data[k][l] = 0.0;
            }
            //rxs_obj.calc_reso_xs_each_case(i, j, ene_data[k], sig_data[k]);
            rxs_obj.set_multi_array_data(i, j);
            rxs_obj.set_multi_array_data_unreso_c(i, j);
            rxs_obj.calc_reso_xs_unreso_c(i, j, ene_data[k], sig_data[k]);
            fout << ene_data[k] << " ";
            for(int l=0; l<4; l++)
            {
              fout << sig_data[k][l] << " ";
            }
            fout << endl;
          }
          fout << endl;
        }
      }
    }
  }
  fout.close();
  // */ 
  
  ifstream fin;
  fin.open("./comp_njoy/calc_result_case_c.dat");
  
  string        line_data;
  Real8         tmp_real;
  vector<Real8> sig_ref;
  sig_ref.resize(4);
  
  //ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  //MF02MT151Converter mf02mt151conv;
  //MF02MT151Parser par_obj;
  //mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);
  //int i_max = static_cast<int>(par_obj.get_NIS());

  reso_data_obj = rxs_obj.get_reso_data_obj();
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      if( i==0 && j==1 )
      {
        int k_max = static_cast<int>(ene_data.size());
        getline(fin, line_data); //Unresolved Resonance (Case C)
        for(int t=0; t<5; t++)
        {
          Real8 temp_val = static_cast<Real8>(t)*100.0;
          rxs_obj.set_temp(temp_val);
          getline(fin, line_data); //TEMP: X.X
          for(int k=0; k<k_max; k++)
          {
            sig_data[k].resize(4);
            for(int l=0; l<4; l++)
            {
              sig_data[k][l] = 0.0;
            }
            //rxs_obj.calc_reso_xs_each_case(i, j, ene_data[k], sig_data[k]);
            rxs_obj.set_multi_array_data(i, j);
            rxs_obj.set_multi_array_data_unreso_c(i, j);
            rxs_obj.calc_reso_xs_unreso_c(i, j, ene_data[k], sig_data[k]);
            fin >> tmp_real >> sig_ref[0] >> sig_ref[1] >> sig_ref[2] >> sig_ref[3];
            for(int l=0; l<4; l++)
            {
              BOOST_CHECK_LE(fabs(sig_data[k][l] - sig_ref[l]),fabs(5.0e-5*sig_ref[l]+1.0e-10));
            }
          }
          getline(fin, line_data);
          getline(fin, line_data);
        }
      }
    }
  }
  fin.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(width_fluctuation_factor_check01)
{
  ResonanceXSCalculator rxs_obj;
  
  vector<Real8> g_val, sig_val;
  g_val.resize(5);
  sig_val.resize(4);
  
  g_val[0] = -2.5;
  g_val[1] = -1.0;
  g_val[2] =  0.0;
  g_val[3] =  1.0;
  g_val[4] =  2.5;
  
  int mu, nu, lamda;
  
  //*
  ofstream fout;
  fout.open("./comp_njoy/calc_width_fluctuation_factor_frendy.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);
  for(int n=0; n<2; n++)
  {
    if( n == 0 )
    {
      mu    = 1;
      nu    = 2;
      lamda = 3;
    }
    else
    {
      mu    = 4;
      nu    = 3;
      lamda = 2;
    }
    
    fout << "Mu, Nu, Lamda:" << mu << " " << nu << " " << lamda << endl; 
    for(int i=0; i<5; i++)
    {
      for(int j=0; j<5; j++)
      {
        for(int k=0; k<5; k++)
        {
          for(int l=0; l<5; l++)
          {
            for(int m=0; m<4; m++)
            {
              sig_val[m] = 0.0;
            }
            rxs_obj.calc_width_fluctuation_factor(g_val[i], g_val[j], g_val[k],
                                                  mu, nu, lamda, sig_val, g_val[l]);
            fout << "  " << g_val[i] << " " << g_val[j] << " " << g_val[k] << " "
                         << sig_val[1] << " " << sig_val[2] << " " << sig_val[3] << endl;
          }
        }
      }
    }
    fout << endl;
    fout << "==================================================" << endl;
    fout << endl;
  }
  fout.close();
  // */ 
  
  
  ifstream fin;
  fin.open("./comp_njoy/calc_width_fluctuation_factor.dat");
  
  string        line_data;
  Real8         tmp_real;
  vector<Real8> sig_ref;
  sig_ref.resize(4);
  
  for(int no=0; no<2; no++)
  {
    if( no == 0 )
    {
      mu    = 1;
      nu    = 2;
      lamda = 3;
    }
    else
    {
      mu    = 4;
      nu    = 3;
      lamda = 2;
    }
    getline(fin, line_data); //Mu, Nu, Lamda
    for(int i=0; i<5; i++)
    {
      for(int j=0; j<5; j++)
      {
        for(int k=0; k<5; k++)
        {
          for(int l=0; l<5; l++)
          {
            for(int m=0; m<4; m++)
            {
              sig_val[m] = 0.0;
            }
            rxs_obj.calc_width_fluctuation_factor(g_val[i], g_val[j], g_val[k],
                                                  mu, nu, lamda, sig_val, g_val[l]);
            fin >> tmp_real >> tmp_real >> tmp_real >> sig_ref[1] >> sig_ref[2] >> sig_ref[3];
            for(int m=1; m<4; m++)
            {
              BOOST_CHECK_LE(fabs(sig_val[m] - sig_ref[m]),fabs(5.0e-5*sig_ref[m]+1.0e-10));
            }
          }
        }
      }
    }
    getline(fin, line_data);
    getline(fin, line_data);
    getline(fin, line_data); // ==================================================
    getline(fin, line_data);
  }
  fin.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

//*
BOOST_AUTO_TEST_CASE(Unresolved_CaseA_check01)
{
  string file_name = "./for_test/test_mf02mt151_lru02_case_a.dat";
  
  Endf6ParserNoCov cp_obj;
  cp_obj.set_file_name(file_name);

  Endf6Converter    conv_obj;
  NuclearDataObject nucl_data_obj;
  conv_obj.convert_endf_format_to_frendy(cp_obj, nucl_data_obj);

  ResonanceEnergyGridLinearizer sre_obj;
  sre_obj.set_nucl_data_obj(nucl_data_obj);

  vector<vector<vector<Real8> > > sre_grid = sre_obj.get_resonance_grid();

  ResonanceXSCalculator rxs_obj;
  rxs_obj.set_temp(100.0);
  rxs_obj.set_err(1.0e-2);
  rxs_obj.set_err_max(2.0e-1);
  rxs_obj.set_err_int(5.0e-3);
  
  vector<vector<vector<Real8> > >          reso_grid;
  vector<vector<vector<vector<Real8> > > > reso_xs;

  BOOST_CHECK_LE(fabs((rxs_obj.get_temp()    - 100.0 )/100.0 ), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err()     - 1.0e-2)/1.0e-2), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err_max() - 2.0e-1)/2.0e-1), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err_int() - 5.0e-3)/5.0e-3), 1.0e-5);

  rxs_obj.set_temp(0.0);
  rxs_obj.set_nucl_data_obj(nucl_data_obj);
  rxs_obj.set_resonance_grid(sre_grid);

  ofstream fout;
  fout.open("./comp_njoy/calc_data_case_a.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);

  vector<Real8>          ene_data;
  vector<vector<Real8> > sig_data;
  ene_data.resize(5);
  sig_data.resize(5);
  ene_data[0] = 1.0e+0;
  ene_data[1] = 5.0e+0;
  ene_data[2] = 1.0e+1;
  ene_data[3] = 5.0e+1;
  ene_data[4] = 1.0e+2;

  ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  MF02MT151Converter mf02mt151conv;
  MF02MT151Parser par_obj;
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  int i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      if( i==0 && j==1 )
      {
        fout << static_cast<int>(ene_data.size()) << endl;
        for(int k=0; k<static_cast<int>(ene_data.size()); k++)
        {
          fout << ene_data[k] << " ";
        }
        fout << endl;
        fout << rxs_obj.get_temp() << endl;
        fout << par_obj.get_LFW()[i]     << " " << par_obj.get_NRO()[i][j] << " "
             << par_obj.get_NAPS()[i][j] << endl;
        fout << par_obj.get_AWRI()[i][j][0] << " " << par_obj.get_SPI()[i][j] << " " 
             << par_obj.get_AP()[i][j] << endl;
        fout << par_obj.get_NLS()[i][j] << endl;
        
        for(int k=0; k<static_cast<int>(par_obj.get_NLS()[i][j]); k++)
        {
          fout << par_obj.get_NJS_LRU02()[i][j][k] << " " << par_obj.get_L_LRU02()[i][j][k] << endl;

          for(int l=0; l<static_cast<int>(par_obj.get_NJS_LRU02()[i][j][k]); l++)
          {
            fout << par_obj.get_D_LRU02_A()[i][j][k][l]    << " " << par_obj.get_AJ_LRU02_A()[i][j][k][l]  << " "
                 << par_obj.get_AMUN_LRU02_A()[i][j][k][l] << " " << par_obj.get_GNO_LRU02_A()[i][j][k][l] << " "
                 << par_obj.get_GG_LRU02_A()[i][j][k][l]   << endl;
          }
        }
        fout << endl;
        fout << "======================================================================" << endl;
        fout << endl;
      }
    }
  }
  fout.close();
}
// */ 

////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(Unresolved_CaseA_check02)
{
  string file_name = "./for_test/test_mf02mt151_lru02_case_a.dat";
  
  Endf6ParserNoCov cp_obj;
  cp_obj.set_file_name(file_name);

  Endf6Converter    conv_obj;
  NuclearDataObject nucl_data_obj;
  conv_obj.convert_endf_format_to_frendy(cp_obj, nucl_data_obj);
  
  ResonanceXSCalculator rxs_obj;
  rxs_obj.set_nucl_data_obj(nucl_data_obj);
  
  vector<Real8>          ene_data;
  vector<vector<Real8> > sig_data;
  ene_data.resize(5);
  sig_data.resize(5);
  ene_data[0] = 1.0e+0;
  ene_data[1] = 5.0e+0;
  ene_data[2] = 1.0e+1;
  ene_data[3] = 5.0e+1;
  ene_data[4] = 1.0e+2;

  //*
  ofstream fout;
  fout.open("./comp_njoy/calc_result_case_a_frendy.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);

  ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  MF02MT151Converter mf02mt151conv;
  MF02MT151Parser par_obj;
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  int i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      if( i==0 && j==1 )
      {
        int k_max = static_cast<int>(ene_data.size());
        fout << "Unresolved Resonance (Case A)" << endl;
        for(int t=0; t<5; t++)
        {
          Real8 temp_val = static_cast<Real8>(t)*100.0;
          fout << "Temp: " << temp_val << "[K]" << endl;
          rxs_obj.set_temp(temp_val);
          for(int k=0; k<k_max; k++)
          {
            sig_data[k].resize(4);
            for(int l=0; l<4; l++)
            {
              sig_data[k][l] = 0.0;
            }
            //rxs_obj.calc_reso_xs_each_case(i, j, ene_data[k], sig_data[k]);
            rxs_obj.set_multi_array_data(i, j);
            rxs_obj.set_multi_array_data_unreso_a(i, j);
            rxs_obj.calc_reso_xs_unreso_a(i, j, ene_data[k], sig_data[k]);
            fout << ene_data[k] << " ";
            for(int l=0; l<4; l++)
            {
              fout << sig_data[k][l] << " ";
            }
            fout << endl;
          }
          fout << endl;
        }
      }
    }
  }
  fout.close();
  // */ 
  
  ifstream fin;
  fin.open("./comp_njoy/calc_result_case_a.dat");
  
  string        line_data;
  Real8         tmp_real;
  vector<Real8> sig_ref;
  sig_ref.resize(4);
  
  //ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  //MF02MT151Converter mf02mt151conv;
  //MF02MT151Parser par_obj;
  //mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);
  //int i_max = static_cast<int>(par_obj.get_NIS());

  reso_data_obj = rxs_obj.get_reso_data_obj();
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      if( i==0 && j==1 )
      {
        int k_max = static_cast<int>(ene_data.size());
        getline(fin, line_data); //Unresolved Resonance (Case A)
        for(int t=0; t<5; t++)
        {
          Real8 temp_val = static_cast<Real8>(t)*100.0;
          rxs_obj.set_temp(temp_val);
          getline(fin, line_data); // TEMP: X.X
          for(int k=0; k<k_max; k++)
          {
            sig_data[k].resize(4);
            for(int l=0; l<4; l++)
            {
              sig_data[k][l] = 0.0;
            }
            //rxs_obj.calc_reso_xs_each_case(i, j, ene_data[k], sig_data[k]);
            rxs_obj.set_multi_array_data(i, j);
            rxs_obj.set_multi_array_data_unreso_a(i, j);
            rxs_obj.calc_reso_xs_unreso_a(i, j, ene_data[k], sig_data[k]);
            fin >> tmp_real >> sig_ref[0] >> sig_ref[1] >> sig_ref[2] >> sig_ref[3];
            for(int l=0; l<4; l++)
            {
              BOOST_CHECK_LE(fabs(sig_data[k][l] - sig_ref[l]),fabs(5.0e-5*sig_ref[l]+1.0e-10));
            }
          }
          getline(fin, line_data);
          getline(fin, line_data);
        }
      }
    }
  }
  fin.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//*
BOOST_AUTO_TEST_CASE(Unresolved_CaseB_check01)
{
  string file_name = "./for_test/test_mf02mt151_lru02_case_b.dat";
  
  Endf6ParserNoCov cp_obj;
  cp_obj.set_file_name(file_name);

  Endf6Converter    conv_obj;
  NuclearDataObject nucl_data_obj;
  conv_obj.convert_endf_format_to_frendy(cp_obj, nucl_data_obj);

  ResonanceEnergyGridLinearizer sre_obj;
  sre_obj.set_nucl_data_obj(nucl_data_obj);

  vector<vector<vector<Real8> > > sre_grid = sre_obj.get_resonance_grid();

  ResonanceXSCalculator rxs_obj;
  rxs_obj.set_temp(100.0);
  rxs_obj.set_err(1.0e-2);
  rxs_obj.set_err_max(2.0e-1);
  rxs_obj.set_err_int(5.0e-3);
  
  vector<vector<vector<Real8> > >          reso_grid;
  vector<vector<vector<vector<Real8> > > > reso_xs;

  BOOST_CHECK_LE(fabs((rxs_obj.get_temp()    - 100.0 )/100.0 ), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err()     - 1.0e-2)/1.0e-2), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err_max() - 2.0e-1)/2.0e-1), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err_int() - 5.0e-3)/5.0e-3), 1.0e-5);

  rxs_obj.set_temp(0.0);
  rxs_obj.set_nucl_data_obj(nucl_data_obj);
  rxs_obj.set_resonance_grid(sre_grid);

  ofstream fout;
  fout.open("./comp_njoy/calc_data_case_b.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);

  vector<Real8>          ene_data;
  vector<vector<Real8> > sig_data;
  ene_data.resize(5);
  sig_data.resize(5);
  ene_data[0] = 1.0e+0;
  ene_data[1] = 5.0e+0;
  ene_data[2] = 1.0e+1;
  ene_data[3] = 5.0e+1;
  ene_data[4] = 1.0e+2;

  ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  MF02MT151Converter mf02mt151conv;
  MF02MT151Parser par_obj;
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  int i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      if( i==0 && j==1 )
      {
        fout << static_cast<int>(ene_data.size()) << endl;
        for(int k=0; k<static_cast<int>(ene_data.size()); k++)
        {
          fout << ene_data[k] << " ";
        }
        fout << endl;
        fout << rxs_obj.get_temp() << endl;
        fout << par_obj.get_LFW()[i]     << " " << par_obj.get_NRO()[i][j] << " "
             << par_obj.get_NAPS()[i][j] << endl;
        fout << par_obj.get_AWRI()[i][j][0] << " " << par_obj.get_SPI()[i][j] << " " 
             << par_obj.get_AP()[i][j]  << endl;
        fout << par_obj.get_NLS()[i][j] << " "  << par_obj.get_NE_LRU02_B()[i][j] << endl;
        for(int k=0; k<static_cast<int>(par_obj.get_NE_LRU02_B()[i][j]); k++)
        {
          fout << par_obj.get_ES_LRU02_B()[i][j][k] << " ";
        }
        fout << endl;
        
        for(int k=0; k<static_cast<int>(par_obj.get_NLS()[i][j]); k++)
        {
          fout << par_obj.get_NJS_LRU02()[i][j][k] << " " << par_obj.get_L_LRU02()[i][j][k] << endl;

          for(int l=0; l<static_cast<int>(par_obj.get_NJS_LRU02()[i][j][k]); l++)
          {
            fout << par_obj.get_MUF_LRU02_B()[i][j][k][l]  << " " << par_obj.get_NE_LRU02_B()[i][j]+6 << endl;
            fout << par_obj.get_D_LRU02_B()[i][j][k][l]    << " " << par_obj.get_AJ_LRU02_B()[i][j][k][l]  << " "
                 << par_obj.get_AMUN_LRU02_B()[i][j][k][l] << " " << par_obj.get_GNO_LRU02_B()[i][j][k][l] << " "
                 << par_obj.get_GG_LRU02_B()[i][j][k][l]   << endl;
            for(int m=0; m<static_cast<int>(par_obj.get_NE_LRU02_B()[i][j]); m++)
            {
              fout << par_obj.get_GF_LRU02_B()[i][j][k][l][m] << " ";
            }
            fout << endl;
          }
        }
        fout << endl;
        fout << "======================================================================" << endl;
        fout << endl;
      }
    }
  }
  fout.close();
}
// */ 

////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(Unresolved_CaseB_check02)
{
  string file_name = "./for_test/test_mf02mt151_lru02_case_b.dat";
  
  Endf6ParserNoCov cp_obj;
  cp_obj.set_file_name(file_name);

  Endf6Converter    conv_obj;
  NuclearDataObject nucl_data_obj;
  conv_obj.convert_endf_format_to_frendy(cp_obj, nucl_data_obj);
  
  ResonanceXSCalculator rxs_obj;
  rxs_obj.set_nucl_data_obj(nucl_data_obj);
  
  vector<Real8>          ene_data;
  vector<vector<Real8> > sig_data;
  ene_data.resize(5);
  sig_data.resize(5);
  ene_data[0] = 1.0e+0;
  ene_data[1] = 5.0e+0;
  ene_data[2] = 1.0e+1;
  ene_data[3] = 5.0e+1;
  ene_data[4] = 1.0e+2;

  //*
  ofstream fout;
  fout.open("./comp_njoy/calc_result_case_b_frendy.dat");
  fout.precision(15);
  fout.setf(ios::scientific);
  fout.setf(ios::showpos);
  fout.setf(ios::showpoint);

  ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  MF02MT151Converter mf02mt151conv;
  MF02MT151Parser par_obj;
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  int i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      if( i==0 && j==1 )
      {
        int k_max = static_cast<int>(ene_data.size());
        fout << "Unresolved Resonance (Case B)" << endl;
        for(int t=0; t<5; t++)
        {
          Real8 temp_val = static_cast<Real8>(t)*100.0;
          fout << "Temp: " << temp_val << "[K]" << endl;
          rxs_obj.set_temp(temp_val);
          for(int k=0; k<k_max; k++)
          {
            sig_data[k].resize(4);
            for(int l=0; l<4; l++)
            {
              sig_data[k][l] = 0.0;
            }
            //rxs_obj.calc_reso_xs_each_case(i, j, ene_data[k], sig_data[k]);
            rxs_obj.set_multi_array_data(i, j);
            rxs_obj.set_multi_array_data_unreso_a(i, j);
            rxs_obj.calc_reso_xs_unreso_a(i, j, ene_data[k], sig_data[k]);
            fout << ene_data[k] << " ";
            for(int l=0; l<4; l++)
            {
              fout << sig_data[k][l] << " ";
            }
            fout << endl;
          }
          fout << endl;
        }
      }
    }
  }
  fout.close();
  // */ 
  
  
  ifstream fin;
  fin.open("./comp_njoy/calc_result_case_b.dat");
  
  string        line_data;
  Real8         tmp_real;
  vector<Real8> sig_ref;
  sig_ref.resize(4);
  
  //ResonanceDataContainer reso_data_obj = rxs_obj.get_reso_data_obj();
  //MF02MT151Converter mf02mt151conv;
  //MF02MT151Parser par_obj;
  //mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);
  //int i_max = static_cast<int>(par_obj.get_NIS());

  reso_data_obj = rxs_obj.get_reso_data_obj();
  mf02mt151conv.convert_frendy_to_endf_format(reso_data_obj, par_obj);

  i_max = static_cast<int>(par_obj.get_NIS());
  for(int i=0; i<i_max; i++)
  {
    int j_max = static_cast<int>(par_obj.get_NER()[i]);
    for(int j=0; j<j_max; j++)
    {
      if( i==0 && j==1 )
      {
        int k_max = static_cast<int>(ene_data.size());
        getline(fin, line_data); //Unresolved Resonance (Case B)
        for(int t=0; t<5; t++)
        {
          Real8 temp_val = static_cast<Real8>(t)*100.0;
          rxs_obj.set_temp(temp_val);
          getline(fin, line_data); //TEMP: X.X
          for(int k=0; k<k_max; k++)
          {
            sig_data[k].resize(4);
            for(int l=0; l<4; l++)
            {
              sig_data[k][l] = 0.0;
            }
            //rxs_obj.calc_reso_xs_each_case(i, j, ene_data[k], sig_data[k]);
            rxs_obj.set_multi_array_data(i, j);
            rxs_obj.set_multi_array_data_unreso_a(i, j);
            rxs_obj.calc_reso_xs_unreso_a(i, j, ene_data[k], sig_data[k]);
            fin >> tmp_real >> sig_ref[0] >> sig_ref[1] >> sig_ref[2] >> sig_ref[3];
            for(int l=0; l<4; l++)
            {
              BOOST_CHECK_LE(fabs(sig_data[k][l] - sig_ref[l]),fabs(5.0e-5*sig_ref[l]+1.0e-10));
            }
          }
          getline(fin, line_data);
          getline(fin, line_data);
        }
      }
    }
  }
  fin.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(unify_energy_grid)
{
  vector<vector<vector<Real8> > >          ene_ori, ene_cal, ene_ref;
  vector<vector<vector<vector<Real8> > > > sig_ori, sig_cal, sig_ref;
  vector<Real8> sig_tmp;
  
  int         nis = 2;
  vector<int> ner;
  ner.resize(nis);
  ene_ori.resize(nis);
  sig_ori.resize(nis);
  ene_ref.resize(nis);
  sig_ref.resize(nis);
  for(int i=0; i<nis; i++)
  {
    ner[i] = i+2;
    ene_ori[i].resize(ner[i]);
    sig_ori[i].resize(ner[i]);
    ene_ref[i].resize(ner[i]);
    sig_ref[i].resize(ner[i]);
  }
  
  //Set dummy data for MF02MT151Parser class
  MF02MT151Parser mf02mt151;
  Real            za_val  = 92235.0;
  Real            awr_val = 92.0;
  Integer         nis_val = static_cast<Integer>(nis);

  vector<Real>    zai_val, abn_val;
  vector<Integer> ner_val, lfw_val;
  vector<vector<Real> >    el_val, eh_val, spi_val, ap_val;
  vector<vector<Integer> > reso_region_flg_val, xs_formula_flg_val, nro_val, naps_val, nls_val;
  vector<vector<vector<Real> > > awri_val;
  
  zai_val.resize(nis);
  abn_val.resize(nis);
  ner_val.resize(nis);
  lfw_val.resize(nis);
  el_val.resize(nis);
  eh_val.resize(nis);
  spi_val.resize(nis);
  ap_val.resize(nis);
  reso_region_flg_val.resize(nis);
  xs_formula_flg_val.resize(nis);
  nro_val.resize(nis);
  naps_val.resize(nis);
  nls_val.resize(nis);
  awri_val.resize(nis);
  for(int i=0; i<nis; i++)
  {
    zai_val[i] = 92.0 + static_cast<Real>(i+1);
    abn_val[i] = 0.5;
    ner_val[i] = static_cast<Integer>(ner[i]);
    lfw_val[i] = 0;
    
    el_val[i].resize(ner[i]);
    eh_val[i].resize(ner[i]);
    spi_val[i].resize(ner[i]);
    ap_val[i].resize(ner[i]);
    reso_region_flg_val[i].resize(ner[i]);
    xs_formula_flg_val[i].resize(ner[i]);
    nro_val[i].resize(ner[i]);
    naps_val[i].resize(ner[i]);
    nls_val[i].resize(ner[i]);
    awri_val[i].resize(ner[i]);
    for(int j=0; j<ner[i]; j++)
    {
      spi_val[i][j] = 0.5 + static_cast<Real>(i+1) + static_cast<Real>(j+1);
      ap_val[i][j]  = 1.0 + static_cast<Real>(i)   + static_cast<Real>(j);
     
      reso_region_flg_val[i][j] = 0;
      xs_formula_flg_val[i][j]  = 0;
      nro_val[i][j]             = 0;
      naps_val[i][j]            = 0;
      nls_val[i][j]             = 0;
      awri_val[i][j].push_back(1.0); 
    }
  }
  
  //Set original xs and energy grid data 
  //i=0, j=0
  ene_ori[0][0].push_back(1.10000E+01);
  ene_ori[0][0].push_back(1.40000E+01);
  ene_ori[0][0].push_back(1.70000E+01);
  ene_ori[0][0].push_back(2.00000E+01);
  ene_ori[0][0].push_back(2.30000E+01);
  el_val[0][0] = 1.10000E+01;
  eh_val[0][0] = 2.30000E+01;

  sig_tmp.push_back(1.00000E+02);
  sig_tmp.push_back(2.00000E+02);
  sig_tmp.push_back(3.00000E+02);
  sig_tmp.push_back(4.00000E+02);
  sig_tmp.push_back(5.00000E+02);
  int i_max = static_cast<int>(sig_tmp.size());
  sig_ori[0][0].resize(i_max);
  for(int i=0; i<i_max; i++)
  {
    Real8 multi = 1.0;
    for(int j=0; j<4; j++)
    {
      sig_ori[0][0][i].push_back(sig_tmp[i]*multi);
      multi *= 0.1;
    }
  }
  sig_tmp.clear();
  
  //i=0, j=1
  ene_ori[0][1].push_back(1.20000E+01);
  ene_ori[0][1].push_back(1.30000E+01);
  ene_ori[0][1].push_back(1.50000E+01);
  ene_ori[0][1].push_back(1.80000E+01);
  ene_ori[0][1].push_back(2.20000E+01);
  ene_ori[0][1].push_back(2.70000E+01);
  el_val[0][1] = 1.20000E+01;
  eh_val[0][1] = 2.70000E+01;

  sig_tmp.push_back(1.00000E+03);
  sig_tmp.push_back(2.00000E+03);
  sig_tmp.push_back(3.00000E+03);
  sig_tmp.push_back(4.00000E+03);
  sig_tmp.push_back(5.00000E+03);
  sig_tmp.push_back(6.00000E+03);
  i_max = static_cast<int>(sig_tmp.size());
  sig_ori[0][1].resize(i_max);
  for(int i=0; i<i_max; i++)
  {
    Real8 multi = 1.0;
    for(int j=0; j<4; j++)
    {
      sig_ori[0][1][i].push_back(sig_tmp[i]*multi);
      multi *= 0.1;
    }
  }
  sig_tmp.clear();
  
  //i=1, j=0
  ene_ori[1][0].push_back(1.30000E+01);
  ene_ori[1][0].push_back(1.40000E+01);
  ene_ori[1][0].push_back(1.60000E+01);
  ene_ori[1][0].push_back(1.70000E+01);
  ene_ori[1][0].push_back(1.90000E+01);
  ene_ori[1][0].push_back(2.00000E+01);
  ene_ori[1][0].push_back(2.10000E+01);
  ene_ori[1][0].push_back(2.40000E+01);
  ene_ori[1][0].push_back(2.50000E+01);
  ene_ori[1][0].push_back(2.60000E+01);
  el_val[1][0] = 1.30000E+01;
  eh_val[1][0] = 2.60000E+01;

  sig_tmp.push_back(1.00000E+04);
  sig_tmp.push_back(2.00000E+04);
  sig_tmp.push_back(3.00000E+04);
  sig_tmp.push_back(4.00000E+04);
  sig_tmp.push_back(5.00000E+04);
  sig_tmp.push_back(6.00000E+04);
  sig_tmp.push_back(7.00000E+04);
  sig_tmp.push_back(8.00000E+04);
  sig_tmp.push_back(9.00000E+04);
  sig_tmp.push_back(1.00000E+05);
  i_max = static_cast<int>(sig_tmp.size());
  sig_ori[1][0].resize(i_max);
  for(int i=0; i<i_max; i++)
  {
    Real8 multi = 1.0;
    for(int j=0; j<4; j++)
    {
      sig_ori[1][0][i].push_back(sig_tmp[i]*multi);
      multi *= 0.1;
    }
  }
  sig_tmp.clear();
  
  //i=1, j=1
  ene_ori[1][1].push_back(2.00000E+01);
  ene_ori[1][1].push_back(2.20000E+01);
  ene_ori[1][1].push_back(2.40000E+01);
  ene_ori[1][1].push_back(2.60000E+01);
  ene_ori[1][1].push_back(2.80000E+01);
  ene_ori[1][1].push_back(3.00000E+01);
  ene_ori[1][1].push_back(3.20000E+01);
  el_val[1][1] = 2.00000E+01;
  eh_val[1][1] = 3.20000E+01;

  sig_tmp.push_back(1.00000E+05);
  sig_tmp.push_back(2.00000E+05);
  sig_tmp.push_back(3.00000E+05);
  sig_tmp.push_back(4.00000E+05);
  sig_tmp.push_back(5.00000E+05);
  sig_tmp.push_back(6.00000E+05);
  sig_tmp.push_back(7.00000E+05);
  i_max = static_cast<int>(sig_tmp.size());
  sig_ori[1][1].resize(i_max);
  for(int i=0; i<i_max; i++)
  {
    Real8 multi = 1.0;
    for(int j=0; j<4; j++)
    {
      sig_ori[1][1][i].push_back(sig_tmp[i]*multi);
      multi *= 0.1;
    }
  }
  sig_tmp.clear();
  
  //i=1, j=2
  ene_ori[1][2].push_back(2.10000E+01);
  ene_ori[1][2].push_back(2.20000E+01);
  ene_ori[1][2].push_back(2.30000E+01);
  ene_ori[1][2].push_back(2.40000E+01);
  ene_ori[1][2].push_back(2.50000E+01);
  ene_ori[1][2].push_back(2.60000E+01);
  ene_ori[1][2].push_back(2.70000E+01);
  ene_ori[1][2].push_back(2.80000E+01);
  ene_ori[1][2].push_back(2.90000E+01);
  ene_ori[1][2].push_back(3.00000E+01);
  ene_ori[1][2].push_back(3.10000E+01);
  ene_ori[1][2].push_back(3.20000E+01);
  ene_ori[1][2].push_back(3.30000E+01);
  el_val[1][2] = 2.10000E+01;
  eh_val[1][2] = 3.30000E+01;

  sig_tmp.push_back(1.00000E+06);
  sig_tmp.push_back(2.00000E+06);
  sig_tmp.push_back(3.00000E+06);
  sig_tmp.push_back(4.00000E+06);
  sig_tmp.push_back(5.00000E+06);
  sig_tmp.push_back(6.00000E+06);
  sig_tmp.push_back(7.00000E+06);
  sig_tmp.push_back(8.00000E+06);
  sig_tmp.push_back(9.00000E+06);
  sig_tmp.push_back(1.00000E+07);
  sig_tmp.push_back(1.10000E+07);
  sig_tmp.push_back(1.20000E+07);
  sig_tmp.push_back(1.30000E+07);
  i_max = static_cast<int>(sig_tmp.size());
  sig_ori[1][2].resize(i_max);
  for(int i=0; i<i_max; i++)
  {
    Real8 multi = 1.0;
    for(int j=0; j<4; j++)
    {
      sig_ori[1][2][i].push_back(sig_tmp[i]*multi);
      multi *= 0.1;
    }
  }
  sig_tmp.clear();

  //Copy data to MF02MT151Parser class 
  mf02mt151.set_ZAR(za_val);
  mf02mt151.set_AWR(awr_val);
  mf02mt151.set_NIS(nis_val);
  
  mf02mt151.set_ZAI(zai_val);
  mf02mt151.set_ABN(abn_val);
  mf02mt151.set_NER(ner_val);
  mf02mt151.set_LFW(lfw_val);
  
  mf02mt151.set_AWRI(awri_val);
  mf02mt151.set_EL(el_val);
  mf02mt151.set_EH(eh_val);
  mf02mt151.set_SPI(spi_val);
  mf02mt151.set_AP(ap_val);
  mf02mt151.set_LRU(reso_region_flg_val);
  mf02mt151.set_LRF(xs_formula_flg_val);
  mf02mt151.set_NRO(nro_val);
  mf02mt151.set_NAPS(naps_val);
  mf02mt151.set_NLS(nls_val);
  
  //Set calc resonance cross section by ResonanceXSCalculator class  
  Endf6ParserNoCov cp_obj;
  cp_obj.set_mf02_mt151_data(mf02mt151);

  Endf6Converter    conv_obj;
  NuclearDataObject nucl_data_obj;
  conv_obj.convert_endf_format_to_frendy(cp_obj, nucl_data_obj);
  
  ResonanceXSCalculator rxs_obj;
  rxs_obj.set_temp(0.0);
  rxs_obj.set_err(1.0e-2);
  rxs_obj.set_err_max(2.0e-1);
  rxs_obj.set_err_int(5.0e-3);
  rxs_obj.set_nucl_data_obj(nucl_data_obj);

  rxs_obj.ene_array = ene_ori;
  rxs_obj.sig_array = sig_ori;
  rxs_obj.unify_energy_grid();
  
  ene_cal = rxs_obj.ene_array;
  sig_cal = rxs_obj.sig_array;

  //Calc reference cross sections
  //i=0, j=0
  ene_ref[0][0].push_back(1.10000E+01);
  ene_ref[0][0].push_back(1.20000E+01);
  ene_ref[0][0].push_back(1.30000E+01);
  ene_ref[0][0].push_back(1.40000E+01);
  ene_ref[0][0].push_back(1.50000E+01);
  ene_ref[0][0].push_back(1.60000E+01);
  ene_ref[0][0].push_back(1.70000E+01);
  ene_ref[0][0].push_back(1.80000E+01);
  ene_ref[0][0].push_back(1.90000E+01);
  ene_ref[0][0].push_back(2.00000E+01);
  ene_ref[0][0].push_back(2.10000E+01);
  ene_ref[0][0].push_back(2.20000E+01);
  ene_ref[0][0].push_back(2.30000E+01);

  sig_tmp.push_back(1.00000E+02);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(2.00000E+02);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(3.00000E+02);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(4.00000E+02);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(5.00000E+02);
  i_max = static_cast<int>(sig_tmp.size());
  sig_ref[0][0].resize(i_max);
  for(int i=0; i<i_max; i++)
  {
    Real8 multi = 1.0;
    for(int j=0; j<4; j++)
    {
      sig_ref[0][0][i].push_back(sig_tmp[i]*multi);
      multi *= 0.1;
    }
  }
  sig_tmp.clear();
  
  //i=0, j=1
  ene_ref[0][1].push_back(1.20000E+01);
  ene_ref[0][1].push_back(1.30000E+01);
  ene_ref[0][1].push_back(1.40000E+01);
  ene_ref[0][1].push_back(1.50000E+01);
  ene_ref[0][1].push_back(1.60000E+01);
  ene_ref[0][1].push_back(1.70000E+01);
  ene_ref[0][1].push_back(1.80000E+01);
  ene_ref[0][1].push_back(1.90000E+01);
  ene_ref[0][1].push_back(2.00000E+01);
  ene_ref[0][1].push_back(2.10000E+01);
  ene_ref[0][1].push_back(2.20000E+01);
  ene_ref[0][1].push_back(2.30000E+01);
  ene_ref[0][1].push_back(2.40000E+01);
  ene_ref[0][1].push_back(2.50000E+01);
  ene_ref[0][1].push_back(2.60000E+01);
  ene_ref[0][1].push_back(2.70000E+01);

  sig_tmp.push_back(1.00000E+03);
  sig_tmp.push_back(2.00000E+03);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(3.00000E+03);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(4.00000E+03);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(5.00000E+03);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(6.00000E+03);
  i_max = static_cast<int>(sig_tmp.size());
  sig_ref[0][1].resize(i_max);
  for(int i=0; i<i_max; i++)
  {
    Real8 multi = 1.0;
    for(int j=0; j<4; j++)
    {
      sig_ref[0][1][i].push_back(sig_tmp[i]*multi);
      multi *= 0.1;
    }
  }
  sig_tmp.clear();
  
  //i=1, j=0
  ene_ref[1][0].push_back(1.30000E+01);
  ene_ref[1][0].push_back(1.40000E+01);
  ene_ref[1][0].push_back(1.50000E+01);
  ene_ref[1][0].push_back(1.60000E+01);
  ene_ref[1][0].push_back(1.70000E+01);
  ene_ref[1][0].push_back(1.80000E+01);
  ene_ref[1][0].push_back(1.90000E+01);
  ene_ref[1][0].push_back(2.00000E+01);
  ene_ref[1][0].push_back(2.10000E+01);
  ene_ref[1][0].push_back(2.20000E+01);
  ene_ref[1][0].push_back(2.30000E+01);
  ene_ref[1][0].push_back(2.40000E+01);
  ene_ref[1][0].push_back(2.50000E+01);
  ene_ref[1][0].push_back(2.60000E+01);

  sig_tmp.push_back(1.00000E+04);
  sig_tmp.push_back(2.00000E+04);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(3.00000E+04);
  sig_tmp.push_back(4.00000E+04);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(5.00000E+04);
  sig_tmp.push_back(6.00000E+04);
  sig_tmp.push_back(7.00000E+04);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(8.00000E+04);
  sig_tmp.push_back(9.00000E+04);
  sig_tmp.push_back(1.00000E+05);
  i_max = static_cast<int>(sig_tmp.size());
  sig_ref[1][0].resize(i_max);
  for(int i=0; i<i_max; i++)
  {
    Real8 multi = 1.0;
    for(int j=0; j<4; j++)
    {
      sig_ref[1][0][i].push_back(sig_tmp[i]*multi);
      multi *= 0.1;
    }
  }
  sig_tmp.clear();
  
  //i=1, j=1
  ene_ref[1][1].push_back(2.00000E+01);
  ene_ref[1][1].push_back(2.10000E+01);
  ene_ref[1][1].push_back(2.20000E+01);
  ene_ref[1][1].push_back(2.30000E+01);
  ene_ref[1][1].push_back(2.40000E+01);
  ene_ref[1][1].push_back(2.50000E+01);
  ene_ref[1][1].push_back(2.60000E+01);
  ene_ref[1][1].push_back(2.70000E+01);
  ene_ref[1][1].push_back(2.80000E+01);
  ene_ref[1][1].push_back(2.90000E+01);
  ene_ref[1][1].push_back(3.00000E+01);
  ene_ref[1][1].push_back(3.10000E+01);
  ene_ref[1][1].push_back(3.20000E+01);

  sig_tmp.push_back(1.00000E+05);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(2.00000E+05);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(3.00000E+05);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(4.00000E+05);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(5.00000E+05);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(6.00000E+05);
  sig_tmp.push_back(0.00000E+00);
  sig_tmp.push_back(7.00000E+05);
  i_max = static_cast<int>(sig_tmp.size());
  sig_ref[1][1].resize(i_max);
  for(int i=0; i<i_max; i++)
  {
    Real8 multi = 1.0;
    for(int j=0; j<4; j++)
    {
      sig_ref[1][1][i].push_back(sig_tmp[i]*multi);
      multi *= 0.1;
    }
  }
  sig_tmp.clear();
  
  //i=1, j=2
  ene_ref[1][2].push_back(2.10000E+01);
  ene_ref[1][2].push_back(2.20000E+01);
  ene_ref[1][2].push_back(2.30000E+01);
  ene_ref[1][2].push_back(2.40000E+01);
  ene_ref[1][2].push_back(2.50000E+01);
  ene_ref[1][2].push_back(2.60000E+01);
  ene_ref[1][2].push_back(2.70000E+01);
  ene_ref[1][2].push_back(2.80000E+01);
  ene_ref[1][2].push_back(2.90000E+01);
  ene_ref[1][2].push_back(3.00000E+01);
  ene_ref[1][2].push_back(3.10000E+01);
  ene_ref[1][2].push_back(3.20000E+01);
  ene_ref[1][2].push_back(3.30000E+01);

  sig_tmp.push_back(1.00000E+06);
  sig_tmp.push_back(2.00000E+06);
  sig_tmp.push_back(3.00000E+06);
  sig_tmp.push_back(4.00000E+06);
  sig_tmp.push_back(5.00000E+06);
  sig_tmp.push_back(6.00000E+06);
  sig_tmp.push_back(7.00000E+06);
  sig_tmp.push_back(8.00000E+06);
  sig_tmp.push_back(9.00000E+06);
  sig_tmp.push_back(1.00000E+07);
  sig_tmp.push_back(1.10000E+07);
  sig_tmp.push_back(1.20000E+07);
  sig_tmp.push_back(1.30000E+07);
  i_max = static_cast<int>(sig_tmp.size());
  sig_ref[1][2].resize(i_max);
  for(int i=0; i<i_max; i++)
  {
    Real8 multi = 1.0;
    for(int j=0; j<4; j++)
    {
      sig_ref[1][2][i].push_back(sig_tmp[i]*multi);
      multi *= 0.1;
    }
  }
  sig_tmp.clear();
 
  //Comp calc result and ref data 
  for(int i=0; i<nis; i++)
  {
    for(int j=0; j<ner[i]; j++)
    {
      int k_max = static_cast<int>(ene_cal[i][j].size());
      BOOST_CHECK_EQUAL(ene_cal[i][j].size(), ene_ref[i][j].size());
      BOOST_CHECK_EQUAL(sig_cal[i][j].size(), sig_ref[i][j].size());
      for(int k=0; k<k_max; k++)
      {
        BOOST_CHECK_EQUAL(sig_cal[i][j][k].size(), sig_ref[i][j][k].size());
        BOOST_CHECK_LE(fabs((ene_cal[i][j][k] - ene_ref[i][j][k])/ene_ref[i][j][k]), 1.0e-10);
        
        for(int l=0; l<static_cast<int>(sig_cal[i][j][k].size()); l++)
        {
          if( sig_cal[i][j][k][l] > 1.0e-15 )
          {
            BOOST_CHECK_LE(fabs((sig_cal[i][j][k][l] - sig_ref[i][j][k][l])/sig_cal[i][j][k][l]), 1.0e-10);
          }
          else
          {
            BOOST_CHECK_EQUAL(sig_cal[i][j][k][l], sig_ref[i][j][k][l]);
          }
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(Test_calc_unreso_xs)
{
  Real8  err_val   = 100.0;
  //string file_name = "../EndfUtils/MFxxMTyyyParser/for_test/n_4331_43-Tc-099.dat";
  string file_name = "../EndfUtils/MFxxMTyyyParser/for_test/n_9434_94-Pu-238.dat";
  
  Endf6ParserNoCov       cp_obj;
  ResonanceXSCalculator  res_xs_obj;
  VectorClearer          clr_obj;

  int i_max, j_max;

  vector<vector<vector<Real8> > >          sre_grid;
  vector<vector<vector<Real8> > >          reso_grid;
  vector<vector<vector<vector<Real8> > > > reso_xs;

  //Read endf file
  cp_obj.set_file_name(file_name);

  Endf6Converter    conv_obj;
  NuclearDataObject nucl_data_obj;
  conv_obj.convert_endf_format_to_frendy(cp_obj, nucl_data_obj);
  
  int nis = 1;
  int ner = 2;
  
  sre_grid.resize(nis);
  for(int i=0; i<nis; i++)
  {
    sre_grid[i].resize(ner);
    for(int j=0; j<ner; j++)
    {
      if( j == 0 ) //el:1.0e-5, eh:5.0e2
      {
        sre_grid[i][j].push_back(1.0e-5);
        sre_grid[i][j].push_back(5.0e-5);
        sre_grid[i][j].push_back(1.0e-4);
        sre_grid[i][j].push_back(5.0e-4);
        sre_grid[i][j].push_back(1.0e-3);
        sre_grid[i][j].push_back(5.0e-3);
        sre_grid[i][j].push_back(1.0e-2);
        sre_grid[i][j].push_back(5.0e-2);
        sre_grid[i][j].push_back(1.0e-1);
        sre_grid[i][j].push_back(5.0e-1);
        sre_grid[i][j].push_back(1.0e+0);
        sre_grid[i][j].push_back(5.0e+0);
        sre_grid[i][j].push_back(1.0e+1);
        sre_grid[i][j].push_back(5.0e+1);
        sre_grid[i][j].push_back(1.0e+2);
        sre_grid[i][j].push_back(5.0e+2);
      }
      else //el:5.0e2, eh:6.0e4
      {
        sre_grid[i][j].push_back(5.0e+2);
        sre_grid[i][j].push_back(1.0e+3);
        sre_grid[i][j].push_back(5.0e+3);
        sre_grid[i][j].push_back(1.0e+4);
      }
    }
  }

  //Calc resonance cross section and modified energy grid
  res_xs_obj.set_nucl_data_obj(nucl_data_obj);
  res_xs_obj.set_temp(0.0);
  res_xs_obj.set_err(err_val);
  res_xs_obj.set_resonance_grid(sre_grid);
  reso_grid = res_xs_obj.get_resonance_grid();
  reso_xs   = res_xs_obj.get_resonance_xs();
 
  vector<Real8>          ref_grid = reso_grid[0][1];
  vector<vector<Real8> > ref_xs;
  i_max = static_cast<int>(ref_grid.size());
  ref_xs.resize(i_max);
  for(int i=0; i<i_max; i++)
  {
    j_max = 4;
    ref_xs[i].resize(j_max);
    for(int j=0; j<j_max; j++)
    {
      ref_xs[i][j] = 0.0;
    }
    
    int ii = 0; //nis
    int jj = 1; //ner
    res_xs_obj.set_multi_array_data(ii, jj);
    res_xs_obj.set_multi_array_data_unreso_c(ii, jj);
    res_xs_obj.calc_reso_xs_unreso_c(ii, jj, ref_grid[i], ref_xs[i]);
  }
  
  BOOST_CHECK_EQUAL(reso_grid[0][1].size(), ref_grid.size());
  BOOST_CHECK_EQUAL(reso_xs[0][1].size(),   ref_xs.size());
  i_max = static_cast<int>(ref_grid.size());
  for(int i=0; i<i_max; i++)
  {
    BOOST_CHECK_EQUAL(reso_xs[0][1][i].size(), ref_xs[i].size());
    BOOST_CHECK_LE(fabs((reso_grid[0][1][i] - ref_grid[i])/ref_grid[i]), 1.0e-5);
    j_max = static_cast<int>(ref_xs[i].size());
    for(int j=0; j<j_max; j++)
    {
      BOOST_CHECK_LE(fabs((reso_xs[0][1][i][j] - ref_xs[i][j])/ref_xs[i][j]), 1.0e-5);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(Rmatrix_check01)
{
  string file_name = "../EndfUtils/MFxxMTyyyParser/for_test/n_1725_17-Cl-35.dat";

  Endf6ParserNoCov cp_obj;
  cp_obj.set_file_name(file_name);

  Endf6Converter    conv_obj;
  NuclearDataObject nucl_data_obj;
  conv_obj.convert_endf_format_to_frendy(cp_obj, nucl_data_obj);

  ResonanceEnergyGridLinearizer sre_obj;
  sre_obj.set_nucl_data_obj(nucl_data_obj);

  vector<vector<vector<Real8> > > sre_grid  = sre_obj.get_resonance_grid();

  ResonanceXSCalculator rxs_obj;
  rxs_obj.set_temp(100.0);
  rxs_obj.set_err(1.0e-3);
  rxs_obj.set_err_max(1.0e-2);
  rxs_obj.set_err_int(5.0e-4);

  vector<vector<vector<Real8> > >          reso_grid, reso_q_val;
  vector<vector<vector<vector<Real8> > > > reso_xs;
  vector<vector<vector<Integer> > >        reso_react_list;

  BOOST_CHECK_LE(fabs((rxs_obj.get_temp()    - 100.0 )/100.0 ), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err()     - 1.0e-3)/1.0e-3), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err_max() - 1.0e-2)/1.0e-2), 1.0e-5);
  BOOST_CHECK_LE(fabs((rxs_obj.get_err_int() - 5.0e-4)/5.0e-4), 1.0e-5);

  rxs_obj.set_temp(0.0);
  rxs_obj.set_nucl_data_obj(nucl_data_obj);
  rxs_obj.set_resonance_grid(sre_grid);

  rxs_obj.calc_resonance_xs();
  reso_grid       = rxs_obj.get_resonance_grid();
  reso_xs         = rxs_obj.get_resonance_xs();
  reso_q_val      = rxs_obj.get_resonance_q_val();
  reso_react_list = rxs_obj.get_resonance_react_type_list();


  ifstream fin;
  fin.open("RmatrixLimited/run/for_test/reso_xs_r_matrix_amur.dat");

  vector<Real>          ene_vec;
  vector<vector<Real> > xs_ref;
  vector<Integer>       react_type_ref;
  int ene_no, mt_no;
  fin >> ene_no >> mt_no;
  ene_vec.resize(ene_no);
  xs_ref.resize(ene_no);
  react_type_ref.resize(mt_no);

  for(int i=0; i<mt_no; i++)
  {
    fin >> react_type_ref[i];
  }

  for(int i=0; i<ene_no; i++)
  {
    fin >> ene_vec[i];

    xs_ref[i].resize(mt_no);
    for(int j=0; j<mt_no; j++)
    {
      fin >> xs_ref[i][j];
    }
  }
  fin.close();

  BOOST_CHECK_EQUAL( mt_no, static_cast<int>(reso_react_list[0][0].size()) );
  BOOST_CHECK_EQUAL( mt_no, static_cast<int>(reso_q_val[0][0].size()) );
  for(int i=0; i<mt_no; i++)
  {
    BOOST_CHECK_EQUAL(reso_react_list[0][0][i], react_type_ref[i]);
    if( react_type_ref[i] != 600 )
    {
      BOOST_CHECK_LE(fabs(reso_q_val[0][0][i]), min_value);
    }
    else
    {
      BOOST_CHECK_EQUAL(static_cast<int>(reso_q_val[0][0][i]), 615220);
    }
  }

  vector<vector<Real8> > reso_xs_mod;
  reso_xs_mod.resize(mt_no);
  for(int i=0; i<mt_no; i++)
  {
    int j_max = static_cast<int>(reso_grid[0][0].size());
    reso_xs_mod[i].resize(j_max);
    for(int j=0; j<j_max; j++)
    {
      reso_xs_mod[i][j] = reso_xs[0][0][j][i];
    }
  }

  TabInterpolator ti_obj;
  for(int i=0; i<ene_no; i++)
  {

    for(int j=0; j<mt_no; j++)
    {
      Real xs_cal = 0.0;
      ti_obj.interpolation_1d(int_lin_lin, ene_vec[i], xs_cal, reso_grid[0][0], reso_xs_mod[j]);

      BOOST_CHECK_LE(fabs((xs_cal - xs_ref[i][j])/xs_ref[i][j]), 3.0E-2);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////


BOOST_AUTO_TEST_CASE(check_react_type_list)
{
  vector<vector<vector<Integer> > >        react_type_list;
  vector<vector<vector<Real8> > >          q_array;
  vector<vector<vector<vector<Real8> > > > sig_array;

  int i_max = 2;
  int j_max = 4;
  int k_max = 5;
  int l_max = 10;

  react_type_list.resize(i_max);
  q_array.resize(i_max);
  sig_array.resize(i_max);
  for(int i=0; i<i_max; i++)
  {
    react_type_list[i].resize(j_max);
    q_array[i].resize(j_max);
    sig_array[i].resize(j_max);

    for(int j=0; j<j_max; j++)
    {
      react_type_list[i][j].resize(k_max);
      q_array[i][j].resize(k_max);
      sig_array[i][j].resize(k_max);

      for(int k=0; k<k_max; k++)
      {
        q_array[i][j][k]   = static_cast<Real>(k+1);

 
        sig_array[i][j][k].resize(l_max);
        for(int l=0; l<l_max; l++)
        {
          sig_array[i][j][k][l] = static_cast<Real>(10*(k+1) + l+1);
        }
      }

      if( j==0 )
      {
        react_type_list[i][j][0] =   1;
        react_type_list[i][j][1] =   2;
        react_type_list[i][j][2] =  18;
        react_type_list[i][j][3] = 102;
        react_type_list[i][j][4] =  19;
      }
      else if( j==1 )
      {
        react_type_list[i][j][0] =   1;
        react_type_list[i][j][1] =   2;
        react_type_list[i][j][2] = 102;
        react_type_list[i][j][3] = 600;
        react_type_list[i][j][4] = 800;
      }
      else if( j==2 )
      {
        react_type_list[i][j][0] =   1;
        react_type_list[i][j][1] =  18;
        react_type_list[i][j][2] = 102;
        react_type_list[i][j][3] = 600;
        react_type_list[i][j][4] = 800;
      }
      else if( j==3 )
      {
        react_type_list[i][j][0] =   1;
        react_type_list[i][j][1] = 102;
        react_type_list[i][j][2] =  19;
        react_type_list[i][j][3] = 600;
        react_type_list[i][j][4] = 800;
      }
    }
  }

  ResonanceXSCalculator rxs_obj;
  rxs_obj.react_type_list = react_type_list;
  rxs_obj.q_array         = q_array;
  rxs_obj.sig_array       = sig_array;

  for(int i=0; i<i_max; i++)
  {
    for(int j=0; j<j_max; j++)
    {
      rxs_obj.check_react_type_list(i, j);

      for( int t=0; t<3; t++)
      {
        if( j==0 || j==1 )
        {
          BOOST_CHECK_EQUAL(k_max, static_cast<int>(rxs_obj.react_type_list[i][j].size()));
          BOOST_CHECK_EQUAL(k_max, static_cast<int>(rxs_obj.q_array[i][j].size()));
          BOOST_CHECK_EQUAL(k_max, static_cast<int>(rxs_obj.sig_array[i][j].size()));
        }
        else
        {
          BOOST_CHECK_EQUAL(k_max+1, static_cast<int>(rxs_obj.react_type_list[i][j].size()));
          BOOST_CHECK_EQUAL(k_max+1, static_cast<int>(rxs_obj.q_array[i][j].size()));
          BOOST_CHECK_EQUAL(k_max+1, static_cast<int>(rxs_obj.sig_array[i][j].size()));
        }
  
        for(int k=0; k<k_max; k++)
        {
          BOOST_CHECK_EQUAL(l_max, static_cast<int>(rxs_obj.sig_array[i][j][k].size()));
  
          BOOST_CHECK_EQUAL(rxs_obj.react_type_list[i][j][k], react_type_list[i][j][k]);
          BOOST_CHECK_LE(fabs((rxs_obj.q_array[i][j][k] - q_array[i][j][k]) / q_array[i][j][k]), 1.0E-8);
  
          for(int l=0; l<l_max; l++)
          {
            BOOST_CHECK_LE(fabs((rxs_obj.sig_array[i][j][k][l] - sig_array[i][j][k][l])
                             / sig_array[i][j][k][l]), 1.0E-8);
          }
        }

        if( j==2 )
        {
          BOOST_CHECK_EQUAL(rxs_obj.react_type_list[i][j][k_max], 19);
          BOOST_CHECK_LE(fabs((rxs_obj.q_array[i][j][k_max] - q_array[i][j][1]) / q_array[i][j][1]), 1.0E-8);

          for(int l=0; l<l_max; l++)
          {
            BOOST_CHECK_LE(fabs((rxs_obj.sig_array[i][j][k_max][l] - sig_array[i][j][1][l])
                             / sig_array[i][j][1][l]), 1.0E-8);
          }
        }
        if( j==3 )
        {
          BOOST_CHECK_EQUAL(rxs_obj.react_type_list[i][j][k_max], 18);
          BOOST_CHECK_LE(fabs((rxs_obj.q_array[i][j][k_max] - q_array[i][j][2]) / q_array[i][j][2]), 1.0E-8);

          for(int l=0; l<l_max; l++)
          {
            BOOST_CHECK_LE(fabs((rxs_obj.sig_array[i][j][k_max][l] - sig_array[i][j][2][l])
                             / sig_array[i][j][2][l]), 1.0E-8);
          }
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

