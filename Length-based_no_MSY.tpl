//  *****************************************************************************************************************
// INTEGRATED SPATIALLY-STRUCTURED SIZE-BASED AND SEX-SPECIFIC POPULATION DYNAMICS MODEL

// ASSUMPTIONS:
//		Catches from other fisheries combined to TDGDLF as most catch from this fishery
//		Selectivity is time and space variant based on changes in proportional effort of 6.5 and 7 inch meshes
//              Model is conditioned on catch, which is used to calculate fishing mortality by solving the Baranov
//              catch equation

//note:         Estimation of  Movement Patterns from Conventional and acoustic tagging is conditioned on recaptures.


// TO RUN MCMC
//              Declare variables as sd_report_X
//              ADMB tab /  Run with Args
//                               1) -mcmc number of iterations -mcsave thinning (e.g 10)
//                               2) -mceval (to extract the results from the binary file)
//                               3) read in .txt file into R for chain analysis and posterios
// 
// TO RUN LIKELIHOOD PROFILE
//              Declare variabl as likeprof_X  (stored in a .plt file)
//              ADMB tab /  Run with Args
//                              -lprof

//  *****************************************************************************************************************

GLOBALS_SECTION
  #include <iostream>
  #include <fstream>
  #undef REPORT
  #define REPORT(object) report<<#object "\n" << object << endl;
  #include <time.h>
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;

     

DATA_SECTION
   // model dimensions
  init_int syr
  init_int nyr
  int N_yrs                       //years with catch and cpue observations
  !!N_yrs=1+nyr-syr;      
  init_int nzone
  init_number Yrs_future             //number of years including future catches
  init_int N_inc                     // number of size classes
  int N_sex                         // number of sexes
  !! N_sex=2;
  
  // reproduction
  init_vector Mat(1,N_inc)	      //size at maturity                 
  init_number Sex_ratio                //pup sex ratio
  
  // sizes
  init_number Len_at_age_1
  init_number SD_len_rec
  init_vector Len_bin_mdpt(1,N_inc)   //mid point of size class
  init_vector Len_bin_lbnd(1,N_inc)   //lower bound
  init_vector Len_bin_ubnd(1,N_inc)   //upper bound

  // max age for avoiding plus group
  init_int Max_age
  
  // weight
  init_vector TWT_F(1,N_inc)	      //mean female total weight for  size class l
  init_vector TWT_M(1,N_inc)	      //mean male total weight for  size class l
  
  // natural mortality
  init_vector M(1,N_inc)	      //mean natural mortality  for  size class l
  
  //gear selectivity
  init_matrix SEL(1,N_inc,1,Yrs_future)            //overall selectivity at time                                      
  init_matrix SEL_west(1,N_inc,1,Yrs_future)	      //mean gillnet selectivity  for  size class l per zone and year     
  init_matrix SEL_zn1(1,N_inc,1,Yrs_future)
  init_matrix SEL_zn2(1,N_inc,1,Yrs_future)

  init_vector SEL_6_5(1,N_inc)	      //empirical selectivity for 6.5 inch
  init_vector SEL_7(1,N_inc)	      //empirical selectivity for 7 inch
  
  // steepness 
  init_number steepness
  
  //catch data
  init_ivector iyr(1,N_yrs)                  //the years of catch data   
  init_matrix ct_F(1,Yrs_future,1,nzone)     //female catch (in tonnes) by year and zone
  init_matrix ct_M(1,Yrs_future,1,nzone)     //male catch  (in tonnes)
    
  //size composition data
      //number of observations
         //6.5 inch
  init_matrix N_siz_comp(1,N_yrs,1,nzone)         
  init_matrix N_siz_comp_M(1,N_yrs,1,nzone)
         //7 inch
  init_matrix N_siz_comp_7(1,N_yrs,1,nzone)       
  init_matrix N_siz_comp_M_7(1,N_yrs,1,nzone)
  
      // index of data        
  init_number First_size_obs                               //first size class with data  
  init_number Last_size_obs                                //last size class with data
  init_number First_size_obs_7
  init_number Last_size_obs_7
  init_number First_yr_size_obs                               //first year with data  
  init_number Last_yr_size_obs                                //last year with data


  number g_Dirichlet                                       // number of size classes with observation
  !! g_Dirichlet=1+Last_size_obs-First_size_obs;  
  number g_Dirichlet_7                                       
  !! g_Dirichlet_7=1+Last_size_obs_7-First_size_obs_7;  
  
     //data
          //6.5 inch 
  init_3darray siz_comp_F(1,nzone,1,N_yrs,1,N_inc)      //female  
  init_3darray siz_comp_M(1,nzone,1,N_yrs,1,N_inc)      //male
  
         //7 inch   
  init_3darray siz_comp_F_7(1,nzone,1,N_yrs,1,N_inc)      //female                          
  init_3darray siz_comp_M_7(1,nzone,1,N_yrs,1,N_inc)      //male  
  init_number rho     //weight of size comp likelihood  
  init_int Type_of_Size_like   //0, multinomial; 1:Dirichlet
  
   //CPUE
  init_int do_cpue
  init_int effective_cpue
  init_int syr_cpue            //first year of cpue data
  init_vector CPUE_eff(1,N_yrs)
  init_matrix CPUE(1,N_yrs,1,nzone)
  init_matrix CPUE_SD(1,N_yrs,1,nzone)
 
    
  //Tagging data
 
    //conventional       
  init_number nzones;
  init_number nrec_c;
  init_matrix Rec_dat_c(1,nrec_c,1,3);

  init_number move_min_size;              // Smallest size of moving shark

   
    //acoustic
  init_number nrec_a;
  init_number nATag;
  init_ivector nATag_obs(1,nATag);
  init_ivector TagID_start(1,nATag);
  init_ivector TagID_end(1,nATag);
  init_3darray Rec_dat_a(1,nATag,TagID_start,TagID_end,1,4); //ragged array of recaptures by tagID, (1st dim), recapture events (2nd dim) and release, recapture info (3rd dim)
  
  
 
  //Recruitment error for future projections 
  init_vector Rec_error(1,Yrs_future)               


  //assumed change in catchability
  init_int q_change

  // year monthly returns changed to daily logbooks
  init_int q_daily

  //age and growth data
  init_int nobs
  init_int nobs_M
  init_matrix Growdat(1,nobs,1,2)
  init_matrix Growdat_M(1,nobs_M,1,2)
  init_number rho2      //weighting of growth likelihood
  init_number rho3      //weighting of cpue likelihood
   
  vector L(1,nobs)
  vector AgE(1,nobs)
  !!AgE=column(Growdat,1);
  !!L=column(Growdat,2);
  
  vector L_M(1,nobs_M)
  vector AgE_M(1,nobs_M)
  !!AgE_M=column(Growdat_M,1);
  !!L_M=column(Growdat_M,2);

  //maximum F
  init_number MaxF
 
  // F_init prior
  init_int add_Finit_prior
  init_number mu_Init_F
  init_number SD_Init_F  //in log space
  
  //Phases                                       
  init_int Phase_dummy                           //1
  init_int Phase_lnR_zero
  init_int Phase_lnR_prop_west
  init_int Phase_lnR_prop_zn1
  init_int Phase_lnq                             //5
  init_int Phase_lnq2
  // init_int Phase_lnq_west               
  // init_int Phase_lnq_zn1                 
  // init_int Phase_lnq_zn2
  // init_int Phase_lnq2_west                                       
  // init_int Phase_lnq2_zn1               
  // init_int Phase_lnq2_zn2
  init_int Phase_log_Qdaily
  // init_int Phase_log_Qdaily_west               
  // init_int Phase_log_Qdaily_zn1                        
  // init_int Phase_log_Qdaily_zn2   
  
  init_int Phase_ln_Init_F               
  init_int Phase_log_tau                 
  init_int Phase_k                                //10                 
  init_int Phase_lnLinf                           
  init_int Phase_k_M                   
  init_int Phase_lnLinf_M
  init_int Phase_sd_growth
  init_int Phase_log_p11                            //15      
  init_int Phase_log_p22       
  init_int Phase_log_p21      
  init_int Phase_log_p33

  init_int Do_move

  init_int eof;
  !! if(eof != 999){cout<<"Error reading data\n"; exit(1);}

  //checking data inputs
  //!!cout<<"L_M\n"<<L_M<<endl;exit(1);

  //Declare loop indices
  int from //STM
  int to
  int t    //years or number of tags
  int i    // sizes
  int j    // sizes
  int a    // ages
  int s    // sex
  int y    //days at liberty
  int n    //number of observations per tagID
  int z    //zones
  

  // Newton's method for estimating F
  number r
  number FStart;
  number DerivInt;

  //MCMC
  int niter;	//number of ierations for MCMC
  !! niter=0;

  
PARAMETER_SECTION
  init_number Dummy(Phase_dummy);                   //dummy for controling phases
  
  //Recruitment pars
  init_bounded_number lnR_zero(4,15,Phase_lnR_zero);                // initial recruitment for init. popn. at equil. (thousands)
  init_bounded_number R_prop_west(0.05,0.9,Phase_lnR_prop_west);           // proportion of initial recruitment in West
  init_bounded_number R_prop_zn1(0.05,0.9,Phase_lnR_prop_zn1);             // proportion of initial recruitment in Zone1
  
  //Q pars
  init_bounded_number lnq(-15,1,Phase_lnq);                        // monthly catchability
  init_bounded_number lnq2(-15,1,Phase_lnq2);                      // monthly catchability second period
  init_bounded_number lnq_daily(-15,1,Phase_log_Qdaily);          // daily catchability
 
  // init_bounded_number lnq_west(-15,1,Phase_lnq_west);               // catchability first period by zone
  // init_bounded_number lnq_zn1(-15,1,Phase_lnq_zn1);                 
  // init_bounded_number lnq_zn2(-15,1,Phase_lnq_zn2);                
  
  // init_bounded_number lnq2_west(-15,1,Phase_lnq2_west);             // catchability second period by zone
  // init_bounded_number lnq2_zn1(-15,1,Phase_lnq2_zn1);               
  // init_bounded_number lnq2_zn2(-15,1,Phase_lnq2_zn2);               

  // init_bounded_number lnq_daily_west(-15,1,Phase_log_Qdaily_west);   // daily catchability by zone
  // init_bounded_number lnq_daily_zn1(-15,1,Phase_log_Qdaily_zn1);   
  // init_bounded_number lnq_daily_zn2(-15,1,Phase_log_Qdaily_zn2);   

  //Other pars
  init_bounded_number ln_Init_F(-20,-0.7,Phase_ln_Init_F);                 // F before start of time series
  init_bounded_number ln_tau(-20,0,Phase_log_tau);	               // observation error

  //Growth pars
  init_bounded_number k(0.1,0.6,Phase_k);                       // growth pars females
  init_bounded_number lnLinf(4,6,Phase_lnLinf);
  init_bounded_number k_M(0.1,0.6,Phase_k_M);                  // growth pars males
  init_bounded_number lnLinf_M(4,6,Phase_lnLinf_M);
  init_bounded_number sd_growth(0.1,50,Phase_sd_growth);  

     //movement  pars
  init_bounded_number p11(0,1,Phase_log_p11);          //area 1
  init_bounded_number p22(0,1,Phase_log_p22);          //area 2
  init_bounded_number p21(0,1,Phase_log_p21);
  init_bounded_number p33(0,1,Phase_log_p33);          //area 3


  objective_function_value f

  //likeprof_number lp_ln_Init_F;    //likelihood profile for Fo
  

   //Declare objects for back transforming estimable pars
   number R_zero
   vector R_zero_prop(1,nzone)
   matrix q(1,N_yrs,1,nzone)
   number Init_F
   number tau
   number Linf
   number Linf_M
   
   //Declare objects used in Prelim Calcs
   matrix M_sex(1,N_sex,1,N_inc)
   matrix TwT(1,N_sex,1,N_inc)
   3darray TC(1,N_sex,1,Yrs_future,1,nzone) 
   3darray N_comp(1,N_sex,1,N_yrs,1,nzone)
   3darray SEL_zone_yr(1,nzone,1,N_inc,1,Yrs_future)        
   4darray siz_comp(1,N_sex,1,nzone,1,N_yrs,1,N_inc)
   3darray N_comp_7(1,N_sex,1,N_yrs,1,nzone)                
   4darray siz_comp_7(1,N_sex,1,nzone,1,N_yrs,1,N_inc)      
   
   //Declare objects used in Age_and_growth()
   vector Lpred(1,nobs)
   vector Lpred_M(1,nobs_M)

   //Declare objects used in Move_rate()
   matrix Mov_mat(1,nzones,1,nzones);
   matrix Move_c(1,nzones,1,nzones);
   matrix Move_a(1,nzones,1,nzones);    
    
   vector Days_libery_c(1,nrec_c);
   vector Rels_c(1,nrec_c);
   vector Recs_c(1,nrec_c);

   number Tag_NLL_c;
   number Tag_NLL_a;
   number Tag_like;
   number fpen_tag;
  
   //Declare objects used in Size_transition_matrix()
   3darray STM(1,N_sex,1,N_inc,1,N_inc) // size transition matrix

   // Declare objects used in Prob_size_new_recruits()
   vector Prob_at_len_rec(1,N_inc)
   number x1
   number x2
   number Tot_prob

   //Declare objects used in Calc_fished_and_unfished_spawn_biom_per_recruit()
   number F
     
   //Declare objects used in Calc_len_distn_for_init_popn()
   matrix Surv_Len(1,N_sex,1,N_inc)
   3darray Surv_L(1,N_sex,0,1,1,N_inc)    // first element, 0=unfished, 1=fished
   matrix N_Len(1,N_sex,1,N_inc)   
   3darray N_L(1,N_sex,0,1,1,N_inc)
   3darray Age_Comp(1,N_sex,1,Max_age,1,Yrs_future)
   matrix PerRec_Exp_Age_Comp(1,N_sex,1,Max_age)  //age_length composition of per recruits

   //Declare objects used in Calc_init_spawn_biom_per_recruit()
   matrix SB_LEN_REC(0,1,1,N_inc)
   vector Spawn_biom_per_recruit(0,1)

   //Declare objects used in Calc_stock_rec_params_and_init_recruitment()
   //number Virgin_Spawn_Biom
   sdreport_number Virgin_Spawn_Biom
   number aSRR
   number bSRR
   //vector Virgin_Recruitment(1,nzone)
   sdreport_vector Virgin_Recruitment(1,nzone)
   vector Init_Recruitment(1,nzone)
   number R_prop_pen
   number Init_rec_pen

   //Declare objects used in Virgin_conditions()
   3darray Virgin_Total_biom_sex(1,N_sex,1,1,1,nzone) 
   //3darray Virgin_Vul_biom_sex(1,N_sex,1,1,1,nzone)   //don't calculate because time changing selectivity affects this
   //matrix Virgin_Total_biom(1,1,1,nzone) 
   //matrix Virgin_Vul_biom(1,1,1,nzone)
   sdreport_matrix Virgin_Total_biom(1,1,1,nzone) 
   //sdreport_matrix Virgin_Vul_biom(1,1,1,nzone)
   //vector Virgin_Vul_biom_temp(1,1)
   4darray Virgin_Num_in_len_class(1,N_sex,1,1,1,nzone,1,N_inc)
   

   //Declare objects used in Population_dynamics()
   3darray Total_biom_sex(1,N_sex,1,Yrs_future,1,nzone) 
   3darray Vul_biom_sex(1,N_sex,1,Yrs_future,1,nzone)
   //matrix Total_biom(1,Yrs_future,1,nzone) 
   //matrix Vul_biom(1,Yrs_future,1,nzone)
   sdreport_matrix Total_biom(1,Yrs_future,1,nzone)
   sdreport_matrix Total_biom_rel(1,Yrs_future,1,nzone) 
   sdreport_matrix Vul_biom(1,Yrs_future,1,nzone)
   //sdreport_matrix Vul_biom_rel(1,Yrs_future,1,nzone)
  // matrix depletion(1,1,1,nzone)
  // matrix depletion_vul(1,1,1,nzone)
  // matrix depletion_spawn(1,1,1,nzone)
   sdreport_matrix depletion(1,1,1,nzone)
   //sdreport_matrix depletion_vul(1,1,1,nzone)
   sdreport_matrix depletion_spawn(1,1,1,nzone)
   vector Vul_biom_temp(1,N_yrs)
   number Rec_temp
   //matrix Fem_spawn_biom(1,Yrs_future,1,nzone)
   sdreport_matrix Fem_spawn_biom(1,Yrs_future,1,nzone)
   sdreport_matrix Fem_spawn_biom_rel(1,Yrs_future,1,nzone) 
   4darray Num_in_len_class(1,N_sex,1,Yrs_future+1,1,nzone,1,N_inc)
   4darray Num_in_len_class_temp(1,N_sex,1,Yrs_future+1,1,nzone,1,N_inc)   
   3darray Est_catch(1,N_sex,1,Yrs_future,1,nzone)
   3darray Annual_F_sex(1,N_sex,1,Yrs_future,1,nzone)
   //matrix Annual_F_female(1,Yrs_future,1,nzone)
   //matrix Annual_F_male(1,Yrs_future,1,nzone)
   sdreport_matrix Annual_F_female(1,Yrs_future,1,nzone)
   sdreport_matrix Annual_F_male(1,Yrs_future,1,nzone)
   //matrix Annual_rec(2,Yrs_future+1,1,nzone)
   sdreport_matrix Annual_rec(2,Yrs_future+1,1,nzone)
   4darray Surv_in_len_class(1,N_sex,1,Yrs_future,1,nzone,1,N_inc)   
   number F_in_len_class
   number F_in_len_class_6_5                              
   number F_in_len_class_7                                        
   number Z_in_len_class
   number Z_in_len_class_6_5                             
   number Z_in_len_class_7  
   4darray Annual_Z(1,N_sex,1,Yrs_future,1,nzone,1,N_inc)
   4darray Est_catch_num_in_len_class(1,N_sex,1,Yrs_future,1,nzone,1,N_inc)    // catch in numbers in each zone, year, and len class
   4darray Est_catch_num_in_len_class_6_5(1,N_sex,1,Yrs_future,1,nzone,1,N_inc)  
   4darray Est_catch_num_in_len_class_7(1,N_sex,1,Yrs_future,1,nzone,1,N_inc)   
   4darray Est_catch_in_len_class(1,N_sex,1,Yrs_future,1,nzone,1,N_inc)        // catch in weight in each zone, year, and len class                             
   vector TC_out(1,Yrs_future)                                                 //store all females and males catches
   vector Est_TC_out(1,Yrs_future)
   matrix Mov_mat_annual(1,nzones,1,nzones);
   
   //Declare objects used in Calc_annual_F_values_using_Newtons_method()
   number CatchPen 
   number tempF
   number tempZ
   number Deriv
   number NewF
   number FMort1
   number FMort2
   number tempCatch
   number tempEstCatch
   number tempEstCatchNum
   number ExpCatch1
   number ExpCatch2

   //Declare dummies
   vector test(1,N_sex)
   vector test2(1,N_sex)
   number test3      

   //Declare objects used in Calc_marginal_distn_for_len_in_each_yr()
   3darray Exp_catch_num_6_5(1,N_sex,1,N_yrs,1,nzone)            
   4darray Est_Prop_at_len_6_5(1,N_sex,1,N_yrs,1,nzone,1,N_inc)   
   3darray Exp_catch_num_7(1,N_sex,1,N_yrs,1,nzone)            
   4darray Est_Prop_at_len_7(1,N_sex,1,N_yrs,1,nzone,1,N_inc)
   4darray Est_Prop_at_len_6_5_out(1,N_sex,1,nzone,1,N_yrs,1,N_inc)     
   4darray Est_Prop_at_len_7_out(1,N_sex,1,nzone,1,N_yrs,1,N_inc)


   //Declare objects used in Calc_NLL_for_Age_and_growth()
   vector epsln_fem(1,nobs)    
   vector epsln_male(1,nobs_M)
   number Growth_NLL

   //Declare objects used in Calc_NLL_for_CPUE()
   matrix tau1(1,N_yrs,1,nzone)
   matrix Est_CPUE(1,N_yrs,1,nzone)
   vector Est_CPUE_eff(1,N_yrs)
   
   matrix CPUE_out(syr_cpue,N_yrs,1,nzone)
   matrix CPUE_SD_out(syr_cpue,N_yrs,1,nzone)
   
   matrix Est_CPUE_out(syr_cpue,N_yrs,1,nzone)
   sdreport_matrix log_Est_CPUE_out(syr_cpue,N_yrs,1,nzone)
   matrix tau1_out(syr_cpue,N_yrs,1,nzone)
   
   //vector Sq_res_CPUE(1,N_yrs)
   number CPUE_NLL
   //number Sum_sq_CPUE
   //number St_dev_CPUE
   //vector epsilon(1,N_yrs)	//residuals
  
   //Declare objects used in Calc_NLL_for_len_comp_data()
   vector Len_comp_NLL_sex(1,N_sex)
   vector Len_comp_NLL_sex_7(1,N_sex)                          
   number Len_comp_NLL

      //Dirichlet objects
   number sum_Pi_ln_pi_yi
   number Effective_n
   number LL_temp
         //6.5 inch   
   vector npi(First_size_obs,Last_size_obs)
   vector ln_gamma_npi(First_size_obs,Last_size_obs)
   vector lnyi(First_size_obs,Last_size_obs)
         //7 inch   
   vector npi_7(First_size_obs_7,Last_size_obs_7)
   vector ln_gamma_npi_7(First_size_obs_7,Last_size_obs_7)
   vector lnyi_7(First_size_obs_7,Last_size_obs_7)
   
   //Declare objects used in Calc_obj_function()
   number Penalties


PRELIMINARY_CALCS_SECTION
                
   //fill female (dimension 1) and male (dimension 2) objects
     //natural mortality      
   for(int s=1;  s<=N_sex;s++)
   {
    for(int i=1; i<=N_inc;i++)		
     {
       M_sex(s,i)=M(i);   
     }
   }
   
     //mean weight at size
   TwT(1)=TWT_F;
   TwT(2)=TWT_M;

    //Selectivity      
     // time varying 
   if(nzone==1)
   {
     SEL_zone_yr(1)=SEL;
   }else
   {
     // time and space varying 
     SEL_zone_yr(1)=SEL_west;                 
     SEL_zone_yr(2)=SEL_zn1;
     SEL_zone_yr(3)=SEL_zn2;
   }
   
     
     //catch
   TC(1)=ct_F;
   TC(2)=ct_M;

     //catch size composition
        //number of observations
          //6.5
   N_comp(1)=N_siz_comp;
   N_comp(2)=N_siz_comp_M;

          //7
   N_comp_7(1)=N_siz_comp_7;                      
   N_comp_7(2)=N_siz_comp_M_7;   

        //observations                         
            //6.5
   siz_comp(1)=siz_comp_F;
   siz_comp(2)=siz_comp_M;

          //7
   siz_comp_7(1)=siz_comp_F_7;                      
   siz_comp_7(2)=siz_comp_M_7;   

  
   //Conventional tagging data
   if(Do_move==1)
   {
     Days_libery_c=column(Rec_dat_c,1);
     Rels_c=column(Rec_dat_c,2);
     Recs_c=column(Rec_dat_c,3);
   }

  
   //required for using Newton's method to estimate fishing mortality
   FStart = 0.1;
   DerivInt = 0.000001;

   //output relevant cpues
   for (t=syr_cpue; t<=N_yrs;t++)
   {
     for(z=1; z<=nzone;z++)
     {
      CPUE_out(t,z)=CPUE(t,z);
      CPUE_SD_out(t,z)=CPUE_SD(t,z);
     }
   }
   

PROCEDURE_SECTION
   Age_and_growth();
   if(Do_move==1) Move_rate();
   Size_transition_matrix();
   Prob_size_new_recruits();
   Calc_fished_and_unfished_spawn_biom_per_recruit();
   Calc_stock_rec_params_and_init_recruitment();
   Virgin_conditions();
   Population_dynamics();
   Calc_marginal_distn_for_len_in_each_yr();
   Calc_NLL_for_Age_and_growth();
   Calc_NLL_for_CPUE();
   Calc_NLL_for_len_comp_data();
   Calc_obj_function();
   if(mceval_phase()) write_mcmc_report();	//write outputs using the sampled posteriors 


FUNCTION Age_and_growth
   Linf=mfexp(lnLinf);
   Linf_M=mfexp(lnLinf_M);

      //2 par vonB
   Lpred=Len_at_age_1+(Linf-Len_at_age_1)*(1-mfexp(-k*AgE));  
   Lpred_M=Len_at_age_1+(Linf_M-Len_at_age_1)*(1-mfexp(-k_M*AgE_M));
  
      //traditional vonB
   //Lpred=Linf*(1-mfexp(-k*(AgE-t0)));   
   //Lpred_M=Linf_M*(1-mfexp(-k_M*(AgE_M-t0_M)));
   
   
FUNCTION Move_rate
  Tag_NLL_c.initialize();                         
  Tag_NLL_a.initialize();  
  fpen_tag.initialize();
  
   //1. fill in daily Movement matrix
  Mov_mat(1,1)=p11;
  Mov_mat(1,2)=1-p11;
  Mov_mat(1,3)=0;

  Mov_mat(2,1)=p21;
  Mov_mat(2,2)=p22;
  Mov_mat(2,3)=1-(p22+p21);

  Mov_mat(3,1)=0;
  Mov_mat(3,2)=1-p33;
  Mov_mat(3,3)=p33;

   //penalty to keep all row elements positive
  Mov_mat(1,2)=posfun(Mov_mat(1,2),0.000001,fpen_tag);
  fpen_tag+=fpen_tag;
  fpen_tag += (1.0 - sum(Mov_mat(1))) * (1.0 - sum(Mov_mat(1)));

  Mov_mat(2,3)=posfun(Mov_mat(2,3),0.000001,fpen_tag);
  fpen_tag+=fpen_tag;
  fpen_tag += (1.0 - sum(Mov_mat(2))) * (1.0 - sum(Mov_mat(2)));
 
  Mov_mat(3,2)=posfun(Mov_mat(3,2),0.000001,fpen_tag);
  fpen_tag+=fpen_tag;
  fpen_tag += (1.0 - sum(Mov_mat(3))) * (1.0 - sum(Mov_mat(3)));

   //penalty to keep all diagonal elements <=1
  Mov_mat(1,1)=1-posfun((1-Mov_mat(1,1)),0.001,fpen_tag);
  fpen_tag+=fpen_tag;
 
  Mov_mat(2,2)=1-posfun((1-Mov_mat(2,2)),0.0001,fpen_tag);
  fpen_tag+=fpen_tag;
 
  Mov_mat(3,3)=1-posfun((1-Mov_mat(3,3)),0.0001,fpen_tag);
  fpen_tag+=fpen_tag;

  //Normalise matrix                              
  for(i=1; i<=nzones;i++) Mov_mat(i)=Mov_mat(i)/sum(Mov_mat(i));


 
  //2. Predict location of recaptured shark

    //2.1 Conventional tags 
  for(t=1; t<=nrec_c;t++)
  {
    Move_c=Mov_mat;
    if(Days_libery_c(t)>1)
    {    
     for(y=2; y<=Days_libery_c(t);y++) Move_c=Move_c*Mov_mat;      
    }
    int from=value(Rels_c(t));          //C++ syntax: need value() to convet to int
    int to=value(Recs_c(t));
    dvariable Pred_Prob_c=Move_c(from,to);     
    Tag_NLL_c += - log(Pred_Prob_c);  // negative log likelihood of observations
    
  }  //end t convetional
      
    //2.2 Acoustic tags    
  for(t=1; t<=nATag;t++)                          // loop over each tag ID
  {
    //extract release and recapture info for each TagID
   dvector Days_libery_a=column(Rec_dat_a(t),1);
   dvector Rels_a=column(Rec_dat_a(t),2);
   dvector Recs_a=column(Rec_dat_a(t),3);
    
   Tag_like.initialize();
   
    //loop over observations for each tag ID
   for(n=TagID_start(t); n<=TagID_end(t);n++)                 
   {
      //get position and time per observation
     int N_days_t_a=Days_libery_a(n);
     int REL_a=Rels_a(n);
     int REC_a=Recs_a(n);
      
      //re set matrix
     Move_a=Mov_mat;

      //calculate position after N_days
     if(N_days_t_a>1)
     {
       for(y=2; y<=N_days_t_a;y++) Move_a=Move_a*Mov_mat;
     }
     dvariable Pred_Prob_a=Move_a(REL_a,REC_a);
     Tag_like += - log(Pred_Prob_a);      //negative log likelihood of observations     
   }

     //Calculate likelihood of each tag ID
    Tag_NLL_a+=(Tag_like/nATag_obs(t));   //weight by number of observations per tagID
    // Tag_NLL_a+=Tag_like; 
    
  }  //end t acoustic


     

FUNCTION Size_transition_matrix
   for(from=1;from<=N_inc;from++)
   { 
     for(to=1;to<=N_inc;to++)
     {
       STM(1,from,to) =  mfexp(-pow(Len_bin_mdpt(to)-(Linf*(1-mfexp(-k))+(Len_bin_mdpt(from)*mfexp(-k))),2)/(2*pow(sd_growth,2)));
       STM(2,from,to) =  mfexp(-pow(Len_bin_mdpt(to)-(Linf_M*(1-mfexp(-k_M))+(Len_bin_mdpt(from)*mfexp(-k_M))),2)/(2*pow(sd_growth,2)));

       //truncate to avoid size shrinking
     //  if(to<from) STM(1,from,to)=0;  
     //  if(to<from) STM(2,from,to)=0;
     }
     //normalise
     STM(1,from)=STM(1,from)/sum(STM(1,from));
     STM(2,from)=STM(2,from)/sum(STM(2,from));
    }
    
     //transpose to have columns sum to 1
     STM(1)=trans(STM(1));  
     STM(2)=trans(STM(2));

       
FUNCTION Prob_size_new_recruits
  // size distribution of new recruits, assumed to come from a N(Mean.neonate, SD.neonate)            
  for (i=1; i<=N_inc;i++)
  {
    x1 = (Len_bin_lbnd(i) - Len_at_age_1) / SD_len_rec; // z score, i.e. (x-mu)/sd
    x2 = (Len_bin_ubnd(i) - Len_at_age_1) / SD_len_rec;
    Prob_at_len_rec(i) = cumd_norm(x2) - cumd_norm(x1);
  }
    //normalise probs
    Tot_prob=sum(Prob_at_len_rec);
    Prob_at_len_rec/=Tot_prob;
    
    

FUNCTION Calc_fished_and_unfished_spawn_biom_per_recruit
  //note: Initialise model accounting for impact of fishing before start of time series
  //       the model considers that individuals recruit at age 1.
  //       Loop through a large number of age classes to build up an expected length distribution

  //steps: Calculate the length composition for
  //         1) initial unfished population at equilibrium  to allow  S_zero calculation
  //         2) repeat for population with estimated initial (low level) of fishing.
            
  N_L.initialize();                  
                                     
  //Step 1. unfished spawning biomass per recruit
  F = 0;
  t = 0;
  Calc_len_distn_for_init_popn();
  Calc_init_spawn_biom_per_recruit();

  //Step 2. initial fished spawning biomass per recruit
  Init_F=mfexp(ln_Init_F);
  //lp_ln_Init_F=ln_Init_F;    //likeprofile
  F = Init_F;
  t = 1;
  Calc_len_distn_for_init_popn();
  Calc_init_spawn_biom_per_recruit();
  

FUNCTION Calc_len_distn_for_init_popn
  Surv_Len.initialize();                                  
  N_Len.initialize();

  //sex loop
  for (s=1; s<=N_sex;s++)
  {
     //age loop
  for (a=1; a<=Max_age;a++)              
  { 
    //1. calculate number of animals surviving per size class
    for (i=1; i<=N_inc;i++)               //size loop          
    {
       if (a == 1)
       {
        Surv_Len(s,i) = Prob_at_len_rec(i)*Sex_ratio;
       }else if(a==Max_age)     //plus group
       {
        Surv_Len(s,i) = N_Len(s,i)*mfexp(-(M_sex(s,i) + (SEL(i,1) * F)))/(1-mfexp(-(M_sex(s,i) + (SEL(i,1) * F))));
       }else
       {
        Surv_Len(s,i) = N_Len(s,i)*mfexp(-(M_sex(s,i) + (SEL(i,1) * F)));
       }
       
       // sum survival in each length class over all ages
       Surv_L(s,t, i) += Surv_Len(s,i);         //accumulate survivors
     } // i
     
    //2. calculate numbers at size considering growth of each size class and survival
    //note: grow animals, to determine numbers surviving in each length
    //           class after mortality and growth
    
    N_Len = 0;   // set numbers to zero before multiplying by growth transition matrix
    
    if (a == 1) // recruits - don't grow as recruits already have sizes
    {  
      N_Len(s) = Surv_Len(s);
    }else      // grow other size classes
    {
     for (i=1; i<=N_inc;i++)                             
     {  //from 
       for (j=1; j<=N_inc;j++)
       { //to
         N_Len(s,j) += (Surv_Len(s,i) * STM(s,j, i));    //columns sums to 1
       } // j     
     } // i
    }
     
    // sum numbers in each length class over all ages and matrix of age-length composition
    for (j=1; j<=N_inc;j++)
    {
       N_L(s,t,j) += N_Len(s,j);
    } // j
     
    // keep track of age composition for the initial exploited popn. at equil. 
    if (t == 1)
    {
       PerRec_Exp_Age_Comp(s,a) = sum(N_Len(s));
    }
        
   } // a   
  }//s
   

FUNCTION Calc_init_spawn_biom_per_recruit
  // Calculate the spawning biomass for the initial population at equilibrium
  Spawn_biom_per_recruit(t)=0;

  for (i=1; i<=N_inc;i++)
  {
   // calculate the spawning biomass in each length class, in kg                      
   SB_LEN_REC(t, i) = N_L(1,t, i)*TwT(1,i)*Mat(i); 
 
    // sum for the ages
    Spawn_biom_per_recruit(t) += SB_LEN_REC(t, i);                                    
  }
    

FUNCTION Calc_stock_rec_params_and_init_recruitment

   R_prop_pen.initialize();
   
   //back transform estimable pars
   R_zero=mfexp(lnR_zero);                //Total recruitment
   if(nzone==1)
   {
     R_zero_prop(1)=1;
   }else
   {
     R_zero_prop(1)=R_prop_west;   //proportion of total recruitment per zone
     R_zero_prop(2)=R_prop_zn1;
     R_zero_prop(3)=1-(R_zero_prop(1)+R_zero_prop(2));
     R_zero_prop(3)=posfun(R_zero_prop(3),0.000001,R_prop_pen);   //penalty to keep sum of Rec prop =1
     R_prop_pen += (1.0 - sum(R_zero_prop)) * (1.0 - sum(R_zero_prop));
     R_zero_prop=R_zero_prop/sum(R_zero_prop);  //normalise           
   }
   
             
  // calculate the virgin recruitment
  Virgin_Spawn_Biom = R_zero * Spawn_biom_per_recruit(0);

  
  // Bev-Holt stock recruitment parameters
  aSRR = (Virgin_Spawn_Biom / R_zero) * ((1 - steepness) / (4 * steepness));
  bSRR = (steepness - 0.2) / (0.8 * steepness * R_zero);

  //test if recovering the recruits (test3 should be equal to R_zero)
   //test3=Virgin_Spawn_Biom/(aSRR+bSRR*Virgin_Spawn_Biom);  
   //cout<<"R_zero "<<R_zero<<"  test3 "<<test3<<endl;exit(1);
  
  // calculate the initial recruitment in each zone
  Init_rec_pen=0;
  for(z=1; z<=nzone;z++)
  {
   //Virgin conditions
    Virgin_Recruitment(z) =R_zero_prop(z)*R_zero;

   //Finit conditions
    Init_Recruitment(z) =R_zero_prop(z)*( (Spawn_biom_per_recruit(1) - aSRR) / (bSRR * Spawn_biom_per_recruit(1)));
    Init_Recruitment(z)=1+posfun(Init_Recruitment(z)-1,0.001,Init_rec_pen);   //penalty to have at least 1 (in 1000s) individuals recruiting
  }

  
FUNCTION Virgin_conditions
  Virgin_Total_biom_sex.initialize();
  //Virgin_Vul_biom_sex.initialize();
  Virgin_Num_in_len_class.initialize();

  for (s=1; s<=N_sex;s++)                        //loop over each sex
  {
   //Iinitial numbers in population @ equilbrium, given initial recruitment and estimated numbers per recruit
   for(z=1; z<=nzone;z++)                        //loop over zone
   {
     for (i=1; i<=N_inc;i++)                     //loop over size bins
     {
       Virgin_Num_in_len_class(s,1,z,i) = N_L(s,1,i) * Virgin_Recruitment(z);  
     } //i  
   }//z
   
   for (t=1; t<=1;t++)                    //time
   {
     for(z=1; z<=nzone;z++)                        //loop over zone
     {
      //1. calculate biomasses @ time 1
      for (i=1; i<=N_inc;i++)                     //loop over size bins
      {
          //total biomass (already in tonnes)
       Virgin_Total_biom_sex(s,t,z)+=Virgin_Num_in_len_class(s,t,z,i)*TwT(s,i);
       
          //vulnerable biomass (already in tonnes)
       //Virgin_Vul_biom_sex(s,t,z)+=Virgin_Num_in_len_class(s,t,z,i)*TwT(s,i)*SEL_zone_yr(z,i,t) ;
               
      } //i    
    } // z
   } //t
  }//s

  //Sum female and male total and vulnerable biomass
  for(z=1; z<=nzone;z++)
  {
   Virgin_Total_biom.colfill(z,column(Virgin_Total_biom_sex(1),z)+column(Virgin_Total_biom_sex(2),z));
   //Virgin_Vul_biom.colfill(z,column(Virgin_Vul_biom_sex(1),z)+column(Virgin_Vul_biom_sex(2),z));
  }

   
FUNCTION Population_dynamics
  Total_biom_sex.initialize();
  Vul_biom_sex.initialize();
  Fem_spawn_biom.initialize();
  Num_in_len_class.initialize();
  Num_in_len_class_temp.initialize();
  Est_catch.initialize();

  for (s=1; s<=N_sex;s++)                         //loop over each sex
  {
   //Iinitial numbers in population @ equilbrium, given initial recruitment and estimated numbers per recruit
   for(z=1; z<=nzone;z++)                        //loop over zone
   {
     for (i=1; i<=N_inc;i++)                     //loop over size bins
     {
       Num_in_len_class(s,1,z,i) = N_L(s,1,i) * Init_Recruitment(z);  
     } //i  
   }//z
               
   //  Age_Comp(s,1,1) = sum(Num_in_len_class(s,1));
 
   // time series calculations 
   for (t=1; t<=Yrs_future;t++)                    //loop over time
   {
     for(z=1; z<=nzone;z++)                        //loop over zone
     {

      //1. calculate biomasses @ time t
      for (i=1; i<=N_inc;i++)                     //loop over size bins
      {
          //total biomass (already in tonnes)
       Total_biom_sex(s,t,z)+=Num_in_len_class(s,t,z,i)*TwT(s,i);
       
          //vulnerable biomass (already in tonnes)
       Vul_biom_sex(s,t,z)+=Num_in_len_class(s,t,z,i)*TwT(s,i)*SEL_zone_yr(z,i,t) ;
      

          //spawning biomass (already in tonnes)
       if(s==1) Fem_spawn_biom(t,z)+=Num_in_len_class(s,t,z,i)*TwT(s,i)*Mat(i);          
      } //i
     
    
      //2. calculate  fully-selected fishing mortality @ time t
      //note:  Fishing mortality applied @ beginning of year.
      //       Init_F used to set up initial equilibrium conditions
    
        //2.1 get fully selected annual F
      if (TC(s,t,z) >0)
      {
         Calc_annual_F_values_using_Newtons_method();  //numerically solve Baranov
      }
      else
      {
         Annual_F_sex(s,t,z)=0;    // set F to 0 if no catch 
      }
   
 
        //2.2 apply M and F by size to calculate population numbers    
      for (i=1; i<=N_inc;i++)
      {
        // calculate the length-specific fishing mortality accounting for selectivity
        F_in_len_class = SEL_zone_yr(z,i,t) * Annual_F_sex(s,t,z);
        F_in_len_class_6_5 = SEL_6_5(i) * Annual_F_sex(s,t,z);            
        F_in_len_class_7 = SEL_7(i) * Annual_F_sex(s,t,z); 
    
        // calculate the total mortality
        Z_in_len_class = F_in_len_class + M_sex(s,i);
        Z_in_len_class_6_5 = F_in_len_class_6_5 + M_sex(s,i);            
        Z_in_len_class_7 = F_in_len_class_7 + M_sex(s,i);   

        // calculate predicted numbers per size class that will survive to the beginning of next year
        Surv_in_len_class(s,t,z,i) = Num_in_len_class(s,t,z,i) * mfexp(-Z_in_len_class);
    
        // calculate predicted catch in numbers per size class, taken thru current year
        Est_catch_num_in_len_class(s,t,z,i)=Num_in_len_class(s,t,z,i)*(F_in_len_class/Z_in_len_class)*(1-mfexp(-Z_in_len_class));
        Est_catch_num_in_len_class_6_5(s,t,z,i)=Num_in_len_class(s,t,z,i)*(F_in_len_class_6_5/Z_in_len_class_6_5)*(1-mfexp(-Z_in_len_class_6_5)); 
        Est_catch_num_in_len_class_7(s,t,z,i)=Num_in_len_class(s,t,z,i)*(F_in_len_class_7/Z_in_len_class_7)*(1-mfexp(-Z_in_len_class_7));         
    
        // calculate predicted catch in weight (tonnes)
        Est_catch_in_len_class(s,t,z,i) = Est_catch_num_in_len_class(s,t,z,i) * TwT(s,i);     
        
        // sum predicted catch (tonnes) over all length classes
        Est_catch(s,t,z) += Est_catch_in_len_class(s,t,z,i);
               
        // grow  animals and place them into the next time step        
        for (j=1; j<=N_inc;j++)
        {
         Num_in_len_class(s,t+1,z,j) += Surv_in_len_class(s,t,z,i) * STM(s,j,i);      
        } // j
        
        // check if creating fish when growing 
        //if(s==2 & t==10 & z==2) cout<<"Surv "<<sum(Surv_in_len_class(s,t,z))<<"  Num_in_len_class "<<sum(Num_in_len_class(s,t+1,z))<<endl;
       
        //store total mortality per sex, zone, time and size bin
        Annual_Z(s,t,z,i) = Z_in_len_class;

       } // i              
    } // z


    //3. Calculate the expected recruitment by zone, for following year
    Rec_temp=0;
    if(s==1)
    {
       Rec_temp=sum(Fem_spawn_biom(t)) / (aSRR + (bSRR * sum(Fem_spawn_biom(t))));
       for(z=1; z<=nzone;z++)
       {
          Annual_rec(t+1,z) = Rec_temp*R_zero_prop(z);
        
        //add recruitment variability to projections into the future
         if(t>N_yrs)
         {
          Annual_rec(t+1,z)=Annual_rec(t+1,z)*Rec_error(t);
          if(Annual_rec(t+1,z)> Virgin_Recruitment(z)) Annual_rec(t+1,z)=Virgin_Recruitment(z); 
         }


       }        
     }


    //4. Add expected recruits to the population, for next year
    for(z=1; z<=nzone;z++) for (i=1; i<=N_inc;i++)  Num_in_len_class(s,t+1,z,i) += Prob_at_len_rec(i) * Annual_rec(t+1,z)*Sex_ratio;
     
   
    //5. Move animals among zones
    if(Do_move==1)
    {
           //derive annual movement matrix from daily movement matrix                          
      Mov_mat_annual=Mov_mat;
      for(y=2; y<=365;y++) Mov_mat_annual=Mov_mat_annual*Mov_mat;
                                                                                           
   
       // first calculate numbers after movement for all size classes
      for(from=1;from<=nzone;from++) for(to=1;to<=nzone;to++) Num_in_len_class_temp(s,t,to)+=Num_in_len_class(s,t,from)* Mov_mat_annual(from,to) ;   
        
      // second replace numbers for size assumed to be able to move
      for(to=1;to<=nzone;to++) Num_in_len_class(s,t,to)(move_min_size,N_inc)=Num_in_len_class_temp(s,t,to)(move_min_size,N_inc);
    }
    
   } //t
  }//s


  //6. Sum female and male total and vulnerable biomass
  for(z=1; z<=nzone;z++)
  {
   Total_biom.colfill(z,column(Total_biom_sex(1),z)+column(Total_biom_sex(2),z));
   Vul_biom.colfill(z,column(Vul_biom_sex(1),z)+column(Vul_biom_sex(2),z));
  }
  
    //relative biomasses      
  for (t=1; t<=Yrs_future;t++) for(z=1; z<=nzone;z++)
  {
     Total_biom_rel(t,z)=Total_biom(t,z)/Virgin_Total_biom(1,z);
     //Vul_biom_rel(t,z)=Vul_biom(t,z)/Virgin_Vul_biom(1,z);
     Fem_spawn_biom_rel(t,z)=Fem_spawn_biom(t,z)/Virgin_Spawn_Biom*R_zero_prop(z);
  }
  
  //depletion      
  for(z=1; z<=nzone;z++)
  {
     depletion(1,z)=Total_biom(N_yrs,z)/Virgin_Total_biom(1,z);
     //depletion_vul(1,z)=Vul_biom(N_yrs,z)/Virgin_Vul_biom(1,z);
     depletion_spawn(1,z)=Fem_spawn_biom(N_yrs,z)/Virgin_Spawn_Biom*R_zero_prop(z);
  }
  

  // get predicted and observed catches to check they are equal (as model is conditioned on catch)
  for (t=1; t<=Yrs_future;t++)    
  {
   TC_out(t)=sum(TC(1,t))+sum(TC(2,t));
   Est_TC_out(t)=sum(Est_catch(1,t))+sum(Est_catch(2,t));
  }


FUNCTION Calc_annual_F_values_using_Newtons_method
    CatchPen=0.0;
    FMort1 = FStart;
    tempCatch = TC(s,t,z);   //observed catch

    // penalize objective function if F becomes too large
    // note: this sets the observed catch as posfun of the catch for a maximum allowable F value
    //       minus the actual observed catch 
    tempEstCatch = 0;
    for (i=1; i<=N_inc;i++)
    {
       tempF = MaxF* SEL_zone_yr(z,i,t);
       tempZ = tempF + M_sex(s,i);
       tempEstCatchNum = (tempF/tempZ)*(1-mfexp(-tempZ))*Num_in_len_class(s,t,z,i); 
       tempEstCatch += tempEstCatchNum * TwT(s,i); 
    }       
       tempCatch = tempEstCatch - posfun(tempEstCatch - tempCatch,0.01,CatchPen);
    
    for (r=1; r<=4;r++)
    {
      ExpCatch1 = 0;
      ExpCatch2 = 0;

      // calculate the survival and estimated catches,
      // given the value of F
      for (i=1; i<=N_inc;i++)
      {

       F_in_len_class = SEL_zone_yr(z,i,t) * FMort1;
       Z_in_len_class = F_in_len_class + M_sex(s,i);
       Est_catch_num_in_len_class(s,t,z,i)=Num_in_len_class(s,t,z,i)*(F_in_len_class/Z_in_len_class)*(1-mfexp(-Z_in_len_class));
       ExpCatch1 += Est_catch_num_in_len_class(s,t,z,i) * TwT(s,i);
      } // i
      
      FMort2 = FMort1 + DerivInt;

      for (i=1; i<=N_inc;i++)
      {
        F_in_len_class = SEL_zone_yr(z,i,t) * FMort2;
        Z_in_len_class = F_in_len_class + M_sex(s,i);
        Est_catch_num_in_len_class(s,t,z,i)=Num_in_len_class(s,t,z,i)*(F_in_len_class/Z_in_len_class)*(1-mfexp(-Z_in_len_class));
        ExpCatch2 += Est_catch_num_in_len_class(s,t,z,i) * TwT(s,i);
      } // i

      // calculate the derivative
      Deriv = (ExpCatch2 - ExpCatch1) / DerivInt;

      // calculate the new value for F
      NewF = FMort1 - ((ExpCatch1 - tempCatch) / Deriv);

      FMort1 = NewF;
      
    }//r
    
    Annual_F_sex(s,t,z) = NewF;


FUNCTION Calc_marginal_distn_for_len_in_each_yr
  Exp_catch_num_6_5.initialize();                                 
  Exp_catch_num_7.initialize();    
  Est_Prop_at_len_6_5.initialize();      
  Est_Prop_at_len_7.initialize();    

  for(s=1; s<=N_sex;s++)                                          //loop over sexes
  {
     // calculate total expected catch numbers
     for (t=1; t<=N_yrs;t++)                                      //loop over years
     {
      for(z=1; z<=nzone;z++)                                      //loop over zones
      {
       for (i=1; i<=N_inc;i++)                                    //loop over size bins
       {
         Exp_catch_num_6_5(s,t,z) += Est_catch_num_in_len_class_6_5(s,t,z,i);
          Exp_catch_num_7(s,t,z) += Est_catch_num_in_len_class_7(s,t,z,i);
       }  //i
      } //z
     } //t

     // calculate expected catch proportions
     for (t=1; t<=N_yrs;t++)                                      //loop over years
     {
      for(z=1; z<=nzone;z++)                                      //loop over zones
      {
       for (i=1; i<=N_inc;i++)                                    //loop over size bins
       {
         Est_Prop_at_len_6_5(s,t,z,i) = Est_catch_num_in_len_class_6_5(s,t,z,i) / Exp_catch_num_6_5(s,t,z);
         Est_Prop_at_len_7(s,t,z,i) = Est_catch_num_in_len_class_7(s,t,z,i) / Exp_catch_num_7(s,t,z);

           //output estimated prop at size in nice format   
          Est_Prop_at_len_6_5_out(s,z,t,i)=Est_Prop_at_len_6_5(s,t,z,i);
          Est_Prop_at_len_7_out(s,z,t,i)=Est_Prop_at_len_7(s,t,z,i);
          
       } //i
      } //z
     } //t
  }//s

 
FUNCTION Calc_NLL_for_Age_and_growth
   Growth_NLL = 0;
   
      //Robust likelihood
   //note: manipulating sd_growth, which is used in STM
   //Growth_NLL=robust_regression(L,Lpred,sd_growth);   
   //Growth_NLL+=robust_regression(L_M,Lpred_M,sd_growth);

      //Normal likelihood
    epsln_fem=L-Lpred;
    Growth_NLL=nobs*log(sd_growth)+sum(square(epsln_fem))/(2.*sd_growth*sd_growth);
    epsln_male=L_M-Lpred_M;
    Growth_NLL+=nobs_M*log(sd_growth)+sum(square(epsln_male))/(2.*sd_growth*sd_growth);  



FUNCTION Calc_NLL_for_CPUE
  Est_CPUE.initialize();
  Est_CPUE_eff.initialize();
  tau1.initialize();
  CPUE_NLL = 0;
   
  if(do_cpue!=0)  
  {
   tau=exp(ln_tau);               //model error
     //backtransform q
      //note: use different catchabilities according q type, period and zones 
    if(effective_cpue==1)
    {
     for(t=1; t<=N_yrs;t++)
     {
      if(q_change>0)
      {
        if(t<=q_change)
        {
           q(t)=mfexp(lnq);               
        }else
        {
           q(t)=mfexp(lnq2);           
        }
      }else
      {
       q(t)=mfexp(lnq);
      }
     } //t
    }else
    {
     //one zone
     if(nzone==1)
     {
      for(t=1; t<=N_yrs;t++)
      {
        if(q_change>0)
        {
          if(t<=q_change)
          { 
             q(t,1)=mfexp(lnq);               
          }else
          {
            if(t<q_daily)
            {
              q(t,1)=mfexp(lnq2);
              if (Phase_lnq2<0) q(t,1)=mfexp(lnq);
            }else
            {
              q(t,1)=mfexp(lnq_daily);
            }
          }
        }else
        {
          if(t<q_daily)
          {
            q(t,1)=mfexp(lnq);
          }else
          {
            q(t,1)=mfexp(lnq_daily);
         }
       }
     } //t  
    }else
     { 
       for(t=1; t<=N_yrs;t++)
       {         
         if(q_change>0)
         {
           if(t<=q_change)
           {
             q(t,1)=mfexp(lnq);
             q(t,2)=mfexp(lnq);
             q(t,3)=mfexp(lnq);
             // q(t,1)=mfexp(lnq_west);
             // q(t,2)=mfexp(lnq_zn1);
             // q(t,3)=mfexp(lnq_zn2);
           }else
           {
             if(t<q_daily)
             {
               q(t,1)=mfexp(lnq2);
               q(t,2)=mfexp(lnq2);
               q(t,3)=mfexp(lnq2);
               if (Phase_lnq2<0)
               {
                q(t,1)=mfexp(lnq);
                q(t,2)=mfexp(lnq);
                q(t,3)=mfexp(lnq);
               }

               //q(t,1)=mfexp(lnq2_west);
               //q(t,2)=mfexp(lnq2_zn1);
               //q(t,3)=mfexp(lnq2_zn2);
             }else
             {
               q(t,1)=mfexp(lnq_daily);
               q(t,2)=mfexp(lnq_daily);
               q(t,3)=mfexp(lnq_daily);
               //q(t,1)=mfexp(lnq_daily_west);
               //q(t,2)=mfexp(lnq_daily_zn1);
               //q(t,3)=mfexp(lnq_daily_zn2);
             }
           }
         }else
         {
           if(t<q_daily)
           {
             q(t,1)=mfexp(lnq);
             q(t,2)=mfexp(lnq);
             q(t,3)=mfexp(lnq);
             //q(t,1)=mfexp(lnq_west);
             //q(t,2)=mfexp(lnq_zn1);
             //q(t,3)=mfexp(lnq_zn2);
           }else
           {
             q(t,1)=mfexp(lnq_daily);
             q(t,2)=mfexp(lnq_daily);
             q(t,3)=mfexp(lnq_daily);
             //q(t,1)=mfexp(lnq_daily_west);
             //q(t,2)=mfexp(lnq_daily_zn1);
             //q(t,3)=mfexp(lnq_daily_zn2);
           }
         }
       } //t    
     }
   }
   

   //calculate predicted cpue
   if(effective_cpue==1)
   {
    for(t=1; t<=N_yrs;t++)
    {
      Vul_biom_temp(t)=sum(Vul_biom(t));  //combine zones
      Est_CPUE_eff(t)=q(t,1)*Vul_biom_temp(t);     //same q across zones so use first zone
    }       
   }else
   {
     for (t=1; t<=N_yrs;t++) for(z=1; z<=nzone;z++) Est_CPUE(t,z)=q(t,z)*Vul_biom(t,z);     
   }
  
    // cout<<"R_zero_prop= "<<R_zero_prop<<endl;    
    // cout<<"bSRR\n"<<bSRR<<endl;
    //  cout<<"Virgin_Spawn_Biom\n"<<Virgin_Spawn_Biom<<endl;
    //  cout<<"Fem_spawn_biom\n"<<Fem_spawn_biom<<endl;
    //   cout<<"Virgin_Recruitment\n"<<Virgin_Recruitment<<endl;
    //   cout<<"Init_Recruitment\n"<<Init_Recruitment<<endl;  
    //   cout<<"Annual_rec\n"<<Annual_rec<<endl;
    //   cout<<"Q\n"<<q<<endl;
    //   cout<<"Vul_biom\n"<<Vul_biom<<endl;
     //  cout<<"CPUE\n"<<CPUE<<endl;
     //   cout<<"Est_CPUE\n"<<Est_CPUE<<endl;
     //    cout<<"TC_out\n"<<TC_out<<endl;
     //     cout<<"Est_TC_out\n"<<Est_TC_out<<endl;
     //     cout<<"depletion\n"<<depletion<<endl;
     //     cout<<"Total_biom\n"<<Total_biom<<endl;
     //exit(1);   
   
  //calculate neg log like
  
    //--Francis 2011 T2.5 negloglike
   if(effective_cpue==1)
   {
    for (t=syr_cpue; t<=N_yrs;t++)
    {
      if(CPUE_eff(t)>0) CPUE_NLL+=log(1)+0.5*(square((log(CPUE_eff(t)/Est_CPUE_eff(t)))/1));    
    }  //t
   }else
   {
    for (t=syr_cpue; t<=N_yrs;t++)
    {
     for(z=1; z<=nzone;z++)
     {
      if(CPUE(t,z)>0)   // only use CPUE with data (missing cpue is set to <0 in R)
      {
       //combine model error and observation error from cpue standardisation    
       tau1(t,z)=sqrt(square(tau)+square(CPUE_SD(t,z)));        //total SD is (sum of variances)^1/2
       CPUE_NLL+=log(tau1(t,z))+0.5*(square((log(CPUE(t,z)/Est_CPUE(t,z)))/tau1(t,z)));       
      }
     }  //z
    }  //t
   }
  

    //--Simple sums of squares
  // epsilon=log(CPUE)-log(Est_CPUE);
  //  int n=size_count(epsilon);
  // CPUE_NLL+=0.5*n*log(tau)+(0.5/square(tau))*norm2(epsilon);


     //--Alex's approach
  // Sq_res_CPUE.initialize();
  // Sum_sq_CPUE = 0;
  // for (t=1; t<=N_yrs;t++)
  // {      
  //    // Calculate LL associated with CPUE, kg / fishing day
  //   Est_CPUE(t) = q * Vul_biom(t);
       
  //   // Calculate the squared residuals
  //   Sq_res_CPUE(t) = (CPUE(t) - Est_CPUE(t)) * (CPUE(t) - Est_CPUE(t));
  //   Sum_sq_CPUE += Sq_res_CPUE(t);
  // } //t
  // St_dev_CPUE = sqrt(Sum_sq_CPUE / N_yrs);
  // CPUE_NLL = (N_yrs / 2.) * (log(2 * M_PI) + (2 * log(St_dev_CPUE) + 1));


      //--Simon's
  // SigmaSrd = ss/ncnt;
  //   Catch_Like(ir) += Weighting for that likelihood  * (ss/(2*SigmaSrd))+0.5*number of obs*log(2.0*pi+SigmaSrd);

  }//end do_cpue!=0

     //output relevant cpues
   for (t=syr_cpue; t<=N_yrs;t++)
   {
     for(z=1; z<=nzone;z++)
     {
      Est_CPUE_out(t,z)=Est_CPUE(t,z);
      log_Est_CPUE_out(t,z)=log(Est_CPUE(t,z));
      tau1_out(t,z)=tau1(t,z);
     }
   }
  
   

FUNCTION Calc_NLL_for_len_comp_data
  Len_comp_NLL = 0;                               
  
  for (s=1; s<=N_sex;s++)
  {
    Len_comp_NLL_sex(s)=0;
    Len_comp_NLL_sex_7(s)=0;
    
    for (t=First_yr_size_obs; t<=Last_yr_size_obs;t++)
    {         
     for(z=1; z<=nzone;z++)
     {
       //6.5 inch
       sum_Pi_ln_pi_yi.initialize();
       Effective_n.initialize();
       npi.initialize();
       ln_gamma_npi.initialize();
       lnyi.initialize();
       LL_temp.initialize();    
       if (N_comp(s,t,z) > 0)       //do calculation if there are observations for that year
       {    
         //---Francis 2011 T2.6 (multinomial)
         if(Type_of_Size_like==0)
         {
           for (i=First_size_obs; i<=Last_size_obs;i++)
           {
             if (siz_comp(s,z,t,i)> 0)Len_comp_NLL_sex(s) += -1.0*N_comp(s,t,z) * siz_comp(s,z,t,i) * log(Est_Prop_at_len_6_5(s,t,z,i));
           }           
         }
         
         //---Dirichlet (Schnute and Haigh (2007))
           //note: This distribution is 'self-weighting'so no need to weight the length compo data outside the model.
           // This distribution does not allow zeros so, as recommended by Francis (2014), 
           // distribution tails were compressed (i.e. remove small and large sizes for which 
           //  no data had been recorded in any year), and a small constant added (0.0001) to
           // observed proportions for any remaining length categories with zero values in a particular year.
         if(Type_of_Size_like==1)
         {
           for (i=First_size_obs; i<=Last_size_obs;i++)
           {
            sum_Pi_ln_pi_yi+=Est_Prop_at_len_6_5(s,t,z,i)* log((Est_Prop_at_len_6_5(s,t,z,i)/siz_comp(s,z,t,i)));            
           }
           Effective_n=((g_Dirichlet-1)/2)*(1/sum_Pi_ln_pi_yi);
           
           for (i=First_size_obs; i<=Last_size_obs;i++)
           {
             npi(i)=Effective_n*Est_Prop_at_len_6_5(s,t,z,i);
             ln_gamma_npi(i)=gammln(npi(i));
             lnyi(i)=log(siz_comp(s,z,t,i));
             LL_temp+=ln_gamma_npi(i)-(npi(i)*lnyi(i));         
           }
           Len_comp_NLL_sex(s)+=LL_temp-gammln(Effective_n);
         }
      }   // end   if (N_comp(s,t,z) > 0)

       //7 inch
       sum_Pi_ln_pi_yi.initialize();
       Effective_n.initialize();
       npi_7.initialize();
       ln_gamma_npi_7.initialize();
       lnyi_7.initialize();
       LL_temp.initialize();      
       if (N_comp_7(s,t,z) > 0)
       {    
         //---Francis 2011 T2.6 (multinomial)
         if(Type_of_Size_like==0)
         {
           for (i=First_size_obs_7; i<=Last_size_obs_7;i++)
           {
            if (siz_comp_7(s,z,t,i)> 0) Len_comp_NLL_sex_7(s) += -1.0*N_comp_7(s,t,z) * siz_comp_7(s,z,t,i) * log(Est_Prop_at_len_7(s,t,z,i));
           }             
         }

         //---Dirichlet (Schnute and Haigh (2007)) 
         if(Type_of_Size_like==1)
         {
           for (i=First_size_obs_7; i<=Last_size_obs_7;i++)
           {
            sum_Pi_ln_pi_yi+=Est_Prop_at_len_7(s,t,z,i)* log((Est_Prop_at_len_7(s,t,z,i)/siz_comp_7(s,z,t,i)));
           }
           Effective_n=((g_Dirichlet_7-1)/2)*(1/sum_Pi_ln_pi_yi);
           for (i=First_size_obs_7; i<=Last_size_obs_7;i++)
           {
             npi_7(i)=Effective_n*Est_Prop_at_len_7(s,t,z,i);
             ln_gamma_npi_7(i)=gammln(npi_7(i));
             lnyi_7(i)=log(siz_comp_7(s,z,t,i));
             LL_temp+=ln_gamma_npi_7(i)-(npi_7(i)*lnyi_7(i));
           }   
           Len_comp_NLL_sex_7(s)+=LL_temp-gammln(Effective_n);
         }           
      }    // end   if (N_comp_7(s,t,z) > 0)
      
     }  //z     
    } //t
  }//s
  
  //Combine male and female from 6.5 and 7 inch likelihoods
  Len_comp_NLL=Len_comp_NLL_sex(1)+Len_comp_NLL_sex(2)+Len_comp_NLL_sex_7(1)+Len_comp_NLL_sex_7(2);

  //extract female and male fishing mortality
  if(nzone==1)
  {
     Annual_F_female=Annual_F_sex(1);
     Annual_F_male=Annual_F_sex(2);
  }else
  {
     Annual_F_female=1;
     Annual_F_male=1;
  }


FUNCTION Calc_obj_function
   Penalties=0;
   //if(!last_phase())
   //{
     Penalties=1e1*CatchPen+1e5*R_prop_pen+1e5*fpen_tag+1e2*Init_rec_pen;     
   //}
   
   f = rho2 * Growth_NLL + rho3 * CPUE_NLL + rho * Len_comp_NLL + Penalties + Tag_NLL_a + Tag_NLL_c;
   cout<<"Growth_NLL="<<Growth_NLL<<" CPUE_NLL="<<CPUE_NLL<<" Len_comp_NLL="<<Len_comp_NLL<<" Penalties="<<Penalties<<" Tag_NLL_a="<<Tag_NLL_a<<" Tag_NLL_c="<<Tag_NLL_c<<endl;
   // cout<<"R_zero_prop= "<<R_zero_prop<<"  R_prop_pen= "<<R_prop_pen<<endl;
   // exit(1);
   
   //add F_init prior
   if(add_Finit_prior==1)f+=dlnorm(Init_F,log(mu_Init_F),SD_Init_F);
     

FUNCTION write_mcmc_report
   if(niter==0)
   {
     ofstream ofs("PARS.mcmc");
      ofs<<"R_zero\t tau\t k\t Linf\t k_M\t Linf_M\t sd_growth"<<endl;
     ofstream ofs1("Fem_spawn_biom.mcmc");
      ofs1<<"Fem_spawn_biom"<<endl;  
     ofstream ofs2("Vul_biom.mcmc");
      ofs2<<"Vul_biom"<<endl;
     ofstream ofs3("Total_biom.mcmc");
      ofs3<<"Total_biom"<<endl;
     ofstream ofs4("Fem_spawn_biom_rel.mcmc");
      ofs4<<"Fem_spawn_biom_rel"<<endl;
     ofstream ofs6("Total_biom_rel.mcmc");
      ofs6<<"Total_biom_rel"<<endl;      
      ofstream ofs7("Annual_rec.mcmc");
      ofs7<<"Annual_rec"<<endl;  
     ofstream ofs8("Init_Recruitment.mcmc");
      ofs8<<"Init_Recruitment"<<endl;
     ofstream ofs9("Annual_F_M.mcmc");    
      ofs9<<"Annual_F_M"<<endl;
     ofstream ofs10("Annual_F_F.mcmc");    
      ofs10<<"Annual_F_F"<<endl;
     ofstream ofs11("Virgin_Spawn_Biom.mcmc");    
      ofs11<<"Virgin_Spawn_Biom"<<endl;
     ofstream ofs13("Virgin_Total_biom.mcmc");    
      ofs13<<"Virgin_Total_biom"<<endl;
    }		//a loop to avoid writing all the staff again
	
    ofstream ofs("PARS.mcmc",ios::app);	// append stuff (iop=input output)
     ofs<<R_zero<<"\t"<<tau<<"\t"<<k<<"\t"<<Linf<<"\t"<<k_M<<"\t"<<Linf_M<<"\t"<<sd_growth<<endl;  
    ofstream ofs1("Fem_spawn_biom.mcmc",ios::app);	
     ofs1<<Fem_spawn_biom<<endl;
    ofstream ofs2("Vul_biom.mcmc",ios::app);	
     ofs2<<Vul_biom<<endl;
    ofstream ofs3("Total_biom.mcmc",ios::app);	
     ofs3<<Total_biom<<endl;
    ofstream ofs4("Fem_spawn_biom_rel.mcmc",ios::app);	
     ofs4<<Fem_spawn_biom_rel<<endl;
    ofstream ofs6("Total_biom_rel.mcmc",ios::app);	
     ofs6<<Total_biom_rel<<endl;     
    ofstream ofs7("Annual_rec.mcmc",ios::app);	
     ofs7<<Annual_rec<<endl;
    ofstream ofs8("Init_Recruitment.mcmc",ios::app);	
     ofs8<<Init_Recruitment<<endl;
    ofstream ofs9("Annual_F_M.mcmc",ios::app);	
     ofs9<<Annual_F_sex(2)<<endl;
    ofstream ofs10("Annual_F_F.mcmc",ios::app);	
     ofs10<<Annual_F_sex(1)<<endl;
    ofstream ofs11("Virgin_Spawn_Biom.mcmc",ios::app);	
     ofs11<<Virgin_Spawn_Biom<<endl;
    ofstream ofs13("Virgin_Total_biom.mcmc",ios::app);	
     ofs13<<Virgin_Total_biom<<endl;

    niter++;	// to make niter go to next


REPORT_SECTION
   REPORT(iyr);                  //years in data time series
   REPORT(TC_out);               //observed catch
   REPORT(Est_TC_out);           //predicted catch
  
   
   if(effective_cpue==0) REPORT(CPUE_out);                 //observed cpue
   if(effective_cpue==0) REPORT(Est_CPUE_out);             //predicted cpue
   if(effective_cpue==0) REPORT(CPUE_SD_out);              //SD from observed CPUE
   if(effective_cpue==0) REPORT(tau1_out);                 //total error for cpue (model plus cpue standardisation)
   if(effective_cpue==1) REPORT(CPUE_eff);                 //effective cpue
   if(effective_cpue==1) REPORT(Est_CPUE_eff);             //predicted effective cpue 

   REPORT(Prob_at_len_rec);      //distribution of recruits
   REPORT(Init_Recruitment);     //initial recruitment
   REPORT(Annual_rec);           //annual recruitment
   REPORT(Spawn_biom_per_recruit);
   REPORT(Virgin_Spawn_Biom);
   REPORT(PerRec_Exp_Age_Comp);
   REPORT(N_L(1));
   REPORT(N_L(2));

   REPORT(Fem_spawn_biom);       //female mature biomass
   REPORT(Total_biom);           //total biomass
   REPORT(Vul_biom);             //vulnerable biomass
   REPORT(Annual_F_sex);         //annual F by sex
   
   REPORT(AgE);                  //growth outputs
   REPORT(L);        
   REPORT(Lpred);        
   REPORT(sd_growth);        
   REPORT(AgE_M);        
   REPORT(L_M);
   REPORT(Lpred_M);
   REPORT(aSRR);
   REPORT(bSRR);
   REPORT(STM(1));               //female size transition matrix
   REPORT(STM(2));               //male size transition matrix
   
   REPORT(Len_bin_mdpt);         //mid point size bin 
   for(z=1; z<=nzone;z++)
   {
    REPORT(siz_comp_F(z));    //observed female size composition 6.5 inch
    REPORT(siz_comp_M(z));    
    REPORT(siz_comp_F_7(z));  //observed female size composition 7 inch
    REPORT(siz_comp_M_7(z));
   }          
   REPORT(N_siz_comp);           //size composition sample size
   REPORT(N_siz_comp_M);
   REPORT(N_siz_comp_7);           
   REPORT(N_siz_comp_M_7);

   for(z=1; z<=nzone;z++)
   {
    REPORT(Est_Prop_at_len_6_5_out(1,z));    //predicted female size composition 6.5 inch
    REPORT(Est_Prop_at_len_6_5_out(2,z));    //male
    REPORT(Est_Prop_at_len_7_out(1,z));      //predicted female size composition 7 inch
    REPORT(Est_Prop_at_len_7_out(2,z));      //male
   }
 
   REPORT(Num_in_len_class(1));  //females numbers at size for each year
   REPORT(Num_in_len_class(2));  //males

   REPORT(Mov_mat);              //daily movement
   REPORT(Mov_mat_annual);       //annual movement


 

TOP_OF_MAIN_SECTION
	time(&start);
	// arrmblsize = 50000000;
	// gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	// gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	// gradient_structure::set_MAX_NVAR_OFFSET(5000);
	// gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
 arrmblsize=8000000;
 gradient_structure::set_NUM_DEPENDENT_VARIABLES(1800);
FINAL_SECTION
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<endl<<endl<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"*******************************************"<<endl;
 
