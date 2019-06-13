//  *****************************************************************************************************************
// SIMPFENDORFER ET AL.'S 2000 AGE AND SEX STRUCTURED POPULATION DYNAMICS MODEL
//
//  NOTES:  
//		Observation error only
//		Model is conditioned commercial catch, used to estimate fishing mortality
//  		Created by  Matias Braccini on 2015-09-23
//  		Copyright (c) 2016. All rights reserved.
//              Fitted to nominal cpue only
// ASSUMPTIONS:
//		Catches from other fisheries combined to TDGDLF as most catch from this fishery
//		Selectivty is known for commercial fishery (used 6.5 inch mesh).
//  *****************************************************************************************************************

GLOBALS_SECTION
  #include <iostream>
  #include <fstream>
  #undef REPORT
  #define REPORT(object) report<<#object "\n" << object << endl;

DATA_SECTION
   // model dimensions
  init_int syr
  init_int nyr	
  init_int nzone
  init_int future_yrs            //end of future projections
  init_int nage_F
  init_int nage_M

  // reproduction
  init_vector fec(1,nage_F)	      //mean number of pups at age class a
  init_number sx_ratio                //pup sex ratio
  init_number breed                   //breeding frequency
  init_number age_mat                 // age at  maturity
  !!fec=fec*breed;                    //annual fecundity
  
  // natural mortality
  init_vector M(1,nage_F)	      //natural mortality  for  age class a
  vector M_M(1,nage_F)	              //define vector for male and female M
  vector M_F(1,nage_F)	              

  !!for(int i=1; i<=nage_F;i++)		//filling male and female M
  !!{
  !!	M_F(i)=M(i);
  !!	M_M(i)=M(i);
  !!}
  

  //gear selectivity
  init_vector SEL_F(1,nage_F)	      //mean gillnet selectivity for females in age class a
  init_vector SEL_M(1,nage_F)	      //mean gillnet selectivity for males in  age class a
  
  //age classes
  init_vector age_F(1,nage_F)	      
  init_vector age_M(1,nage_M)	      

  // weight
  init_vector TWT_F(1,nage_F)	      //mean female total weight for age class a
  init_vector TWT_M(1,nage_F)	      //mean male total weight for  age class a

  //catch data
  init_ivector iyr(syr,nyr)               //years of catch data
  init_vector ct_F(syr,future_yrs)        //female catch
  init_vector ct_M(syr,future_yrs)        //male catch
  vector TC_out(syr,future_yrs)           //total catch (in kg)
  !!TC_out=(ct_F+ct_M);
  
  //cpue
  init_int effective_cpue
  init_int syr_cpue                   //first year of cpue data
  init_vector CPUE(syr,nyr)

  //Recruitment error for future projections      //not used in this scenario
  init_vector Rec_error(syr,future_yrs)            

  //assumed change in catchability
  init_int q_change
  init_int q_daily   //daily cpue
  
  //Phases                                       
  init_int Phase_dummy
  init_int Phase_RSTAR                 
  init_int Phase_z
  init_int Phase_Q1  
  init_int Phase_Q2
  init_int Phase_Qdaily
  init_int Phase_Fo                 
  

  //Number of parameters to be estimated
  init_int n_par

    //cpue variance
  init_int Do_var      
  init_number VR1
  init_number VR2

  
  //error trap to ensure you read all the data.
  init_int eof;
  !! if(eof != 999){cout<<"Error reading data"<<endl; exit(1);};
  

PARAMETER_SECTION
  init_number Dummy(Phase_dummy);                   //dummy for controling phases
  init_number RSTAR(Phase_RSTAR);
  init_number Z(Phase_z);
  init_bounded_number Q1(1e-6,1e2,Phase_Q1);  //first period catchability
  init_bounded_number Q2(1e-6,1e2,Phase_Q2);  //second period catchability
  init_bounded_number Qdaily(1e-6,1e2,Phase_Qdaily)    // daily cpue q
  init_number Fo(Phase_Fo);                   // pre 1975 fishing mortality
  
 
  objective_function_value f
  //likeprof_number lp_ln_Init_Fo;    //likelihood profile for Fo

  //Declare objects used in virgin_cond()
  number rstar;
  vector No_M(1,nage_F);
  vector No_F(1,nage_F);
  vector egg_o(1,nage_F);
  number Egg_o;
  number minF;

  //Declare objects used in init_model()
  number z;
  number Fpen;
  number Zpen;
  number z_max;
  vector Per_recruit(1,nage_F);
  vector egg_per_rec(1,nage_F);
  number Egg_per_rec;
  number Ro;
  number a;
  number b;
  vector N_1975_M(1,nage_F);
  vector N_1975_F(1,nage_F);
  vector Bo_M(1,nage_F);
  vector Bo_F(1,nage_F);
  number B_1975;
  vector bo_mature(1,nage_F);
  

  //Declare objects used in dynamic_model()
  number Bpen;
  matrix N_F(syr,future_yrs,1,nage_F);
  matrix N_M(syr,future_yrs,1,nage_F);
  vector Annual_F(syr,future_yrs);
  matrix c_F(syr,future_yrs,1,nage_F);
  matrix c_M(syr,future_yrs,1,nage_F);
  matrix c_F_wgt(syr,future_yrs,1,nage_F);
  matrix c_M_wgt(syr,future_yrs,1,nage_F);
  vector Est_TC_out(syr,future_yrs);
  vector EGGS(syr,future_yrs);
  vector Annual_rec(syr,future_yrs);
  matrix B_M(syr,future_yrs,1,nage_F);
  matrix B_F(syr,future_yrs,1,nage_F);
  matrix b_mature(syr,future_yrs,1,nage_F);
  
  // vector Total_biom(syr,future_yrs);      //turn on later
  // vector Fem_spawn_biom(syr,future_yrs);
  // vector Vul_biom(syr,future_yrs);
  // number depletion;
  
  sdreport_vector Total_biom(syr,future_yrs);      //turn on later
  sdreport_vector Fem_spawn_biom(syr,future_yrs);
  sdreport_vector Total_biom_rel(syr,future_yrs);      
  sdreport_vector Fem_spawn_biom_rel(syr,future_yrs);
  sdreport_vector Vul_biom(syr,future_yrs);
  sdreport_number depletion;

  sdreport_number Virgin_Total_biom;    
  sdreport_number Virgin_Vul_biom;      
  sdreport_number Virgin_Spawn_Biom;    



  
  //Declare objects used in observation_model()
  vector Q(syr,nyr);
  vector Est_CPUE(syr,nyr);         // predicted cpue 
  vector epsilon(syr_cpue,nyr);	// residuals
  number Penalties;
  //vector tau1(syr,nyr);
  number var1;
  number var2;

  vector Est_CPUE_out(syr_cpue,nyr);
  sdreport_vector log_Est_CPUE_out(syr_cpue,nyr);
  vector CPUE_out(syr_cpue,nyr);
    
  //Declare objects used in calc_obj_function()
  number sums_squares;
  number deg_free;
  number resid_var;
  number num_obs;


PRELIMINARY_CALCS_SECTION




PROCEDURE_SECTION
  virgin_cond();
  init_model();
  dynamic_model();
  observation_model();
  calc_obj_function();


FUNCTION virgin_cond
     int i;

     //virgin recruitment
     //rstar=mfexp(ln_RSTAR);
     rstar=RSTAR*1e4;
     
     //Numbers at age     
     No_M(1)=rstar*sx_ratio;
     No_F(1)=rstar*sx_ratio;
     for(i=2;i<=nage_F;i++)
     {
       No_F(i)=No_F(i-1)*exp(- M_F(i));
       if(age_F(i)==(nage_F-1))  No_F(i)= No_F(i-1)*exp(- M_F(i))/M_F(i);   //plus group used by Simpfendorfer et al 2000

       No_M(i)=No_M(i-1)*exp(- M_M(i));
       if(age_F(i)==(nage_M-1))  No_M(i)= No_M(i-1)*exp(- M_M(i))/M_M(i);
       if(age_F(i)>(nage_M-1))  No_M(i)=0;         
     }
     
     //eggs
     egg_o=elem_prod(No_F,fec);   //numbers of embryos at age
     Egg_o=sum(egg_o);            //total number of embryos 
     

FUNCTION init_model
    //note: this accounts for impact of fishing between 1940 and 1974
    int i;
    Fpen=0.0;
    Zpen=0.0;

    //backtransform z
    z=cos(Z);    
    z=z*z;
    z=-0.2*(pow(z,0.5)-2.15);

    //lp_ln_Init_Fo=Fo;    //likeprofile
    
    //avoid Fo<0
    Fo=pow(pow(Fo,2),0.5);   
      

    //numbers per recruit
    Per_recruit(1)=1;
    for(i=2;i<=nage_F;i++)
    {
      Per_recruit(i)=Per_recruit(i-1)*exp(-(M_F(i)+Fo));
      if(age_F(i)==(nage_F-1))  Per_recruit(i)=Per_recruit(i)/(1-exp(-(M_F(i)+Fo)));  //correct plus group
      if(age_F(i)==(nage_F-1))  Per_recruit(i)=Per_recruit(i-1)*exp(-(M_F(i)+Fo))/(M_F(i)+Fo);  //plus group used by Simpfendorfer
    }   
    
    //female eggs per recruit
    egg_per_rec=elem_prod(Per_recruit,fec*sx_ratio);             
    Egg_per_rec=sum(egg_per_rec);           
       
    //penalty for z>Z.max (Simpfendorfer 2000)
    //z_max=Egg_o/(4*rstar+Egg_o);          //maximimum possible z
    //z=z_max-posfun(z_max-z,0.1,Zpen);
   
    //penalty for z<Z.min (Simpfendorfer 2000)
    //z=posfun(z,0.205,Zpen);
    


    //initial recruitment
    a=Egg_o/rstar*(1-(z-0.2)/(0.8*z));  // R.star sets the scale of initial recruitment 
    b=(z-0.2)/(0.8*z*rstar);
    Ro=(Egg_per_rec-a)/(b*Egg_per_rec);   

  
    //numbers 1975
      // recruits
    N_1975_M(1)=Ro*sx_ratio;
    N_1975_F(1)=Ro*sx_ratio;
    
      // other age classes
    for(i=2;i<=nage_F;i++)
    {
       //female
      N_1975_F(i)=N_1975_F(i-1)*exp(-(M_F(i)+Fo));
      if(age_F(i)==(nage_F-1)) N_1975_F(i)=N_1975_F(i-1)*exp(-(M_F(i)+Fo))/(M_F(i)+Fo);  //plus group used by Simpfendorer et al 2000

       //male
      N_1975_M(i)=N_1975_M(i-1)*exp(-(M_M(i)+Fo));
      if(age_F(i)==(nage_M-1)) N_1975_M(i)=N_1975_M(i-1)*exp(-(M_M(i)+Fo))/(M_M(i)+Fo);
      if(age_F(i)>(nage_M-1))  N_1975_M(i)=0;         
    }   
    
    // Biomasses
    Bo_M=elem_prod(TWT_M,No_M);
    Bo_F=elem_prod(TWT_F,No_F);
    Virgin_Total_biom=sum(Bo_M+Bo_F);
    
    Virgin_Vul_biom=sum(elem_prod(Bo_M,SEL_M)+elem_prod(Bo_F,SEL_F));
    B_1975=sum(elem_prod(TWT_M,N_1975_M)+elem_prod(TWT_F,N_1975_F));
    for(i=1;i<=nage_F;i++)
    {
      if(age_F(i)<age_mat) bo_mature(i)=0;
      if(age_F(i)>=age_mat) bo_mature(i)=Bo_F(i);
    }
    Virgin_Spawn_Biom=sum(bo_mature);
    

FUNCTION dynamic_model
    int i,j;
    Bpen=0.0;
    
    for(i=syr;i<=future_yrs;i++)
    {
      
      if(i==syr)
      {
           N_F(i)(1,nage_F)=N_1975_F;
           N_M(i)(1,nage_F)=N_1975_M;           
      }else
      {
        
        for(j=2;j<=nage_F;j++)
        {
          //female
         if(age_F(j)<(nage_F-1))
         {
           N_F(i,j)=(N_F(i-1,j-1)-c_F(i-1,j-1))*exp(-M_F(j));
          }else
          {
           N_F(i,j)=(N_F(i-1,j-1)-c_F(i-1,j-1)+N_F(i-1,j)-c_F(i-1,j))*exp(-M_F(j)); //plus group
          }
          
            //male
          if(age_F(j)<(nage_M-1)) N_M(i,j)=(N_M(i-1,j-1)-c_M(i-1,j-1))*exp(-M_M(j));
          if(age_F(j)==(nage_M-1)) N_M(i,j)=(N_M(i-1,j-1)-c_M(i-1,j-1)+N_M(i-1,j)-c_M(i-1,j))*exp(-M_M(j));
          if(age_F(j)>(nage_M-1)) N_M(i,j)=0;
        }
      }
      
      //eggs
      EGGS(i)=sum(elem_prod(row(N_F,i)(2,nage_F),fec(2,nage_F)));  
            
      //recruitment
      Annual_rec(i)=EGGS(i)/(a+b*EGGS(i));
      if(i>nyr)                                          
      {
        Annual_rec(i)=Annual_rec(i)*Rec_error(i);   //add variability in recruitment for future projections
      }
      if(i<=future_yrs)
      {
        N_M(i,1)=sx_ratio*Annual_rec(i);      
        N_F(i,1)=sx_ratio*Annual_rec(i);        
      }


      //total biomass
      B_M(i)=elem_prod(TWT_M,N_M(i));
      B_F(i)=elem_prod(TWT_F,N_F(i));
      Total_biom(i)=sum(B_M(i)+B_F(i));
 
      
      //female mature biomass
      for(j=1;j<=(nage_F);j++)
      {
        if(age_F(j)<age_mat) b_mature(i,j)=0;
        if(age_F(j)>=age_mat) b_mature(i,j)=B_F(i,j);
      }
      Fem_spawn_biom(i)=sum(b_mature(i));
      
      //exploitable biomass
      Vul_biom(i)=sum(elem_prod(B_M(i),SEL_M)+elem_prod(B_F(i),SEL_F));
      //Vul_biom(i)=posfun(Vul_biom(i),0.01,Bpen);
      
      //fishing mortality      
      Annual_F(i)=TC_out(i)/Vul_biom(i);
      if(Annual_F(i)>1) Annual_F(i)=1;
      
      //predicted catch numbers
      c_F(i)(1,nage_F)=elem_prod(N_F(i)*Annual_F(i),SEL_F);
      c_M(i)(1,nage_F)=elem_prod(N_M(i)*Annual_F(i),SEL_M);

 
      //predicted catch weight
      c_F_wgt(i)=elem_prod(c_F(i),TWT_F);
      c_M_wgt(i)=elem_prod(c_M(i),TWT_M);
      Est_TC_out(i)=sum(c_F_wgt(i)+c_M_wgt(i));
      
    } //end i
    
    //depletion level
    depletion=Total_biom(nyr)/Virgin_Total_biom;
    for(i=syr;i<=future_yrs;i++)
    {
     Total_biom_rel(i)=Total_biom(i)/Virgin_Total_biom;
     Fem_spawn_biom_rel(i)=Fem_spawn_biom(i)/Virgin_Spawn_Biom;
    }
    
    
    
    //cout<<"N_F\n"<<column(N_F,1)<<endl;
    //cout<<"B_M\n"<<row(B_M,1975)<<endl;
    //exit(1);

     
FUNCTION observation_model
     int i;
     //use different catchabilities according to period
      for(i=syr;i<=nyr;i++)                      
     {
       if(effective_cpue==1)
       {
         if(q_change>0)
         {
           if(i<=q_change)
           {
              Q(i)=Q1*1e-6;
              //Q(i)=mfexp(Q1);
           }else
           {
              Q(i)=Q2*1e-6;
              //Q(i)=mfexp(Q2);
           }
         }else
         {
           Q(i)=Q1*1e-6;
           //Q(i)=mfexp(Q1);
         }
       
       }else
       {
         if(q_change>0)
         {
          if(i<=q_change)
          {
            Q(i)=Q1*1e-6;
            //Q(i)=mfexp(Q1);

          //Haddon 2001 increase in efficiency
        // if(i<=1993)Q(i)=mfexp(log_Q);
        // if(i>1993 & i<2005)Q(i)=mfexp(log_Q+(i-1993)*log(1.1));
        // if(i>=2005) Q(i)=mfexp(log_Q+(i-1993)*log(1.05));
        //cout<<"i "<<i<<"  log_Q "<<log_Q<<"  (i-1993) "<<(i-1993)<<" Q(i) "<<Q(i)<<endl; 
          }else
          {
            if(i<q_daily)
            {
               Q(i)=Q2*1e-6;
               //Q(i)=mfexp(Q2);
            }else
            {
               Q(i)=Qdaily*1e-6;
               //Q(i)=mfexp(Qdaily);
            }
          }    
         }else
         {
            if(i<q_daily)
            {
               Q(i)=Q1*1e-6;
               //Q(i)=mfexp(Q1);
            }else
            {
               Q(i)=Qdaily*1e-6;
               //Q(i)=mfexp(Qdaily);
            }
         }
       }    
       
      }//i
     
      
     //estimated cpue
     Est_CPUE=elem_prod(Q,Vul_biom(syr,nyr));			
     
     //residuals
     for(i=syr;i<=nyr;i++) 
     {
     if(CPUE(i)>0) epsilon(i)=log(CPUE(i)+1e-6)-log(Est_CPUE(i)+1e-6);
     }

          
        //output quantities
    Est_CPUE_out=Est_CPUE(syr_cpue,nyr);
    log_Est_CPUE_out=log(Est_CPUE(syr_cpue,nyr));    
    CPUE_out=CPUE(syr_cpue,nyr);
 

   
     //residual variance
     if(Do_var==1)
     {
      if(q_change>0)
      {
        var1=sum(pow((epsilon(syr,q_change)-mean(epsilon(syr,q_change))),2))/((q_change-syr));
        var2=sum(pow((epsilon((q_change+1),nyr)-mean(epsilon((q_change+1),nyr))),2))/((nyr-q_change-1));
      }else
      {
       var1=sum(pow((epsilon-mean(epsilon)),2))/((nyr-syr));
       var2=var1;
      }
     }
     if(Do_var==0)     //fixed var as done in Simpfendorfer et al 2001
     {
        var1=VR1;
        var2=VR2;
     }

     
     //corrected residuals
     if(q_change>0)
     {
       epsilon(syr,q_change)=epsilon(syr,q_change)/(pow(var1,0.5));
       epsilon((q_change+1),nyr)=epsilon((q_change+1),nyr)/(pow(var2,0.5));
     }


 
FUNCTION calc_obj_function
     //dvariable tau2=exp(log_tau2);
     dvariable CPUE_like;
     dvariable Catch_like;

     Penalties=0.0;
     //Penalties=100*Fpen+1e7*Zpen+1e7*Bpen;
           
      
     //cpue likelihood
     int n=size_count(epsilon);
     CPUE_like=0;

     //--Sums of squares (Simpfendorfer's approach)
     sums_squares=sum(pow(epsilon,2));
     num_obs=nyr-syr+1;
     deg_free=((num_obs-n_par-1));
     resid_var=sums_squares/deg_free;
     CPUE_like=-((sums_squares/(2*resid_var))-(num_obs*log(pow(resid_var*2*3.1417,0.5))));

  
     //--Simple lognormal likelihood
     //CPUE_like+=0.5*n*log(tau2)+(0.5/square(tau2))*norm2(epsilon);
     
     //--Francis 2011 T2.5 negloglike     
     // int i;     
     // for (i=syr; i<=nyr;i++)
     // {
     //  //combine model error and observation error from cpue standardisation
     //  tau1(i)=sqrt(square(tau2)+square(CPUE_SD(i)));        //total SD is (sum of variances)^1/2
     //  CPUE_like+=log(tau1(i))+0.5*(square((log(CPUE(i)/Est_CPUE(i)))/tau1(i)));
     //  }

     
     //catch penalty (keep similar to observed catch)
     // Catch_like=0.0;
     //Catch_like+=sum(square(TC_out-Est_TC_out));
 
     //combine loglikelihoods
     f=CPUE_like+Penalties;
     //cout<<"rstar="<<rstar<<" z="<<z<<" Fo="<<Fo<<"  Q1="<<Q1<<" Q2="<<Q2<<endl;
     // cout<<"100*Fpen="<<100*Fpen<<"  1e7*Zpen="<<1e7*Zpen<<"  1e7*Bpen="<<1e7*Bpen<<"  CPUE_like="<<CPUE_like<<endl;
     // cout<<"Est_CPUE\n"<<Est_CPUE<<endl;
     // cout<<"CPUE\n"<<CPUE<<endl;
     //   cout<<"epsilon "<<epsilon<<endl;
     //   cout<<"sums_squares "<<sums_squares<<endl;
     //  cout<<"resid_var "<<resid_var<<endl;
     //  cout<<"CPUE_like "<<CPUE_like<<endl;
      //cout<<" "<<RSTAR<<" "<<Z<<" "<<Q1<<" "<<Q2<<endl;
      //cout<<"Penalties= "<<Penalties<<"  z="<<z<<"  zmax="<<z_max<<"  CPUE_like "<<CPUE_like<<endl;
     // exit(1);


REPORT_SECTION
    REPORT(Fpen);
    REPORT(Zpen);
    REPORT(Bpen);
    REPORT(CPUE_out);                 //observed cpue
    REPORT(Est_CPUE_out);             //predicted cpue
    REPORT(Est_TC_out);
    REPORT(TC_out);
    REPORT(Annual_F);
    REPORT(depletion);
    REPORT(Total_biom);
    REPORT(Fem_spawn_biom);
    REPORT(Vul_biom);  
    REPORT(epsilon);
    REPORT(fec);
    REPORT(SEL_F);
    REPORT(SEL_M);
    REPORT(TWT_F);
    REPORT(TWT_M);
    REPORT(ct_F);
    REPORT(ct_M);
    REPORT(c_F_wgt);
    REPORT(c_M_wgt);
    REPORT(Annual_rec);
 
 

// TOP_OF_MAIN_SECTION
// 	time(&start);
// 	arrmblsize = 50000000;
// 	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
// 	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
// 	gradient_structure::set_MAX_NVAR_OFFSET(5000);
// 	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
 

// GLOBALS_SECTION
// //anything declared in the GLOBAL section is global throughout!!
// 	/**
// 	\def REPORT(object)
// 	Prints name and value of \a object on ADMB report %ofstream file.
// 	*/
// 	#undef REPORT
// 	#define REPORT(object) report << #object "\n" << object << endl;	// Report macro. object is all the stuff I declared in REPORT()

// 	#include <admodel.h>
// 	#include <time.h>
// 	time_t start,finish;
// 	long hour,minute,second;
// 	double elapsed_time;
	
// FINAL_SECTION
// 	time(&finish);
// 	elapsed_time=difftime(finish,start);
// 	hour=long(elapsed_time)/3600;
// 	minute=long(elapsed_time)%3600/60;
// 	second=(long(elapsed_time)%3600)%60;
// 	cout<<endl<<endl<<"*******************************************"<<endl;
// 	cout<<"--Start time: "<<ctime(&start)<<endl;
// 	cout<<"--Finish time: "<<ctime(&finish)<<endl;
// 	cout<<"--Runtime: ";
// 	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
// 	cout<<"*******************************************"<<endl;

