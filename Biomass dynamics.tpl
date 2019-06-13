
//  Surplus production model


DATA_SECTION
  init_int syr                    //start modelling year  
  init_int nyr	                  //end modelling eyar
  init_int nzone                  //number of zones
  init_int future_yrs            //end of future projections	
  init_ivector iyr(syr,nyr)          //years of catch data
  init_vector TC_out(syr,future_yrs)
  
  init_int effective_cpue         //use effective cpue (1) or standardise cpue (0)
  init_int syr_cpue               //start year of cpue data
  
  init_vector CPUE(syr,nyr)
  init_vector CPUE_SD(syr,nyr)
  
  init_int q_change               //first change in q
  init_int q_daily                //change in q due to daily reporting
  
  init_number p0                 //population proportion at start of time series  
  
  // r prior
  init_number r_max
  init_int add_r_prior
  init_number mu_r                 
  init_number SD_r
  
  //Phases                                      
  init_int Phase_dummy
  init_int Phase_r                 
  init_int Phase_log_k
  init_int Phase_Q1  
  init_int Phase_Q2  
  init_int Phase_Qdaily                 
  init_int Phase_log_tau2	
  
  init_int eof
  //error trap to ensure you read all the data.
  !! if(eof != 999){cout<<"Error reading data"<<endl; exit(1);};
  int niter;
  !! niter=0;
       

PARAMETER_SECTION
  init_number Dummy(Phase_dummy);                   //dummy for controling phases                                 
  init_bounded_number r(0.01,1.0,Phase_r)	
  init_bounded_number log_k(1,20.0,Phase_log_k)	//min bound set to max observed catch
  init_bounded_number Q1(0,1.0,Phase_Q1)    //Qs bounded between 0 and 1 because working with proportions
  init_bounded_number Q2(0,1.0,Phase_Q2)
  init_bounded_number Qdaily(0,1,Phase_Qdaily)
  init_bounded_number log_tau2(-10.0,1.0,Phase_log_tau2)	// the TOTAL variance (observation and process)

  objective_function_value f;
	
  number k;				// carrying capacity
  number bpen;                           // penalty to keep biomass positive
  vector pt(syr,future_yrs);		//population size relative to k
  vector Total_biom_temp(syr,future_yrs);      //population biomass  
  vector U(syr,future_yrs);		        // Exploitation rate
  vector Est_CPUE(syr,nyr);	// predicted cpue index
  vector epsilon(syr,nyr);	// predicted observation error residuals
  vector q(syr,nyr);
  number Penalties;
  number CPUE_like;
  number FPen;      //penaltity for initial fishing mortality set to initial value provided
  number rpen;   //penalty to keep r below max possible for sharks
  number Upen;   //penalty to keep U below 1
  number Fo;

  sdreport_vector Total_biom(syr+1,future_yrs);	// output biomass (+1 because not estimating first year)
  sdreport_vector Total_biom_rel(syr+1,future_yrs);	
  sdreport_vector F(syr+1,future_yrs);		// Fishing mortality (+1 because not estimating first year)
  sdreport_number depletion; 	//depletion level (i.e. current B/k) an MCMC. This gets copied to the .rep file
  sdreport_number fmsy;
  sdreport_number msy;
  sdreport_number bmsy;
  sdreport_number Virgin_Total_biom;
  vector tau1(syr,nyr);          //observation error  

  vector CPUE_out(syr,nyr);
  vector Est_CPUE_out(syr,nyr);
  sdreport_vector log_Est_CPUE_out(syr,nyr);
  vector CPUE_SD_out(syr,nyr);
  

	
PROCEDURE_SECTION 
  population_model();		
  observation_model();		
  calc_objective_function();	
  if(mceval_phase()) write_mcmc_report();	
  
  fmsy=0.5*r;	 
  msy=k*r/4.0;
  bmsy=0.5*k;
  
  

	
FUNCTION population_model
  int i;			
  k=exp(log_k);
  
  rpen=0.0;   //penalty to avoid unrealistic r values
  r=r_max-posfun(r_max-r,0.1,rpen);  //r_max, max value for sharks (blue shark)

  bpen=0.0;    // penalty to negative biomass
  pt(syr)=p0;
  for(i=syr+1;i<=future_yrs;i++)
  {
   pt(i)=pt(i-1)+r*pt(i-1)*(1.0-pt(i-1))-TC_out(i-1)/k; 
   pt(i)=posfun(pt(i),0.001,bpen);	//keep biomass positive   
  }
  
  Total_biom_temp=pt*k;
  for(i=syr+1;i<=future_yrs;i++)Total_biom(i)=Total_biom_temp(i);  //extract sdreport biomass

  //exploitation rate
  U=elem_div(TC_out,Total_biom_temp);
  Upen=0.0;                       
  U=1-posfun(1-U,0.1,Upen);    // penalty to keep U<=1
  
  //fishing mortality
  Fo=-log(1-U(syr));  
  for(i=syr+1;i<=future_yrs;i++) F(i)=-log(1-U(i));
  FPen=0.0;
  //FPen=(F(syr+1)-Fo)*(F(syr+1)-Fo);   //set initial F to Fo provided in .dat file
  
  //if(min(Total_biom_temp)<0)exit(1);
  depletion=Total_biom_temp(nyr)/k;
  Virgin_Total_biom=k;
  for(i=syr+1;i<=future_yrs;i++)  Total_biom_rel(i)= Total_biom(i)/k;     

	
FUNCTION observation_model
   int i;
   for(i=syr;i<=nyr;i++)
     {
       if(effective_cpue==1)
       {
         if(q_change>0)
         {
           if(i<=q_change)
           {
             q(i)=Q1;               
           }else
           {
              q(i)=Q2;
           }
         }else
         {
          q(i)=Q1;
         }
       }else
       {
         if(q_change>0)
         {
           if(i<=q_change)
           {
             q(i)=Q1;               
           }else
           {
              if(i<q_daily)
              {
                q(i)=Q2;
              }else
              {
                q(i)=Qdaily;
              }
            }
         }else
         {
           if(i<q_daily)
           {
              q(i)=Q1;
           }else
           {
              q(i)=Qdaily;
           }
          }
         }
     }// end i
     	
    Est_CPUE=elem_prod(q,Total_biom_temp(syr,nyr));   //predicted cpue index				

    for(i=syr_cpue;i<=nyr;i++)
    {
     if(CPUE(i)>0)  epsilon+=log(CPUE(i))-log(Est_CPUE(i));		//epsilon vector
    }
    

      // cout<<"Q\n"<<q<<endl;
      // cout<<"pt\n"<<pt<<endl;
      // cout<<"Total_biom\n"<<Total_biom<<endl;
      // cout<<"Est_CPUE\n"<<Est_CPUE<<endl;
      // cout<<"CPUE\n"<<CPUE<<endl;exit(1);

      //output quantities
     CPUE_out=CPUE(syr_cpue,nyr);
     Est_CPUE_out=Est_CPUE(syr_cpue,nyr);
     log_Est_CPUE_out=log(Est_CPUE(syr_cpue,nyr));
     CPUE_SD_out=CPUE_SD(syr_cpue,nyr);


FUNCTION calc_objective_function
     dvariable tau2=exp(log_tau2);
     
     Penalties=0.0;
     //if(!last_phase()) Penalties=FPen+1e6*rpen+1e-2*Upen;
     Penalties=1e3*bpen+1e5*rpen+1e-8*Upen;
    
     //cpue likelihood
     int i; 
     int n=size_count(epsilon);
     CPUE_like=0;
     
       
         //--sums of squares
      // for (i=syr_cpue; i<=nyr;i++)
      // {
      //   CPUE_like+= (log(CPUE(i))-log(Est_CPUE(i))) * (log(CPUE(i))-log(Est_CPUE(i)));
      // }
       
          
        //--simple lognormal
     //CPUE_like+=0.5*n*log(tau2)+(0.5/square(tau2))*norm2(epsilon);
     //CPUE_like=0.5*n*log(tau2)+0.5*norm2(epsilon)/tau2;

      //--Francis 2011 T2.5 negloglike
      for (i=syr_cpue; i<=nyr;i++)
      {
         if(CPUE(i)>0)
         {
           //combine model error and observation error from cpue standardisation
           tau1(i)=sqrt(square(tau2)+square(CPUE_SD(i)));        //total SD is (sum of variances)^1/2
           CPUE_like+=log(tau1(i))+0.5*(square((log(CPUE(i)/Est_CPUE(i)))/tau1(i)));
         }
       }
      cout<<"CPUE_like="<<CPUE_like<<"  1e3*bpen="<<1e3*bpen<<"  1e5*rpen="<<1e5*rpen<<"  1e-8*Upen="<<1e-8*Upen<<endl;
      f=CPUE_like+Penalties;
     
     //add r prior
     if(add_r_prior==1)f+=dlnorm(r,mu_r,SD_r)*1e2;

 
      

REPORT_SECTION
        REPORT(CPUE_out);                 //observed cpue
        REPORT(Est_CPUE_out);             //predicted cpue
        REPORT(CPUE_SD_out);
	REPORT(Total_biom);
        REPORT(F);
        REPORT(Fo);
        REPORT(depletion);
        REPORT(fmsy);
        REPORT(msy);
        REPORT(bmsy);
        REPORT(TC_out);
        REPORT(Penalties);
        REPORT(CPUE_like);
        

FUNCTION write_mcmc_report
	 //double fmsy=value(0.5*r);	//open up the file and append the MSY, FMSY, BMSY calculaitons
	// double msy=value(k*r/4.0);
	 //double bmsy=value(0.5*k);
	
	if(niter==0)
        {
	  ofstream ofs("pm.mcmc");
	  ofs<<"r\t k\t fmsy\t msy\t bmsy\t depletion"<<endl;
          ofstream ofs1("Total_biom_temp.mcmc");
          ofs1<<"Total_biom_temp"<<endl;            
	}		//a loop to avoid writing all the staff again
	
	ofstream ofs("pm.mcmc",ios::app);	//open file pm.mcmc and make sure that I will append stuff (iop=input output)
	ofs<<r<<"\t"<<k<<"\t"<<fmsy<<"\t"<<msy<<"\t"<<bmsy<<"\t"<<depletion<<endl;	//appending all this stuff in ofs
        ofstream ofs1("Total_biom_temp.mcmc",ios::app);	
        ofs1<<Total_biom_temp<<endl;
	niter++;	// to make niter go to next
	
	
TOP_OF_MAIN_SECTION
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
 

GLOBALS_SECTION
//anything declared in the GLOBAL section is global throughout!!
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;	// Report macro. object is all the stuff I declared in REPORT()

	#include <admodel.h>
	#include <time.h>
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	
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

