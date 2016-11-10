// MODEL NAME: O2pc_MPC
// SHORT DESCRIPTION: Modeling Total O2 concentration as a function of t and x in pc.
// Accounts for Myoglobin binding and consumption.

import nsrunit;   unit conversion on;

math O2pc_MPC {

//%START TIME_DOMAIN
realDomain t sec;t.min=0;t.max=100;t.delta=0.1;
//%END TIME_DOMAIN

//%START SPATIALDOMAIN
realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid; // BTEX space domain 
//%END SPATIALDOMAIN


//%START PARAMS
//real Fpc ml/(min*g);           // Flow of pc. Assume zero.
real Vpc ml/g;                 // Volume of region
real Wpc dimensionless;	       // fractional water content of pc
real VWpc ml/g;	               // volume of water content in pc
real O2Dpc = 1e-4 cm^2/sec;    // diffusion coefficient for O2 in pc 
real O2PSpc ml/(g*min);        // PS for O2 across PC Membrane
//%END PARAMS

//%START SOLUBILITY_PARAM
real alphaO2(t,x) M/mmHg;      // Solubility of O2 in pc
//%END SOLUBILITY_PARAM

//%START CONSUMPTION_VARS
real Vmax_pc = 5.0e-6 mol/min/g; // Vmax for O2 consumpt by cytochrome oxidase in PCs ;
real Km_pc = 7.0e-7 M;         // Km for O2 on cytochrome oxidase in PCs;
real O2Gpc(t,x) ml/min/g;      // gulosity for O2 consumption in PCs
real MRO2pc(t,x) mol/min/g;    // O2 consumption rate in PCs 
//%END CONSUMPTION_VARS

//%START PC_Mb_RATE_CONSTS
// pc MbRateConsts:
real  CMb    = 0.4e-3 M,	// CONC of Mb in PCs = 5/(16800*FCV) M
      CMbpc  = CMb/Wpc,		// CONC of Mb in the water space of PCs
      P50Mb  = 2.4 mmHg;	// P50 for MbO2 saturation
// pc END MbRateConsts
//%END PC_Mb_RATE_CONSTS

//%START INPUT
real O2_isf_in(t,x) M;     // Free O2 in ISF
//%END INPUT

//%START MbO2_VARS
real KMbO2(t,x) M^(-1);	// Hill's coefficient for MbO2 saturation
real SMbO2(t,x);	// MbO2 saturation (unitless)
real CMbO2(t,x) M;      // concentration of MbO2 in RBCs
//%END MbO2_VARS

//%START VARIABLES
real O2_pc(t,x) M;      // O2 in pc
real TO2_pc(t,x) M;     // total O2 in pc
//%END VARIABLES

//%START VAR_OUTPUT
real pO2_pc(t,x) mmHg;   // convert to partial press
     pO2_pc = O2_pc/alphaO2;
//%END VAR_OUTPUT

//%START VAR_INITS
real TO2Cpct0(x) M;
real KMbO2t0(x) M^(-1);
real SMbO2t0(x);
real pO2pct0 mmHg;
real O2Cpct0 M;
real CMbO2t0(x) M;
 
//%END VAR_INITS


 when (t=t.min) {	// O2 PDE INITIAL CONDITIONS
//%START O2_PC_IC
      TO2_pc = TO2Cpct0;  // O2 PDE INITIAL CONDITION
//%END O2_PC_IC
 } 



 when (x=x.min) {	// TO2 RBC PDE LEFT BOUNDARY CONDITIONS
//%START O2_PC_LBC
     O2Dpc*TO2_pc:x =0;
//%END O2_PC_LBC  
 }   
 when (x=x.max) {	// TO2 RBC PDE RIGHT BOUNDARY CONDITIONS
//%START O2_PC_RBC
     O2Dpc*TO2_pc:x = 0;
     
//%END O2_PC_RBC
 }  


//%START O2_PC_PDES
// O2 PDE:
 TO2_pc:t = (O2PSpc/VWpc)*(O2_isf_in-O2_pc)
           + O2Dpc*(TO2_pc:x:x)
           - (O2Gpc/VWpc)*O2_pc;  // Consumption term             
//%END O2_PC_PDES

//%START MbO2_BINDING
// MbO2 binding:
 O2_pc = if (TO2_pc <= 0) 0
         else (((1+KMbO2*(CMbpc-TO2_pc))^2 + 4*KMbO2*TO2_pc)^0.5
             - (1+KMbO2*(CMbpc-TO2_pc)))/(2*KMbO2);

// HillEqs for Mb:
      KMbO2  = (alphaO2*P50Mb)^(-1);
      SMbO2  = KMbO2*O2_pc/(1+KMbO2*O2_pc);  // Myoglobin saturation in pc
      CMbO2  = CMbpc*SMbO2; 
//%END MbO2_BINDING

//%START O2_CONSUMPTION
     O2Gpc = Vmax_pc/(Km_pc+O2_pc);
     MRO2pc = O2Gpc*O2_pc;     // <--- get an estimate of O2 consumed 
//%END O2_CONSUMPTION

// To make stand-alone model:
real alphaO20 = 1.1148E-6;     // Initial solubility
     alphaO2 = 1.1148E-6 ;     // Solubility of O2 in pc

real pO2pct0 = 90 mmHg;

     
     O2Cpct0 =  pO2pct0 * alphaO20; 
     KMbO2t0  = (alphaO20*P50Mb)^(-1);
     SMbO2t0  = KMbO2t0*O2Cpct0/(1+KMbO2t0*O2Cpct0);
     CMbO2t0  = CMbpc*SMbO2t0;
     TO2Cpct0 = O2Cpct0+CMbO2t0;

     O2_isf_in = 0;  // M
     O2PSpc = 0;

     Vpc = 0.01;      // ml/g
     Wpc = 0.72;
     VWpc = Vpc*Wpc;

}
