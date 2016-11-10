// MODEL NAME: O2isf_MPC
// SHORT DESCRIPTION: Modeling Total O2 concentration as a function of t and x in isf.

import nsrunit;   unit conversion on;

math O2isf_MPC {

//%START TIME_DOMAIN
realDomain t sec;t.min=0;t.max=100;t.delta=0.1;
//%END TIME_DOMAIN

//%START SPATIALDOMAIN
realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid; // BTEX space domain 
//%END SPATIALDOMAIN


//%START PARAMS
real Fisf ml/(min*g);           // Flow of isf. Assume zero.
real Visf ml/g;                 // Volume of region
real Wisf dimensionless;	       // fractional water content of isf
real VWisf ml/g;	               // volume of water content in isf
real O2Disf = 1e-4 cm^2/sec;    // diffusion coefficient for O2 in RBCs 
real O2PSpc ml/(g*min);       // PS for O2 across PC Membrane
real O2PScap ml/(g*min);       // PS for O2 across Capillary Membrane
//%END PARAMS

//%START SOLUBILITY_PARAM
real alphaO2(t,x) M/mmHg;      // Solubility of O2 in isf
//%END SOLUBILITY_PARAM

//%START INPUT
real O2_pl_in(t,x) M;     // Free O2 in plasma
real O2_pc_in(t,x) M;     // Free O2 in pc
//%END INPUT

//%START VARIABLES
real O2_isf(t,x) M;      // O2 in isf
//%END VARIABLES

//%START VAR_OUTPUT
real pO2_isf(t,x) mmHg;   // convert to partial press
     pO2_isf = O2_isf/alphaO2;
//%END VAR_OUTPUT

//%START VAR_INITS
real O2Cisft0 M;
//%END VAR_INITS


 when (t=t.min) {	// O2 PDE INITIAL CONDITIONS
//%START O2_ISF_IC
      O2_isf = O2Cisft0;  // O2 PDE INITIAL CONDITION
//%END O2_ISF_IC
 } 



 when (x=x.min) {	// TO2 RBC PDE LEFT BOUNDARY CONDITIONS
//%START O2_ISF_LBC
     O2Disf*O2_isf:x =0;
//%END O2_ISF_LBC  
 }   
 when (x=x.max) {	// TO2 RBC PDE RIGHT BOUNDARY CONDITIONS
//%START O2_ISF_RBC
     O2Disf*O2_isf:x = 0;
//%END O2_ISF_RBC
 }  


//%START O2_ISF_PDES
// O2 PDE:
 O2_isf:t = + (O2PScap/VWisf)*(O2_pl_in-O2_isf)
           - (O2PSpc/VWisf)*(O2_isf-O2_pc_in)
           + O2Disf*(O2_isf:x:x);             
//%END O2_ISF_PDES



// To make stand-alone model:
real alphaO20 = 1.1148E-6 M/mmHg; // Solubility of O2 in isf
     alphaO2 = 1.1148E-6;
real pO2isft0 = 90 mmHg;
  
     O2Cisft0 =  pO2isft0 * alphaO20; //.01459;
     O2_pl_in = 0;  // M
     O2_pc_in = 0.0;
     O2PSpc = 0;
     O2PScap = 0;
     Fisf =1;
     Visf = 0.01;      // ml/g
     Wisf = 0.72;
     VWisf = Visf*Wisf;

}
