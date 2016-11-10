// MODEL NAME: O2plasma_MPC
// SHORT DESCRIPTION: Modeling Total O2 concentration as a function of t and x in plasma.

import nsrunit;   unit conversion on;

math O2plasma_MPC {

//%START TIME_DOMAIN
realDomain t sec;t.min=0;t.max=100;t.delta=0.1;
//%END TIME_DOMAIN

//%START SPATIALDOMAIN
realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid; // BTEX space domain 
//%END SPATIALDOMAIN


//%START PARAMS
real Fpl ml/(min*g);           // Flow of pls.
real Vpl ml/g;                 // Volume of region
real Wpl dimensionless;	       // fractional water content of plasma
real VWpl ml/g;	               // volume of water content in plasma
real O2Dpl = 1e-4 cm^2/sec;    // diffusion coefficient for O2 in RBCs 
real O2PSrbc ml/(g*min);       // PS for O2 across RBC Membrane
real O2PScap ml/(g*min);       // PS for O2 across Capillary Membrane
//%END PARAMS

//%START FLOW_INPUT
real O2_pl_Fin(t) M;   // O2 coming into the plasma region (x.min)
//%END FLOW_INPUT
//%START INPUT
real O2_rbc_in(t,x) M;     // Free O2 in RBCs
real O2_isf_in(t,x) M;     // Free O2 in plasma
//%END INPUT

//%START VARIABLES
real O2_pl(t,x) M;      // O2 in plasma
//%END VARIABLES

//%START VAR_OUTPUT
real O2_plOut(t) M;

//%END VAR_OUTPUT

//%START VAR_INITS
real O2Cplt0 M;
//%END VAR_INITS


 when (t=t.min) {	// O2 PDE INITIAL CONDITIONS
//%START O2_PLASMA_IC
      O2_pl = O2Cplt0;  // O2 PDE INITIAL CONDITION
//%END O2_PLASMA_IC
 } 



 when (x=x.min) {	// TO2 RBC PDE LEFT BOUNDARY CONDITIONS
//%START O2_PLASMA_LBC
     (-Fpl*L/Vpl)*(O2_pl-O2_pl_Fin)+O2Dpl*O2_pl:x =0;
//%END O2_PLASMA_LBC  
 }   
 when (x=x.max) {	// TO2 RBC PDE RIGHT BOUNDARY CONDITIONS
//%START O2_PLASMA_RBC
     O2Dpl*O2_pl:x = 0;
     O2_plOut = O2_pl;
//%END O2_PLASMA_RBC
 }  


//%START O2_PLASMA_PDES
// O2 PDE:
 O2_pl:t = - (Fpl/Vpl)*L*(O2_pl:x) 
           + (O2PSrbc/VWpl)*(O2_rbc_in-O2_pl)
           - (O2PScap/VWpl)*(O2_pl-O2_isf_in)
           + O2Dpl*(O2_pl:x:x);             
//%END O2_PLASMA_PDES



// To make stand-alone model:
 real alphaO2 = 1.1148E-6 M/mmHg; // Solubility of O2 in plasma
 real pO2_pl(t,x) mmHg;
 real pO2plt0 = 90 mmHg;
     pO2_pl = O2_pl/alphaO2;
      
     O2Cplt0 =  pO2plt0 * alphaO2; //.01459;
     O2_rbc_in = 0;  // M
     O2_isf_in = 0.0;
     O2PSrbc = 0;
     O2PScap = 0;
     O2_pl_Fin = pO2plt0 * alphaO2;
     Fpl =1;
     Vpl = 0.01;      // ml/g
     Wpl = 0.72;
     VWpl = Vpl*Wpl;

}
