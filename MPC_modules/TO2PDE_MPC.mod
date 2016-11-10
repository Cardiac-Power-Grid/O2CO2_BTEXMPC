// MODEL NAME: TO2PDE_MPC
// SHORT DESCRIPTION: Modeling Total O2 concentration as a function of t and x.

import nsrunit;   unit conversion on;

math TO2ODE_MPC {

//%START TIME_DOMAIN
realDomain t sec;t.min=0;t.max=100;t.delta=0.1;
//%END TIME_DOMAIN

//%START SPATIALDOMAIN
realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid; // BTEX space domain 
//%END SPATIALDOMAIN


//%START PARAMS
real Frbc ml/(min*g);           // Flow of RBCs.
real Vrbc ml/g;                 // Volume of region
real Wrbc dimensionless;	// fractional water content of RBCs
real VWrbc ml/g;	        // volume of water content in RBCs
real O2Drbc = 1e-4 cm^2/sec;   // diffusion coefficient for O2 in RBCs 
real O2PSrbc ml/(g*min);       // PS for O2 across RBC Membrane
//%END PARAMS

//%START FLOW_INPUT
real TO2_rbc_Fin(t) M;   // TO2 coming into the rbc region (x.min)
//%END FLOW_INPUT

//%START INPUT
real O2_rbc_in(t,x) M;     // Free O2 in RBCs
real O2_pl_in(t,x) M;        // Free O2 in plasma
//%END INPUT

//%START VARIABLES
real TO2_rbc(t,x) M;     // Total O2 in RBCs
//%END VARIABLES

//%START VAR_OUTPUT
real TO2_rbcOut(t) M;

//%END VAR_OUTPUT

//%START VAR_INITS
real TO2Crbct0 M;
//%END VAR_INITS


 when (t=t.min) {	// TO2 PDE INITIAL CONDITIONS
//%START TO2_PDE_IC
      TO2_rbc = TO2Crbct0;  // TO2 PDE INITIAL CONDITION
//%END TO2_PDE_IC
     } 



 when (x=x.min) {	// TO2 RBC PDE LEFT BOUNDARY CONDITIONS
//%START TO2_PDE_LBC
     (-Frbc*L/Vrbc)*(TO2_rbc-TO2_rbc_Fin)+O2Drbc*TO2_rbc:x =0;
//%END TO2_PDE_LBC  
 }   
 when (x=x.max) {	// TO2 RBC PDE RIGHT BOUNDARY CONDITIONS
//%START TO2_PDE_RBC
     O2Drbc*TO2_rbc:x = 0;
     TO2_rbcOut = TO2_rbc;
//%END TO2_PDE_RBC
 }  


//%START TO2_PDES
// TO2 PDE:
 TO2_rbc:t = - (Frbc/Vrbc)*L*(TO2_rbc:x) 
                   - (O2PSrbc/VWrbc)*(O2_rbc_in-O2_pl_in)
                   + O2Drbc*(TO2_rbc:x:x);             
//%END TO2_PDES



// To make stand-alone model:
     TO2Crbct0 =.01459;
     O2_rbc_in = .0023;  // M
     O2_pl_in = 0023;
     O2PSrbc = 0;
     TO2_rbc_Fin = 0.02;
     Frbc =1;
     Vrbc = 0.01;      // ml/g
     Wrbc = 0.72;
     VWrbc = Vrbc*Wrbc;

}
