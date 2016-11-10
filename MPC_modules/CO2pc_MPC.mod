// MODEL NAME: CO2pc_MPC
// SHORT DESCRIPTION: Modeling Total CO2 concentration as a function of t and x in pc region.


import nsrunit;   unit conversion on;

math CO2pc_MPC {

//%START TIME_DOMAIN
realDomain t sec;t.min=0;t.max=100;t.delta=0.1;
//%END TIME_DOMAIN

//%START SPATIALDOMAIN
realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid; // BTEX space domain 
//%END SPATIALDOMAIN

//%START RXN_PARAMETERS
real CFpc = 10000;		// catalytic factor in pc due to CA
real BCpc = 45 mM;		// buffering capacity in pc
real kp1    = 0.12 1/sec,	// forward rate constant in CO2+H2O reaction
     km1    = 89 1/sec;         // backward rate constant in CO2+H2O reaction
real KH2CO3 = 5.5e-4 M;	        // EQUIL constant in H2CO3 ionization
//%END RXN_PARAMETERS

//%START PARAMS
//real Fpc ml/(min*g);           // Flow of pc. Assume zero.
real Vpc ml/g;                 // Volume of region
real Wpc dimensionless;	       // fractional water content of pc
real VWpc ml/g;	               // volume of water content in pc
real CO2Dpc = 1e-4 cm^2/sec;   // diffusion coefficient for CO2 in pc 
real CO2PSpc ml/(g*min);       // PS for CO2 across pc Membrane
//%END PARAMS

//%START SOLUBILITY_PARAM
real alphaCO2(t,x) M/mmHg;      // Solubility of CO2 in pc
//%END SOLUBILITY_PARAM

//%START CONSUMPTION_VARS
real O2Gpc(t,x) ml/min/g;      // gulosity for O2 consumption in PCs
real RQ = 0.80 dimensionless;  // respiratory quotient (unitless)
//%END CONSUMPTION_VARS

//%START HCO3_H_INPUT
real HCO3m_pc_in(t,x) M;    // HCO3- in pc
real Hp_pc_in(t,x) M;       // H+ conc in pc
//%END HCO3_H_INPUT

//%START O2_INPUT
real O2_pc_in(t,x) M;       // Needed to calculate CO2 production
//%END O2_INPUT

//%START INPUT
real CO2_isf_in(t,x) M;     // Free CO2 from isf
//%END INPUT

//%START VARIABLES
real CO2_pc(t,x) M;      // CO2 in pc
//%END VARIABLES

//%START VAR_OUTPUT
real pCO2_pc(t,x) mmHg;   // convert to partial press
     pCO2_pc = CO2_pc/alphaCO2;
//%END VAR_OUTPUT

//%START VAR_INITS
real CO2Cpct0 M;
//%END VAR_INITS


 when (t=t.min) {	// CO2 ISF PDE INITIAL CONDITIONS
//%START CO2_PDE_IC
      CO2_pc = CO2Cpct0;
//%END CO2_PDE_IC
 } 

 when (x=x.min) {	// CO2 ISF PDE LEFT BOUNDARY CONDITIONS
//%START CO2_PDE_LBC
     CO2Dpc*CO2_pc:x =0;
//%END CO2_PDE_LBC
 }   //  END CO2 PDE BC

 when (x=x.max) {	// CO2 ISF PDE RIGHT BOUNDARY CONDITIONS
//%START CO2_PDE_RBC
     CO2Dpc*CO2_pc:x = 0;
//%END CO2_PDE_RBC
 }  


//%START CO2_PDES
// PDEs:
 CO2_pc:t =   (CO2PSpc/VWpc)*(CO2_isf_in-CO2_pc)
             + RQ*(O2Gpc/VWpc)*O2_pc_in      // CO2 production calc.          
             + CO2Dpc*(CO2_pc:x:x) 
             - CFpc*(kp1*CO2_pc-(km1/KH2CO3)*HCO3m_pc_in*Hp_pc_in);
//%END CO2_PDES



// To make stand-alone model:
real alphaCO20 = 2.8472E-5 ;  // Solubility of CO2 in pc
     alphaCO2 = 2.8472E-5 ;
 real pCO2pct0 = 40 mmHg;
  
     CO2Cpct0 = pCO2pct0 * alphaCO20; // .00239;
     CO2_isf_in = 0.0;
     CO2PSpc = 0;
     O2_pc_in = 1;
     Hp_pc_in = 5.7E-8;
     HCO3m_pc_in = 0.015;
     Vpc = 0.01;      // ml/g
     Wpc = 0.72;
     VWpc = Vpc*Wpc;

real Vmax_pc = 5.0e-6 mol/min/g; // Vmax for O2 consumpt by cytochrome oxidase in PCs ;
real Km_pc = 7.0e-7 M;           // Km for O2 on cytochrome oxidase in PCs;
     O2Gpc = Vmax_pc/(Km_pc+O2_pc_in);
}
