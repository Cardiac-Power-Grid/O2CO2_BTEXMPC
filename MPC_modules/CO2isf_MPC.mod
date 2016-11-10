// MODEL NAME: CO2isf_MPC
// SHORT DESCRIPTION: Modeling Total CO2 concentration as a function of t and x in isf region.


import nsrunit;   unit conversion on;

math CO2isf_MPC {

//%START TIME_DOMAIN
realDomain t sec;t.min=0;t.max=100;t.delta=0.1;
//%END TIME_DOMAIN

//%START SPATIALDOMAIN
realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid; // BTEX space domain 
//%END SPATIALDOMAIN

//%START RXN_PARAMETERS
real CFisf = 5000;		// catalytic factor in isf due to CA
real BCisf = 24 mM;		// buffering capacity in isf
real kp1    = 0.12 1/sec,	// forward rate constant in CO2+H2O reaction
     km1    = 89 1/sec;         // backward rate constant in CO2+H2O reaction
real KH2CO3 = 5.5e-4 M;	        // EQUIL constant in H2CO3 ionization
//%END RXN_PARAMETERS

//%START PARAMS
real Fisf ml/(min*g);           // Flow of isf. Assume zero.
real Visf ml/g;                 // Volume of region
real Wisf dimensionless;	       // fractional water content of isf
real VWisf ml/g;	               // volume of water content in isf
real CO2Disf = 1e-4 cm^2/sec;   // diffusion coefficient for CO2 in isf 
real CO2PSpc ml/(g*min);       // PS for CO2 across pc Membrane
real CO2PScap ml/(g*min);       // PS for CO2 across Capillary Membrane
//%END PARAMS

//%START SOLUBILITY_PARAM
real alphaCO2(t,x) M/mmHg;      // Solubility of CO2 in isf
//%END SOLUBILITY_PARAM

//%START HCO3_HP_INPUT
real HCO3m_isf_in(t,x) M;    // HCO3m in isf
real Hp_isf_in(t,x) M;       // H+ conc in isf
//%END HCO3_HP_INPUT

//%START INPUT
real CO2_pl_in(t,x) M;     // Free CO2 from pl
real CO2_pc_in(t,x) M;     // Free CO2 from pc
//%END INPUT

//%START VARIABLES
real CO2_isf(t,x) M;      // Free CO2 in isf
//%END VARIABLES

//%START VAR_OUTPUT
real pCO2_isf(t,x) mmHg;   // convert to partial press
     pCO2_isf = CO2_isf/alphaCO2;
//%END VAR_OUTPUT

//%START VAR_INITS
real CO2Cisft0 M;
//%END VAR_INITS


 when (t=t.min) {	// CO2 ISF PDE INITIAL CONDITIONS
//%START CO2_PDE_IC
      CO2_isf = CO2Cisft0;
//%END CO2_PDE_IC
 } 

 when (x=x.min) {	// CO2 ISF PDE LEFT BOUNDARY CONDITIONS
//%START CO2_PDE_LBC
     CO2Disf*CO2_isf:x =0;
//%END CO2_PDE_LBC
 }   //  END CO2 PDE BC

 when (x=x.max) {	// CO2 ISF PDE RIGHT BOUNDARY CONDITIONS
//%START CO2_PDE_RBC
     CO2Disf*CO2_isf:x = 0;
//%END CO2_PDE_RBC
 }  


//%START CO2_PDES
// PDEs:
 CO2_isf:t =  - (CO2PScap/VWisf)*(CO2_isf-CO2_pl_in)
              - (CO2PSpc/VWisf)*(CO2_isf-CO2_pc_in) 
              + CO2Disf*(CO2_isf:x:x) 
              - CFisf*(kp1*CO2_isf-(km1/KH2CO3)*HCO3m_isf_in*Hp_isf_in);
//%END CO2_PDES



// To make stand-alone model:
 real alphaCO20 = 2.8472E-5 M/mmHg;  // Solubility of CO2 in isf
     alphaCO2 = 2.8472E-5;
 real pCO2isft0 = 40 mmHg;
     CO2Cisft0 = pCO2isft0 * alphaCO20; // .00239;
     CO2_pl_in = 0.0;
     CO2_pc_in = 0.0;
     CO2PSpc = 0;
     CO2PScap = 0;
     Hp_isf_in = 5.7E-8;
     HCO3m_isf_in = 0.015;
     Visf = 0.01;      // ml/g
     Wisf = 0.72;
     VWisf = Visf*Wisf;
     Fisf = 1;

}
