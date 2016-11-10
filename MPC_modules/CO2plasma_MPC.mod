// MODEL NAME: CO2plasma_MPC
// SHORT DESCRIPTION: Modeling Total CO2 concentration as a function of t and x in plasma region..


import nsrunit;   unit conversion on;

math CO2plasma_MPC {

//%START TIME_DOMAIN
realDomain t sec;t.min=0;t.max=100;t.delta=0.1;
//%END TIME_DOMAIN

//%START SPATIALDOMAIN
realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid; // BTEX space domain 
//%END SPATIALDOMAIN

//%START RXN_PARAMETERS
real CFpl = 100;		// catalytic factor in plasma due to CA
real BCpl = 6 mM;		// buffering capacity in plasma
real kp1    = 0.12 1/sec,	// forward rate constant in CO2+H2O reaction
     km1    = 89 1/sec;         // backward rate constant in CO2+H2O reaction
real KH2CO3 = 5.5e-4 M;	        // EQUIL constant in H2CO3 ionization
//%END RXN_PARAMETERS

//%START PARAMS
real Fpl ml/(min*g);           // Flow of pls.
real Vpl ml/g;                 // Volume of region
real Wpl dimensionless;	       // fractional water content of plasma
real VWpl ml/g;	               // volume of water content in plasma
real CO2Dpl = 1e-4 cm^2/sec;   // diffusion coefficient for CO2 in plasma 
real CO2PSrbc ml/(g*min);       // PS for CO2 across RBC Membrane
real CO2PScap ml/(g*min);       // PS for CO2 across Capillary Membrane
//%END PARAMS

//%START FLOW_INPUT
real CO2_pl_Fin(t) M;   // CO2 coming into the plasma region (x.min)
//%END FLOW_INPUT

//%START INPUT
real CO2_rbc_in(t,x) M;     // Free CO2 from RBCs
real CO2_isf_in(t,x) M;      // CO2 in isf
//%END INPUT

//%START HCO3_HP_INPUT
real HCO3m_pl_in(t,x) M;     // HCO3m in plasma
real Hp_pl_in(t,x) M;        // H+ in plasma
//%END HCO3_HP_INPUT

//%START VARIABLES
real CO2_pl(t,x) M;      // CO2 in plasma
//%END VARIABLES

//%START VAR_OUTPUT
real CO2_plOut(t) M;    // CO2 coming out of plasma region at x.max

//%END VAR_OUTPUT

//%START VAR_INITS
real CO2Cplt0 M;
//%END VAR_INITS


 when (t=t.min) {	// CO2 PDE INITIAL CONDITIONS
//%START CO2_PDE_IC
      CO2_pl = CO2Cplt0;
//%END CO2_PDE_IC
 } 

 when (x=x.min) {	// CO2 PDE LEFT BOUNDARY CONDITIONS
//%START CO2_PDE_LBC
     (-Fpl*L/Vpl)*(CO2_pl-CO2_pl_Fin)+CO2Dpl*CO2_pl:x =0;
//%END CO2_PDE_LBC
 }   //  END CO2 PDE BC

 when (x=x.max) {	// CO2 PDE RIGHT BOUNDARY CONDITIONS
//%START CO2_PDE_RBC
     CO2Dpl*CO2_pl:x = 0;
     CO2_plOut = CO2_pl;
//%END CO2_PDE_RBC
 }  


//%START CO2_PDES
// PDEs:
 CO2_pl:t = - (Fpl/Vpl)*L*(CO2_pl:x) 
                   - (CO2PSrbc/VWpl)*(CO2_pl-CO2_rbc_in)
                   - (CO2PScap/VWpl)*(CO2_pl-CO2_isf_in) 
                   + CO2Dpl*(CO2_pl:x:x) 
                   - CFpl*(kp1*CO2_pl-(km1/KH2CO3)*HCO3m_pl_in*Hp_pl_in);
//%END CO2_PDES



// To make stand-alone model:
 real alphaCO2 = 2.8472E-5 M/mmHg;  // Solubility of CO2 in plasma
 real pCO2_pl(t,x) mmHg;
 real pCO2plt0 = 40 mmHg;
     pCO2_pl = CO2_pl/alphaCO2;
     CO2Cplt0 = pCO2plt0 * alphaCO2; // .00239;
     CO2_rbc_in = 0.0;
     CO2_isf_in = 0.0;
     CO2PSrbc = 0;
     CO2PScap = 0;
     Hp_pl_in = 5.7E-8;
     HCO3m_pl_in = 0.015;
     CO2_pl_Fin = CO2Cplt0;
     Vpl = 0.01;      // ml/g
     Wpl = 0.72;
     VWpl = Vpl*Wpl;
     Fpl = 1;

}
