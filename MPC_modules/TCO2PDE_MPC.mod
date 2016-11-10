// MODEL NAME: TCO2PDE_MPC
// SHORT DESCRIPTION: Modeling Total CO2 concentration as a function of t and x.


import nsrunit;   unit conversion on;

math TCO2PDE_MPC {

//%START TIME_DOMAIN
realDomain t sec;t.min=0;t.max=10;t.delta=0.1;
//%END TIME_DOMAIN

//%START SPATIALDOMAIN
realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid; // BTEX space domain 
//%END SPATIALDOMAIN

//%START RXN_PARAMETERS
real CFrbc = 13000;		// catalytic factor in RBC due to CA
real BCrbc = 54 mM;		// buffering capacity in RBC
real kp1    = 0.12 1/sec,	// forward rate constant in CO2+H2O reaction
     km1    = 89 1/sec;         // backward rate constant in CO2+H2O reaction
real KH2CO3 = 5.5e-4 M;	        // EQUIL constant in H2CO3 ionization
//%END RXN_PARAMETERS

//%START PARAMS
real Frbc ml/(min*g);           // Flow of RBCs.
real Vrbc ml/g;                 // Volume of region
real Wrbc dimensionless;	// fractional water content of RBCs
real VWrbc ml/g;	        // volume of water content in RBCs
real CO2Drbc = 1e-4 cm^2/sec;   // diffusion coefficient for CO2 in RBCs 
real CO2PSrbc ml/(g*min);       // PS for CO2 across RBC Membrane
//%END PARAMS

//%START FLOW_INPUT
real TCO2_rbc_Fin(t) M;   // TCO2 coming into the rbc region (x.min)
//%END FLOW_INPUT

//%START INPUT
real CO2_rbc_in(t,x) M;     // Free CO2 from RBCs
real CO2_pl_in(t,x) M;      // Free CO2 in plasma
real HCO3m_rbc_in(t,x) M;   // HCO3m in RBCs
real Hp_rbc_in(t,x) M;      // H+ conc in RBCs
//%END INPUT

//%START VARIABLES
real TCO2_rbc(t,x) M;    // Total CO2 in RBCs
//%END VARIABLES

//%START VAR_OUTPUT
real TCO2_rbcOut(t) M;

//%END VAR_OUTPUT

//%START VAR_INITS
real TCO2Crbct0 M;
//%END VAR_INITS


 when (t=t.min) {	// TCO2 PDE INITIAL CONDITIONS
//%START TCO2_PDE_IC
      TCO2_rbc = TCO2Crbct0;
//%END TCO2_PDE_IC
 } 



 when (x=x.min) {	// TCO2 PDE LEFT BOUNDARY CONDITIONS
//%START TCO2_PDE_LBC
     (-Frbc*L/Vrbc)*(TCO2_rbc-TCO2_rbc_Fin)+CO2Drbc*TCO2_rbc:x =0;
//%END TCO2_PDE_LBC
 }   //  END TCO2 PDE BC

 when (x=x.max) {	// TCO2 PDE RIGHT BOUNDARY CONDITIONS
//%START TCO2_PDE_RBC
     CO2Drbc*TCO2_rbc:x = 0;
     TCO2_rbcOut = TCO2_rbc;
//%END TCO2_PDE_RBC
 }  


//%START TCO2_PDES
// PDEs:
 TCO2_rbc:t = - (Frbc/Vrbc)*L*(TCO2_rbc:x) 
                   - (CO2PSrbc/VWrbc)*(CO2_rbc_in-CO2_pl_in)
                   + CO2Drbc*(TCO2_rbc:x:x) 
                   - CFrbc*(kp1*CO2_rbc_in-(km1/KH2CO3)*HCO3m_rbc_in*Hp_rbc_in);
//%END TCO2_PDES



// To make stand-alone model:
     TCO2Crbct0 =.00239;
real Hb_rbc = .0052 M;   // Hb conc in rbc
real SHbCO2(t,x) = 0.131;
     CO2_rbc_in = TCO2_rbc-4*Hb_rbc*SHbCO2;  // M  
     CO2_pl_in = 0.0;
     CO2PSrbc = 0;
     Hp_rbc_in = 5.7E-8;
     HCO3m_rbc_in = 0.015;
     TCO2_rbc_Fin = TCO2Crbct0;
     Vrbc = 0.01;      // ml/g
     Wrbc = 0.72;
     VWrbc = Vrbc*Wrbc;
     Frbc = 1;

}
