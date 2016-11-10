// MODEL NAME: HCO3PDE_MPC
// SHORT DESCRIPTION: Modeling HCO3 concentration as a function of CO2 concentration.

import nsrunit;   unit conversion on;

math HCO3PDE_MPC {

//%START TIME_DOMAIN
realDomain t sec;t.min=0;t.max=100;t.delta=0.1;
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
real K1_calc = (kp1/km1)*KH2CO3; // EQUIL constant in overall CO2+H2O reaction
//%END RXN_PARAMETERS

//%START PARAMS
real Frbc ml/(min*g);           // Flow of RBCs.
real Vrbc ml/g;                 // Volume of region
real Wrbc dimensionless;	// fractional water content of RBCs
real VWrbc ml/g;	        // volume of water content in RBCs
real HCO3mDrbc = 1e-4 cm^2/sec; // diffusion coefficient for HCO3- in RBCs
real HCO3mPSrbc ml/(min*g);     // PS for HCO3- across RBC Membrane
real HpDrbc = 1e-4 cm^2/sec;    // diffusion coefficient for H+ in RBCs
real HpPSrbc ml/(min*g);        // PS for H+ across RBC Membrane
real Rrbc dimensionless;               // Gibbs-Donnan ratio [H+]pl/[H+]rbc
//%END PARAMS

//%START RBC_INPUT
real CO2_rbc_in(t,x) M;     // Free CO2 from RBCs
//%END RBC_INPUT

//%START INPUT
real HCO3m_pl_in(t,x) M;     // HCO3m in plasma
real Hp_pl_in(t,x) M;        // H+ in plasma 
//%END INPUT

//%START FLOW_INPUT
real HCO3m_rbc_Fin(t) M;   // HCO3- coming into the rbc region (x=x.min)
real Hp_rbc_Fin(t) M;      // H+ coming into the rbc region (x=x.min)
//%END FLOW_INPUT

//%START VARIABLES
real Hp_rbc(t,x) M;       // H+ conc in RBCs
real pH_rbc(t,x) dimensionless; // pH in RBCs
     pH_rbc = -log(Hp_rbc/(1 M));
real HCO3m_rbc(t,x) M;    // HCO3- in RBCs
//%END VARIABLES

//%START VAR_OUTPUTS
real HCO3m_rbcOut(t) M;  // HCO3- leaving RBCs at x=x.max
real Hp_rbcOut(t) M;     // H+ leaving RBCs at x=x.max
//%END VAR_OUTPUTS

//%START VAR_INITS
// Inputs needed for HCO3PDE_MPC
real HpCrbct0 M;
real pH_rbct0 dimensionless;   // Initial pH in rbc
real CO2Crbct0 M;             // Initial free CO2 conc in rbc
real HCO3mCrbct0 M;           //  Initial HCO3m conc in rbc
     HpCrbct0 = 10^(-pH_rbct0)* (1 M);
     HCO3mCrbct0 = K1_calc*CO2Crbct0/HpCrbct0;
 
//%END VAR_INITS


 when (t=t.min) {	// INITIAL CONDITIONS
//%START PDE_IC
      HCO3m_rbc = HCO3mCrbct0;
      Hp_rbc = HpCrbct0;
//%END PDE_IC
  }



 when (x=x.min) {	// LEFT BOUNDARY CONDITIONS
//%START PDE_LBC
      (-Frbc*L/Vrbc)*(HCO3m_rbc-HCO3m_rbc_Fin)+HCO3mDrbc*HCO3m_rbc:x =0;
      (-Frbc*L/Vrbc)*(Hp_rbc-Hp_rbc_Fin)+HpDrbc*Hp_rbc:x =0;
//%END PDE_LBC
 }
 when (x=x.max) {	// RIGHT BOUNDARY CONDITIONS
//%START PDE_RBC
      HCO3mDrbc*HCO3m_rbc:x = 0;
      HpDrbc*Hp_rbc:x = 0;
      HCO3m_rbcOut = HCO3m_rbc;
      Hp_rbcOut = Hp_rbc; 
//%END PDE_RBC
 }  


//%START PDES
// PDEs:
 HCO3m_rbc:t = - (Frbc/Vrbc)*L*(HCO3m_rbc:x) + HCO3mDrbc*(HCO3m_rbc:x:x)
               - (HCO3mPSrbc/VWrbc)*(HCO3m_rbc-Rrbc*HCO3m_pl_in)
               + CFrbc*(kp1*CO2_rbc_in-(km1/KH2CO3)*HCO3m_rbc*Hp_rbc);
 Hp_rbc:t = - (Frbc/Vrbc)*L*(Hp_rbc:x) + HpDrbc*(Hp_rbc:x:x) 
            - (HpPSrbc/VWrbc)*(Rrbc*Hp_rbc-Hp_pl_in)
            + (2.303/BCrbc)*Hp_rbc * CFrbc*(kp1*CO2_rbc_in-(km1/KH2CO3)*HCO3m_rbc*Hp_rbc);
//%END PDES

//%START PARAMETER_ASSIGNMENT

//%END PARAMETER_ASSIGNMENT

// To make stand-alone model:
    pH_rbct0 = 7.24;   // Input initial condition
    CO2_rbc_in = .0023;  // M
    CO2Crbct0 = CO2_rbc_in(t.min,x.min);
    Vrbc = 0.01;      // ml/g
    Frbc = 1;
    Wrbc = 0.72;
    VWrbc = Vrbc*Wrbc;
    HCO3m_rbc_Fin = .01459;
    Hp_rbc_Fin = 10^(-pH_rbct0)* (1 M);
    HCO3mPSrbc =0;
    HpPSrbc = 0;
    HCO3m_pl_in = 0;
    Hp_pl_in = 0;
    Rrbc = 0.69;

}
