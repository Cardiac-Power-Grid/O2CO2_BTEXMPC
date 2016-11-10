// MODEL NAME: HCO3plasma_MPC
// SHORT DESCRIPTION: Modeling HCO3 concentration as a function of CO2 concentration in the plasma region.

import nsrunit;   unit conversion on;

math HCO3plasma_MPC {

//%START TIME_DOMAIN
realDomain t sec;t.min=0;t.max=10;t.delta=0.1;
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
real K1_calc = (kp1/km1)*KH2CO3; // EQUIL constant in overall CO2+H2O reaction
//%END RXN_PARAMETERS

//%START PARAMS
real Fpl ml/(min*g);           // Flow of pls.
real Vpl ml/g;                 // Volume of region
real Wpl dimensionless;	       // fractional water content of plasma
real VWpl ml/g;	               // volume of water content in plasma
real HCO3mDpl = 1e-4 cm^2/sec; // diffusion coefficient for HCO3- in plasma
real HCO3mPScap ml/(min*g);    // PS for HCO3- across Capillary Membrane
real HCO3mPSrbc ml/(min*g);    // PS for HCO3- across RBC Membrane
real HpDpl = 1e-4 cm^2/sec;    // diffusion coefficient for H+ in plasma
real HpPScap ml/(min*g);       // PS for H+ across Capillary Membrane
real HpPSrbc ml/(min*g);       // PS for H+ across RBC Membrane
real Rcap dimensionless;       // Gibbs-Donnan ratio [H+]pl/[H+]isf
real Rrbc dimensionless;       // Gibbs-Donnan ratio [H+]pl/[H+]rbc
//%END PARAMS

//%START CO2_PLASMA_INPUT
real CO2_pl_in(t,x) M;     // Free CO2 from pl
//%END CO2_PLASMA_INPUT

//%START INPUT
real HCO3m_rbc_in(t,x) M;   // HCO3m in RBCs
real Hp_rbc_in(t,x) M;      // H+ conc in RBCs
real HCO3m_isf_in(t,x) M;    // HCO3m in isf
real Hp_isf_in(t,x) M;       // H+ conc in isf
//%END INPUT 

//%START FLOW_INPUT
real HCO3m_pl_Fin(t) M;   // HCO3m coming into the pl region (x=x.min)
real Hp_pl_Fin(t) M;      // H+ coming into the pl region (x=x.min)

//%END FLOW_INPUT

//%START VARIABLES
real Hp_pl(t,x) M;       // H+ conc in plasma
real pH_pl(t,x) dimensionless; // pH in plasma
     pH_pl = -log(Hp_pl/(1 M));
real HCO3m_pl(t,x) M;    // HCO3m in plasma
//%END VARIABLES

//%START VAR_OUTPUTS
real HCO3m_plOut(t) M;  // HCO3m from HCO3plasma_MPC
real Hp_plOut(t) M;     // H+ from HCO3plasma_MPC
//%END VAR_OUTPUTS

//%START VAR_INITS
// Inputs needed for HCO3plasma_MPC
real HpCplt0 M;
real pH_plt0 dimensionless;   // Initial pH in pl
real CO2Cplt0 M;             // Initial free CO2 conc in pl
real HCO3mCplt0 M;           // Initial HCO3m conc in pl
     HpCplt0 = 10^(-pH_plt0)* (1 M);
     HCO3mCplt0 = K1_calc*CO2Cplt0/HpCplt0;

//%END VAR_INITS


 when (t=t.min) {	// INITIAL CONDITIONS
//%START PDE_IC
      HCO3m_pl = HCO3mCplt0;
      Hp_pl = HpCplt0;
//%END PDE_IC
  }



 when (x=x.min) {	// LEFT BOUNDARY CONDITIONS
//%START PDE_LBC
      (-Fpl*L/Vpl)*(HCO3m_pl-HCO3m_pl_Fin)+HCO3mDpl*HCO3m_pl:x =0;
      (-Fpl*L/Vpl)*(Hp_pl-Hp_pl_Fin)+HpDpl*Hp_pl:x =0;
//%END PDE_LBC
 }
 when (x=x.max) {	// RIGHT BOUNDARY CONDITIONS
//%START PDE_RBC
      HCO3mDpl*HCO3m_pl:x = 0;
      HpDpl*Hp_pl:x = 0;
      HCO3m_plOut = HCO3m_pl;
      Hp_plOut = Hp_pl; 
//%END PDE_RBC
 }  


//%START PDES
// PDEs:
 HCO3m_pl:t = - (Fpl/Vpl)*L*(HCO3m_pl:x) + HCO3mDpl*(HCO3m_pl:x:x)
               + (HCO3mPSrbc/VWpl)*(HCO3m_rbc_in-Rrbc*HCO3m_pl)
               - (HCO3mPScap/VWpl)*(Rcap*HCO3m_pl-HCO3m_isf_in)
               + CFpl*(kp1*CO2_pl_in-(km1/KH2CO3)*HCO3m_pl*Hp_pl);
 Hp_pl:t = - (Fpl/Vpl)*L*(Hp_pl:x) + HpDpl*(Hp_pl:x:x) 
            + (HpPSrbc/VWpl)*(Rrbc*Hp_rbc_in-Hp_pl)
            - (HpPScap/VWpl)*(Hp_pl-Rcap*Hp_isf_in)
            + (2.303/BCpl)*Hp_pl * CFpl*(kp1*CO2_pl_in-(km1/KH2CO3)*HCO3m_pl*Hp_pl);
//%END PDES

// To make stand-alone model:
    pH_plt0 = 7.24;   // Input initial condition
    CO2_pl_in = .0023;  // M
    CO2Cplt0 = CO2_pl_in(t.min,x.min);
    Vpl = 0.01;      // ml/g
    Fpl = 1;
    Wpl = 0.72;
    VWpl = Vpl*Wpl;
    HCO3m_pl_Fin = .01459;
    Hp_pl_Fin = 10^(-pH_plt0)* (1 M);
    HCO3mPSrbc =0;
    HpPSrbc = 0;
    HCO3mPScap =0;
    HpPScap = 0;
    HCO3m_rbc_in = 0;
    Hp_rbc_in = 0;
    HCO3m_isf_in = 0;
    Hp_isf_in = 0;
    Rcap = 0.63;
    Rrbc = 0.69;

}
