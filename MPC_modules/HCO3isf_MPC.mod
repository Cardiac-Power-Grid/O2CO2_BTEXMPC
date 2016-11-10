// MODEL NAME: HCO3isf_MPC
// SHORT DESCRIPTION: Modeling HCO3 concentration as a function of CO2 concentration in the isf region.

import nsrunit;   unit conversion on;

math HCO3isf_MPC {

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
real K1_calc = (kp1/km1)*KH2CO3; // EQUIL constant in overall CO2+H2O reaction
//%END RXN_PARAMETERS

//%START PARAMS
real Fisf ml/(min*g);           // Flow of isf. Assume zero.
real Visf ml/g;                 // Volume of region
real Wisf dimensionless;	       // fractional water content of isf
real VWisf ml/g;	               // volume of water content in isf
real HCO3mDisf = 1e-4 cm^2/sec; // diffusion coefficient for HCO3- in isf
real HCO3mPScap ml/(min*g);    // PS for HCO3- across Capillary Membrane
real HCO3mPSpc ml/(min*g);    // PS for HCO3- across pc Membrane
real HpDisf = 1e-4 cm^2/sec;    // diffusion coefficient for H+ in isf
real HpPScap ml/(min*g);       // PS for H+ across Capillary Membrane
real HpPSpc ml/(min*g);       // PS for H+ across pc Membrane
real Rpc dimensionless;       // Gibbs-Donnan ratio [H+]isf/[H+]pc
real Rcap dimensionless;       // Gibbs-Donnan ratio [H+]isf/[H+]pl
//%END PARAMS

//%START CO2_INPUT
real CO2_isf_in(t,x) M;      // Free CO2 in isf
//%END CO2_INPUT

//%START INPUT
real HCO3m_pl_in(t,x) M;     // HCO3m in plasma
real Hp_pl_in(t,x) M;        // H+ in plasma
real HCO3m_pc_in(t,x) M;     // HCO3m in parenchymal cell
real Hp_pc_in(t,x) M;        // H+ in parenchymal cell
//%END INPUT

//%START VARIABLES
real Hp_isf(t,x) M;       // H+ conc in isf
real pH_isf(t,x) dimensionless; // pH in isf
     pH_isf = -log(Hp_isf/(1 M));
real HCO3m_isf(t,x) M;    // HCO3m in isf
//%END VARIABLES

//%START VAR_OUTPUTS
real HCO3m_isfOut(t) M;  // HCO3m from HCO3isf_MPC at x=x.max
real Hp_isfOut(t) M;     // H+ out from HCO3isf_MPC at x=x.max
//%END VAR_OUTPUTS

//%START VAR_INITS
// Inputs needed for HCO3isf_MPC
real HpCisft0 M;
real pH_isft0 dimensionless;   // Initial pH in isf
real CO2Cisft0 M;             // Initial free CO2 conc in isf
real HCO3mCisft0 M;           // Initial HCO3m conc in isf
     HpCisft0 = 10^(-pH_isft0)* (1 M);
     HCO3mCisft0 = K1_calc*CO2Cisft0/HpCisft0;
//%END VAR_INITS


 when (t=t.min) {	// INITIAL CONDITIONS
//%START PDE_IC
      HCO3m_isf = HCO3mCisft0;
      Hp_isf = HpCisft0;
//%END PDE_IC
  }



 when (x=x.min) {	// LEFT BOUNDARY CONDITIONS
//%START PDE_LBC
      HCO3mDisf*HCO3m_isf:x =0;
      HpDisf*Hp_isf:x =0;
//%END PDE_LBC
 }
 when (x=x.max) {	// RIGHT BOUNDARY CONDITIONS
//%START PDE_RBC
      HCO3mDisf*HCO3m_isf:x = 0;
      HpDisf*Hp_isf:x = 0;
      HCO3m_isfOut = HCO3m_isf;
      Hp_isfOut = Hp_isf; 
//%END PDE_RBC
 }  


//%START PDES
// PDEs:
 HCO3m_isf:t = HCO3mDisf*(HCO3m_isf:x:x)
               + (HCO3mPScap/VWisf)*(Rcap*HCO3m_pl_in-HCO3m_isf)
               - (HCO3mPSpc/VWisf)*(Rpc*HCO3m_isf-HCO3m_pc_in)
               + CFisf*(kp1*CO2_isf_in-(km1/KH2CO3)*HCO3m_isf*Hp_isf);
 Hp_isf:t = HpDisf*(Hp_isf:x:x) 
            + (HpPScap/VWisf)*(Hp_pl_in-Rcap*Hp_isf)
            - (HpPSpc/VWisf)*(Hp_isf-Rpc*Hp_pc_in)
            + (2.303/BCisf)*Hp_isf * CFisf*(kp1*CO2_isf_in-(km1/KH2CO3)*HCO3m_isf*Hp_isf);
//%END PDES

// To make stand-alone model:

    Rcap = 0.63;
    Rpc = 0.79;
    Visf = 0.01;      // ml/g
    Fisf = 1;
    Wisf = 0.72;
    VWisf = Visf*Wisf;
    CO2Cisft0 = CO2_isf_in(t.min,x.min);
    pH_isft0 = 7.24;   // Input initial condition
    CO2_isf_in = .0023;  // M
    HCO3mPSpc =0;
    HpPSpc = 0;
    HCO3mPScap =0;
    HpPScap = 0;
    HCO3m_pl_in = 0;
    Hp_pl_in = 0;
    HCO3m_pc_in = 0;
    Hp_pc_in = 0;


}

