// MODEL NAME: HCO3pc_MPC
// SHORT DESCRIPTION: Modeling HCO3 concentration as a function of CO2 concentration in the pc region.

import nsrunit;   unit conversion on;

math HCO3pc_MPC {

//%START TIME_DOMAIN
realDomain t sec;t.min=0;t.max=20;t.delta=0.1;
//%END TIME_DOMAIN

//%START SPATIALDOMAIN
realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid; // BTEX space domain 
//%END SPATIALDOMAIN

//%START RXN_PARAMETERS
real CFpc = 10000;              // catalytic factor in pc due to CA
real BCpc = 45 mM;              // buffering capacity in pc
real kp1    = 0.12 1/sec,	// forward rate constant in CO2+H2O reaction
     km1    = 89 1/sec;         // backward rate constant in CO2+H2O reaction
real KH2CO3 = 5.5e-4 M;	        // EQUIL constant in H2CO3 ionization
real K1_calc = (kp1/km1)*KH2CO3; // EQUIL constant in overall CO2+H2O reaction
//%END RXN_PARAMETERS

//%START PARAMS
//real Fpc ml/(min*g);           // Flow of isf. Assume zero.
real Vpc ml/g;                 // Volume of region
real Wpc dimensionless;	       // fractional water content of pc
real VWpc ml/g;	               // volume of water content in pc
real HCO3mDpc = 1e-4 cm^2/sec; // diffusion coefficient for HCO3- in pc
real HCO3mPSpc ml/(min*g);    // PS for HCO3- across pc Membrane
real HpDpc = 1e-4 cm^2/sec;    // diffusion coefficient for H+ in pc
real HpPSpc ml/(min*g);       // PS for H+ across pc Membrane
real Rpc dimensionless;       // Gibbs-Donnan ratio [H+]isf/[H+]pc
//%END PARAMS

//%START CO2_INPUT
real CO2_pc_in(t,x) M;      // CO2 in pc
//%END CO2_INPUT

//%START INPUT
real HCO3m_isf_in(t,x) M;    // HCO3m in isf
real Hp_isf_in(t,x) M;       // H+ conc in isf
//%END INPUT

//%START VARIABLES
real Hp_pc(t,x) M;       // H+ conc in pc
real pH_pc(t,x) dimensionless; // pH in pc
     pH_pc = -log(Hp_pc/(1 M));
real HCO3m_pc(t,x) M;    // HCO3m in pc
//%END VARIABLES

//%START VAR_OUTPUTS

//%END VAR_OUTPUTS

//%START VAR_INITS
// Inputs needed for HCO3pc_MPC
real HpCpct0 M;
real pH_pct0 dimensionless;   // Initial pH in pc
real CO2Cpct0 M;        // Initial free CO2 conc in pc
real HCO3mCpct0 M;      // Initial HCO3m conc in pc
     HpCpct0 = 10^(-pH_pct0)* (1 M);
     HCO3mCpct0 = K1_calc*CO2Cpct0/HpCpct0; 
//%END VAR_INITS


 when (t=t.min) {	// INITIAL CONDITIONS
//%START PDE_IC
      HCO3m_pc = HCO3mCpct0;
      Hp_pc = HpCpct0;
//%END PDE_IC
  }



 when (x=x.min) {	// LEFT BOUNDARY CONDITIONS
//%START PDE_LBC
      HCO3mDpc*HCO3m_pc:x =0;
      HpDpc*Hp_pc:x =0;
//%END PDE_LBC
 }
 when (x=x.max) {	// RIGHT BOUNDARY CONDITIONS
//%START PDE_RBC
      HCO3mDpc*HCO3m_pc:x = 0;
      HpDpc*Hp_pc:x = 0;
//%END PDE_RBC
 }  


//%START PDES
// PDEs:
 HCO3m_pc:t = HCO3mDpc*(HCO3m_pc:x:x)
               + (HCO3mPSpc/VWpc)*(Rpc*HCO3m_isf_in-HCO3m_pc)
               + CFpc*(kp1*CO2_pc_in-(km1/KH2CO3)*HCO3m_pc*Hp_pc);
 Hp_pc:t = HpDpc*(Hp_pc:x:x) 
            + (HpPSpc/VWpc)*(Hp_isf_in-Rpc*Hp_pc)
            + (2.303/BCpc)*Hp_pc * CFpc*(kp1*CO2_pc_in-(km1/KH2CO3)*HCO3m_pc*Hp_pc);
//%END PDES

// To make stand-alone model:
    pH_pct0 = 7.24;   // Input initial condition
    CO2_pc_in = .0023;  // M
    CO2Cpct0 = CO2_pc_in(t.min,x.min);
    Vpc = 0.01;      // ml/g
    Wpc = 0.72;
    VWpc = Vpc*Wpc;

    HCO3mPSpc =0;
    HpPSpc = 0;
    HCO3m_isf_in = 0;
    Hp_isf_in =0;
 
    Rpc = 0.79;

}

