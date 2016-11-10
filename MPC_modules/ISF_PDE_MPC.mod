// MODEL NAME: ISF_PDE_MPC
// SHORT DESCRIPTION: Simple model of Intersitial region (ISF) that calculates the pO2 and pCO2
// in ISF based on Total O2, CO2, H+, and HCO3. Takes into account temperature
// dependence of O2/CO2 solubility. Uses isf module with renaming.


//%START ISF_MODULE_INFO
// * Interstitial fluid (ISF) Variables *:
// pCO2_isf(t,x) mmHg; // free pCO2 in ISF
// CO2_isf(t,x) M;     // free CO2 in ISF
// O2_isf(t,x) M;      // free O2 in ISF
// pO2_isf(t,x) mmHg;  // free pO2 in ISF
// Hp_isf(t,x) M;      // H+ conc in ISF
// pH_isf(t,x) ;       // pH in ISF
// HCO3m_isf(t,x) M;   // HCO3- in ISF
// TempSys_isf(t,x) K; //  Temperature in isf
// alphaO2Sys_isf(t,x) M/mmHg // O2 Solubility coefficient
// alphaCO2Sys_isf(t,x) M/mmHg // CO2 Solubility coefficient

// * ISF Input variables *:
// CO2_pl(t,x) M;     // Free CO2 in plasma (pl)
// O2_pl(t,x) M;      // Free O2 in pl
// HCO3m_pl(t,x) M;   // HCO3m in pl
// Hp_pl(t,x) M;      // H+ in pl
// TempSys_cap(t,x) K;// Temperature in capillary

// TempSys_pc(t,x) K; // Temp in pc (parenchymal cell)
// CO2_pc(t,x) M;     // Free CO2 in pc 
// O2_pc(t,x) M;      // Free O2 in pc
// HCO3m_pc(t,x) M;   // HCO3m in pc
// Hp_pc(t,x) M;      // H+ in pc 

// * ISF Flow input variables (Conc of substrate coming into pipe at x=x.min):
//   None, assume no flow in ISF.
//%END ISF_MODULE_INFO

import nsrunit;
unit conversion on;
 unit celsius = fundamental;

math ISF_pde_MPC {
realDomain t sec;t.min=0;t.max=100;t.delta=0.1;
realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid; // BTEX space domain 

//%START ISF_INPUTS
real HCO3m_pl(t,x) M;     // HCO3m in plasma
real Hp_pl(t,x) M;        // H+ in plasma
real HCO3m_pc(t,x) M;     // HCO3m in parenchymal cell
real Hp_pc(t,x) M;        // H+ in parenchymal cell
real CO2_pl(t,x) M;     // Free CO2 from pl
real CO2_pc(t,x) M;     // Free CO2 from pc
real O2_pl(t,x) M;     // Free O2 in plasma
real O2_pc(t,x) M;     // Free O2 in pc
//%END ISF_INPUTS
real TempSys_cap(t,x) K;
real TempSys_pc(t,x) K;
real MolV  = 22.414 L/mol;     // Liters of gas per mole at std temp.
real specificHeat = 1 kilocal/(g*K*min);   // Heat required to raise mass one K per min
real ThermCoeff = 0.0001 g/sec;  // Thermal coefficient
real Rconst = 62.36358 L*mmHg/K/mole;  // ideal gas const
real TempExp K;       // Initial temperature
//%START ISF_VARS
real Hp_isf(t,x) M;       // H+ conc in isf
real pH_isf(t,x) dimensionless; // pH in isf
     pH_isf = -log(Hp_isf/(1 M));
real HCO3m_isf(t,x) M;    // HCO3m in isf
real CO2_isf(t,x) M;      // Free CO2 in isf
real O2_isf(t,x) M;      // O2 in isf
real TempSys_isf(t,x) K;      //  Temperature in isf
//%END ISF_VARS
// New solubility calcs..... handle up to 40C, for solubility in plasma(O2)/saline(CO2):
// Coefficients to fit solubility curves:
real  alphaO20calc  = 0.0082 	mL/mL/atm,  
      alphaO201calc  = 0.0331 	mL/mL/atm,
      alphaCO20calc = 1.526	mL/mL/atm,  // saline
      alphaCO201calc = 0.132	mL/mL/atm,
      O2k1	= -0.0061 	1/celsius,
      O2k2      = 0.0292	1/celsius,
      CO2k1	= 0.0385	1/celsius,
      CO2k2	= -0.0105 	1/celsius;
//%START SOLUBILITY_COEFF
real alphaCO2Sys_isf(t,x) M/mmHg;      // Solubility of CO2 in isf
real alphaO2Sys_isf(t,x) M/mmHg;      // Solubility of O2 in isf
real O2solSys_isf(t,x)	mL/mL/atm, // ml O2 per ml plasma per atm
     CO2solSys_isf(t,x)	mL/mL/atm;  
 // ANALYTIC SOLUTION
real TempC_Sys_isf(t,x) celsius;
     TempC_Sys_isf = (TempSys_isf-273.15)* (1 celsius)/(1 K);  // convert kelvin to celsius 
     O2solSys_isf = alphaO20calc*exp(-O2k1*TempC_Sys_isf)+alphaO201calc*exp(-O2k2*TempC_Sys_isf);  
     CO2solSys_isf = alphaCO20calc*exp(-CO2k1*TempC_Sys_isf)+alphaCO201calc*exp(-CO2k2*TempC_Sys_isf);
// convert solubilities to M/mmHg;
     alphaO2Sys_isf =  (1 atm)*O2solSys_isf/Rconst/TempSys_isf;
     alphaCO2Sys_isf = (1 atm)* CO2solSys_isf/Rconst/TempSys_isf;
// **** END of solubility calcs........
//%END SOLUBILITY_COEFF
//%START ISF_REGION
// ISF region:
real CFisf = 5000;		// catalytic factor in isf due to CA
real BCisf = 24 mM;		// buffering capacity in isf
real kp1    = 0.12 1/sec,	// forward rate constant in CO2+H2O reaction
     km1    = 89 1/sec;         // backward rate constant in CO2+H2O reaction
real KH2CO3 = 5.5e-4 M;	        // EQUIL constant in H2CO3 ionization
real K1_calc = (kp1/km1)*KH2CO3; // EQUIL constant in overall CO2+H2O reaction
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
real CO2Disf = 1e-4 cm^2/sec;   // diffusion coefficient for CO2 in isf 
real CO2PSpc ml/(g*min);       // PS for CO2 across pc Membrane
real CO2PScap ml/(g*min);       // PS for CO2 across Capillary Membrane
real O2Disf = 1e-4 cm^2/sec;    // diffusion coefficient for O2 in RBCs 
real O2PSpc ml/(g*min);       // PS for O2 across PC Membrane
real O2PScap ml/(g*min);       // PS for O2 across Capillary Membrane
real TDpc(t,x) s^-1;      // Temp diffusion across pc (muscle) membrane
real TDcap(t,x) s^-1;      // Temp diffusion across cap (capillary) membrane
real TCisf cm^2/sec;       // Thermal conductivity (Diffusion) in isf
     TCisf =  ThermCoeff*Visf/L;
real HCO3m_isfOut(t) M;  // HCO3m from HCO3isf_MPC at x=x.max
real Hp_isfOut(t) M;     // H+ out from HCO3isf_MPC at x=x.max
real pCO2_isf(t,x) mmHg;   // convert to partial press
     pCO2_isf = CO2_isf/alphaCO2Sys_isf;
real pO2_isf(t,x) mmHg;   // convert to partial press
     pO2_isf = O2_isf/alphaO2Sys_isf;
// Inputs needed for HCO3isf_MPC
real HpCisft0 M;
real pH_isft0 dimensionless;   // Initial pH in isf
real CO2Cisft0 M;             // Initial free CO2 conc in isf
real HCO3mCisft0 M;           // Initial HCO3m conc in isf
     HpCisft0 = 10^(-pH_isft0)* (1 M);
     HCO3mCisft0 = K1_calc*CO2Cisft0/HpCisft0;
real CO2Cisft0 M;
real O2Cisft0 M;
 when (t=t.min) {	// ISF INITIAL CONDITIONS
      HCO3m_isf = HCO3mCisft0;
      Hp_isf = HpCisft0;
      CO2_isf = CO2Cisft0;
      O2_isf = O2Cisft0;  // O2 PDE INITIAL CONDITION
  TempSys_isf = TempExp;
 } // end ISF ICs
 when (x=x.min) {	// LEFT ISF BOUNDARY CONDITIONS
      HCO3mDisf*HCO3m_isf:x =0;
      HpDisf*Hp_isf:x =0;
     CO2Disf*CO2_isf:x =0;
     O2Disf*O2_isf:x =0;
  TCisf*TempSys_isf:x = 0;  
 } // end isf left BCs
 when (x=x.max) {	// RIGHT ISF BOUNDARY CONDITIONS
      HCO3mDisf*HCO3m_isf:x = 0;
      HpDisf*Hp_isf:x = 0;
      HCO3m_isfOut = HCO3m_isf;
      Hp_isfOut = Hp_isf; 
     CO2Disf*CO2_isf:x = 0;
     O2Disf*O2_isf:x = 0;
  TempSys_isf:x =0;   
 } // end isf right BCs
// *** ISF PDES:
// PDEs:
 CO2_isf:t =  - (CO2PScap/VWisf)*(CO2_isf-CO2_pl)
              - (CO2PSpc/VWisf)*(CO2_isf-CO2_pc) 
              + CO2Disf*(CO2_isf:x:x) 
              - CFisf*(kp1*CO2_isf-(km1/KH2CO3)*HCO3m_isf*Hp_isf);
 HCO3m_isf:t = HCO3mDisf*(HCO3m_isf:x:x)
               + (HCO3mPScap/VWisf)*(Rcap*HCO3m_pl-HCO3m_isf)
               - (HCO3mPSpc/VWisf)*(Rpc*HCO3m_isf-HCO3m_pc)
               + CFisf*(kp1*CO2_isf-(km1/KH2CO3)*HCO3m_isf*Hp_isf);
 Hp_isf:t = HpDisf*(Hp_isf:x:x) 
            + (HpPScap/VWisf)*(Hp_pl-Rcap*Hp_isf)
            - (HpPSpc/VWisf)*(Hp_isf-Rpc*Hp_pc)
            + (2.303/BCisf)*Hp_isf * CFisf*(kp1*CO2_isf-(km1/KH2CO3)*HCO3m_isf*Hp_isf);
// O2 PDE:
 O2_isf:t = + (O2PScap/VWisf)*(O2_pl-O2_isf)
           - (O2PSpc/VWisf)*(O2_isf-O2_pc)
           + O2Disf*(O2_isf:x:x);             
TempSys_isf:t =  TDpc*(TempSys_pc-TempSys_isf) 
                - TDcap*(TempSys_isf-TempSys_cap) 
                + TCisf*(TempSys_isf):x:x ; 
// ------ End of isf pdes
//%END ISF_REGION
//%START O2_CO2_SOLUBILITY_INIT
real alphaO20 = 1.1148E-6 M/mmHg;   // Initial solubility of O2 in isf
real alphaCO20 = 2.8472E-5 M/mmHg;  // Solubility of CO2 in isf
//%END O2_CO2_SOLUBILITY_INIT
//%START ISF_INIT
     pH_isft0 = 7.24;   // Input initial condition
real pO2isft0 = 90 mmHg;
real pCO2isft0 = 40 mmHg;
     O2Cisft0 =  pO2isft0 * alphaO20;
     CO2Cisft0 = pCO2isft0 * alphaCO20; 
//%END ISF_INIT
//%START PS_VALUES
    O2PSpc = 0;
    O2PScap = 0;
    CO2PSpc = 0;
    CO2PScap = 0;
    HCO3mPSpc =0;
    HCO3mPScap =0;
    HpPSpc = 0;
    HpPScap = 0;
//%END PS_VALUES
//%START ISF_ASSIGN_INPUTS
// Inputs for ISF modules:
   // alphaO2Sys_isf = alphaO20;   // Solubility of O2 in isf
   // alphaCO2Sys_isf = alphaCO20; // Solubility of O2 in isf
    O2_pl =0;  // M
    CO2_pl =0;
    HCO3m_pl =0;
    Hp_pl = 0;  
    O2_pc =0;
    CO2_pc=0;
    Hp_pc =0;
    HCO3m_pc =0;
// End of Inputs for ISF modules.
// ** No isf Input at x.min for PDEs, assume Flow in isf is zero.
//%END ISF_ASSIGN_INPUTS
//%START SYS_PARAMS
// System parameters:
real  H = 0.45 dimensionless;   // large vessel hematocrit (unitless)
real  Fb = 1.0 ml/(min*g);	// total blood flow  (blood flow per gram tissue)
//%END SYS_PARAMS
//%START ISF_SYS_PARAMS
real  Fisf = 0;	     // flow of isf, ml/(min*g): Assume zero, Do not use LSFEA solver if substrate const over all x.
real  Visf = 0.2;	        // anatomical volume of isf
real  Wisf = 0.92;              // fractional water content of isf
      VWisf = Visf*Wisf;
//%END ISF_SYS_PARAMS
//%START CAP_SYS_PARAMS
      Rcap = 0.63;
      TDcap =0;                 // Temp (heat) diffusion across cap (capillary) membrane
      TempSys_cap = 310;
//%END CAP_SYS_PARAMS
//%START PC_SYS_PARAMS
      Rpc = 0.79;               // Gibbs-Donnan ratio [H+]isf/[H+]pc
      TDpc = 0;
      TempSys_pc = 310;
//%END PC_SYS_PARAMS
      TempExp = 310;
 
} // End of isf pde model
// This MML file generated from ISF_PDE_MPC.mpc using MPC v1.01.
