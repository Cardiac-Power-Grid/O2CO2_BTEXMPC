// MODEL NAME: PC_PDE_MPC
// SHORT DESCRIPTION: Simple model of Parenchymal cell region (PC) that calculates the pO2 and pCO2
// in PC based on Total O2, CO2, H+, and HCO3. Takes into account temperature
// dependence of O2/CO2 solubility.


//%START PC_MODULE_INFO
// * Parenchymal cell (PC) Variables *:
// TempSys_pc(t,x) K; // Temp in pc
// CO2_pc(t,x) M;     // Free CO2 in pc 
// TO2_pc(t,x) M;     // total O2 in pc
// O2_pc(t,x) M;      // Free O2 in pc
// HCO3m_pc(t,x) M;   // HCO3m in pc
// Hp_pc(t,x) M;      // H+ in pc 

// * PC Input variables *:
// O2Gpc(t,x) ml/min/g;// gulosity for O2 consumption in PCs
// pCO2_isf(t,x) mmHg; // free pCO2 in ISF
// CO2_isf(t,x) M;     // free CO2 in ISF
// O2_isf(t,x) M;      // free O2 in ISF
// pO2_isf(t,x) mmHg;  // free pO2 in rbc
// Hp_isf(t,x) M;      // H+ conc in RBCs
// HCO3m_isf(t,x) M;   // HCO3- in RBCs
// TempSys_isf(t,x) K; //  Temperature in isf

// * PC Flow input variables (Conc of substrate coming into pipe at x=x.min):
//   None, assume no flow in PC.
//%END PC_MODULE_INFO

import nsrunit;
unit conversion on;
 unit celsius = fundamental;

math PC_pde_MPC {
realDomain t sec;t.min=0;t.max=20;t.delta=0.1;
realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid; // BTEX space domain 

// ideal gas const
real Rconst = 62.36358 L*mmHg/K/mole;  // ideal gas const
//%START PC_VARS
real Hp_pc(t,x) M;       // H+ conc in pc
real HCO3m_pc(t,x) M;    // HCO3m in pc
real CO2_pc(t,x) M;      // CO2 in pc
real O2_pc(t,x) M;      // O2 in pc
real TO2_pc(t,x) M;     // total O2 in pc
real TempSys_pc(t,x) K;      //  Temperature in pc
//%END PC_VARS
//%START SOLUBILITY_COEFF
real alphaCO2Sys_pc(t,x) M/mmHg;      // Solubility of CO2 in pc
real alphaO2Sys_pc(t,x) M/mmHg;      // Solubility of O2 in pc
// ******************************************************************************
// New solubility calcs..... handle up to 40C, for solubility in plasma(O2)/saline(CO2):
// Coefficients to fit solubility curves:
real  alphaO20calc  = 0.0082 	mL/mL/atm,  
      alphaO201calc  = 0.0331 	mL/mL/atm,
      alphaCO20calc = 1.526	mL/mL/atm,  // saline
      alphaCO201calc = 0.132	mL/mL/atm,
      O2k1	= -0.0061 	1/celsius,
      O2k2      = 0.0292	1/celsius,
      CO2k1	= 0.0385	1/celsius,
      CO2k2	= -0.0105 	1/celsius,
     O2solSys(t,x)	mL/mL/atm, // ml O2 per ml plasma per atm
     CO2solSys(t,x)	mL/mL/atm;  
 // ANALYTIC SOLUTION
real TempC_Sys(t,x) celsius;
     TempC_Sys = (TempSys_pc-273.15)* (1 celsius)/(1 K);  // convert kelvin to celsius 
     O2solSys = alphaO20calc*exp(-O2k1*TempC_Sys)+alphaO201calc*exp(-O2k2*TempC_Sys);  
     CO2solSys = alphaCO20calc*exp(-CO2k1*TempC_Sys)+alphaCO201calc*exp(-CO2k2*TempC_Sys);
// convert solubilities to M/mmHg;
     alphaO2Sys_pc =  (1 atm)*O2solSys/Rconst/TempSys_pc;
     alphaCO2Sys_pc = (1 atm)* CO2solSys/Rconst/TempSys_pc;
// **** END of solubility calcs........
//%END SOLUBILITY_COEFF
//%START PC_INPUTS
real HCO3m_isf(t,x) M;    // HCO3m in isf
real Hp_isf(t,x) M;       // H+ conc in isf
real CO2_isf(t,x) M;     // Free CO2 from isf
real O2_isf(t,x) M;     // Free O2 in ISF
//%END PC_INPUTS
real TempSys_isf(t,x) K;     // Input: Temperature in ISF region
real MolV  = 22.414 L/mol;     // Liters of gas per mole at std temp.
real specificHeat = 1 kilocal/(g*K*min);   // Heat required to raise mass one K per min
real ThermCoeff = 0.0001 g/sec;  // Thermal coefficient
real TempExp K;       // Initial temperature
//%START PC_REGION
// Parenchymal cell region:
real CFpc = 10000;              // catalytic factor in pc due to CA
real BCpc = 45 mM;              // buffering capacity in pc
real kp1    = 0.12 1/sec,	// forward rate constant in CO2+H2O reaction
     km1    = 89 1/sec;         // backward rate constant in CO2+H2O reaction
real KH2CO3 = 5.5e-4 M;	        // EQUIL constant in H2CO3 ionization
real K1_calc = (kp1/km1)*KH2CO3; // EQUIL constant in overall CO2+H2O reaction
//real Fpc ml/(min*g);           // Flow of isf. Assume zero.
real Vpc ml/g;                 // Volume of region
real Wpc dimensionless;	       // fractional water content of pc
real VWpc ml/g;	               // volume of water content in pc
real HCO3mDpc = 1e-4 cm^2/sec; // diffusion coefficient for HCO3- in pc
real HCO3mPSpc ml/(min*g);    // PS for HCO3- across pc Membrane
real HpDpc = 1e-4 cm^2/sec;    // diffusion coefficient for H+ in pc
real HpPSpc ml/(min*g);       // PS for H+ across pc Membrane
real Rpc dimensionless;       // Gibbs-Donnan ratio [H+]isf/[H+]pc
//real Fpc ml/(min*g);           // Flow of pc. Assume zero.
real CO2Dpc = 1e-4 cm^2/sec;   // diffusion coefficient for CO2 in pc 
real CO2PSpc ml/(g*min);       // PS for CO2 across pc Membrane
real O2Dpc = 1e-4 cm^2/sec;    // diffusion coefficient for O2 in pc 
real O2PSpc ml/(g*min);        // PS for O2 across PC Membrane
real TDpc(t,x) s^-1;      // Temp diffusion across pc (muscle) membrane
real TCpc cm^2/sec;       // Thermal conductivity (Diffusion) in pc
     TCpc =  ThermCoeff*Vpc/L;
real O2Gpc(t,x) ml/min/g;      // gulosity for O2 consumption in PCs
real RQ = 0.80 dimensionless;  // respiratory quotient (unitless)
real Vmax_pc = 5.0e-6 mol/min/g; // Vmax for O2 consumpt by cytochrome oxidase in PCs ;
real Km_pc = 7.0e-7 M;         // Km for O2 on cytochrome oxidase in PCs;
real MRO2pc(t,x) mol/min/g;    // O2 consumption rate in PCs 
real O2energy_cal kilocal/mol/sec;         // Energy released per mole O2 consumed per sec.
real O2energy_calperliter = 1000 kilocal/L/sec; // Energy released per liter of O2 consumed per sec
     O2energy_cal = O2energy_calperliter*MolV;
real HeatGain(t,x) kilocal/min/g/sec;
     HeatGain = O2energy_cal*MRO2pc;
real TempGain(t,x) K/sec;
     TempGain = HeatGain/specificHeat;
// pc MbRateConsts:
real  CMb    = 0.4e-3 M,	// CONC of Mb in PCs = 5/(16800*FCV) M
      CMbpc  = CMb/Wpc,		// CONC of Mb in the water space of PCs
      P50Mb  = 2.4 mmHg;	// P50 for MbO2 saturation
// pc END MbRateConsts
real KMbO2(t,x) M^(-1);	// Hill's coefficient for MbO2 saturation
real SMbO2(t,x);	// MbO2 saturation (unitless)
real CMbO2(t,x) M;      // concentration of MbO2 in RBCs
real pCO2_pc(t,x) mmHg;   // convert to partial press
     pCO2_pc = CO2_pc/alphaCO2Sys_pc;
real pO2_pc(t,x) mmHg;   // convert to partial press
     pO2_pc = O2_pc/alphaO2Sys_pc;
// Inputs needed for HCO3pc_MPC
real HpCpct0 M;
real pH_pct0 dimensionless;   // Initial pH in pc
real CO2Cpct0 M;        // Initial free CO2 conc in pc
real HCO3mCpct0 M;      // Initial HCO3m conc in pc
     HpCpct0 = 10^(-pH_pct0)* (1 M);
     HCO3mCpct0 = K1_calc*CO2Cpct0/HpCpct0; 
real CO2Cpct0 M;
real TO2Cpct0(x) M;
real KMbO2t0(x) M^(-1);
real SMbO2t0(x);
real pO2pct0 mmHg;
real O2Cpct0 M;
real CMbO2t0(x) M;
 
 when (t=t.min) {	// PC INITIAL CONDITIONS
      HCO3m_pc = HCO3mCpct0;
      Hp_pc = HpCpct0;
      CO2_pc = CO2Cpct0;
      TO2_pc = TO2Cpct0;  // O2 PDE INITIAL CONDITION
  TempSys_pc = TempExp;
 } // end pc ICs
 when (x=x.min) {	// PC LEFT BOUNDARY CONDITIONS
      HCO3mDpc*HCO3m_pc:x =0;
      HpDpc*Hp_pc:x =0;
     CO2Dpc*CO2_pc:x =0;
     O2Dpc*TO2_pc:x =0;
  TCpc*TempSys_pc:x = 0;  
 } // end pc left BCs
 when (x=x.max) {	// PC RIGHT BOUNDARY CONDITIONS
      HCO3mDpc*HCO3m_pc:x = 0;
      HpDpc*Hp_pc:x = 0;
     CO2Dpc*CO2_pc:x = 0;
     O2Dpc*TO2_pc:x = 0;
     
  TempSys_pc:x =0;   
 } // end pc right BCs
// PDEs:
 CO2_pc:t =   (CO2PSpc/VWpc)*(CO2_isf-CO2_pc)
             + RQ*(O2Gpc/VWpc)*O2_pc      // CO2 production calc.          
             + CO2Dpc*(CO2_pc:x:x) 
             - CFpc*(kp1*CO2_pc-(km1/KH2CO3)*HCO3m_pc*Hp_pc);
 HCO3m_pc:t = HCO3mDpc*(HCO3m_pc:x:x)
               + (HCO3mPSpc/VWpc)*(Rpc*HCO3m_isf-HCO3m_pc)
               + CFpc*(kp1*CO2_pc-(km1/KH2CO3)*HCO3m_pc*Hp_pc);
 Hp_pc:t = HpDpc*(Hp_pc:x:x) 
            + (HpPSpc/VWpc)*(Hp_isf-Rpc*Hp_pc)
            + (2.303/BCpc)*Hp_pc * CFpc*(kp1*CO2_pc-(km1/KH2CO3)*HCO3m_pc*Hp_pc);
// O2 PDE:
 TO2_pc:t = (O2PSpc/VWpc)*(O2_isf-O2_pc)
           + O2Dpc*(TO2_pc:x:x)
           - (O2Gpc/VWpc)*O2_pc;  // Consumption term             
TempSys_pc:t = -TDpc*(TempSys_pc-TempSys_isf) 
               +(O2energy_cal/specificHeat)*MRO2pc  // Temp gain
               + TCpc*(TempSys_pc):x:x ; //pc
// ------ End of pc pdes
// MbO2 binding:
 O2_pc = if (TO2_pc <= 0) 0
         else (((1+KMbO2*(CMbpc-TO2_pc))^2 + 4*KMbO2*TO2_pc)^0.5
             - (1+KMbO2*(CMbpc-TO2_pc)))/(2*KMbO2);
// HillEqs for Mb:
      KMbO2  = (alphaO2Sys_pc*P50Mb)^(-1);
      SMbO2  = KMbO2*O2_pc/(1+KMbO2*O2_pc);  // Myoglobin saturation in pc
      CMbO2  = CMbpc*SMbO2; 
     O2Gpc = Vmax_pc/(Km_pc+O2_pc);
     MRO2pc = O2Gpc*O2_pc;     // <--- get an estimate of O2 consumed 
//%END PC_REGION
//%START O2_CO2_SOLUBILITY_INIT
real alphaO20 = 1.1148E-6;     // Initial solubility
real alphaCO20 = 2.8472E-5 M/mmHg;  // Solubility of CO2 in pc
//%END O2_CO2_SOLUBILITY_INIT
//%START PC_INIT
     pH_pct0 = 7.24;   // Input initial condition
real pO2pct0 = 90 mmHg;
real pCO2pct0 = 40 mmHg;
     O2Cpct0 =  pO2pct0 * alphaO20;
     CO2Cpct0 = pCO2pct0 * alphaCO20; 
     KMbO2t0  = (alphaO20*P50Mb)^(-1);
     SMbO2t0  = KMbO2t0*O2Cpct0/(1+KMbO2t0*O2Cpct0);
     CMbO2t0  = CMbpc*SMbO2t0;
     TO2Cpct0 = O2Cpct0+CMbO2t0;
//%END PC_INIT
//%START PC_ASSIGN_INPUTS
// Set inputs needed by PC modules to whole model variables/params
    O2_isf = 0;  // M
    CO2_isf = 0.0;
    HCO3m_isf = 0;
    Hp_isf = 0; 
 //   alphaO2Sys = alphaO20; // Solubility of O2 in pc
 //   alphaCO2Sys = alphaCO20; // Solubility of O2 in pc
// ** pc Fin at x.min =0, assume Flow in pc is zero
//%END PC_ASSIGN_INPUTS
//%START PS_VALUES
    O2PSpc = 0;
    CO2PSpc = 0;
    HCO3mPSpc =0;
    HpPSpc = 0;
//%END PS_VALUES
//%START SYS_PARAMS
// System parameters:
real  H = 0.45 dimensionless;   // large vessel hematocrit (unitless)
real  Fb = 1.0 ml/(min*g);	// total blood flow  (blood flow per gram tissue)
      TempExp = 310;
//%END SYS_PARAMS
//%START PC_SYS_PARAMS
real  Fpc = 0;	   // flow of pc, ml/(min*g):  Assume zero, Do not use LSFEA solver if substrate const over all x.
real  Vpc = 0.01;	        // anatomical volume of pc
real  Wpc = 0.72;               // fractional water content of pc
      VWpc = Vpc*Wpc;
      Rpc = 0.79;               // Gibbs-Donnan ratio [H+]isf/[H+]pc
      TDpc =0;                  // Temp diffusion across pc (muscle) membrane
//%END PC_SYS_PARAMS
   TempSys_isf = 310;
} // End of pc pde model
// This MML file generated from PC_PDE_MPC.mpc using MPC v1.01.
