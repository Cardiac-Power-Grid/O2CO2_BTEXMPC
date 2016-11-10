
// Calculate free CO2 in RBCs based on total CO2.
// Compute free CO2 and O2 (or PCO2 and PO2) from total O2 and CO2 (TO2 and TCO2)
// through numerical inversion method that uses SHbO2CO2_EJAP2016 routine.





import nsrunit;
unit conversion on;
 unit celsius = fundamental;

math TotCO2_CO2InvertMPC {

 realDomain t sec; t.min=0; t.max=300; t.delta=1;   // time domain
 realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid;  
//%START CO2INVERT
// ***  CO2 invert to get pCO2:
real TotCO2_Invert_in(t,x) M;  // Tot CO2 used for inversion to get pCO2  
real pO2_CO2_Invert_in(t,x) mmHg;  // pO2 used for inversion to get pCO2
real pHrbc_in(t,x);    // pHrbc used for inversion to get pO2/CO2
real DPGrbc_in(t,x) M; // DPGrbc used for inversion to get pO2/CO2
real Temp_in(t,x) K;   // Temp used for inversion to get pO2/CO2
real Hbrbc_in(t,x) M;  // Hbrbc used for inversion to get pO2/CO2
real Hct_in;           // Vol fraction of RBCs in blood (Hematocrit), needed?
real alphaCO2Sys(t,x) M/mmHg; // solubility of CO2 in plasma.
real alphaO2Sys(t,x) M/mmHg; // solubility of O2 in water. Calc can be done in SHbO2CO2_EJAP2016.
real pCO2initguess =40 mmHg;   // Init guess.
real pCO2new_CO2_Invert_0(t,x) = pCO2initguess;  
real pHrbc(t,x) dimensionless; real DPGrbc(t,x) M ; // pH and DPG conc in RBCs
real TempSHbO2CO2(t,x) K;          // Physiological temperature, Used in Dash2016SHbO2CO2MPC
real Hct dimensionless ;   // Volume fraction of RBCs in blood (Hematocrit), unitless
real mol2ml = 22400 mol/ml;          // Conversion factor from mol of gas to ml of gas at STP
real AmHbBl = Hct*(150 g)/(0.45 L);  // Amount of hemoglobin in gm per liter of blood, gm/L
real MWHb = 64100 g/mol;             // Molecular weight of hemoglobin, gm/mol
real HbBl M;
     HbBl = AmHbBl/MWHb;             // Concentration of hemoglobin in blood, mol/L or M
// TotCO2_FreeCO2_Inversion: calcs pCO2 given TotCO2, pO2, alphaCO2, 
real pCO2_out(t,x) mmHg;  // free CO2 in rbc
real FCO2_out(t,x) M;     // free CO2 in rbc
real SHbCO2_out(t,x);     // Value from TCO2 inversion calcs
// IterativeVars
real pCO2c_CO2_Invert_1(t,x) mmHg;
real pCO2p_CO2_Invert_1(t,x) mmHg;
real pCO2m_CO2_Invert_1(t,x) mmHg;
real dCO2totCO2_Invert_1c(t,x);
real funcCO2_CO2_Invert_1(t,x);
real dfuncCO2_CO2_Invert_1(t,x);
real pCO2new_CO2_Invert_1(t,x) mmHg;
real pCO2c_CO2_Invert_2(t,x) mmHg;
real pCO2p_CO2_Invert_2(t,x) mmHg;
real pCO2m_CO2_Invert_2(t,x) mmHg;
real dCO2totCO2_Invert_2c(t,x);
real funcCO2_CO2_Invert_2(t,x);
real dfuncCO2_CO2_Invert_2(t,x);
real pCO2new_CO2_Invert_2(t,x) mmHg;
real pCO2c_CO2_Invert_3(t,x) mmHg;
real pCO2p_CO2_Invert_3(t,x) mmHg;
real pCO2m_CO2_Invert_3(t,x) mmHg;
real dCO2totCO2_Invert_3c(t,x);
real funcCO2_CO2_Invert_3(t,x);
real dfuncCO2_CO2_Invert_3(t,x);
real pCO2new_CO2_Invert_3(t,x) mmHg;
real pCO2c_CO2_Invert_4(t,x) mmHg;
real pCO2p_CO2_Invert_4(t,x) mmHg;
real pCO2m_CO2_Invert_4(t,x) mmHg;
real dCO2totCO2_Invert_4c(t,x);
real funcCO2_CO2_Invert_4(t,x);
real dfuncCO2_CO2_Invert_4(t,x);
real pCO2new_CO2_Invert_4(t,x) mmHg;
real pCO2c_CO2_Invert_5(t,x) mmHg;
real pCO2p_CO2_Invert_5(t,x) mmHg;
real pCO2m_CO2_Invert_5(t,x) mmHg;
real dCO2totCO2_Invert_5c(t,x);
real funcCO2_CO2_Invert_5(t,x);
real dfuncCO2_CO2_Invert_5(t,x);
real pCO2new_CO2_Invert_5(t,x) mmHg;
real pCO2c_CO2_Invert_6(t,x) mmHg;
real pCO2p_CO2_Invert_6(t,x) mmHg;
real pCO2m_CO2_Invert_6(t,x) mmHg;
real dCO2totCO2_Invert_6c(t,x);
real funcCO2_CO2_Invert_6(t,x);
real dfuncCO2_CO2_Invert_6(t,x);
real pCO2new_CO2_Invert_6(t,x) mmHg;
// Parameters those are fixed in the model (i.e. water fractions, RBCs hemoglobin
// concentration, equilibrium constants, and Hill coefficient)
real Hbrbc(t,x) M;
     Hbrbc = HbBl/Hct;              // Concentration of hemoglobin in RBCs, M; 
real Wpl = 0.94 dimensionless;      // fractional water space of plasma; unitless
real Wrbc = 0.65 dimensionless;     // fractional water space of RBCs; unitless
real Wbl = (1-Hct)*Wpl+Hct*Wrbc;    // fractional water space of blood; unitless
real K1 = 10^(-6.12);               // CO2 hydration reaction equilibrium constant 
real K2 = 21.5e-6 dimensionless;    // CO2 + HbNH2 equilibrium constant; unitless
real K2dp = 1e-6 M;                 // HbNHCOOH dissociation constant; M
real K2p = K2/K2dp;                 // kf2p/kb2p; 1/M
real K3 = 11.3e-6 dimensionless;                  // CO2 + O2HbNH2 equilibrium constant; unitless
real K3dp = 1e-6 M;                 // O2HbNHCOOH dissociation constant; M
real K3p 1/M;                       // kf3p/kb3p; 
     K3p = K3/K3dp;                
real K5dp = 2.4e-8 M;               // HbNH3+ dissociation constant; M
real K6dp = 1.2e-8 M;               // O2HbNH3+ dissociation constant; M
real Rrbc = 0.69;                   // Gibbs-Donnan ratio across the RBC membrane
// Variables those are fixed in the model with values at standard physiological conditions
// (i.e. pO20, pCO20, pHpl0, pHrbc0, DPGrbc0, Temp0)
real pO20 = 100 mmHg;               // standard O2 partial pressure in blood; mmHg
real pCO20 = 40 mmHg;               // standard CO2 partial pressure in blood; mmHg
real pHrbc0 = 7.24 dimensionless;   // standard pH in RBCs; unitless
real pHpl0 = pHrbc0-log(Rrbc);      // standard pH in plsama; unitless
real DPGrbc0 = 4.65e-3;             // standard 2,3-DPG concentration in RBCs; M
real Temp0 = 310 K;                 // standard temperature in blood; deg K
real P500 = 26.8 mmHg;              // standard pO2 for 50% SHbO2; mmHg
real pO2_CO2_Invert_1(t,x) mmHg;    // partial pressure in rbc
real pO2_CO2_Invert_2(t,x) mmHg;    // partial pressure in rbc
real pO2_CO2_Invert_3(t,x) mmHg;    // partial pressure in rbc
real pO2_CO2_Invert_4(t,x) mmHg;    // partial pressure in rbc
real pO2_CO2_Invert_5(t,x) mmHg;    // partial pressure in rbc
real pO2_CO2_Invert_6(t,x) mmHg;    // partial pressure in rbc
real pCO2_CO2_Invert_1c(t,x) mmHg;   // partial pressure in rbc
real pCO2_CO2_Invert_1p(t,x) mmHg;   // partial pressure in rbc
real pCO2_CO2_Invert_1m(t,x) mmHg;   // partial pressure in rbc
real pCO2_CO2_Invert_2c(t,x) mmHg;   // partial pressure in rbc
real pCO2_CO2_Invert_2p(t,x) mmHg;   // partial pressure in rbc
real pCO2_CO2_Invert_2m(t,x) mmHg;   // partial pressure in rbc
real pCO2_CO2_Invert_3c(t,x) mmHg;   // partial pressure in rbc
real pCO2_CO2_Invert_3p(t,x) mmHg;   // partial pressure in rbc
real pCO2_CO2_Invert_3m(t,x) mmHg;   // partial pressure in rbc
real pCO2_CO2_Invert_4c(t,x) mmHg;   // partial pressure in rbc
real pCO2_CO2_Invert_4p(t,x) mmHg;   // partial pressure in rbc
real pCO2_CO2_Invert_4m(t,x) mmHg;   // partial pressure in rbc
real pCO2_CO2_Invert_5c(t,x) mmHg;   // partial pressure in rbc
real pCO2_CO2_Invert_5p(t,x) mmHg;   // partial pressure in rbc
real pCO2_CO2_Invert_5m(t,x) mmHg;   // partial pressure in rbc
real pCO2_CO2_Invert_6c(t,x) mmHg;   // partial pressure in rbc
real pCO2_CO2_Invert_6p(t,x) mmHg;   // partial pressure in rbc
real pCO2_CO2_Invert_6m(t,x) mmHg;   // partial pressure in rbc
// Assign inputs to SHbO2CO2 calcs to get pCO2 given TCO2 in rbc
 Hct = Hct_in;
 pCO2_CO2_Invert_1c=pCO2c_CO2_Invert_1;
 pCO2_CO2_Invert_1p=pCO2p_CO2_Invert_1;
 pCO2_CO2_Invert_1m=pCO2m_CO2_Invert_1;
 pHrbc=pHrbc_in;
 DPGrbc = DPGrbc_in;
 pO2_CO2_Invert_1= pO2_CO2_Invert_in;
 TempSHbO2CO2 =Temp_in;  // Used in Dash2016SHbO2CO2MPC
 pCO2_CO2_Invert_2c=pCO2c_CO2_Invert_2;
 pCO2_CO2_Invert_2p=pCO2p_CO2_Invert_2;
 pCO2_CO2_Invert_2m=pCO2m_CO2_Invert_2;
 pO2_CO2_Invert_2= pO2_CO2_Invert_in;
 pCO2_CO2_Invert_3c=pCO2c_CO2_Invert_3;
 pCO2_CO2_Invert_3p=pCO2p_CO2_Invert_3;
 pCO2_CO2_Invert_3m=pCO2m_CO2_Invert_3;
 pO2_CO2_Invert_3= pO2_CO2_Invert_in;
 pCO2_CO2_Invert_4c=pCO2c_CO2_Invert_4;
 pCO2_CO2_Invert_4p=pCO2p_CO2_Invert_4;
 pCO2_CO2_Invert_4m=pCO2m_CO2_Invert_4;
 pO2_CO2_Invert_4= pO2_CO2_Invert_in;
 pCO2_CO2_Invert_5c=pCO2c_CO2_Invert_5;
 pCO2_CO2_Invert_5p=pCO2p_CO2_Invert_5;
 pCO2_CO2_Invert_5m=pCO2m_CO2_Invert_5;
 pO2_CO2_Invert_5= pO2_CO2_Invert_in;
 pCO2_CO2_Invert_6c=pCO2c_CO2_Invert_6;
 pCO2_CO2_Invert_6p=pCO2p_CO2_Invert_6;
 pCO2_CO2_Invert_6m=pCO2m_CO2_Invert_6;
 pO2_CO2_Invert_6= pO2_CO2_Invert_in;
// Calculation of intermediate variables in the computations of SHbO2 and SHbCO2
      // Rrbc = Hpl/Hrbc = 10^-(pHpl-pHrbc)
real delpHrbc(t,x) = pHrbc-pHrbc0; 
real delpCO2_CO2_Invert_1c(t,x) = pCO2_CO2_Invert_1c-pCO20; 
real delDPGrbc(t,x) = DPGrbc-DPGrbc0;
real Hrbc(t,x) M;
     Hrbc  = 10^(-pHrbc); 
real delTemp(t,x) = TempSHbO2CO2-Temp0;
real delpCO2_CO2_Invert_1p(t,x) = pCO2_CO2_Invert_1p-pCO20; 
real delpCO2_CO2_Invert_1m(t,x) = pCO2_CO2_Invert_1m-pCO20; 
real delpCO2_CO2_Invert_2c(t,x) = pCO2_CO2_Invert_2c-pCO20; 
real delpCO2_CO2_Invert_2p(t,x) = pCO2_CO2_Invert_2p-pCO20; 
real delpCO2_CO2_Invert_2m(t,x) = pCO2_CO2_Invert_2m-pCO20; 
real delpCO2_CO2_Invert_3c(t,x) = pCO2_CO2_Invert_3c-pCO20; 
real delpCO2_CO2_Invert_3p(t,x) = pCO2_CO2_Invert_3p-pCO20; 
real delpCO2_CO2_Invert_3m(t,x) = pCO2_CO2_Invert_3m-pCO20; 
real delpCO2_CO2_Invert_4c(t,x) = pCO2_CO2_Invert_4c-pCO20; 
real delpCO2_CO2_Invert_4p(t,x) = pCO2_CO2_Invert_4p-pCO20; 
real delpCO2_CO2_Invert_4m(t,x) = pCO2_CO2_Invert_4m-pCO20; 
real delpCO2_CO2_Invert_5c(t,x) = pCO2_CO2_Invert_5c-pCO20; 
real delpCO2_CO2_Invert_5p(t,x) = pCO2_CO2_Invert_5p-pCO20; 
real delpCO2_CO2_Invert_5m(t,x) = pCO2_CO2_Invert_5m-pCO20; 
real delpCO2_CO2_Invert_6c(t,x) = pCO2_CO2_Invert_6c-pCO20; 
real delpCO2_CO2_Invert_6p(t,x) = pCO2_CO2_Invert_6p-pCO20; 
real delpCO2_CO2_Invert_6m(t,x) = pCO2_CO2_Invert_6m-pCO20; 
real Rconst = 62.36358 L*mmHg/K/mole;  // ideal gas const
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
     TempC_Sys = (TempSHbO2CO2-273.15)* (1 celsius)/(1 K);  // convert kelvin to celsius 
     O2solSys = alphaO20calc*exp(-O2k1*TempC_Sys)+alphaO201calc*exp(-O2k2*TempC_Sys);  
     CO2solSys = alphaCO20calc*exp(-CO2k1*TempC_Sys)+alphaCO201calc*exp(-CO2k2*TempC_Sys);
// convert solubilities to M/mmHg;
     alphaO2Sys =  (1 atm)*O2solSys/Rconst/TempSHbO2CO2;
     alphaCO2Sys = (1 atm)* CO2solSys/Rconst/TempSHbO2CO2;
// **** END of solubility calcs........
real O2_CO2_Invert_1c(t,x) M;             // free O2_CO2_Invert_1c
     O2_CO2_Invert_1c = alphaO2Sys*pO2_CO2_Invert_1; 
real CO2_CO2_Invert_1c(t,x) M;            // free CO2_CO2_Invert_1c
     CO2_CO2_Invert_1c = alphaCO2Sys*pCO2_CO2_Invert_1c;  
real O2_CO2_Invert_1p(t,x) M;             // free O2_CO2_Invert_1p
     O2_CO2_Invert_1p = alphaO2Sys*pO2_CO2_Invert_1; 
real CO2_CO2_Invert_1p(t,x) M;            // free CO2_CO2_Invert_1p
     CO2_CO2_Invert_1p = alphaCO2Sys*pCO2_CO2_Invert_1p;  
real O2_CO2_Invert_1m(t,x) M;             // free O2_CO2_Invert_1m
     O2_CO2_Invert_1m = alphaO2Sys*pO2_CO2_Invert_1; 
real CO2_CO2_Invert_1m(t,x) M;            // free CO2_CO2_Invert_1m
     CO2_CO2_Invert_1m = alphaCO2Sys*pCO2_CO2_Invert_1m;  
real O2_CO2_Invert_2c(t,x) M;             // free O2_CO2_Invert_2c
     O2_CO2_Invert_2c = alphaO2Sys*pO2_CO2_Invert_2; 
real CO2_CO2_Invert_2c(t,x) M;            // free CO2_CO2_Invert_2c
     CO2_CO2_Invert_2c = alphaCO2Sys*pCO2_CO2_Invert_2c;  
real O2_CO2_Invert_2p(t,x) M;             // free O2_CO2_Invert_2p
     O2_CO2_Invert_2p = alphaO2Sys*pO2_CO2_Invert_2; 
real CO2_CO2_Invert_2p(t,x) M;            // free CO2_CO2_Invert_2p
     CO2_CO2_Invert_2p = alphaCO2Sys*pCO2_CO2_Invert_2p;  
real O2_CO2_Invert_2m(t,x) M;             // free O2_CO2_Invert_2m
     O2_CO2_Invert_2m = alphaO2Sys*pO2_CO2_Invert_2; 
real CO2_CO2_Invert_2m(t,x) M;            // free CO2_CO2_Invert_2m
     CO2_CO2_Invert_2m = alphaCO2Sys*pCO2_CO2_Invert_2m;  
real O2_CO2_Invert_3c(t,x) M;             // free O2_CO2_Invert_3c
     O2_CO2_Invert_3c = alphaO2Sys*pO2_CO2_Invert_3; 
real CO2_CO2_Invert_3c(t,x) M;            // free CO2_CO2_Invert_3c
     CO2_CO2_Invert_3c = alphaCO2Sys*pCO2_CO2_Invert_3c;  
real O2_CO2_Invert_3p(t,x) M;             // free O2_CO2_Invert_3p
     O2_CO2_Invert_3p = alphaO2Sys*pO2_CO2_Invert_3; 
real CO2_CO2_Invert_3p(t,x) M;            // free CO2_CO2_Invert_3p
     CO2_CO2_Invert_3p = alphaCO2Sys*pCO2_CO2_Invert_3p;  
real O2_CO2_Invert_3m(t,x) M;             // free O2_CO2_Invert_3m
     O2_CO2_Invert_3m = alphaO2Sys*pO2_CO2_Invert_3; 
real CO2_CO2_Invert_3m(t,x) M;            // free CO2_CO2_Invert_3m
     CO2_CO2_Invert_3m = alphaCO2Sys*pCO2_CO2_Invert_3m;  
real O2_CO2_Invert_4c(t,x) M;             // free O2_CO2_Invert_4c
     O2_CO2_Invert_4c = alphaO2Sys*pO2_CO2_Invert_4; 
real CO2_CO2_Invert_4c(t,x) M;            // free CO2_CO2_Invert_4c
     CO2_CO2_Invert_4c = alphaCO2Sys*pCO2_CO2_Invert_4c;  
real O2_CO2_Invert_4p(t,x) M;             // free O2_CO2_Invert_4p
     O2_CO2_Invert_4p = alphaO2Sys*pO2_CO2_Invert_4; 
real CO2_CO2_Invert_4p(t,x) M;            // free CO2_CO2_Invert_4p
     CO2_CO2_Invert_4p = alphaCO2Sys*pCO2_CO2_Invert_4p;  
real O2_CO2_Invert_4m(t,x) M;             // free O2_CO2_Invert_4m
     O2_CO2_Invert_4m = alphaO2Sys*pO2_CO2_Invert_4; 
real CO2_CO2_Invert_4m(t,x) M;            // free CO2_CO2_Invert_4m
     CO2_CO2_Invert_4m = alphaCO2Sys*pCO2_CO2_Invert_4m;  
real O2_CO2_Invert_5c(t,x) M;             // free O2_CO2_Invert_5c
     O2_CO2_Invert_5c = alphaO2Sys*pO2_CO2_Invert_5; 
real CO2_CO2_Invert_5c(t,x) M;            // free CO2_CO2_Invert_5c
     CO2_CO2_Invert_5c = alphaCO2Sys*pCO2_CO2_Invert_5c;  
real O2_CO2_Invert_5p(t,x) M;             // free O2_CO2_Invert_5p
     O2_CO2_Invert_5p = alphaO2Sys*pO2_CO2_Invert_5; 
real CO2_CO2_Invert_5p(t,x) M;            // free CO2_CO2_Invert_5p
     CO2_CO2_Invert_5p = alphaCO2Sys*pCO2_CO2_Invert_5p;  
real O2_CO2_Invert_5m(t,x) M;             // free O2_CO2_Invert_5m
     O2_CO2_Invert_5m = alphaO2Sys*pO2_CO2_Invert_5; 
real CO2_CO2_Invert_5m(t,x) M;            // free CO2_CO2_Invert_5m
     CO2_CO2_Invert_5m = alphaCO2Sys*pCO2_CO2_Invert_5m;  
real O2_CO2_Invert_6c(t,x) M;             // free O2_CO2_Invert_6c
     O2_CO2_Invert_6c = alphaO2Sys*pO2_CO2_Invert_6; 
real CO2_CO2_Invert_6c(t,x) M;            // free CO2_CO2_Invert_6c
     CO2_CO2_Invert_6c = alphaCO2Sys*pCO2_CO2_Invert_6c;  
real O2_CO2_Invert_6p(t,x) M;             // free O2_CO2_Invert_6p
     O2_CO2_Invert_6p = alphaO2Sys*pO2_CO2_Invert_6; 
real CO2_CO2_Invert_6p(t,x) M;            // free CO2_CO2_Invert_6p
     CO2_CO2_Invert_6p = alphaCO2Sys*pCO2_CO2_Invert_6p;  
real O2_CO2_Invert_6m(t,x) M;             // free O2_CO2_Invert_6m
     O2_CO2_Invert_6m = alphaO2Sys*pO2_CO2_Invert_6; 
real CO2_CO2_Invert_6m(t,x) M;            // free CO2_CO2_Invert_6m
     CO2_CO2_Invert_6m = alphaCO2Sys*pCO2_CO2_Invert_6m;  
// p50 adjustments:
real P501_CO2_Invert_1c(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_CO2_Invert_1c(t,x) = P500 + 1.7*(4.28e-2*delpCO2_CO2_Invert_1c + 3.64e-5*(1 mmHg^-1)*delpCO2_CO2_Invert_1c^2); // all standard conditions, except CO2
real P503_CO2_Invert_1c(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_CO2_Invert_1c(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_CO2_Invert_1c(t,x) mmHg;             // Final P50_CO2_Invert_1c value for 50% O2 binding to Hb
     P50_CO2_Invert_1c = P500*(P501_CO2_Invert_1c/P500)*(P502_CO2_Invert_1c/P500)*(P503_CO2_Invert_1c/P500)*(P504_CO2_Invert_1c/P500); 
real C50_CO2_Invert_1c(t,x) M;
     C50_CO2_Invert_1c = alphaO2Sys*P50_CO2_Invert_1c; 
real P501_CO2_Invert_1p(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_CO2_Invert_1p(t,x) = P500 + 1.7*(4.28e-2*delpCO2_CO2_Invert_1p + 3.64e-5*(1 mmHg^-1)*delpCO2_CO2_Invert_1p^2); // all standard conditions, except CO2
real P503_CO2_Invert_1p(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_CO2_Invert_1p(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_CO2_Invert_1p(t,x) mmHg;             // Final P50_CO2_Invert_1p value for 50% O2 binding to Hb
     P50_CO2_Invert_1p = P500*(P501_CO2_Invert_1p/P500)*(P502_CO2_Invert_1p/P500)*(P503_CO2_Invert_1p/P500)*(P504_CO2_Invert_1p/P500); 
real C50_CO2_Invert_1p(t,x) M;
     C50_CO2_Invert_1p = alphaO2Sys*P50_CO2_Invert_1p; 
real P501_CO2_Invert_1m(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_CO2_Invert_1m(t,x) = P500 + 1.7*(4.28e-2*delpCO2_CO2_Invert_1m + 3.64e-5*(1 mmHg^-1)*delpCO2_CO2_Invert_1m^2); // all standard conditions, except CO2
real P503_CO2_Invert_1m(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_CO2_Invert_1m(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_CO2_Invert_1m(t,x) mmHg;             // Final P50_CO2_Invert_1m value for 50% O2 binding to Hb
     P50_CO2_Invert_1m = P500*(P501_CO2_Invert_1m/P500)*(P502_CO2_Invert_1m/P500)*(P503_CO2_Invert_1m/P500)*(P504_CO2_Invert_1m/P500); 
real C50_CO2_Invert_1m(t,x) M;
     C50_CO2_Invert_1m = alphaO2Sys*P50_CO2_Invert_1m; 
real P501_CO2_Invert_2c(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_CO2_Invert_2c(t,x) = P500 + 1.7*(4.28e-2*delpCO2_CO2_Invert_2c + 3.64e-5*(1 mmHg^-1)*delpCO2_CO2_Invert_2c^2); // all standard conditions, except CO2
real P503_CO2_Invert_2c(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_CO2_Invert_2c(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_CO2_Invert_2c(t,x) mmHg;             // Final P50_CO2_Invert_2c value for 50% O2 binding to Hb
     P50_CO2_Invert_2c = P500*(P501_CO2_Invert_2c/P500)*(P502_CO2_Invert_2c/P500)*(P503_CO2_Invert_2c/P500)*(P504_CO2_Invert_2c/P500); 
real C50_CO2_Invert_2c(t,x) M;
     C50_CO2_Invert_2c = alphaO2Sys*P50_CO2_Invert_2c; 
real P501_CO2_Invert_2p(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_CO2_Invert_2p(t,x) = P500 + 1.7*(4.28e-2*delpCO2_CO2_Invert_2p + 3.64e-5*(1 mmHg^-1)*delpCO2_CO2_Invert_2p^2); // all standard conditions, except CO2
real P503_CO2_Invert_2p(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_CO2_Invert_2p(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_CO2_Invert_2p(t,x) mmHg;             // Final P50_CO2_Invert_2p value for 50% O2 binding to Hb
     P50_CO2_Invert_2p = P500*(P501_CO2_Invert_2p/P500)*(P502_CO2_Invert_2p/P500)*(P503_CO2_Invert_2p/P500)*(P504_CO2_Invert_2p/P500); 
real C50_CO2_Invert_2p(t,x) M;
     C50_CO2_Invert_2p = alphaO2Sys*P50_CO2_Invert_2p; 
real P501_CO2_Invert_2m(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_CO2_Invert_2m(t,x) = P500 + 1.7*(4.28e-2*delpCO2_CO2_Invert_2m + 3.64e-5*(1 mmHg^-1)*delpCO2_CO2_Invert_2m^2); // all standard conditions, except CO2
real P503_CO2_Invert_2m(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_CO2_Invert_2m(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_CO2_Invert_2m(t,x) mmHg;             // Final P50_CO2_Invert_2m value for 50% O2 binding to Hb
     P50_CO2_Invert_2m = P500*(P501_CO2_Invert_2m/P500)*(P502_CO2_Invert_2m/P500)*(P503_CO2_Invert_2m/P500)*(P504_CO2_Invert_2m/P500); 
real C50_CO2_Invert_2m(t,x) M;
     C50_CO2_Invert_2m = alphaO2Sys*P50_CO2_Invert_2m; 
real P501_CO2_Invert_3c(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_CO2_Invert_3c(t,x) = P500 + 1.7*(4.28e-2*delpCO2_CO2_Invert_3c + 3.64e-5*(1 mmHg^-1)*delpCO2_CO2_Invert_3c^2); // all standard conditions, except CO2
real P503_CO2_Invert_3c(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_CO2_Invert_3c(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_CO2_Invert_3c(t,x) mmHg;             // Final P50_CO2_Invert_3c value for 50% O2 binding to Hb
     P50_CO2_Invert_3c = P500*(P501_CO2_Invert_3c/P500)*(P502_CO2_Invert_3c/P500)*(P503_CO2_Invert_3c/P500)*(P504_CO2_Invert_3c/P500); 
real C50_CO2_Invert_3c(t,x) M;
     C50_CO2_Invert_3c = alphaO2Sys*P50_CO2_Invert_3c; 
real P501_CO2_Invert_3p(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_CO2_Invert_3p(t,x) = P500 + 1.7*(4.28e-2*delpCO2_CO2_Invert_3p + 3.64e-5*(1 mmHg^-1)*delpCO2_CO2_Invert_3p^2); // all standard conditions, except CO2
real P503_CO2_Invert_3p(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_CO2_Invert_3p(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_CO2_Invert_3p(t,x) mmHg;             // Final P50_CO2_Invert_3p value for 50% O2 binding to Hb
     P50_CO2_Invert_3p = P500*(P501_CO2_Invert_3p/P500)*(P502_CO2_Invert_3p/P500)*(P503_CO2_Invert_3p/P500)*(P504_CO2_Invert_3p/P500); 
real C50_CO2_Invert_3p(t,x) M;
     C50_CO2_Invert_3p = alphaO2Sys*P50_CO2_Invert_3p; 
real P501_CO2_Invert_3m(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_CO2_Invert_3m(t,x) = P500 + 1.7*(4.28e-2*delpCO2_CO2_Invert_3m + 3.64e-5*(1 mmHg^-1)*delpCO2_CO2_Invert_3m^2); // all standard conditions, except CO2
real P503_CO2_Invert_3m(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_CO2_Invert_3m(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_CO2_Invert_3m(t,x) mmHg;             // Final P50_CO2_Invert_3m value for 50% O2 binding to Hb
     P50_CO2_Invert_3m = P500*(P501_CO2_Invert_3m/P500)*(P502_CO2_Invert_3m/P500)*(P503_CO2_Invert_3m/P500)*(P504_CO2_Invert_3m/P500); 
real C50_CO2_Invert_3m(t,x) M;
     C50_CO2_Invert_3m = alphaO2Sys*P50_CO2_Invert_3m; 
real P501_CO2_Invert_4c(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_CO2_Invert_4c(t,x) = P500 + 1.7*(4.28e-2*delpCO2_CO2_Invert_4c + 3.64e-5*(1 mmHg^-1)*delpCO2_CO2_Invert_4c^2); // all standard conditions, except CO2
real P503_CO2_Invert_4c(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_CO2_Invert_4c(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_CO2_Invert_4c(t,x) mmHg;             // Final P50_CO2_Invert_4c value for 50% O2 binding to Hb
     P50_CO2_Invert_4c = P500*(P501_CO2_Invert_4c/P500)*(P502_CO2_Invert_4c/P500)*(P503_CO2_Invert_4c/P500)*(P504_CO2_Invert_4c/P500); 
real C50_CO2_Invert_4c(t,x) M;
     C50_CO2_Invert_4c = alphaO2Sys*P50_CO2_Invert_4c; 
real P501_CO2_Invert_4p(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_CO2_Invert_4p(t,x) = P500 + 1.7*(4.28e-2*delpCO2_CO2_Invert_4p + 3.64e-5*(1 mmHg^-1)*delpCO2_CO2_Invert_4p^2); // all standard conditions, except CO2
real P503_CO2_Invert_4p(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_CO2_Invert_4p(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_CO2_Invert_4p(t,x) mmHg;             // Final P50_CO2_Invert_4p value for 50% O2 binding to Hb
     P50_CO2_Invert_4p = P500*(P501_CO2_Invert_4p/P500)*(P502_CO2_Invert_4p/P500)*(P503_CO2_Invert_4p/P500)*(P504_CO2_Invert_4p/P500); 
real C50_CO2_Invert_4p(t,x) M;
     C50_CO2_Invert_4p = alphaO2Sys*P50_CO2_Invert_4p; 
real P501_CO2_Invert_4m(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_CO2_Invert_4m(t,x) = P500 + 1.7*(4.28e-2*delpCO2_CO2_Invert_4m + 3.64e-5*(1 mmHg^-1)*delpCO2_CO2_Invert_4m^2); // all standard conditions, except CO2
real P503_CO2_Invert_4m(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_CO2_Invert_4m(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_CO2_Invert_4m(t,x) mmHg;             // Final P50_CO2_Invert_4m value for 50% O2 binding to Hb
     P50_CO2_Invert_4m = P500*(P501_CO2_Invert_4m/P500)*(P502_CO2_Invert_4m/P500)*(P503_CO2_Invert_4m/P500)*(P504_CO2_Invert_4m/P500); 
real C50_CO2_Invert_4m(t,x) M;
     C50_CO2_Invert_4m = alphaO2Sys*P50_CO2_Invert_4m; 
real P501_CO2_Invert_5c(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_CO2_Invert_5c(t,x) = P500 + 1.7*(4.28e-2*delpCO2_CO2_Invert_5c + 3.64e-5*(1 mmHg^-1)*delpCO2_CO2_Invert_5c^2); // all standard conditions, except CO2
real P503_CO2_Invert_5c(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_CO2_Invert_5c(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_CO2_Invert_5c(t,x) mmHg;             // Final P50_CO2_Invert_5c value for 50% O2 binding to Hb
     P50_CO2_Invert_5c = P500*(P501_CO2_Invert_5c/P500)*(P502_CO2_Invert_5c/P500)*(P503_CO2_Invert_5c/P500)*(P504_CO2_Invert_5c/P500); 
real C50_CO2_Invert_5c(t,x) M;
     C50_CO2_Invert_5c = alphaO2Sys*P50_CO2_Invert_5c; 
real P501_CO2_Invert_5p(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_CO2_Invert_5p(t,x) = P500 + 1.7*(4.28e-2*delpCO2_CO2_Invert_5p + 3.64e-5*(1 mmHg^-1)*delpCO2_CO2_Invert_5p^2); // all standard conditions, except CO2
real P503_CO2_Invert_5p(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_CO2_Invert_5p(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_CO2_Invert_5p(t,x) mmHg;             // Final P50_CO2_Invert_5p value for 50% O2 binding to Hb
     P50_CO2_Invert_5p = P500*(P501_CO2_Invert_5p/P500)*(P502_CO2_Invert_5p/P500)*(P503_CO2_Invert_5p/P500)*(P504_CO2_Invert_5p/P500); 
real C50_CO2_Invert_5p(t,x) M;
     C50_CO2_Invert_5p = alphaO2Sys*P50_CO2_Invert_5p; 
real P501_CO2_Invert_5m(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_CO2_Invert_5m(t,x) = P500 + 1.7*(4.28e-2*delpCO2_CO2_Invert_5m + 3.64e-5*(1 mmHg^-1)*delpCO2_CO2_Invert_5m^2); // all standard conditions, except CO2
real P503_CO2_Invert_5m(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_CO2_Invert_5m(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_CO2_Invert_5m(t,x) mmHg;             // Final P50_CO2_Invert_5m value for 50% O2 binding to Hb
     P50_CO2_Invert_5m = P500*(P501_CO2_Invert_5m/P500)*(P502_CO2_Invert_5m/P500)*(P503_CO2_Invert_5m/P500)*(P504_CO2_Invert_5m/P500); 
real C50_CO2_Invert_5m(t,x) M;
     C50_CO2_Invert_5m = alphaO2Sys*P50_CO2_Invert_5m; 
real P501_CO2_Invert_6c(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_CO2_Invert_6c(t,x) = P500 + 1.7*(4.28e-2*delpCO2_CO2_Invert_6c + 3.64e-5*(1 mmHg^-1)*delpCO2_CO2_Invert_6c^2); // all standard conditions, except CO2
real P503_CO2_Invert_6c(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_CO2_Invert_6c(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_CO2_Invert_6c(t,x) mmHg;             // Final P50_CO2_Invert_6c value for 50% O2 binding to Hb
     P50_CO2_Invert_6c = P500*(P501_CO2_Invert_6c/P500)*(P502_CO2_Invert_6c/P500)*(P503_CO2_Invert_6c/P500)*(P504_CO2_Invert_6c/P500); 
real C50_CO2_Invert_6c(t,x) M;
     C50_CO2_Invert_6c = alphaO2Sys*P50_CO2_Invert_6c; 
real P501_CO2_Invert_6p(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_CO2_Invert_6p(t,x) = P500 + 1.7*(4.28e-2*delpCO2_CO2_Invert_6p + 3.64e-5*(1 mmHg^-1)*delpCO2_CO2_Invert_6p^2); // all standard conditions, except CO2
real P503_CO2_Invert_6p(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_CO2_Invert_6p(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_CO2_Invert_6p(t,x) mmHg;             // Final P50_CO2_Invert_6p value for 50% O2 binding to Hb
     P50_CO2_Invert_6p = P500*(P501_CO2_Invert_6p/P500)*(P502_CO2_Invert_6p/P500)*(P503_CO2_Invert_6p/P500)*(P504_CO2_Invert_6p/P500); 
real C50_CO2_Invert_6p(t,x) M;
     C50_CO2_Invert_6p = alphaO2Sys*P50_CO2_Invert_6p; 
real P501_CO2_Invert_6m(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_CO2_Invert_6m(t,x) = P500 + 1.7*(4.28e-2*delpCO2_CO2_Invert_6m + 3.64e-5*(1 mmHg^-1)*delpCO2_CO2_Invert_6m^2); // all standard conditions, except CO2
real P503_CO2_Invert_6m(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_CO2_Invert_6m(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_CO2_Invert_6m(t,x) mmHg;             // Final P50_CO2_Invert_6m value for 50% O2 binding to Hb
     P50_CO2_Invert_6m = P500*(P501_CO2_Invert_6m/P500)*(P502_CO2_Invert_6m/P500)*(P503_CO2_Invert_6m/P500)*(P504_CO2_Invert_6m/P500); 
real C50_CO2_Invert_6m(t,x) M;
     C50_CO2_Invert_6m = alphaO2Sys*P50_CO2_Invert_6m; 
// equilibirum consts are not directly tied to CO2 and O2 conc:
real BPH1(t,x) dimensionless;
     BPH1 = 1+K2dp/Hrbc;         // Binding polynomial involving K2dp and Hrbc 
real BPH2(t,x) = 1+K3dp/Hrbc;      // Binding polynomial involving K3dp and Hrbc 
real BPH3(t,x) = 1+Hrbc/K5dp;      // Binding polynomial involving K5dp and Hrbc 
real BPH4(t,x) = 1+Hrbc/K6dp;      // Binding polynomial involving K6dp and Hrbc 
// Hill coefficient; unitless (redefined as a function pO2_CO2_Invert_1)
real alpha = 2.8 dimensionless; real beta = 1.2 dimensionless; real gamma = 29.2 mmHg; // Roughton et al data
real nH_CO2_Invert_1c(t,x) dimensionless;
     nH_CO2_Invert_1c = alpha-beta*10^(-pO2_CO2_Invert_1/gamma); // pO2_CO2_Invert_1 dependent variable nH_CO2_Invert_1c
real nH_CO2_Invert_1p(t,x) dimensionless;
     nH_CO2_Invert_1p = alpha-beta*10^(-pO2_CO2_Invert_1/gamma); // pO2_CO2_Invert_1 dependent variable nH_CO2_Invert_1p
real nH_CO2_Invert_1m(t,x) dimensionless;
     nH_CO2_Invert_1m = alpha-beta*10^(-pO2_CO2_Invert_1/gamma); // pO2_CO2_Invert_1 dependent variable nH_CO2_Invert_1m
// Hill coefficient; unitless (redefined as a function pO2_CO2_Invert_2)
real nH_CO2_Invert_2c(t,x) dimensionless;
     nH_CO2_Invert_2c = alpha-beta*10^(-pO2_CO2_Invert_2/gamma); // pO2_CO2_Invert_2 dependent variable nH_CO2_Invert_2c
real nH_CO2_Invert_2p(t,x) dimensionless;
     nH_CO2_Invert_2p = alpha-beta*10^(-pO2_CO2_Invert_2/gamma); // pO2_CO2_Invert_2 dependent variable nH_CO2_Invert_2p
real nH_CO2_Invert_2m(t,x) dimensionless;
     nH_CO2_Invert_2m = alpha-beta*10^(-pO2_CO2_Invert_2/gamma); // pO2_CO2_Invert_2 dependent variable nH_CO2_Invert_2m
// Hill coefficient; unitless (redefined as a function pO2_CO2_Invert_3)
real nH_CO2_Invert_3c(t,x) dimensionless;
     nH_CO2_Invert_3c = alpha-beta*10^(-pO2_CO2_Invert_3/gamma); // pO2_CO2_Invert_3 dependent variable nH_CO2_Invert_3c
real nH_CO2_Invert_3p(t,x) dimensionless;
     nH_CO2_Invert_3p = alpha-beta*10^(-pO2_CO2_Invert_3/gamma); // pO2_CO2_Invert_3 dependent variable nH_CO2_Invert_3p
real nH_CO2_Invert_3m(t,x) dimensionless;
     nH_CO2_Invert_3m = alpha-beta*10^(-pO2_CO2_Invert_3/gamma); // pO2_CO2_Invert_3 dependent variable nH_CO2_Invert_3m
// Hill coefficient; unitless (redefined as a function pO2_CO2_Invert_4)
real nH_CO2_Invert_4c(t,x) dimensionless;
     nH_CO2_Invert_4c = alpha-beta*10^(-pO2_CO2_Invert_4/gamma); // pO2_CO2_Invert_4 dependent variable nH_CO2_Invert_4c
real nH_CO2_Invert_4p(t,x) dimensionless;
     nH_CO2_Invert_4p = alpha-beta*10^(-pO2_CO2_Invert_4/gamma); // pO2_CO2_Invert_4 dependent variable nH_CO2_Invert_4p
real nH_CO2_Invert_4m(t,x) dimensionless;
     nH_CO2_Invert_4m = alpha-beta*10^(-pO2_CO2_Invert_4/gamma); // pO2_CO2_Invert_4 dependent variable nH_CO2_Invert_4m
// Hill coefficient; unitless (redefined as a function pO2_CO2_Invert_5)
real nH_CO2_Invert_5c(t,x) dimensionless;
     nH_CO2_Invert_5c = alpha-beta*10^(-pO2_CO2_Invert_5/gamma); // pO2_CO2_Invert_5 dependent variable nH_CO2_Invert_5c
real nH_CO2_Invert_5p(t,x) dimensionless;
     nH_CO2_Invert_5p = alpha-beta*10^(-pO2_CO2_Invert_5/gamma); // pO2_CO2_Invert_5 dependent variable nH_CO2_Invert_5p
real nH_CO2_Invert_5m(t,x) dimensionless;
     nH_CO2_Invert_5m = alpha-beta*10^(-pO2_CO2_Invert_5/gamma); // pO2_CO2_Invert_5 dependent variable nH_CO2_Invert_5m
// Hill coefficient; unitless (redefined as a function pO2_CO2_Invert_6)
real nH_CO2_Invert_6c(t,x) dimensionless;
     nH_CO2_Invert_6c = alpha-beta*10^(-pO2_CO2_Invert_6/gamma); // pO2_CO2_Invert_6 dependent variable nH_CO2_Invert_6c
real nH_CO2_Invert_6p(t,x) dimensionless;
     nH_CO2_Invert_6p = alpha-beta*10^(-pO2_CO2_Invert_6/gamma); // pO2_CO2_Invert_6 dependent variable nH_CO2_Invert_6p
real nH_CO2_Invert_6m(t,x) dimensionless;
     nH_CO2_Invert_6m = alpha-beta*10^(-pO2_CO2_Invert_6/gamma); // pO2_CO2_Invert_6 dependent variable nH_CO2_Invert_6m
// Compute the apparent equilibrium constant of Hb with O2_CO2_Invert_1c and CO2_CO2_Invert_1c (KHbO2 and KHbCO2); 
// O2_CO2_Invert_1c and CO2_CO2_Invert_1c saturations of Hb (SHbO2 and SHbCO2); and O2_CO2_Invert_1c and CO2_CO2_Invert_1c contents in blood. 
real K4p_CO2_Invert_1c(t,x) 1/M;  
     K4p_CO2_Invert_1c = (1 M^-1)*((O2_CO2_Invert_1c*(1 M^-1))^(nH_CO2_Invert_1c-1)*(K2p*BPH1*CO2_CO2_Invert_1c+BPH3))/( ((1 M^-1)*C50_CO2_Invert_1c)^nH_CO2_Invert_1c*(K3p*BPH2*CO2_CO2_Invert_1c+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_CO2_Invert_1p and CO2_CO2_Invert_1p (KHbO2 and KHbCO2); 
// O2_CO2_Invert_1p and CO2_CO2_Invert_1p saturations of Hb (SHbO2 and SHbCO2); and O2_CO2_Invert_1p and CO2_CO2_Invert_1p contents in blood. 
real K4p_CO2_Invert_1p(t,x) 1/M;  
     K4p_CO2_Invert_1p = (1 M^-1)*((O2_CO2_Invert_1p*(1 M^-1))^(nH_CO2_Invert_1p-1)*(K2p*BPH1*CO2_CO2_Invert_1p+BPH3))/( ((1 M^-1)*C50_CO2_Invert_1p)^nH_CO2_Invert_1p*(K3p*BPH2*CO2_CO2_Invert_1p+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_CO2_Invert_1m and CO2_CO2_Invert_1m (KHbO2 and KHbCO2); 
// O2_CO2_Invert_1m and CO2_CO2_Invert_1m saturations of Hb (SHbO2 and SHbCO2); and O2_CO2_Invert_1m and CO2_CO2_Invert_1m contents in blood. 
real K4p_CO2_Invert_1m(t,x) 1/M;  
     K4p_CO2_Invert_1m = (1 M^-1)*((O2_CO2_Invert_1m*(1 M^-1))^(nH_CO2_Invert_1m-1)*(K2p*BPH1*CO2_CO2_Invert_1m+BPH3))/( ((1 M^-1)*C50_CO2_Invert_1m)^nH_CO2_Invert_1m*(K3p*BPH2*CO2_CO2_Invert_1m+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_CO2_Invert_2c and CO2_CO2_Invert_2c (KHbO2 and KHbCO2); 
// O2_CO2_Invert_2c and CO2_CO2_Invert_2c saturations of Hb (SHbO2 and SHbCO2); and O2_CO2_Invert_2c and CO2_CO2_Invert_2c contents in blood. 
real K4p_CO2_Invert_2c(t,x) 1/M;  
     K4p_CO2_Invert_2c = (1 M^-1)*((O2_CO2_Invert_2c*(1 M^-1))^(nH_CO2_Invert_2c-1)*(K2p*BPH1*CO2_CO2_Invert_2c+BPH3))/( ((1 M^-1)*C50_CO2_Invert_2c)^nH_CO2_Invert_2c*(K3p*BPH2*CO2_CO2_Invert_2c+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_CO2_Invert_2p and CO2_CO2_Invert_2p (KHbO2 and KHbCO2); 
// O2_CO2_Invert_2p and CO2_CO2_Invert_2p saturations of Hb (SHbO2 and SHbCO2); and O2_CO2_Invert_2p and CO2_CO2_Invert_2p contents in blood. 
real K4p_CO2_Invert_2p(t,x) 1/M;  
     K4p_CO2_Invert_2p = (1 M^-1)*((O2_CO2_Invert_2p*(1 M^-1))^(nH_CO2_Invert_2p-1)*(K2p*BPH1*CO2_CO2_Invert_2p+BPH3))/( ((1 M^-1)*C50_CO2_Invert_2p)^nH_CO2_Invert_2p*(K3p*BPH2*CO2_CO2_Invert_2p+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_CO2_Invert_2m and CO2_CO2_Invert_2m (KHbO2 and KHbCO2); 
// O2_CO2_Invert_2m and CO2_CO2_Invert_2m saturations of Hb (SHbO2 and SHbCO2); and O2_CO2_Invert_2m and CO2_CO2_Invert_2m contents in blood. 
real K4p_CO2_Invert_2m(t,x) 1/M;  
     K4p_CO2_Invert_2m = (1 M^-1)*((O2_CO2_Invert_2m*(1 M^-1))^(nH_CO2_Invert_2m-1)*(K2p*BPH1*CO2_CO2_Invert_2m+BPH3))/( ((1 M^-1)*C50_CO2_Invert_2m)^nH_CO2_Invert_2m*(K3p*BPH2*CO2_CO2_Invert_2m+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_CO2_Invert_3c and CO2_CO2_Invert_3c (KHbO2 and KHbCO2); 
// O2_CO2_Invert_3c and CO2_CO2_Invert_3c saturations of Hb (SHbO2 and SHbCO2); and O2_CO2_Invert_3c and CO2_CO2_Invert_3c contents in blood. 
real K4p_CO2_Invert_3c(t,x) 1/M;  
     K4p_CO2_Invert_3c = (1 M^-1)*((O2_CO2_Invert_3c*(1 M^-1))^(nH_CO2_Invert_3c-1)*(K2p*BPH1*CO2_CO2_Invert_3c+BPH3))/( ((1 M^-1)*C50_CO2_Invert_3c)^nH_CO2_Invert_3c*(K3p*BPH2*CO2_CO2_Invert_3c+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_CO2_Invert_3p and CO2_CO2_Invert_3p (KHbO2 and KHbCO2); 
// O2_CO2_Invert_3p and CO2_CO2_Invert_3p saturations of Hb (SHbO2 and SHbCO2); and O2_CO2_Invert_3p and CO2_CO2_Invert_3p contents in blood. 
real K4p_CO2_Invert_3p(t,x) 1/M;  
     K4p_CO2_Invert_3p = (1 M^-1)*((O2_CO2_Invert_3p*(1 M^-1))^(nH_CO2_Invert_3p-1)*(K2p*BPH1*CO2_CO2_Invert_3p+BPH3))/( ((1 M^-1)*C50_CO2_Invert_3p)^nH_CO2_Invert_3p*(K3p*BPH2*CO2_CO2_Invert_3p+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_CO2_Invert_3m and CO2_CO2_Invert_3m (KHbO2 and KHbCO2); 
// O2_CO2_Invert_3m and CO2_CO2_Invert_3m saturations of Hb (SHbO2 and SHbCO2); and O2_CO2_Invert_3m and CO2_CO2_Invert_3m contents in blood. 
real K4p_CO2_Invert_3m(t,x) 1/M;  
     K4p_CO2_Invert_3m = (1 M^-1)*((O2_CO2_Invert_3m*(1 M^-1))^(nH_CO2_Invert_3m-1)*(K2p*BPH1*CO2_CO2_Invert_3m+BPH3))/( ((1 M^-1)*C50_CO2_Invert_3m)^nH_CO2_Invert_3m*(K3p*BPH2*CO2_CO2_Invert_3m+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_CO2_Invert_4c and CO2_CO2_Invert_4c (KHbO2 and KHbCO2); 
// O2_CO2_Invert_4c and CO2_CO2_Invert_4c saturations of Hb (SHbO2 and SHbCO2); and O2_CO2_Invert_4c and CO2_CO2_Invert_4c contents in blood. 
real K4p_CO2_Invert_4c(t,x) 1/M;  
     K4p_CO2_Invert_4c = (1 M^-1)*((O2_CO2_Invert_4c*(1 M^-1))^(nH_CO2_Invert_4c-1)*(K2p*BPH1*CO2_CO2_Invert_4c+BPH3))/( ((1 M^-1)*C50_CO2_Invert_4c)^nH_CO2_Invert_4c*(K3p*BPH2*CO2_CO2_Invert_4c+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_CO2_Invert_4p and CO2_CO2_Invert_4p (KHbO2 and KHbCO2); 
// O2_CO2_Invert_4p and CO2_CO2_Invert_4p saturations of Hb (SHbO2 and SHbCO2); and O2_CO2_Invert_4p and CO2_CO2_Invert_4p contents in blood. 
real K4p_CO2_Invert_4p(t,x) 1/M;  
     K4p_CO2_Invert_4p = (1 M^-1)*((O2_CO2_Invert_4p*(1 M^-1))^(nH_CO2_Invert_4p-1)*(K2p*BPH1*CO2_CO2_Invert_4p+BPH3))/( ((1 M^-1)*C50_CO2_Invert_4p)^nH_CO2_Invert_4p*(K3p*BPH2*CO2_CO2_Invert_4p+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_CO2_Invert_4m and CO2_CO2_Invert_4m (KHbO2 and KHbCO2); 
// O2_CO2_Invert_4m and CO2_CO2_Invert_4m saturations of Hb (SHbO2 and SHbCO2); and O2_CO2_Invert_4m and CO2_CO2_Invert_4m contents in blood. 
real K4p_CO2_Invert_4m(t,x) 1/M;  
     K4p_CO2_Invert_4m = (1 M^-1)*((O2_CO2_Invert_4m*(1 M^-1))^(nH_CO2_Invert_4m-1)*(K2p*BPH1*CO2_CO2_Invert_4m+BPH3))/( ((1 M^-1)*C50_CO2_Invert_4m)^nH_CO2_Invert_4m*(K3p*BPH2*CO2_CO2_Invert_4m+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_CO2_Invert_5c and CO2_CO2_Invert_5c (KHbO2 and KHbCO2); 
// O2_CO2_Invert_5c and CO2_CO2_Invert_5c saturations of Hb (SHbO2 and SHbCO2); and O2_CO2_Invert_5c and CO2_CO2_Invert_5c contents in blood. 
real K4p_CO2_Invert_5c(t,x) 1/M;  
     K4p_CO2_Invert_5c = (1 M^-1)*((O2_CO2_Invert_5c*(1 M^-1))^(nH_CO2_Invert_5c-1)*(K2p*BPH1*CO2_CO2_Invert_5c+BPH3))/( ((1 M^-1)*C50_CO2_Invert_5c)^nH_CO2_Invert_5c*(K3p*BPH2*CO2_CO2_Invert_5c+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_CO2_Invert_5p and CO2_CO2_Invert_5p (KHbO2 and KHbCO2); 
// O2_CO2_Invert_5p and CO2_CO2_Invert_5p saturations of Hb (SHbO2 and SHbCO2); and O2_CO2_Invert_5p and CO2_CO2_Invert_5p contents in blood. 
real K4p_CO2_Invert_5p(t,x) 1/M;  
     K4p_CO2_Invert_5p = (1 M^-1)*((O2_CO2_Invert_5p*(1 M^-1))^(nH_CO2_Invert_5p-1)*(K2p*BPH1*CO2_CO2_Invert_5p+BPH3))/( ((1 M^-1)*C50_CO2_Invert_5p)^nH_CO2_Invert_5p*(K3p*BPH2*CO2_CO2_Invert_5p+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_CO2_Invert_5m and CO2_CO2_Invert_5m (KHbO2 and KHbCO2); 
// O2_CO2_Invert_5m and CO2_CO2_Invert_5m saturations of Hb (SHbO2 and SHbCO2); and O2_CO2_Invert_5m and CO2_CO2_Invert_5m contents in blood. 
real K4p_CO2_Invert_5m(t,x) 1/M;  
     K4p_CO2_Invert_5m = (1 M^-1)*((O2_CO2_Invert_5m*(1 M^-1))^(nH_CO2_Invert_5m-1)*(K2p*BPH1*CO2_CO2_Invert_5m+BPH3))/( ((1 M^-1)*C50_CO2_Invert_5m)^nH_CO2_Invert_5m*(K3p*BPH2*CO2_CO2_Invert_5m+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_CO2_Invert_6c and CO2_CO2_Invert_6c (KHbO2 and KHbCO2); 
// O2_CO2_Invert_6c and CO2_CO2_Invert_6c saturations of Hb (SHbO2 and SHbCO2); and O2_CO2_Invert_6c and CO2_CO2_Invert_6c contents in blood. 
real K4p_CO2_Invert_6c(t,x) 1/M;  
     K4p_CO2_Invert_6c = (1 M^-1)*((O2_CO2_Invert_6c*(1 M^-1))^(nH_CO2_Invert_6c-1)*(K2p*BPH1*CO2_CO2_Invert_6c+BPH3))/( ((1 M^-1)*C50_CO2_Invert_6c)^nH_CO2_Invert_6c*(K3p*BPH2*CO2_CO2_Invert_6c+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_CO2_Invert_6p and CO2_CO2_Invert_6p (KHbO2 and KHbCO2); 
// O2_CO2_Invert_6p and CO2_CO2_Invert_6p saturations of Hb (SHbO2 and SHbCO2); and O2_CO2_Invert_6p and CO2_CO2_Invert_6p contents in blood. 
real K4p_CO2_Invert_6p(t,x) 1/M;  
     K4p_CO2_Invert_6p = (1 M^-1)*((O2_CO2_Invert_6p*(1 M^-1))^(nH_CO2_Invert_6p-1)*(K2p*BPH1*CO2_CO2_Invert_6p+BPH3))/( ((1 M^-1)*C50_CO2_Invert_6p)^nH_CO2_Invert_6p*(K3p*BPH2*CO2_CO2_Invert_6p+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_CO2_Invert_6m and CO2_CO2_Invert_6m (KHbO2 and KHbCO2); 
// O2_CO2_Invert_6m and CO2_CO2_Invert_6m saturations of Hb (SHbO2 and SHbCO2); and O2_CO2_Invert_6m and CO2_CO2_Invert_6m contents in blood. 
real K4p_CO2_Invert_6m(t,x) 1/M;  
     K4p_CO2_Invert_6m = (1 M^-1)*((O2_CO2_Invert_6m*(1 M^-1))^(nH_CO2_Invert_6m-1)*(K2p*BPH1*CO2_CO2_Invert_6m+BPH3))/( ((1 M^-1)*C50_CO2_Invert_6m)^nH_CO2_Invert_6m*(K3p*BPH2*CO2_CO2_Invert_6m+BPH4));
real KHbCO2_CO2_Invert_1c(t,x) 1/M;              // Apparent equilibrium constants of Hb with CO2_CO2_Invert_1c
     KHbCO2_CO2_Invert_1c = (K2p*BPH1+K3p*K4p_CO2_Invert_1c*BPH2*O2_CO2_Invert_1c)/(BPH3+K4p_CO2_Invert_1c*BPH4*O2_CO2_Invert_1c);
real SHbCO2_CO2_Invert_1c(t,x) dimensionless;
     SHbCO2_CO2_Invert_1c = KHbCO2_CO2_Invert_1c*CO2_CO2_Invert_1c/(1+KHbCO2_CO2_Invert_1c*CO2_CO2_Invert_1c);
real KHbCO2_CO2_Invert_1p(t,x) 1/M;              // Apparent equilibrium constants of Hb with CO2_CO2_Invert_1p
     KHbCO2_CO2_Invert_1p = (K2p*BPH1+K3p*K4p_CO2_Invert_1p*BPH2*O2_CO2_Invert_1p)/(BPH3+K4p_CO2_Invert_1p*BPH4*O2_CO2_Invert_1p);
real SHbCO2_CO2_Invert_1p(t,x) dimensionless;
     SHbCO2_CO2_Invert_1p = KHbCO2_CO2_Invert_1p*CO2_CO2_Invert_1p/(1+KHbCO2_CO2_Invert_1p*CO2_CO2_Invert_1p);
real KHbCO2_CO2_Invert_1m(t,x) 1/M;              // Apparent equilibrium constants of Hb with CO2_CO2_Invert_1m
     KHbCO2_CO2_Invert_1m = (K2p*BPH1+K3p*K4p_CO2_Invert_1m*BPH2*O2_CO2_Invert_1m)/(BPH3+K4p_CO2_Invert_1m*BPH4*O2_CO2_Invert_1m);
real SHbCO2_CO2_Invert_1m(t,x) dimensionless;
     SHbCO2_CO2_Invert_1m = KHbCO2_CO2_Invert_1m*CO2_CO2_Invert_1m/(1+KHbCO2_CO2_Invert_1m*CO2_CO2_Invert_1m);
real KHbCO2_CO2_Invert_2c(t,x) 1/M;              // Apparent equilibrium constants of Hb with CO2_CO2_Invert_2c
     KHbCO2_CO2_Invert_2c = (K2p*BPH1+K3p*K4p_CO2_Invert_2c*BPH2*O2_CO2_Invert_2c)/(BPH3+K4p_CO2_Invert_2c*BPH4*O2_CO2_Invert_2c);
real SHbCO2_CO2_Invert_2c(t,x) dimensionless;
     SHbCO2_CO2_Invert_2c = KHbCO2_CO2_Invert_2c*CO2_CO2_Invert_2c/(1+KHbCO2_CO2_Invert_2c*CO2_CO2_Invert_2c);
real KHbCO2_CO2_Invert_2p(t,x) 1/M;              // Apparent equilibrium constants of Hb with CO2_CO2_Invert_2p
     KHbCO2_CO2_Invert_2p = (K2p*BPH1+K3p*K4p_CO2_Invert_2p*BPH2*O2_CO2_Invert_2p)/(BPH3+K4p_CO2_Invert_2p*BPH4*O2_CO2_Invert_2p);
real SHbCO2_CO2_Invert_2p(t,x) dimensionless;
     SHbCO2_CO2_Invert_2p = KHbCO2_CO2_Invert_2p*CO2_CO2_Invert_2p/(1+KHbCO2_CO2_Invert_2p*CO2_CO2_Invert_2p);
real KHbCO2_CO2_Invert_2m(t,x) 1/M;              // Apparent equilibrium constants of Hb with CO2_CO2_Invert_2m
     KHbCO2_CO2_Invert_2m = (K2p*BPH1+K3p*K4p_CO2_Invert_2m*BPH2*O2_CO2_Invert_2m)/(BPH3+K4p_CO2_Invert_2m*BPH4*O2_CO2_Invert_2m);
real SHbCO2_CO2_Invert_2m(t,x) dimensionless;
     SHbCO2_CO2_Invert_2m = KHbCO2_CO2_Invert_2m*CO2_CO2_Invert_2m/(1+KHbCO2_CO2_Invert_2m*CO2_CO2_Invert_2m);
real KHbCO2_CO2_Invert_3c(t,x) 1/M;              // Apparent equilibrium constants of Hb with CO2_CO2_Invert_3c
     KHbCO2_CO2_Invert_3c = (K2p*BPH1+K3p*K4p_CO2_Invert_3c*BPH2*O2_CO2_Invert_3c)/(BPH3+K4p_CO2_Invert_3c*BPH4*O2_CO2_Invert_3c);
real SHbCO2_CO2_Invert_3c(t,x) dimensionless;
     SHbCO2_CO2_Invert_3c = KHbCO2_CO2_Invert_3c*CO2_CO2_Invert_3c/(1+KHbCO2_CO2_Invert_3c*CO2_CO2_Invert_3c);
real KHbCO2_CO2_Invert_3p(t,x) 1/M;              // Apparent equilibrium constants of Hb with CO2_CO2_Invert_3p
     KHbCO2_CO2_Invert_3p = (K2p*BPH1+K3p*K4p_CO2_Invert_3p*BPH2*O2_CO2_Invert_3p)/(BPH3+K4p_CO2_Invert_3p*BPH4*O2_CO2_Invert_3p);
real SHbCO2_CO2_Invert_3p(t,x) dimensionless;
     SHbCO2_CO2_Invert_3p = KHbCO2_CO2_Invert_3p*CO2_CO2_Invert_3p/(1+KHbCO2_CO2_Invert_3p*CO2_CO2_Invert_3p);
real KHbCO2_CO2_Invert_3m(t,x) 1/M;              // Apparent equilibrium constants of Hb with CO2_CO2_Invert_3m
     KHbCO2_CO2_Invert_3m = (K2p*BPH1+K3p*K4p_CO2_Invert_3m*BPH2*O2_CO2_Invert_3m)/(BPH3+K4p_CO2_Invert_3m*BPH4*O2_CO2_Invert_3m);
real SHbCO2_CO2_Invert_3m(t,x) dimensionless;
     SHbCO2_CO2_Invert_3m = KHbCO2_CO2_Invert_3m*CO2_CO2_Invert_3m/(1+KHbCO2_CO2_Invert_3m*CO2_CO2_Invert_3m);
real KHbCO2_CO2_Invert_4c(t,x) 1/M;              // Apparent equilibrium constants of Hb with CO2_CO2_Invert_4c
     KHbCO2_CO2_Invert_4c = (K2p*BPH1+K3p*K4p_CO2_Invert_4c*BPH2*O2_CO2_Invert_4c)/(BPH3+K4p_CO2_Invert_4c*BPH4*O2_CO2_Invert_4c);
real SHbCO2_CO2_Invert_4c(t,x) dimensionless;
     SHbCO2_CO2_Invert_4c = KHbCO2_CO2_Invert_4c*CO2_CO2_Invert_4c/(1+KHbCO2_CO2_Invert_4c*CO2_CO2_Invert_4c);
real KHbCO2_CO2_Invert_4p(t,x) 1/M;              // Apparent equilibrium constants of Hb with CO2_CO2_Invert_4p
     KHbCO2_CO2_Invert_4p = (K2p*BPH1+K3p*K4p_CO2_Invert_4p*BPH2*O2_CO2_Invert_4p)/(BPH3+K4p_CO2_Invert_4p*BPH4*O2_CO2_Invert_4p);
real SHbCO2_CO2_Invert_4p(t,x) dimensionless;
     SHbCO2_CO2_Invert_4p = KHbCO2_CO2_Invert_4p*CO2_CO2_Invert_4p/(1+KHbCO2_CO2_Invert_4p*CO2_CO2_Invert_4p);
real KHbCO2_CO2_Invert_4m(t,x) 1/M;              // Apparent equilibrium constants of Hb with CO2_CO2_Invert_4m
     KHbCO2_CO2_Invert_4m = (K2p*BPH1+K3p*K4p_CO2_Invert_4m*BPH2*O2_CO2_Invert_4m)/(BPH3+K4p_CO2_Invert_4m*BPH4*O2_CO2_Invert_4m);
real SHbCO2_CO2_Invert_4m(t,x) dimensionless;
     SHbCO2_CO2_Invert_4m = KHbCO2_CO2_Invert_4m*CO2_CO2_Invert_4m/(1+KHbCO2_CO2_Invert_4m*CO2_CO2_Invert_4m);
real KHbCO2_CO2_Invert_5c(t,x) 1/M;              // Apparent equilibrium constants of Hb with CO2_CO2_Invert_5c
     KHbCO2_CO2_Invert_5c = (K2p*BPH1+K3p*K4p_CO2_Invert_5c*BPH2*O2_CO2_Invert_5c)/(BPH3+K4p_CO2_Invert_5c*BPH4*O2_CO2_Invert_5c);
real SHbCO2_CO2_Invert_5c(t,x) dimensionless;
     SHbCO2_CO2_Invert_5c = KHbCO2_CO2_Invert_5c*CO2_CO2_Invert_5c/(1+KHbCO2_CO2_Invert_5c*CO2_CO2_Invert_5c);
real KHbCO2_CO2_Invert_5p(t,x) 1/M;              // Apparent equilibrium constants of Hb with CO2_CO2_Invert_5p
     KHbCO2_CO2_Invert_5p = (K2p*BPH1+K3p*K4p_CO2_Invert_5p*BPH2*O2_CO2_Invert_5p)/(BPH3+K4p_CO2_Invert_5p*BPH4*O2_CO2_Invert_5p);
real SHbCO2_CO2_Invert_5p(t,x) dimensionless;
     SHbCO2_CO2_Invert_5p = KHbCO2_CO2_Invert_5p*CO2_CO2_Invert_5p/(1+KHbCO2_CO2_Invert_5p*CO2_CO2_Invert_5p);
real KHbCO2_CO2_Invert_5m(t,x) 1/M;              // Apparent equilibrium constants of Hb with CO2_CO2_Invert_5m
     KHbCO2_CO2_Invert_5m = (K2p*BPH1+K3p*K4p_CO2_Invert_5m*BPH2*O2_CO2_Invert_5m)/(BPH3+K4p_CO2_Invert_5m*BPH4*O2_CO2_Invert_5m);
real SHbCO2_CO2_Invert_5m(t,x) dimensionless;
     SHbCO2_CO2_Invert_5m = KHbCO2_CO2_Invert_5m*CO2_CO2_Invert_5m/(1+KHbCO2_CO2_Invert_5m*CO2_CO2_Invert_5m);
real KHbCO2_CO2_Invert_6c(t,x) 1/M;              // Apparent equilibrium constants of Hb with CO2_CO2_Invert_6c
     KHbCO2_CO2_Invert_6c = (K2p*BPH1+K3p*K4p_CO2_Invert_6c*BPH2*O2_CO2_Invert_6c)/(BPH3+K4p_CO2_Invert_6c*BPH4*O2_CO2_Invert_6c);
real SHbCO2_CO2_Invert_6c(t,x) dimensionless;
     SHbCO2_CO2_Invert_6c = KHbCO2_CO2_Invert_6c*CO2_CO2_Invert_6c/(1+KHbCO2_CO2_Invert_6c*CO2_CO2_Invert_6c);
real KHbCO2_CO2_Invert_6p(t,x) 1/M;              // Apparent equilibrium constants of Hb with CO2_CO2_Invert_6p
     KHbCO2_CO2_Invert_6p = (K2p*BPH1+K3p*K4p_CO2_Invert_6p*BPH2*O2_CO2_Invert_6p)/(BPH3+K4p_CO2_Invert_6p*BPH4*O2_CO2_Invert_6p);
real SHbCO2_CO2_Invert_6p(t,x) dimensionless;
     SHbCO2_CO2_Invert_6p = KHbCO2_CO2_Invert_6p*CO2_CO2_Invert_6p/(1+KHbCO2_CO2_Invert_6p*CO2_CO2_Invert_6p);
real KHbCO2_CO2_Invert_6m(t,x) 1/M;              // Apparent equilibrium constants of Hb with CO2_CO2_Invert_6m
     KHbCO2_CO2_Invert_6m = (K2p*BPH1+K3p*K4p_CO2_Invert_6m*BPH2*O2_CO2_Invert_6m)/(BPH3+K4p_CO2_Invert_6m*BPH4*O2_CO2_Invert_6m);
real SHbCO2_CO2_Invert_6m(t,x) dimensionless;
     SHbCO2_CO2_Invert_6m = KHbCO2_CO2_Invert_6m*CO2_CO2_Invert_6m/(1+KHbCO2_CO2_Invert_6m*CO2_CO2_Invert_6m);
real CO2free_CO2_Invert_1c(t,x) = Wrbc*CO2_CO2_Invert_1c;  // M (mol CO2_CO2_Invert_1c per L rbc)
real CO2bicarb_CO2_Invert_1c(t,x) = (Wrbc*Rrbc)*(K1*CO2_CO2_Invert_1c/Hrbc); // M (mol CO2_CO2_Invert_1c per L rbc)
real CO2bound_CO2_Invert_1c(t,x) = 4*Hbrbc*SHbCO2_CO2_Invert_1c; // M (mol CO2_CO2_Invert_1c per L rbc)
real CO2tot1_CO2_Invert_1c(t,x) = CO2free_CO2_Invert_1c+CO2bound_CO2_Invert_1c; // M (mol CO2_CO2_Invert_1c per L rbc)
real CO2tot2_CO2_Invert_1c(t,x) = CO2free_CO2_Invert_1c+CO2bicarb_CO2_Invert_1c+CO2bound_CO2_Invert_1c; // M (mol CO2_CO2_Invert_1c per L rbc)
real CO2free_CO2_Invert_1p(t,x) = Wrbc*CO2_CO2_Invert_1p;  // M (mol CO2_CO2_Invert_1p per L rbc)
real CO2bicarb_CO2_Invert_1p(t,x) = (Wrbc*Rrbc)*(K1*CO2_CO2_Invert_1p/Hrbc); // M (mol CO2_CO2_Invert_1p per L rbc)
real CO2bound_CO2_Invert_1p(t,x) = 4*Hbrbc*SHbCO2_CO2_Invert_1p; // M (mol CO2_CO2_Invert_1p per L rbc)
real CO2tot1_CO2_Invert_1p(t,x) = CO2free_CO2_Invert_1p+CO2bound_CO2_Invert_1p; // M (mol CO2_CO2_Invert_1p per L rbc)
real CO2tot2_CO2_Invert_1p(t,x) = CO2free_CO2_Invert_1p+CO2bicarb_CO2_Invert_1p+CO2bound_CO2_Invert_1p; // M (mol CO2_CO2_Invert_1p per L rbc)
real CO2free_CO2_Invert_1m(t,x) = Wrbc*CO2_CO2_Invert_1m;  // M (mol CO2_CO2_Invert_1m per L rbc)
real CO2bicarb_CO2_Invert_1m(t,x) = (Wrbc*Rrbc)*(K1*CO2_CO2_Invert_1m/Hrbc); // M (mol CO2_CO2_Invert_1m per L rbc)
real CO2bound_CO2_Invert_1m(t,x) = 4*Hbrbc*SHbCO2_CO2_Invert_1m; // M (mol CO2_CO2_Invert_1m per L rbc)
real CO2tot1_CO2_Invert_1m(t,x) = CO2free_CO2_Invert_1m+CO2bound_CO2_Invert_1m; // M (mol CO2_CO2_Invert_1m per L rbc)
real CO2tot2_CO2_Invert_1m(t,x) = CO2free_CO2_Invert_1m+CO2bicarb_CO2_Invert_1m+CO2bound_CO2_Invert_1m; // M (mol CO2_CO2_Invert_1m per L rbc)
real CO2free_CO2_Invert_2c(t,x) = Wrbc*CO2_CO2_Invert_2c;  // M (mol CO2_CO2_Invert_2c per L rbc)
real CO2bicarb_CO2_Invert_2c(t,x) = (Wrbc*Rrbc)*(K1*CO2_CO2_Invert_2c/Hrbc); // M (mol CO2_CO2_Invert_2c per L rbc)
real CO2bound_CO2_Invert_2c(t,x) = 4*Hbrbc*SHbCO2_CO2_Invert_2c; // M (mol CO2_CO2_Invert_2c per L rbc)
real CO2tot1_CO2_Invert_2c(t,x) = CO2free_CO2_Invert_2c+CO2bound_CO2_Invert_2c; // M (mol CO2_CO2_Invert_2c per L rbc)
real CO2tot2_CO2_Invert_2c(t,x) = CO2free_CO2_Invert_2c+CO2bicarb_CO2_Invert_2c+CO2bound_CO2_Invert_2c; // M (mol CO2_CO2_Invert_2c per L rbc)
real CO2free_CO2_Invert_2p(t,x) = Wrbc*CO2_CO2_Invert_2p;  // M (mol CO2_CO2_Invert_2p per L rbc)
real CO2bicarb_CO2_Invert_2p(t,x) = (Wrbc*Rrbc)*(K1*CO2_CO2_Invert_2p/Hrbc); // M (mol CO2_CO2_Invert_2p per L rbc)
real CO2bound_CO2_Invert_2p(t,x) = 4*Hbrbc*SHbCO2_CO2_Invert_2p; // M (mol CO2_CO2_Invert_2p per L rbc)
real CO2tot1_CO2_Invert_2p(t,x) = CO2free_CO2_Invert_2p+CO2bound_CO2_Invert_2p; // M (mol CO2_CO2_Invert_2p per L rbc)
real CO2tot2_CO2_Invert_2p(t,x) = CO2free_CO2_Invert_2p+CO2bicarb_CO2_Invert_2p+CO2bound_CO2_Invert_2p; // M (mol CO2_CO2_Invert_2p per L rbc)
real CO2free_CO2_Invert_2m(t,x) = Wrbc*CO2_CO2_Invert_2m;  // M (mol CO2_CO2_Invert_2m per L rbc)
real CO2bicarb_CO2_Invert_2m(t,x) = (Wrbc*Rrbc)*(K1*CO2_CO2_Invert_2m/Hrbc); // M (mol CO2_CO2_Invert_2m per L rbc)
real CO2bound_CO2_Invert_2m(t,x) = 4*Hbrbc*SHbCO2_CO2_Invert_2m; // M (mol CO2_CO2_Invert_2m per L rbc)
real CO2tot1_CO2_Invert_2m(t,x) = CO2free_CO2_Invert_2m+CO2bound_CO2_Invert_2m; // M (mol CO2_CO2_Invert_2m per L rbc)
real CO2tot2_CO2_Invert_2m(t,x) = CO2free_CO2_Invert_2m+CO2bicarb_CO2_Invert_2m+CO2bound_CO2_Invert_2m; // M (mol CO2_CO2_Invert_2m per L rbc)
real CO2free_CO2_Invert_3c(t,x) = Wrbc*CO2_CO2_Invert_3c;  // M (mol CO2_CO2_Invert_3c per L rbc)
real CO2bicarb_CO2_Invert_3c(t,x) = (Wrbc*Rrbc)*(K1*CO2_CO2_Invert_3c/Hrbc); // M (mol CO2_CO2_Invert_3c per L rbc)
real CO2bound_CO2_Invert_3c(t,x) = 4*Hbrbc*SHbCO2_CO2_Invert_3c; // M (mol CO2_CO2_Invert_3c per L rbc)
real CO2tot1_CO2_Invert_3c(t,x) = CO2free_CO2_Invert_3c+CO2bound_CO2_Invert_3c; // M (mol CO2_CO2_Invert_3c per L rbc)
real CO2tot2_CO2_Invert_3c(t,x) = CO2free_CO2_Invert_3c+CO2bicarb_CO2_Invert_3c+CO2bound_CO2_Invert_3c; // M (mol CO2_CO2_Invert_3c per L rbc)
real CO2free_CO2_Invert_3p(t,x) = Wrbc*CO2_CO2_Invert_3p;  // M (mol CO2_CO2_Invert_3p per L rbc)
real CO2bicarb_CO2_Invert_3p(t,x) = (Wrbc*Rrbc)*(K1*CO2_CO2_Invert_3p/Hrbc); // M (mol CO2_CO2_Invert_3p per L rbc)
real CO2bound_CO2_Invert_3p(t,x) = 4*Hbrbc*SHbCO2_CO2_Invert_3p; // M (mol CO2_CO2_Invert_3p per L rbc)
real CO2tot1_CO2_Invert_3p(t,x) = CO2free_CO2_Invert_3p+CO2bound_CO2_Invert_3p; // M (mol CO2_CO2_Invert_3p per L rbc)
real CO2tot2_CO2_Invert_3p(t,x) = CO2free_CO2_Invert_3p+CO2bicarb_CO2_Invert_3p+CO2bound_CO2_Invert_3p; // M (mol CO2_CO2_Invert_3p per L rbc)
real CO2free_CO2_Invert_3m(t,x) = Wrbc*CO2_CO2_Invert_3m;  // M (mol CO2_CO2_Invert_3m per L rbc)
real CO2bicarb_CO2_Invert_3m(t,x) = (Wrbc*Rrbc)*(K1*CO2_CO2_Invert_3m/Hrbc); // M (mol CO2_CO2_Invert_3m per L rbc)
real CO2bound_CO2_Invert_3m(t,x) = 4*Hbrbc*SHbCO2_CO2_Invert_3m; // M (mol CO2_CO2_Invert_3m per L rbc)
real CO2tot1_CO2_Invert_3m(t,x) = CO2free_CO2_Invert_3m+CO2bound_CO2_Invert_3m; // M (mol CO2_CO2_Invert_3m per L rbc)
real CO2tot2_CO2_Invert_3m(t,x) = CO2free_CO2_Invert_3m+CO2bicarb_CO2_Invert_3m+CO2bound_CO2_Invert_3m; // M (mol CO2_CO2_Invert_3m per L rbc)
real CO2free_CO2_Invert_4c(t,x) = Wrbc*CO2_CO2_Invert_4c;  // M (mol CO2_CO2_Invert_4c per L rbc)
real CO2bicarb_CO2_Invert_4c(t,x) = (Wrbc*Rrbc)*(K1*CO2_CO2_Invert_4c/Hrbc); // M (mol CO2_CO2_Invert_4c per L rbc)
real CO2bound_CO2_Invert_4c(t,x) = 4*Hbrbc*SHbCO2_CO2_Invert_4c; // M (mol CO2_CO2_Invert_4c per L rbc)
real CO2tot1_CO2_Invert_4c(t,x) = CO2free_CO2_Invert_4c+CO2bound_CO2_Invert_4c; // M (mol CO2_CO2_Invert_4c per L rbc)
real CO2tot2_CO2_Invert_4c(t,x) = CO2free_CO2_Invert_4c+CO2bicarb_CO2_Invert_4c+CO2bound_CO2_Invert_4c; // M (mol CO2_CO2_Invert_4c per L rbc)
real CO2free_CO2_Invert_4p(t,x) = Wrbc*CO2_CO2_Invert_4p;  // M (mol CO2_CO2_Invert_4p per L rbc)
real CO2bicarb_CO2_Invert_4p(t,x) = (Wrbc*Rrbc)*(K1*CO2_CO2_Invert_4p/Hrbc); // M (mol CO2_CO2_Invert_4p per L rbc)
real CO2bound_CO2_Invert_4p(t,x) = 4*Hbrbc*SHbCO2_CO2_Invert_4p; // M (mol CO2_CO2_Invert_4p per L rbc)
real CO2tot1_CO2_Invert_4p(t,x) = CO2free_CO2_Invert_4p+CO2bound_CO2_Invert_4p; // M (mol CO2_CO2_Invert_4p per L rbc)
real CO2tot2_CO2_Invert_4p(t,x) = CO2free_CO2_Invert_4p+CO2bicarb_CO2_Invert_4p+CO2bound_CO2_Invert_4p; // M (mol CO2_CO2_Invert_4p per L rbc)
real CO2free_CO2_Invert_4m(t,x) = Wrbc*CO2_CO2_Invert_4m;  // M (mol CO2_CO2_Invert_4m per L rbc)
real CO2bicarb_CO2_Invert_4m(t,x) = (Wrbc*Rrbc)*(K1*CO2_CO2_Invert_4m/Hrbc); // M (mol CO2_CO2_Invert_4m per L rbc)
real CO2bound_CO2_Invert_4m(t,x) = 4*Hbrbc*SHbCO2_CO2_Invert_4m; // M (mol CO2_CO2_Invert_4m per L rbc)
real CO2tot1_CO2_Invert_4m(t,x) = CO2free_CO2_Invert_4m+CO2bound_CO2_Invert_4m; // M (mol CO2_CO2_Invert_4m per L rbc)
real CO2tot2_CO2_Invert_4m(t,x) = CO2free_CO2_Invert_4m+CO2bicarb_CO2_Invert_4m+CO2bound_CO2_Invert_4m; // M (mol CO2_CO2_Invert_4m per L rbc)
real CO2free_CO2_Invert_5c(t,x) = Wrbc*CO2_CO2_Invert_5c;  // M (mol CO2_CO2_Invert_5c per L rbc)
real CO2bicarb_CO2_Invert_5c(t,x) = (Wrbc*Rrbc)*(K1*CO2_CO2_Invert_5c/Hrbc); // M (mol CO2_CO2_Invert_5c per L rbc)
real CO2bound_CO2_Invert_5c(t,x) = 4*Hbrbc*SHbCO2_CO2_Invert_5c; // M (mol CO2_CO2_Invert_5c per L rbc)
real CO2tot1_CO2_Invert_5c(t,x) = CO2free_CO2_Invert_5c+CO2bound_CO2_Invert_5c; // M (mol CO2_CO2_Invert_5c per L rbc)
real CO2tot2_CO2_Invert_5c(t,x) = CO2free_CO2_Invert_5c+CO2bicarb_CO2_Invert_5c+CO2bound_CO2_Invert_5c; // M (mol CO2_CO2_Invert_5c per L rbc)
real CO2free_CO2_Invert_5p(t,x) = Wrbc*CO2_CO2_Invert_5p;  // M (mol CO2_CO2_Invert_5p per L rbc)
real CO2bicarb_CO2_Invert_5p(t,x) = (Wrbc*Rrbc)*(K1*CO2_CO2_Invert_5p/Hrbc); // M (mol CO2_CO2_Invert_5p per L rbc)
real CO2bound_CO2_Invert_5p(t,x) = 4*Hbrbc*SHbCO2_CO2_Invert_5p; // M (mol CO2_CO2_Invert_5p per L rbc)
real CO2tot1_CO2_Invert_5p(t,x) = CO2free_CO2_Invert_5p+CO2bound_CO2_Invert_5p; // M (mol CO2_CO2_Invert_5p per L rbc)
real CO2tot2_CO2_Invert_5p(t,x) = CO2free_CO2_Invert_5p+CO2bicarb_CO2_Invert_5p+CO2bound_CO2_Invert_5p; // M (mol CO2_CO2_Invert_5p per L rbc)
real CO2free_CO2_Invert_5m(t,x) = Wrbc*CO2_CO2_Invert_5m;  // M (mol CO2_CO2_Invert_5m per L rbc)
real CO2bicarb_CO2_Invert_5m(t,x) = (Wrbc*Rrbc)*(K1*CO2_CO2_Invert_5m/Hrbc); // M (mol CO2_CO2_Invert_5m per L rbc)
real CO2bound_CO2_Invert_5m(t,x) = 4*Hbrbc*SHbCO2_CO2_Invert_5m; // M (mol CO2_CO2_Invert_5m per L rbc)
real CO2tot1_CO2_Invert_5m(t,x) = CO2free_CO2_Invert_5m+CO2bound_CO2_Invert_5m; // M (mol CO2_CO2_Invert_5m per L rbc)
real CO2tot2_CO2_Invert_5m(t,x) = CO2free_CO2_Invert_5m+CO2bicarb_CO2_Invert_5m+CO2bound_CO2_Invert_5m; // M (mol CO2_CO2_Invert_5m per L rbc)
real CO2free_CO2_Invert_6c(t,x) = Wrbc*CO2_CO2_Invert_6c;  // M (mol CO2_CO2_Invert_6c per L rbc)
real CO2bicarb_CO2_Invert_6c(t,x) = (Wrbc*Rrbc)*(K1*CO2_CO2_Invert_6c/Hrbc); // M (mol CO2_CO2_Invert_6c per L rbc)
real CO2bound_CO2_Invert_6c(t,x) = 4*Hbrbc*SHbCO2_CO2_Invert_6c; // M (mol CO2_CO2_Invert_6c per L rbc)
real CO2tot1_CO2_Invert_6c(t,x) = CO2free_CO2_Invert_6c+CO2bound_CO2_Invert_6c; // M (mol CO2_CO2_Invert_6c per L rbc)
real CO2tot2_CO2_Invert_6c(t,x) = CO2free_CO2_Invert_6c+CO2bicarb_CO2_Invert_6c+CO2bound_CO2_Invert_6c; // M (mol CO2_CO2_Invert_6c per L rbc)
real CO2free_CO2_Invert_6p(t,x) = Wrbc*CO2_CO2_Invert_6p;  // M (mol CO2_CO2_Invert_6p per L rbc)
real CO2bicarb_CO2_Invert_6p(t,x) = (Wrbc*Rrbc)*(K1*CO2_CO2_Invert_6p/Hrbc); // M (mol CO2_CO2_Invert_6p per L rbc)
real CO2bound_CO2_Invert_6p(t,x) = 4*Hbrbc*SHbCO2_CO2_Invert_6p; // M (mol CO2_CO2_Invert_6p per L rbc)
real CO2tot1_CO2_Invert_6p(t,x) = CO2free_CO2_Invert_6p+CO2bound_CO2_Invert_6p; // M (mol CO2_CO2_Invert_6p per L rbc)
real CO2tot2_CO2_Invert_6p(t,x) = CO2free_CO2_Invert_6p+CO2bicarb_CO2_Invert_6p+CO2bound_CO2_Invert_6p; // M (mol CO2_CO2_Invert_6p per L rbc)
real CO2free_CO2_Invert_6m(t,x) = Wrbc*CO2_CO2_Invert_6m;  // M (mol CO2_CO2_Invert_6m per L rbc)
real CO2bicarb_CO2_Invert_6m(t,x) = (Wrbc*Rrbc)*(K1*CO2_CO2_Invert_6m/Hrbc); // M (mol CO2_CO2_Invert_6m per L rbc)
real CO2bound_CO2_Invert_6m(t,x) = 4*Hbrbc*SHbCO2_CO2_Invert_6m; // M (mol CO2_CO2_Invert_6m per L rbc)
real CO2tot1_CO2_Invert_6m(t,x) = CO2free_CO2_Invert_6m+CO2bound_CO2_Invert_6m; // M (mol CO2_CO2_Invert_6m per L rbc)
real CO2tot2_CO2_Invert_6m(t,x) = CO2free_CO2_Invert_6m+CO2bicarb_CO2_Invert_6m+CO2bound_CO2_Invert_6m; // M (mol CO2_CO2_Invert_6m per L rbc)
real pCO2old_CO2_Invert_1(t,x) = pCO2new_CO2_Invert_0; // Initial guess for pCO2, otherwise it was previous iteration's value (pCO2new_CO2_Invert_1)
// Newton-Raphson's iterative method for computation of pCO2 from TotCO2
// iter:
  pCO2c_CO2_Invert_1 = pCO2old_CO2_Invert_1;
  pCO2p_CO2_Invert_1 = pCO2old_CO2_Invert_1 + 1e-2*pCO2old_CO2_Invert_1;
  pCO2m_CO2_Invert_1 = pCO2old_CO2_Invert_1 - 1e-2*pCO2old_CO2_Invert_1;
//    Inputc = [pCO2c1,pO2,pHrbc,DPGrbc,Temp,Hbrbc,Hct];
//    [Outputc] = SHbO2CO2_EJAP2016(Inputc); CO2totc1 = Outputc{2}(1);
  dCO2totCO2_Invert_1c = (CO2tot1_CO2_Invert_1p-CO2tot1_CO2_Invert_1m)/(pCO2p_CO2_Invert_1-pCO2m_CO2_Invert_1);  // Derivative by central difference
  funcCO2_CO2_Invert_1 = TotCO2_Invert_in - CO2tot1_CO2_Invert_1c;   // Function to be solved (find pCO2 such that funcO2 = 0) 
  dfuncCO2_CO2_Invert_1 = -dCO2totCO2_Invert_1c;        // As TotCO2 is an input (constant), dTotCO2 = 0
  pCO2new_CO2_Invert_1 = (pCO2old_CO2_Invert_1 - funcCO2_CO2_Invert_1/dfuncCO2_CO2_Invert_1) ;// Newton-Raphson iterative formula
real errCO2_CO2_Invert_1(t,x) = abs(pCO2new_CO2_Invert_1 - pCO2old_CO2_Invert_1)/pCO2new_CO2_Invert_1;   // use this as a check.
real pCO2old_CO2_Invert_2(t,x) = pCO2new_CO2_Invert_1; // Initial guess for pCO2, otherwise it was previous iteration's value (pCO2new_CO2_Invert_2)
  pCO2c_CO2_Invert_2 = pCO2old_CO2_Invert_2;
  pCO2p_CO2_Invert_2 = pCO2old_CO2_Invert_2 + 1e-2*pCO2old_CO2_Invert_2;
  pCO2m_CO2_Invert_2 = pCO2old_CO2_Invert_2 - 1e-2*pCO2old_CO2_Invert_2;
  dCO2totCO2_Invert_2c = (CO2tot1_CO2_Invert_2p-CO2tot1_CO2_Invert_2m)/(pCO2p_CO2_Invert_2-pCO2m_CO2_Invert_2);  // Derivative by central difference
  funcCO2_CO2_Invert_2 = TotCO2_Invert_in - CO2tot1_CO2_Invert_2c;   // Function to be solved (find pCO2 such that funcO2 = 0) 
  dfuncCO2_CO2_Invert_2 = -dCO2totCO2_Invert_2c;        // As TotCO2 is an input (constant), dTotCO2 = 0
  pCO2new_CO2_Invert_2 = (pCO2old_CO2_Invert_2 - funcCO2_CO2_Invert_2/dfuncCO2_CO2_Invert_2) ;// Newton-Raphson iterative formula
real errCO2_CO2_Invert_2(t,x) = abs(pCO2new_CO2_Invert_2 - pCO2old_CO2_Invert_2)/pCO2new_CO2_Invert_2;   // use this as a check.
real pCO2old_CO2_Invert_3(t,x) = pCO2new_CO2_Invert_2; // Initial guess for pCO2, otherwise it was previous iteration's value (pCO2new_CO2_Invert_3)
  pCO2c_CO2_Invert_3 = pCO2old_CO2_Invert_3;
  pCO2p_CO2_Invert_3 = pCO2old_CO2_Invert_3 + 1e-2*pCO2old_CO2_Invert_3;
  pCO2m_CO2_Invert_3 = pCO2old_CO2_Invert_3 - 1e-2*pCO2old_CO2_Invert_3;
  dCO2totCO2_Invert_3c = (CO2tot1_CO2_Invert_3p-CO2tot1_CO2_Invert_3m)/(pCO2p_CO2_Invert_3-pCO2m_CO2_Invert_3);  // Derivative by central difference
  funcCO2_CO2_Invert_3 = TotCO2_Invert_in - CO2tot1_CO2_Invert_3c;   // Function to be solved (find pCO2 such that funcO2 = 0) 
  dfuncCO2_CO2_Invert_3 = -dCO2totCO2_Invert_3c;        // As TotCO2 is an input (constant), dTotCO2 = 0
  pCO2new_CO2_Invert_3 = (pCO2old_CO2_Invert_3 - funcCO2_CO2_Invert_3/dfuncCO2_CO2_Invert_3) ;// Newton-Raphson iterative formula
real errCO2_CO2_Invert_3(t,x) = abs(pCO2new_CO2_Invert_3 - pCO2old_CO2_Invert_3)/pCO2new_CO2_Invert_3;   // use this as a check.
real pCO2old_CO2_Invert_4(t,x) = pCO2new_CO2_Invert_3; // Initial guess for pCO2, otherwise it was previous iteration's value (pCO2new_CO2_Invert_4)
  pCO2c_CO2_Invert_4 = pCO2old_CO2_Invert_4;
  pCO2p_CO2_Invert_4 = pCO2old_CO2_Invert_4 + 1e-2*pCO2old_CO2_Invert_4;
  pCO2m_CO2_Invert_4 = pCO2old_CO2_Invert_4 - 1e-2*pCO2old_CO2_Invert_4;
  dCO2totCO2_Invert_4c = (CO2tot1_CO2_Invert_4p-CO2tot1_CO2_Invert_4m)/(pCO2p_CO2_Invert_4-pCO2m_CO2_Invert_4);  // Derivative by central difference
  funcCO2_CO2_Invert_4 = TotCO2_Invert_in - CO2tot1_CO2_Invert_4c;   // Function to be solved (find pCO2 such that funcO2 = 0) 
  dfuncCO2_CO2_Invert_4 = -dCO2totCO2_Invert_4c;        // As TotCO2 is an input (constant), dTotCO2 = 0
  pCO2new_CO2_Invert_4 = (pCO2old_CO2_Invert_4 - funcCO2_CO2_Invert_4/dfuncCO2_CO2_Invert_4) ;// Newton-Raphson iterative formula
real errCO2_CO2_Invert_4(t,x) = abs(pCO2new_CO2_Invert_4 - pCO2old_CO2_Invert_4)/pCO2new_CO2_Invert_4;   // use this as a check.
real pCO2old_CO2_Invert_5(t,x) = pCO2new_CO2_Invert_4; // Initial guess for pCO2, otherwise it was previous iteration's value (pCO2new_CO2_Invert_5)
  pCO2c_CO2_Invert_5 = pCO2old_CO2_Invert_5;
  pCO2p_CO2_Invert_5 = pCO2old_CO2_Invert_5 + 1e-2*pCO2old_CO2_Invert_5;
  pCO2m_CO2_Invert_5 = pCO2old_CO2_Invert_5 - 1e-2*pCO2old_CO2_Invert_5;
  dCO2totCO2_Invert_5c = (CO2tot1_CO2_Invert_5p-CO2tot1_CO2_Invert_5m)/(pCO2p_CO2_Invert_5-pCO2m_CO2_Invert_5);  // Derivative by central difference
  funcCO2_CO2_Invert_5 = TotCO2_Invert_in - CO2tot1_CO2_Invert_5c;   // Function to be solved (find pCO2 such that funcO2 = 0) 
  dfuncCO2_CO2_Invert_5 = -dCO2totCO2_Invert_5c;        // As TotCO2 is an input (constant), dTotCO2 = 0
  pCO2new_CO2_Invert_5 = (pCO2old_CO2_Invert_5 - funcCO2_CO2_Invert_5/dfuncCO2_CO2_Invert_5) ;// Newton-Raphson iterative formula
real errCO2_CO2_Invert_5(t,x) = abs(pCO2new_CO2_Invert_5 - pCO2old_CO2_Invert_5)/pCO2new_CO2_Invert_5;   // use this as a check.
real pCO2old_CO2_Invert_6(t,x) = pCO2new_CO2_Invert_5; // Initial guess for pCO2, otherwise it was previous iteration's value (pCO2new_CO2_Invert_6)
  pCO2c_CO2_Invert_6 = pCO2old_CO2_Invert_6;
  pCO2p_CO2_Invert_6 = pCO2old_CO2_Invert_6 + 1e-2*pCO2old_CO2_Invert_6;
  pCO2m_CO2_Invert_6 = pCO2old_CO2_Invert_6 - 1e-2*pCO2old_CO2_Invert_6;
  dCO2totCO2_Invert_6c = (CO2tot1_CO2_Invert_6p-CO2tot1_CO2_Invert_6m)/(pCO2p_CO2_Invert_6-pCO2m_CO2_Invert_6);  // Derivative by central difference
  funcCO2_CO2_Invert_6 = TotCO2_Invert_in - CO2tot1_CO2_Invert_6c;   // Function to be solved (find pCO2 such that funcO2 = 0) 
  dfuncCO2_CO2_Invert_6 = -dCO2totCO2_Invert_6c;        // As TotCO2 is an input (constant), dTotCO2 = 0
  pCO2new_CO2_Invert_6 = (pCO2old_CO2_Invert_6 - funcCO2_CO2_Invert_6/dfuncCO2_CO2_Invert_6) ;// Newton-Raphson iterative formula
real errCO2_CO2_Invert_6(t,x) = abs(pCO2new_CO2_Invert_6 - pCO2old_CO2_Invert_6)/pCO2new_CO2_Invert_6;   // use this as a check.
pCO2_out = pCO2new_CO2_Invert_6;              // Final Output of inversion
FCO2_out = alphaCO2Sys*pCO2_out; // 
SHbCO2_out = SHbCO2_CO2_Invert_6c;             // Final SHbCO2_CO2_Invert_6c from inversion
//%END CO2INVERT
// Set initial values:
Hct_in =0.45;
Temp_in=310;  // Physiological temperature (function of t,x)
TotCO2_Invert_in = 0.02827;
pO2_CO2_Invert_in = 90;    // Needed to determine pCO2
pHrbc_in = 7.238;
DPGrbc_in =.00465; 
Hbrbc_in = .0052;
//alphaCO2_in = 2.8472E-5;
}
// This MML file generated from TotCO2freeCO2Invert.mpc using MPC v1.01.
  
