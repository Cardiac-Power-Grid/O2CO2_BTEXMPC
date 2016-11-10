
// Calculate free O2 in blood based on total O2.
// Compute free O2 and CO2 (or PO2 and PCO2) from total O2 and CO2 (TO2 and TCO2)
// through numerical inversion method that uses SHbO2CO2_EJAP2016 routine.



import nsrunit;
unit conversion on;
 unit celsius = fundamental;

math TotO2_O2InvertMPC {

 realDomain t sec; t.min=0; t.max=300; t.delta=1;   // time domain
 realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid;  
//%START O2INVERT
real TotO2_Invert_in(t,x) M;   // Tot O2 used for inversion to get pO2 
real pCO2_Invert_in(t,x) mmHg; // pCO2 used for inversion to get pO2
real pHrbc_in(t,x);    // pHrbc used for inversion to get pO2/CO2
real DPGrbc_in(t,x) M; // DPGrbc used for inversion to get pO2/CO2
real Temp_in(t,x) K;   // Temp used for inversion to get pO2/CO2
real Hbrbc_in(t,x) M;  // Hbrbc used for inversion to get pO2/CO2
real Hct_in;           // Vol fraction of RBCs in blood (Hematocrit), needed?
real alphaO2Sys(t,x) M/mmHg; // solubility of O2 in water. Calc can be done in SHbO2CO2_EJAP2016.
real P50_in mmHg;   // Init pO2 guess (P50 for O2) for inversion to get pO2.
real pO2new_0(t,x) = P50_in;
real pCO2t0_in(x) mmHg;  // Initial pCO2 at t=0 for all x
real alphaCO2Sys(t,x) M/mmHg; // solubility of CO2 in plasma.
real pHrbc(t,x) dimensionless; real DPGrbc(t,x) M ; // pH and DPG conc in RBCs
real TempSHbO2CO2(t,x) K;          // Physiological temperature, Used in Dash2016SHbO2CO2MPC
real Hct dimensionless ;   // Volume fraction of RBCs in blood (Hematocrit), unitless
real mol2ml = 22400 mol/ml;          // Conversion factor from mol of gas to ml of gas at STP
real AmHbBl = Hct*(150 g)/(0.45 L);  // Amount of hemoglobin in gm per liter of blood, gm/L
real MWHb = 64100 g/mol;             // Molecular weight of hemoglobin, gm/mol
real HbBl M;
     HbBl = AmHbBl/MWHb;             // Concentration of hemoglobin in blood, mol/L or M
// TotO2_FreeO2_Inversion: calcs pO2 given TotO2, pCO2, alphaO2, 
real pO2_out(t,x) mmHg;  // free O2 in rbc
real FO2_out(t,x) M;     // free O2 in rbc
real SHbO2_out(t,x);    
// IterativeVars
real pO2c_1(t,x) mmHg;
real pO2p_1(t,x) mmHg;
real pO2m_1(t,x) mmHg;
real dO2tot1c(t,x);
real funcO2_1(t,x);
real dfuncO2_1(t,x);
real pO2new_1(t,x) mmHg;
real pO2c_2(t,x) mmHg;
real pO2p_2(t,x) mmHg;
real pO2m_2(t,x) mmHg;
real dO2tot2c(t,x);
real funcO2_2(t,x);
real dfuncO2_2(t,x);
real pO2new_2(t,x) mmHg;
real pO2c_3(t,x) mmHg;
real pO2p_3(t,x) mmHg;
real pO2m_3(t,x) mmHg;
real dO2tot3c(t,x);
real funcO2_3(t,x);
real dfuncO2_3(t,x);
real pO2new_3(t,x) mmHg;
real pO2c_4(t,x) mmHg;
real pO2p_4(t,x) mmHg;
real pO2m_4(t,x) mmHg;
real dO2tot4c(t,x);
real funcO2_4(t,x);
real dfuncO2_4(t,x);
real pO2new_4(t,x) mmHg;
real pO2c_5(t,x) mmHg;
real pO2p_5(t,x) mmHg;
real pO2m_5(t,x) mmHg;
real dO2tot5c(t,x);
real funcO2_5(t,x);
real dfuncO2_5(t,x);
real pO2new_5(t,x) mmHg;
real pO2c_6(t,x) mmHg;
real pO2p_6(t,x) mmHg;
real pO2m_6(t,x) mmHg;
real dO2tot6c(t,x);
real funcO2_6(t,x);
real dfuncO2_6(t,x);
real pO2new_6(t,x) mmHg;
real pO2c_7(t,x) mmHg;
real pO2p_7(t,x) mmHg;
real pO2m_7(t,x) mmHg;
real dO2tot7c(t,x);
real funcO2_7(t,x);
real dfuncO2_7(t,x);
real pO2new_7(t,x) mmHg;
real pO2c_8(t,x) mmHg;
real pO2p_8(t,x) mmHg;
real pO2m_8(t,x) mmHg;
real dO2tot8c(t,x);
real funcO2_8(t,x);
real dfuncO2_8(t,x);
real pO2new_8(t,x) mmHg;
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
real pO2_1c(t,x) mmHg;    // partial pressure in rbc
real pO2_1p(t,x) mmHg;    // partial pressure in rbc
real pO2_1m(t,x) mmHg;    // partial pressure in rbc
real pO2_2c(t,x) mmHg;    // partial pressure in rbc
real pO2_2p(t,x) mmHg;    // partial pressure in rbc
real pO2_2m(t,x) mmHg;    // partial pressure in rbc
real pO2_3c(t,x) mmHg;    // partial pressure in rbc
real pO2_3p(t,x) mmHg;    // partial pressure in rbc
real pO2_3m(t,x) mmHg;    // partial pressure in rbc
real pO2_4c(t,x) mmHg;    // partial pressure in rbc
real pO2_4p(t,x) mmHg;    // partial pressure in rbc
real pO2_4m(t,x) mmHg;    // partial pressure in rbc
real pO2_5c(t,x) mmHg;    // partial pressure in rbc
real pO2_5p(t,x) mmHg;    // partial pressure in rbc
real pO2_5m(t,x) mmHg;    // partial pressure in rbc
real pO2_6c(t,x) mmHg;    // partial pressure in rbc
real pO2_6p(t,x) mmHg;    // partial pressure in rbc
real pO2_6m(t,x) mmHg;    // partial pressure in rbc
real pO2_7c(t,x) mmHg;    // partial pressure in rbc
real pO2_7p(t,x) mmHg;    // partial pressure in rbc
real pO2_7m(t,x) mmHg;    // partial pressure in rbc
real pO2_8c(t,x) mmHg;    // partial pressure in rbc
real pO2_8p(t,x) mmHg;    // partial pressure in rbc
real pO2_8m(t,x) mmHg;    // partial pressure in rbc
// Assign an declare variable inputs to SHbO2CO2 calcs to get pO2 given TO2 in rbc
 Hct = Hct_in;
 pO2_1c=pO2c_1;
 pO2_1p=pO2p_1;
 pO2_1m=pO2m_1;
 pHrbc=pHrbc_in;
 DPGrbc = DPGrbc_in;
 TempSHbO2CO2 =Temp_in;  // Used in Dash2016SHbO2CO2MPC
realState pCO2_1(t,x);
  when (t=t.min) pCO2_1 = pCO2t0_in;
  event(t>t.min) {pCO2_1= pCO2_Invert_in;} // wait untill pCO2rbc calculated
 
 pO2_2c=pO2c_2;
 pO2_2p=pO2p_2;
 pO2_2m=pO2m_2;
realState pCO2_2(t,x);
  when (t=t.min) pCO2_2 = pCO2t0_in;
  event(t>t.min) {pCO2_2= pCO2_Invert_in;} // wait untill pCO2rbc calculated
 pO2_3c=pO2c_3;
 pO2_3p=pO2p_3;
 pO2_3m=pO2m_3;
realState pCO2_3(t,x);
  when (t=t.min) pCO2_3 = pCO2t0_in;
  event(t>t.min) {pCO2_3= pCO2_Invert_in;} // wait untill pCO2rbc calculated
 pO2_4c=pO2c_4;
 pO2_4p=pO2p_4;
 pO2_4m=pO2m_4;
realState pCO2_4(t,x);
  when (t=t.min) pCO2_4 = pCO2t0_in;
  event(t>t.min) {pCO2_4= pCO2_Invert_in;} // wait untill pCO2rbc calculated
 pO2_5c=pO2c_5;
 pO2_5p=pO2p_5;
 pO2_5m=pO2m_5;
realState pCO2_5(t,x);
  when (t=t.min) pCO2_5 = pCO2t0_in;
  event(t>t.min) {pCO2_5= pCO2_Invert_in;} // wait untill pCO2rbc calculated
 pO2_6c=pO2c_6;
 pO2_6p=pO2p_6;
 pO2_6m=pO2m_6;
realState pCO2_6(t,x);
  when (t=t.min) pCO2_6 = pCO2t0_in;
  event(t>t.min) {pCO2_6= pCO2_Invert_in;} // wait untill pCO2rbc calculated
 pO2_7c=pO2c_7;
 pO2_7p=pO2p_7;
 pO2_7m=pO2m_7;
realState pCO2_7(t,x);
  when (t=t.min) pCO2_7 = pCO2t0_in;
  event(t>t.min) {pCO2_7= pCO2_Invert_in;} // wait untill pCO2rbc calculated
 pO2_8c=pO2c_8;
 pO2_8p=pO2p_8;
 pO2_8m=pO2m_8;
realState pCO2_8(t,x);
  when (t=t.min) pCO2_8 = pCO2t0_in;
  event(t>t.min) {pCO2_8= pCO2_Invert_in;} // wait untill pCO2rbc calculated
// delpCO2_ is the same for all iterations as pCO2 is same for all when inverting TO2 to get pO2.
// Calculation of intermediate variables in the computations of SHbO2 and SHbCO2
      // Rrbc = Hpl/Hrbc = 10^-(pHpl-pHrbc)
real delpHrbc(t,x) = pHrbc-pHrbc0; 
real delpCO2_1c(t,x) = pCO2_1-pCO20; 
real delDPGrbc(t,x) = DPGrbc-DPGrbc0;
real Hrbc(t,x) M;
     Hrbc  = 10^(-pHrbc); 
real delTemp(t,x) = TempSHbO2CO2-Temp0;
real delpCO2_1p(t,x) = pCO2_1-pCO20; 
real delpCO2_1m(t,x) = pCO2_1-pCO20; 
real delpCO2_2c(t,x) = pCO2_2-pCO20; 
real delpCO2_2p(t,x) = pCO2_2-pCO20; 
real delpCO2_2m(t,x) = pCO2_2-pCO20; 
real delpCO2_3c(t,x) = pCO2_3-pCO20; 
real delpCO2_3p(t,x) = pCO2_3-pCO20; 
real delpCO2_3m(t,x) = pCO2_3-pCO20; 
real delpCO2_4c(t,x) = pCO2_4-pCO20; 
real delpCO2_4p(t,x) = pCO2_4-pCO20; 
real delpCO2_4m(t,x) = pCO2_4-pCO20; 
real delpCO2_5c(t,x) = pCO2_5-pCO20; 
real delpCO2_5p(t,x) = pCO2_5-pCO20; 
real delpCO2_5m(t,x) = pCO2_5-pCO20; 
real delpCO2_6c(t,x) = pCO2_6-pCO20; 
real delpCO2_6p(t,x) = pCO2_6-pCO20; 
real delpCO2_6m(t,x) = pCO2_6-pCO20; 
real delpCO2_7c(t,x) = pCO2_7-pCO20; 
real delpCO2_7p(t,x) = pCO2_7-pCO20; 
real delpCO2_7m(t,x) = pCO2_7-pCO20; 
real delpCO2_8c(t,x) = pCO2_8-pCO20; 
real delpCO2_8p(t,x) = pCO2_8-pCO20; 
real delpCO2_8m(t,x) = pCO2_8-pCO20; 
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
real O2_1c(t,x) M;             // free O2_1c
     O2_1c = alphaO2Sys*pO2_1c; 
real CO2_1c(t,x) M;            // free CO2_1c
     CO2_1c = alphaCO2Sys*pCO2_1;  
real O2_1p(t,x) M;             // free O2_1p
     O2_1p = alphaO2Sys*pO2_1p; 
real CO2_1p(t,x) M;            // free CO2_1p
     CO2_1p = alphaCO2Sys*pCO2_1;  
real O2_1m(t,x) M;             // free O2_1m
     O2_1m = alphaO2Sys*pO2_1m; 
real CO2_1m(t,x) M;            // free CO2_1m
     CO2_1m = alphaCO2Sys*pCO2_1;  
real O2_2c(t,x) M;             // free O2_2c
     O2_2c = alphaO2Sys*pO2_2c; 
real CO2_2c(t,x) M;            // free CO2_2c
     CO2_2c = alphaCO2Sys*pCO2_2;  
real O2_2p(t,x) M;             // free O2_2p
     O2_2p = alphaO2Sys*pO2_2p; 
real CO2_2p(t,x) M;            // free CO2_2p
     CO2_2p = alphaCO2Sys*pCO2_2;  
real O2_2m(t,x) M;             // free O2_2m
     O2_2m = alphaO2Sys*pO2_2m; 
real CO2_2m(t,x) M;            // free CO2_2m
     CO2_2m = alphaCO2Sys*pCO2_2;  
real O2_3c(t,x) M;             // free O2_3c
     O2_3c = alphaO2Sys*pO2_3c; 
real CO2_3c(t,x) M;            // free CO2_3c
     CO2_3c = alphaCO2Sys*pCO2_3;  
real O2_3p(t,x) M;             // free O2_3p
     O2_3p = alphaO2Sys*pO2_3p; 
real CO2_3p(t,x) M;            // free CO2_3p
     CO2_3p = alphaCO2Sys*pCO2_3;  
real O2_3m(t,x) M;             // free O2_3m
     O2_3m = alphaO2Sys*pO2_3m; 
real CO2_3m(t,x) M;            // free CO2_3m
     CO2_3m = alphaCO2Sys*pCO2_3;  
real O2_4c(t,x) M;             // free O2_4c
     O2_4c = alphaO2Sys*pO2_4c; 
real CO2_4c(t,x) M;            // free CO2_4c
     CO2_4c = alphaCO2Sys*pCO2_4;  
real O2_4p(t,x) M;             // free O2_4p
     O2_4p = alphaO2Sys*pO2_4p; 
real CO2_4p(t,x) M;            // free CO2_4p
     CO2_4p = alphaCO2Sys*pCO2_4;  
real O2_4m(t,x) M;             // free O2_4m
     O2_4m = alphaO2Sys*pO2_4m; 
real CO2_4m(t,x) M;            // free CO2_4m
     CO2_4m = alphaCO2Sys*pCO2_4;  
real O2_5c(t,x) M;             // free O2_5c
     O2_5c = alphaO2Sys*pO2_5c; 
real CO2_5c(t,x) M;            // free CO2_5c
     CO2_5c = alphaCO2Sys*pCO2_5;  
real O2_5p(t,x) M;             // free O2_5p
     O2_5p = alphaO2Sys*pO2_5p; 
real CO2_5p(t,x) M;            // free CO2_5p
     CO2_5p = alphaCO2Sys*pCO2_5;  
real O2_5m(t,x) M;             // free O2_5m
     O2_5m = alphaO2Sys*pO2_5m; 
real CO2_5m(t,x) M;            // free CO2_5m
     CO2_5m = alphaCO2Sys*pCO2_5;  
real O2_6c(t,x) M;             // free O2_6c
     O2_6c = alphaO2Sys*pO2_6c; 
real CO2_6c(t,x) M;            // free CO2_6c
     CO2_6c = alphaCO2Sys*pCO2_6;  
real O2_6p(t,x) M;             // free O2_6p
     O2_6p = alphaO2Sys*pO2_6p; 
real CO2_6p(t,x) M;            // free CO2_6p
     CO2_6p = alphaCO2Sys*pCO2_6;  
real O2_6m(t,x) M;             // free O2_6m
     O2_6m = alphaO2Sys*pO2_6m; 
real CO2_6m(t,x) M;            // free CO2_6m
     CO2_6m = alphaCO2Sys*pCO2_6;  
real O2_7c(t,x) M;             // free O2_7c
     O2_7c = alphaO2Sys*pO2_7c; 
real CO2_7c(t,x) M;            // free CO2_7c
     CO2_7c = alphaCO2Sys*pCO2_7;  
real O2_7p(t,x) M;             // free O2_7p
     O2_7p = alphaO2Sys*pO2_7p; 
real CO2_7p(t,x) M;            // free CO2_7p
     CO2_7p = alphaCO2Sys*pCO2_7;  
real O2_7m(t,x) M;             // free O2_7m
     O2_7m = alphaO2Sys*pO2_7m; 
real CO2_7m(t,x) M;            // free CO2_7m
     CO2_7m = alphaCO2Sys*pCO2_7;  
real O2_8c(t,x) M;             // free O2_8c
     O2_8c = alphaO2Sys*pO2_8c; 
real CO2_8c(t,x) M;            // free CO2_8c
     CO2_8c = alphaCO2Sys*pCO2_8;  
real O2_8p(t,x) M;             // free O2_8p
     O2_8p = alphaO2Sys*pO2_8p; 
real CO2_8p(t,x) M;            // free CO2_8p
     CO2_8p = alphaCO2Sys*pCO2_8;  
real O2_8m(t,x) M;             // free O2_8m
     O2_8m = alphaO2Sys*pO2_8m; 
real CO2_8m(t,x) M;            // free CO2_8m
     CO2_8m = alphaCO2Sys*pCO2_8;  
// p50 adjustments:
real P501_1c(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_1c(t,x) = P500 + 1.7*(4.28e-2*delpCO2_1c + 3.64e-5*(1 mmHg^-1)*delpCO2_1c^2); // all standard conditions, except CO2
real P503_1c(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_1c(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_1c(t,x) mmHg;             // Final P50_1c value for 50% O2 binding to Hb
     P50_1c = P500*(P501_1c/P500)*(P502_1c/P500)*(P503_1c/P500)*(P504_1c/P500); 
real C50_1c(t,x) M;
     C50_1c = alphaO2Sys*P50_1c; 
real P501_1p(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_1p(t,x) = P500 + 1.7*(4.28e-2*delpCO2_1p + 3.64e-5*(1 mmHg^-1)*delpCO2_1p^2); // all standard conditions, except CO2
real P503_1p(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_1p(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_1p(t,x) mmHg;             // Final P50_1p value for 50% O2 binding to Hb
     P50_1p = P500*(P501_1p/P500)*(P502_1p/P500)*(P503_1p/P500)*(P504_1p/P500); 
real C50_1p(t,x) M;
     C50_1p = alphaO2Sys*P50_1p; 
real P501_1m(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_1m(t,x) = P500 + 1.7*(4.28e-2*delpCO2_1m + 3.64e-5*(1 mmHg^-1)*delpCO2_1m^2); // all standard conditions, except CO2
real P503_1m(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_1m(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_1m(t,x) mmHg;             // Final P50_1m value for 50% O2 binding to Hb
     P50_1m = P500*(P501_1m/P500)*(P502_1m/P500)*(P503_1m/P500)*(P504_1m/P500); 
real C50_1m(t,x) M;
     C50_1m = alphaO2Sys*P50_1m; 
real P501_2c(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_2c(t,x) = P500 + 1.7*(4.28e-2*delpCO2_2c + 3.64e-5*(1 mmHg^-1)*delpCO2_2c^2); // all standard conditions, except CO2
real P503_2c(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_2c(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_2c(t,x) mmHg;             // Final P50_2c value for 50% O2 binding to Hb
     P50_2c = P500*(P501_2c/P500)*(P502_2c/P500)*(P503_2c/P500)*(P504_2c/P500); 
real C50_2c(t,x) M;
     C50_2c = alphaO2Sys*P50_2c; 
real P501_2p(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_2p(t,x) = P500 + 1.7*(4.28e-2*delpCO2_2p + 3.64e-5*(1 mmHg^-1)*delpCO2_2p^2); // all standard conditions, except CO2
real P503_2p(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_2p(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_2p(t,x) mmHg;             // Final P50_2p value for 50% O2 binding to Hb
     P50_2p = P500*(P501_2p/P500)*(P502_2p/P500)*(P503_2p/P500)*(P504_2p/P500); 
real C50_2p(t,x) M;
     C50_2p = alphaO2Sys*P50_2p; 
real P501_2m(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_2m(t,x) = P500 + 1.7*(4.28e-2*delpCO2_2m + 3.64e-5*(1 mmHg^-1)*delpCO2_2m^2); // all standard conditions, except CO2
real P503_2m(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_2m(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_2m(t,x) mmHg;             // Final P50_2m value for 50% O2 binding to Hb
     P50_2m = P500*(P501_2m/P500)*(P502_2m/P500)*(P503_2m/P500)*(P504_2m/P500); 
real C50_2m(t,x) M;
     C50_2m = alphaO2Sys*P50_2m; 
real P501_3c(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_3c(t,x) = P500 + 1.7*(4.28e-2*delpCO2_3c + 3.64e-5*(1 mmHg^-1)*delpCO2_3c^2); // all standard conditions, except CO2
real P503_3c(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_3c(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_3c(t,x) mmHg;             // Final P50_3c value for 50% O2 binding to Hb
     P50_3c = P500*(P501_3c/P500)*(P502_3c/P500)*(P503_3c/P500)*(P504_3c/P500); 
real C50_3c(t,x) M;
     C50_3c = alphaO2Sys*P50_3c; 
real P501_3p(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_3p(t,x) = P500 + 1.7*(4.28e-2*delpCO2_3p + 3.64e-5*(1 mmHg^-1)*delpCO2_3p^2); // all standard conditions, except CO2
real P503_3p(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_3p(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_3p(t,x) mmHg;             // Final P50_3p value for 50% O2 binding to Hb
     P50_3p = P500*(P501_3p/P500)*(P502_3p/P500)*(P503_3p/P500)*(P504_3p/P500); 
real C50_3p(t,x) M;
     C50_3p = alphaO2Sys*P50_3p; 
real P501_3m(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_3m(t,x) = P500 + 1.7*(4.28e-2*delpCO2_3m + 3.64e-5*(1 mmHg^-1)*delpCO2_3m^2); // all standard conditions, except CO2
real P503_3m(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_3m(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_3m(t,x) mmHg;             // Final P50_3m value for 50% O2 binding to Hb
     P50_3m = P500*(P501_3m/P500)*(P502_3m/P500)*(P503_3m/P500)*(P504_3m/P500); 
real C50_3m(t,x) M;
     C50_3m = alphaO2Sys*P50_3m; 
real P501_4c(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_4c(t,x) = P500 + 1.7*(4.28e-2*delpCO2_4c + 3.64e-5*(1 mmHg^-1)*delpCO2_4c^2); // all standard conditions, except CO2
real P503_4c(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_4c(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_4c(t,x) mmHg;             // Final P50_4c value for 50% O2 binding to Hb
     P50_4c = P500*(P501_4c/P500)*(P502_4c/P500)*(P503_4c/P500)*(P504_4c/P500); 
real C50_4c(t,x) M;
     C50_4c = alphaO2Sys*P50_4c; 
real P501_4p(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_4p(t,x) = P500 + 1.7*(4.28e-2*delpCO2_4p + 3.64e-5*(1 mmHg^-1)*delpCO2_4p^2); // all standard conditions, except CO2
real P503_4p(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_4p(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_4p(t,x) mmHg;             // Final P50_4p value for 50% O2 binding to Hb
     P50_4p = P500*(P501_4p/P500)*(P502_4p/P500)*(P503_4p/P500)*(P504_4p/P500); 
real C50_4p(t,x) M;
     C50_4p = alphaO2Sys*P50_4p; 
real P501_4m(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_4m(t,x) = P500 + 1.7*(4.28e-2*delpCO2_4m + 3.64e-5*(1 mmHg^-1)*delpCO2_4m^2); // all standard conditions, except CO2
real P503_4m(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_4m(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_4m(t,x) mmHg;             // Final P50_4m value for 50% O2 binding to Hb
     P50_4m = P500*(P501_4m/P500)*(P502_4m/P500)*(P503_4m/P500)*(P504_4m/P500); 
real C50_4m(t,x) M;
     C50_4m = alphaO2Sys*P50_4m; 
real P501_5c(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_5c(t,x) = P500 + 1.7*(4.28e-2*delpCO2_5c + 3.64e-5*(1 mmHg^-1)*delpCO2_5c^2); // all standard conditions, except CO2
real P503_5c(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_5c(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_5c(t,x) mmHg;             // Final P50_5c value for 50% O2 binding to Hb
     P50_5c = P500*(P501_5c/P500)*(P502_5c/P500)*(P503_5c/P500)*(P504_5c/P500); 
real C50_5c(t,x) M;
     C50_5c = alphaO2Sys*P50_5c; 
real P501_5p(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_5p(t,x) = P500 + 1.7*(4.28e-2*delpCO2_5p + 3.64e-5*(1 mmHg^-1)*delpCO2_5p^2); // all standard conditions, except CO2
real P503_5p(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_5p(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_5p(t,x) mmHg;             // Final P50_5p value for 50% O2 binding to Hb
     P50_5p = P500*(P501_5p/P500)*(P502_5p/P500)*(P503_5p/P500)*(P504_5p/P500); 
real C50_5p(t,x) M;
     C50_5p = alphaO2Sys*P50_5p; 
real P501_5m(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_5m(t,x) = P500 + 1.7*(4.28e-2*delpCO2_5m + 3.64e-5*(1 mmHg^-1)*delpCO2_5m^2); // all standard conditions, except CO2
real P503_5m(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_5m(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_5m(t,x) mmHg;             // Final P50_5m value for 50% O2 binding to Hb
     P50_5m = P500*(P501_5m/P500)*(P502_5m/P500)*(P503_5m/P500)*(P504_5m/P500); 
real C50_5m(t,x) M;
     C50_5m = alphaO2Sys*P50_5m; 
real P501_6c(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_6c(t,x) = P500 + 1.7*(4.28e-2*delpCO2_6c + 3.64e-5*(1 mmHg^-1)*delpCO2_6c^2); // all standard conditions, except CO2
real P503_6c(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_6c(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_6c(t,x) mmHg;             // Final P50_6c value for 50% O2 binding to Hb
     P50_6c = P500*(P501_6c/P500)*(P502_6c/P500)*(P503_6c/P500)*(P504_6c/P500); 
real C50_6c(t,x) M;
     C50_6c = alphaO2Sys*P50_6c; 
real P501_6p(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_6p(t,x) = P500 + 1.7*(4.28e-2*delpCO2_6p + 3.64e-5*(1 mmHg^-1)*delpCO2_6p^2); // all standard conditions, except CO2
real P503_6p(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_6p(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_6p(t,x) mmHg;             // Final P50_6p value for 50% O2 binding to Hb
     P50_6p = P500*(P501_6p/P500)*(P502_6p/P500)*(P503_6p/P500)*(P504_6p/P500); 
real C50_6p(t,x) M;
     C50_6p = alphaO2Sys*P50_6p; 
real P501_6m(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_6m(t,x) = P500 + 1.7*(4.28e-2*delpCO2_6m + 3.64e-5*(1 mmHg^-1)*delpCO2_6m^2); // all standard conditions, except CO2
real P503_6m(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_6m(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_6m(t,x) mmHg;             // Final P50_6m value for 50% O2 binding to Hb
     P50_6m = P500*(P501_6m/P500)*(P502_6m/P500)*(P503_6m/P500)*(P504_6m/P500); 
real C50_6m(t,x) M;
     C50_6m = alphaO2Sys*P50_6m; 
real P501_7c(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_7c(t,x) = P500 + 1.7*(4.28e-2*delpCO2_7c + 3.64e-5*(1 mmHg^-1)*delpCO2_7c^2); // all standard conditions, except CO2
real P503_7c(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_7c(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_7c(t,x) mmHg;             // Final P50_7c value for 50% O2 binding to Hb
     P50_7c = P500*(P501_7c/P500)*(P502_7c/P500)*(P503_7c/P500)*(P504_7c/P500); 
real C50_7c(t,x) M;
     C50_7c = alphaO2Sys*P50_7c; 
real P501_7p(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_7p(t,x) = P500 + 1.7*(4.28e-2*delpCO2_7p + 3.64e-5*(1 mmHg^-1)*delpCO2_7p^2); // all standard conditions, except CO2
real P503_7p(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_7p(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_7p(t,x) mmHg;             // Final P50_7p value for 50% O2 binding to Hb
     P50_7p = P500*(P501_7p/P500)*(P502_7p/P500)*(P503_7p/P500)*(P504_7p/P500); 
real C50_7p(t,x) M;
     C50_7p = alphaO2Sys*P50_7p; 
real P501_7m(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_7m(t,x) = P500 + 1.7*(4.28e-2*delpCO2_7m + 3.64e-5*(1 mmHg^-1)*delpCO2_7m^2); // all standard conditions, except CO2
real P503_7m(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_7m(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_7m(t,x) mmHg;             // Final P50_7m value for 50% O2 binding to Hb
     P50_7m = P500*(P501_7m/P500)*(P502_7m/P500)*(P503_7m/P500)*(P504_7m/P500); 
real C50_7m(t,x) M;
     C50_7m = alphaO2Sys*P50_7m; 
real P501_8c(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_8c(t,x) = P500 + 1.7*(4.28e-2*delpCO2_8c + 3.64e-5*(1 mmHg^-1)*delpCO2_8c^2); // all standard conditions, except CO2
real P503_8c(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_8c(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_8c(t,x) mmHg;             // Final P50_8c value for 50% O2 binding to Hb
     P50_8c = P500*(P501_8c/P500)*(P502_8c/P500)*(P503_8c/P500)*(P504_8c/P500); 
real C50_8c(t,x) M;
     C50_8c = alphaO2Sys*P50_8c; 
real P501_8p(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_8p(t,x) = P500 + 1.7*(4.28e-2*delpCO2_8p + 3.64e-5*(1 mmHg^-1)*delpCO2_8p^2); // all standard conditions, except CO2
real P503_8p(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_8p(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_8p(t,x) mmHg;             // Final P50_8p value for 50% O2 binding to Hb
     P50_8p = P500*(P501_8p/P500)*(P502_8p/P500)*(P503_8p/P500)*(P504_8p/P500); 
real C50_8p(t,x) M;
     C50_8p = alphaO2Sys*P50_8p; 
real P501_8m(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502_8m(t,x) = P500 + 1.7*(4.28e-2*delpCO2_8m + 3.64e-5*(1 mmHg^-1)*delpCO2_8m^2); // all standard conditions, except CO2
real P503_8m(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504_8m(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50_8m(t,x) mmHg;             // Final P50_8m value for 50% O2 binding to Hb
     P50_8m = P500*(P501_8m/P500)*(P502_8m/P500)*(P503_8m/P500)*(P504_8m/P500); 
real C50_8m(t,x) M;
     C50_8m = alphaO2Sys*P50_8m; 
// equilibirum consts are not directly tied to CO2 and O2 conc:
real BPH1(t,x) dimensionless;
     BPH1 = 1+K2dp/Hrbc;         // Binding polynomial involving K2dp and Hrbc 
real BPH2(t,x) = 1+K3dp/Hrbc;      // Binding polynomial involving K3dp and Hrbc 
real BPH3(t,x) = 1+Hrbc/K5dp;      // Binding polynomial involving K5dp and Hrbc 
real BPH4(t,x) = 1+Hrbc/K6dp;      // Binding polynomial involving K6dp and Hrbc 
// Hill coefficient; unitless (redefined as a function pO2_1c)
real alpha = 2.8 dimensionless; real beta = 1.2 dimensionless; real gamma = 29.2 mmHg; // Roughton et al data
real nH_1c(t,x) dimensionless;
     nH_1c = alpha-beta*10^(-pO2_1c/gamma); // pO2_1c dependent variable nH_1c
// Hill coefficient; unitless (redefined as a function pO2_1p)
real nH_1p(t,x) dimensionless;
     nH_1p = alpha-beta*10^(-pO2_1p/gamma); // pO2_1p dependent variable nH_1p
// Hill coefficient; unitless (redefined as a function pO2_1m)
real nH_1m(t,x) dimensionless;
     nH_1m = alpha-beta*10^(-pO2_1m/gamma); // pO2_1m dependent variable nH_1m
// Hill coefficient; unitless (redefined as a function pO2_2c)
real nH_2c(t,x) dimensionless;
     nH_2c = alpha-beta*10^(-pO2_2c/gamma); // pO2_2c dependent variable nH_2c
// Hill coefficient; unitless (redefined as a function pO2_2p)
real nH_2p(t,x) dimensionless;
     nH_2p = alpha-beta*10^(-pO2_2p/gamma); // pO2_2p dependent variable nH_2p
// Hill coefficient; unitless (redefined as a function pO2_2m)
real nH_2m(t,x) dimensionless;
     nH_2m = alpha-beta*10^(-pO2_2m/gamma); // pO2_2m dependent variable nH_2m
// Hill coefficient; unitless (redefined as a function pO2_3c)
real nH_3c(t,x) dimensionless;
     nH_3c = alpha-beta*10^(-pO2_3c/gamma); // pO2_3c dependent variable nH_3c
// Hill coefficient; unitless (redefined as a function pO2_3p)
real nH_3p(t,x) dimensionless;
     nH_3p = alpha-beta*10^(-pO2_3p/gamma); // pO2_3p dependent variable nH_3p
// Hill coefficient; unitless (redefined as a function pO2_3m)
real nH_3m(t,x) dimensionless;
     nH_3m = alpha-beta*10^(-pO2_3m/gamma); // pO2_3m dependent variable nH_3m
// Hill coefficient; unitless (redefined as a function pO2_4c)
real nH_4c(t,x) dimensionless;
     nH_4c = alpha-beta*10^(-pO2_4c/gamma); // pO2_4c dependent variable nH_4c
// Hill coefficient; unitless (redefined as a function pO2_4p)
real nH_4p(t,x) dimensionless;
     nH_4p = alpha-beta*10^(-pO2_4p/gamma); // pO2_4p dependent variable nH_4p
// Hill coefficient; unitless (redefined as a function pO2_4m)
real nH_4m(t,x) dimensionless;
     nH_4m = alpha-beta*10^(-pO2_4m/gamma); // pO2_4m dependent variable nH_4m
// Hill coefficient; unitless (redefined as a function pO2_5c)
real nH_5c(t,x) dimensionless;
     nH_5c = alpha-beta*10^(-pO2_5c/gamma); // pO2_5c dependent variable nH_5c
// Hill coefficient; unitless (redefined as a function pO2_5p)
real nH_5p(t,x) dimensionless;
     nH_5p = alpha-beta*10^(-pO2_5p/gamma); // pO2_5p dependent variable nH_5p
// Hill coefficient; unitless (redefined as a function pO2_5m)
real nH_5m(t,x) dimensionless;
     nH_5m = alpha-beta*10^(-pO2_5m/gamma); // pO2_5m dependent variable nH_5m
// Hill coefficient; unitless (redefined as a function pO2_6c)
real nH_6c(t,x) dimensionless;
     nH_6c = alpha-beta*10^(-pO2_6c/gamma); // pO2_6c dependent variable nH_6c
// Hill coefficient; unitless (redefined as a function pO2_6p)
real nH_6p(t,x) dimensionless;
     nH_6p = alpha-beta*10^(-pO2_6p/gamma); // pO2_6p dependent variable nH_6p
// Hill coefficient; unitless (redefined as a function pO2_6m)
real nH_6m(t,x) dimensionless;
     nH_6m = alpha-beta*10^(-pO2_6m/gamma); // pO2_6m dependent variable nH_6m
// Hill coefficient; unitless (redefined as a function pO2_7c)
real nH_7c(t,x) dimensionless;
     nH_7c = alpha-beta*10^(-pO2_7c/gamma); // pO2_7c dependent variable nH_7c
// Hill coefficient; unitless (redefined as a function pO2_7p)
real nH_7p(t,x) dimensionless;
     nH_7p = alpha-beta*10^(-pO2_7p/gamma); // pO2_7p dependent variable nH_7p
// Hill coefficient; unitless (redefined as a function pO2_7m)
real nH_7m(t,x) dimensionless;
     nH_7m = alpha-beta*10^(-pO2_7m/gamma); // pO2_7m dependent variable nH_7m
// Hill coefficient; unitless (redefined as a function pO2_8c)
real nH_8c(t,x) dimensionless;
     nH_8c = alpha-beta*10^(-pO2_8c/gamma); // pO2_8c dependent variable nH_8c
// Hill coefficient; unitless (redefined as a function pO2_8p)
real nH_8p(t,x) dimensionless;
     nH_8p = alpha-beta*10^(-pO2_8p/gamma); // pO2_8p dependent variable nH_8p
// Hill coefficient; unitless (redefined as a function pO2_8m)
real nH_8m(t,x) dimensionless;
     nH_8m = alpha-beta*10^(-pO2_8m/gamma); // pO2_8m dependent variable nH_8m
// Compute the apparent equilibrium constant of Hb with O2_1c and CO2_1c (KHbO2 and KHbCO2); 
// O2_1c and CO2_1c saturations of Hb (SHbO2 and SHbCO2); and O2_1c and CO2_1c contents in blood. 
real K4p_1c(t,x) 1/M;  
     K4p_1c = (1 M^-1)*((O2_1c*(1 M^-1))^(nH_1c-1)*(K2p*BPH1*CO2_1c+BPH3))/( ((1 M^-1)*C50_1c)^nH_1c*(K3p*BPH2*CO2_1c+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_1p and CO2_1p (KHbO2 and KHbCO2); 
// O2_1p and CO2_1p saturations of Hb (SHbO2 and SHbCO2); and O2_1p and CO2_1p contents in blood. 
real K4p_1p(t,x) 1/M;  
     K4p_1p = (1 M^-1)*((O2_1p*(1 M^-1))^(nH_1p-1)*(K2p*BPH1*CO2_1p+BPH3))/( ((1 M^-1)*C50_1p)^nH_1p*(K3p*BPH2*CO2_1p+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_1m and CO2_1m (KHbO2 and KHbCO2); 
// O2_1m and CO2_1m saturations of Hb (SHbO2 and SHbCO2); and O2_1m and CO2_1m contents in blood. 
real K4p_1m(t,x) 1/M;  
     K4p_1m = (1 M^-1)*((O2_1m*(1 M^-1))^(nH_1m-1)*(K2p*BPH1*CO2_1m+BPH3))/( ((1 M^-1)*C50_1m)^nH_1m*(K3p*BPH2*CO2_1m+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_2c and CO2_2c (KHbO2 and KHbCO2); 
// O2_2c and CO2_2c saturations of Hb (SHbO2 and SHbCO2); and O2_2c and CO2_2c contents in blood. 
real K4p_2c(t,x) 1/M;  
     K4p_2c = (1 M^-1)*((O2_2c*(1 M^-1))^(nH_2c-1)*(K2p*BPH1*CO2_2c+BPH3))/( ((1 M^-1)*C50_2c)^nH_2c*(K3p*BPH2*CO2_2c+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_2p and CO2_2p (KHbO2 and KHbCO2); 
// O2_2p and CO2_2p saturations of Hb (SHbO2 and SHbCO2); and O2_2p and CO2_2p contents in blood. 
real K4p_2p(t,x) 1/M;  
     K4p_2p = (1 M^-1)*((O2_2p*(1 M^-1))^(nH_2p-1)*(K2p*BPH1*CO2_2p+BPH3))/( ((1 M^-1)*C50_2p)^nH_2p*(K3p*BPH2*CO2_2p+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_2m and CO2_2m (KHbO2 and KHbCO2); 
// O2_2m and CO2_2m saturations of Hb (SHbO2 and SHbCO2); and O2_2m and CO2_2m contents in blood. 
real K4p_2m(t,x) 1/M;  
     K4p_2m = (1 M^-1)*((O2_2m*(1 M^-1))^(nH_2m-1)*(K2p*BPH1*CO2_2m+BPH3))/( ((1 M^-1)*C50_2m)^nH_2m*(K3p*BPH2*CO2_2m+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_3c and CO2_3c (KHbO2 and KHbCO2); 
// O2_3c and CO2_3c saturations of Hb (SHbO2 and SHbCO2); and O2_3c and CO2_3c contents in blood. 
real K4p_3c(t,x) 1/M;  
     K4p_3c = (1 M^-1)*((O2_3c*(1 M^-1))^(nH_3c-1)*(K2p*BPH1*CO2_3c+BPH3))/( ((1 M^-1)*C50_3c)^nH_3c*(K3p*BPH2*CO2_3c+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_3p and CO2_3p (KHbO2 and KHbCO2); 
// O2_3p and CO2_3p saturations of Hb (SHbO2 and SHbCO2); and O2_3p and CO2_3p contents in blood. 
real K4p_3p(t,x) 1/M;  
     K4p_3p = (1 M^-1)*((O2_3p*(1 M^-1))^(nH_3p-1)*(K2p*BPH1*CO2_3p+BPH3))/( ((1 M^-1)*C50_3p)^nH_3p*(K3p*BPH2*CO2_3p+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_3m and CO2_3m (KHbO2 and KHbCO2); 
// O2_3m and CO2_3m saturations of Hb (SHbO2 and SHbCO2); and O2_3m and CO2_3m contents in blood. 
real K4p_3m(t,x) 1/M;  
     K4p_3m = (1 M^-1)*((O2_3m*(1 M^-1))^(nH_3m-1)*(K2p*BPH1*CO2_3m+BPH3))/( ((1 M^-1)*C50_3m)^nH_3m*(K3p*BPH2*CO2_3m+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_4c and CO2_4c (KHbO2 and KHbCO2); 
// O2_4c and CO2_4c saturations of Hb (SHbO2 and SHbCO2); and O2_4c and CO2_4c contents in blood. 
real K4p_4c(t,x) 1/M;  
     K4p_4c = (1 M^-1)*((O2_4c*(1 M^-1))^(nH_4c-1)*(K2p*BPH1*CO2_4c+BPH3))/( ((1 M^-1)*C50_4c)^nH_4c*(K3p*BPH2*CO2_4c+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_4p and CO2_4p (KHbO2 and KHbCO2); 
// O2_4p and CO2_4p saturations of Hb (SHbO2 and SHbCO2); and O2_4p and CO2_4p contents in blood. 
real K4p_4p(t,x) 1/M;  
     K4p_4p = (1 M^-1)*((O2_4p*(1 M^-1))^(nH_4p-1)*(K2p*BPH1*CO2_4p+BPH3))/( ((1 M^-1)*C50_4p)^nH_4p*(K3p*BPH2*CO2_4p+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_4m and CO2_4m (KHbO2 and KHbCO2); 
// O2_4m and CO2_4m saturations of Hb (SHbO2 and SHbCO2); and O2_4m and CO2_4m contents in blood. 
real K4p_4m(t,x) 1/M;  
     K4p_4m = (1 M^-1)*((O2_4m*(1 M^-1))^(nH_4m-1)*(K2p*BPH1*CO2_4m+BPH3))/( ((1 M^-1)*C50_4m)^nH_4m*(K3p*BPH2*CO2_4m+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_5c and CO2_5c (KHbO2 and KHbCO2); 
// O2_5c and CO2_5c saturations of Hb (SHbO2 and SHbCO2); and O2_5c and CO2_5c contents in blood. 
real K4p_5c(t,x) 1/M;  
     K4p_5c = (1 M^-1)*((O2_5c*(1 M^-1))^(nH_5c-1)*(K2p*BPH1*CO2_5c+BPH3))/( ((1 M^-1)*C50_5c)^nH_5c*(K3p*BPH2*CO2_5c+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_5p and CO2_5p (KHbO2 and KHbCO2); 
// O2_5p and CO2_5p saturations of Hb (SHbO2 and SHbCO2); and O2_5p and CO2_5p contents in blood. 
real K4p_5p(t,x) 1/M;  
     K4p_5p = (1 M^-1)*((O2_5p*(1 M^-1))^(nH_5p-1)*(K2p*BPH1*CO2_5p+BPH3))/( ((1 M^-1)*C50_5p)^nH_5p*(K3p*BPH2*CO2_5p+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_5m and CO2_5m (KHbO2 and KHbCO2); 
// O2_5m and CO2_5m saturations of Hb (SHbO2 and SHbCO2); and O2_5m and CO2_5m contents in blood. 
real K4p_5m(t,x) 1/M;  
     K4p_5m = (1 M^-1)*((O2_5m*(1 M^-1))^(nH_5m-1)*(K2p*BPH1*CO2_5m+BPH3))/( ((1 M^-1)*C50_5m)^nH_5m*(K3p*BPH2*CO2_5m+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_6c and CO2_6c (KHbO2 and KHbCO2); 
// O2_6c and CO2_6c saturations of Hb (SHbO2 and SHbCO2); and O2_6c and CO2_6c contents in blood. 
real K4p_6c(t,x) 1/M;  
     K4p_6c = (1 M^-1)*((O2_6c*(1 M^-1))^(nH_6c-1)*(K2p*BPH1*CO2_6c+BPH3))/( ((1 M^-1)*C50_6c)^nH_6c*(K3p*BPH2*CO2_6c+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_6p and CO2_6p (KHbO2 and KHbCO2); 
// O2_6p and CO2_6p saturations of Hb (SHbO2 and SHbCO2); and O2_6p and CO2_6p contents in blood. 
real K4p_6p(t,x) 1/M;  
     K4p_6p = (1 M^-1)*((O2_6p*(1 M^-1))^(nH_6p-1)*(K2p*BPH1*CO2_6p+BPH3))/( ((1 M^-1)*C50_6p)^nH_6p*(K3p*BPH2*CO2_6p+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_6m and CO2_6m (KHbO2 and KHbCO2); 
// O2_6m and CO2_6m saturations of Hb (SHbO2 and SHbCO2); and O2_6m and CO2_6m contents in blood. 
real K4p_6m(t,x) 1/M;  
     K4p_6m = (1 M^-1)*((O2_6m*(1 M^-1))^(nH_6m-1)*(K2p*BPH1*CO2_6m+BPH3))/( ((1 M^-1)*C50_6m)^nH_6m*(K3p*BPH2*CO2_6m+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_7c and CO2_7c (KHbO2 and KHbCO2); 
// O2_7c and CO2_7c saturations of Hb (SHbO2 and SHbCO2); and O2_7c and CO2_7c contents in blood. 
real K4p_7c(t,x) 1/M;  
     K4p_7c = (1 M^-1)*((O2_7c*(1 M^-1))^(nH_7c-1)*(K2p*BPH1*CO2_7c+BPH3))/( ((1 M^-1)*C50_7c)^nH_7c*(K3p*BPH2*CO2_7c+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_7p and CO2_7p (KHbO2 and KHbCO2); 
// O2_7p and CO2_7p saturations of Hb (SHbO2 and SHbCO2); and O2_7p and CO2_7p contents in blood. 
real K4p_7p(t,x) 1/M;  
     K4p_7p = (1 M^-1)*((O2_7p*(1 M^-1))^(nH_7p-1)*(K2p*BPH1*CO2_7p+BPH3))/( ((1 M^-1)*C50_7p)^nH_7p*(K3p*BPH2*CO2_7p+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_7m and CO2_7m (KHbO2 and KHbCO2); 
// O2_7m and CO2_7m saturations of Hb (SHbO2 and SHbCO2); and O2_7m and CO2_7m contents in blood. 
real K4p_7m(t,x) 1/M;  
     K4p_7m = (1 M^-1)*((O2_7m*(1 M^-1))^(nH_7m-1)*(K2p*BPH1*CO2_7m+BPH3))/( ((1 M^-1)*C50_7m)^nH_7m*(K3p*BPH2*CO2_7m+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_8c and CO2_8c (KHbO2 and KHbCO2); 
// O2_8c and CO2_8c saturations of Hb (SHbO2 and SHbCO2); and O2_8c and CO2_8c contents in blood. 
real K4p_8c(t,x) 1/M;  
     K4p_8c = (1 M^-1)*((O2_8c*(1 M^-1))^(nH_8c-1)*(K2p*BPH1*CO2_8c+BPH3))/( ((1 M^-1)*C50_8c)^nH_8c*(K3p*BPH2*CO2_8c+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_8p and CO2_8p (KHbO2 and KHbCO2); 
// O2_8p and CO2_8p saturations of Hb (SHbO2 and SHbCO2); and O2_8p and CO2_8p contents in blood. 
real K4p_8p(t,x) 1/M;  
     K4p_8p = (1 M^-1)*((O2_8p*(1 M^-1))^(nH_8p-1)*(K2p*BPH1*CO2_8p+BPH3))/( ((1 M^-1)*C50_8p)^nH_8p*(K3p*BPH2*CO2_8p+BPH4));
// Compute the apparent equilibrium constant of Hb with O2_8m and CO2_8m (KHbO2 and KHbCO2); 
// O2_8m and CO2_8m saturations of Hb (SHbO2 and SHbCO2); and O2_8m and CO2_8m contents in blood. 
real K4p_8m(t,x) 1/M;  
     K4p_8m = (1 M^-1)*((O2_8m*(1 M^-1))^(nH_8m-1)*(K2p*BPH1*CO2_8m+BPH3))/( ((1 M^-1)*C50_8m)^nH_8m*(K3p*BPH2*CO2_8m+BPH4));
real KHbO2_1c(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_1c
     KHbO2_1c = K4p_1c*(K3p*BPH2*CO2_1c+BPH4)/(K2p*BPH1*CO2_1c+BPH3); 
real SHbO2_1c(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_1c = KHbO2_1c*O2_1c/(1+KHbO2_1c*O2_1c);
real KHbO2_1p(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_1p
     KHbO2_1p = K4p_1p*(K3p*BPH2*CO2_1p+BPH4)/(K2p*BPH1*CO2_1p+BPH3); 
real SHbO2_1p(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_1p = KHbO2_1p*O2_1p/(1+KHbO2_1p*O2_1p);
real KHbO2_1m(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_1m
     KHbO2_1m = K4p_1m*(K3p*BPH2*CO2_1m+BPH4)/(K2p*BPH1*CO2_1m+BPH3); 
real SHbO2_1m(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_1m = KHbO2_1m*O2_1m/(1+KHbO2_1m*O2_1m);
real KHbO2_2c(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_2c
     KHbO2_2c = K4p_2c*(K3p*BPH2*CO2_2c+BPH4)/(K2p*BPH1*CO2_2c+BPH3); 
real SHbO2_2c(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_2c = KHbO2_2c*O2_2c/(1+KHbO2_2c*O2_2c);
real KHbO2_2p(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_2p
     KHbO2_2p = K4p_2p*(K3p*BPH2*CO2_2p+BPH4)/(K2p*BPH1*CO2_2p+BPH3); 
real SHbO2_2p(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_2p = KHbO2_2p*O2_2p/(1+KHbO2_2p*O2_2p);
real KHbO2_2m(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_2m
     KHbO2_2m = K4p_2m*(K3p*BPH2*CO2_2m+BPH4)/(K2p*BPH1*CO2_2m+BPH3); 
real SHbO2_2m(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_2m = KHbO2_2m*O2_2m/(1+KHbO2_2m*O2_2m);
real KHbO2_3c(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_3c
     KHbO2_3c = K4p_3c*(K3p*BPH2*CO2_3c+BPH4)/(K2p*BPH1*CO2_3c+BPH3); 
real SHbO2_3c(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_3c = KHbO2_3c*O2_3c/(1+KHbO2_3c*O2_3c);
real KHbO2_3p(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_3p
     KHbO2_3p = K4p_3p*(K3p*BPH2*CO2_3p+BPH4)/(K2p*BPH1*CO2_3p+BPH3); 
real SHbO2_3p(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_3p = KHbO2_3p*O2_3p/(1+KHbO2_3p*O2_3p);
real KHbO2_3m(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_3m
     KHbO2_3m = K4p_3m*(K3p*BPH2*CO2_3m+BPH4)/(K2p*BPH1*CO2_3m+BPH3); 
real SHbO2_3m(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_3m = KHbO2_3m*O2_3m/(1+KHbO2_3m*O2_3m);
real KHbO2_4c(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_4c
     KHbO2_4c = K4p_4c*(K3p*BPH2*CO2_4c+BPH4)/(K2p*BPH1*CO2_4c+BPH3); 
real SHbO2_4c(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_4c = KHbO2_4c*O2_4c/(1+KHbO2_4c*O2_4c);
real KHbO2_4p(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_4p
     KHbO2_4p = K4p_4p*(K3p*BPH2*CO2_4p+BPH4)/(K2p*BPH1*CO2_4p+BPH3); 
real SHbO2_4p(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_4p = KHbO2_4p*O2_4p/(1+KHbO2_4p*O2_4p);
real KHbO2_4m(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_4m
     KHbO2_4m = K4p_4m*(K3p*BPH2*CO2_4m+BPH4)/(K2p*BPH1*CO2_4m+BPH3); 
real SHbO2_4m(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_4m = KHbO2_4m*O2_4m/(1+KHbO2_4m*O2_4m);
real KHbO2_5c(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_5c
     KHbO2_5c = K4p_5c*(K3p*BPH2*CO2_5c+BPH4)/(K2p*BPH1*CO2_5c+BPH3); 
real SHbO2_5c(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_5c = KHbO2_5c*O2_5c/(1+KHbO2_5c*O2_5c);
real KHbO2_5p(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_5p
     KHbO2_5p = K4p_5p*(K3p*BPH2*CO2_5p+BPH4)/(K2p*BPH1*CO2_5p+BPH3); 
real SHbO2_5p(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_5p = KHbO2_5p*O2_5p/(1+KHbO2_5p*O2_5p);
real KHbO2_5m(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_5m
     KHbO2_5m = K4p_5m*(K3p*BPH2*CO2_5m+BPH4)/(K2p*BPH1*CO2_5m+BPH3); 
real SHbO2_5m(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_5m = KHbO2_5m*O2_5m/(1+KHbO2_5m*O2_5m);
real KHbO2_6c(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_6c
     KHbO2_6c = K4p_6c*(K3p*BPH2*CO2_6c+BPH4)/(K2p*BPH1*CO2_6c+BPH3); 
real SHbO2_6c(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_6c = KHbO2_6c*O2_6c/(1+KHbO2_6c*O2_6c);
real KHbO2_6p(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_6p
     KHbO2_6p = K4p_6p*(K3p*BPH2*CO2_6p+BPH4)/(K2p*BPH1*CO2_6p+BPH3); 
real SHbO2_6p(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_6p = KHbO2_6p*O2_6p/(1+KHbO2_6p*O2_6p);
real KHbO2_6m(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_6m
     KHbO2_6m = K4p_6m*(K3p*BPH2*CO2_6m+BPH4)/(K2p*BPH1*CO2_6m+BPH3); 
real SHbO2_6m(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_6m = KHbO2_6m*O2_6m/(1+KHbO2_6m*O2_6m);
real KHbO2_7c(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_7c
     KHbO2_7c = K4p_7c*(K3p*BPH2*CO2_7c+BPH4)/(K2p*BPH1*CO2_7c+BPH3); 
real SHbO2_7c(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_7c = KHbO2_7c*O2_7c/(1+KHbO2_7c*O2_7c);
real KHbO2_7p(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_7p
     KHbO2_7p = K4p_7p*(K3p*BPH2*CO2_7p+BPH4)/(K2p*BPH1*CO2_7p+BPH3); 
real SHbO2_7p(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_7p = KHbO2_7p*O2_7p/(1+KHbO2_7p*O2_7p);
real KHbO2_7m(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_7m
     KHbO2_7m = K4p_7m*(K3p*BPH2*CO2_7m+BPH4)/(K2p*BPH1*CO2_7m+BPH3); 
real SHbO2_7m(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_7m = KHbO2_7m*O2_7m/(1+KHbO2_7m*O2_7m);
real KHbO2_8c(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_8c
     KHbO2_8c = K4p_8c*(K3p*BPH2*CO2_8c+BPH4)/(K2p*BPH1*CO2_8c+BPH3); 
real SHbO2_8c(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_8c = KHbO2_8c*O2_8c/(1+KHbO2_8c*O2_8c);
real KHbO2_8p(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_8p
     KHbO2_8p = K4p_8p*(K3p*BPH2*CO2_8p+BPH4)/(K2p*BPH1*CO2_8p+BPH3); 
real SHbO2_8p(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_8p = KHbO2_8p*O2_8p/(1+KHbO2_8p*O2_8p);
real KHbO2_8m(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2_8m
     KHbO2_8m = K4p_8m*(K3p*BPH2*CO2_8m+BPH4)/(K2p*BPH1*CO2_8m+BPH3); 
real SHbO2_8m(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2_8m = KHbO2_8m*O2_8m/(1+KHbO2_8m*O2_8m);
real O2free_1c(t,x) = Wrbc*O2_1c; // M (mol O2_1c per L rbc)
real O2bound_1c(t,x) = 4*Hbrbc*SHbO2_1c; // M (mol O2_1c per L rbc)
real O2tot_1c(t,x) = O2free_1c+O2bound_1c;   // M (mol O2_1c per L rbc)
real O2free_1p(t,x) = Wrbc*O2_1p; // M (mol O2_1p per L rbc)
real O2bound_1p(t,x) = 4*Hbrbc*SHbO2_1p; // M (mol O2_1p per L rbc)
real O2tot_1p(t,x) = O2free_1p+O2bound_1p;   // M (mol O2_1p per L rbc)
real O2free_1m(t,x) = Wrbc*O2_1m; // M (mol O2_1m per L rbc)
real O2bound_1m(t,x) = 4*Hbrbc*SHbO2_1m; // M (mol O2_1m per L rbc)
real O2tot_1m(t,x) = O2free_1m+O2bound_1m;   // M (mol O2_1m per L rbc)
real O2free_2c(t,x) = Wrbc*O2_2c; // M (mol O2_2c per L rbc)
real O2bound_2c(t,x) = 4*Hbrbc*SHbO2_2c; // M (mol O2_2c per L rbc)
real O2tot_2c(t,x) = O2free_2c+O2bound_2c;   // M (mol O2_2c per L rbc)
real O2free_2p(t,x) = Wrbc*O2_2p; // M (mol O2_2p per L rbc)
real O2bound_2p(t,x) = 4*Hbrbc*SHbO2_2p; // M (mol O2_2p per L rbc)
real O2tot_2p(t,x) = O2free_2p+O2bound_2p;   // M (mol O2_2p per L rbc)
real O2free_2m(t,x) = Wrbc*O2_2m; // M (mol O2_2m per L rbc)
real O2bound_2m(t,x) = 4*Hbrbc*SHbO2_2m; // M (mol O2_2m per L rbc)
real O2tot_2m(t,x) = O2free_2m+O2bound_2m;   // M (mol O2_2m per L rbc)
real O2free_3c(t,x) = Wrbc*O2_3c; // M (mol O2_3c per L rbc)
real O2bound_3c(t,x) = 4*Hbrbc*SHbO2_3c; // M (mol O2_3c per L rbc)
real O2tot_3c(t,x) = O2free_3c+O2bound_3c;   // M (mol O2_3c per L rbc)
real O2free_3p(t,x) = Wrbc*O2_3p; // M (mol O2_3p per L rbc)
real O2bound_3p(t,x) = 4*Hbrbc*SHbO2_3p; // M (mol O2_3p per L rbc)
real O2tot_3p(t,x) = O2free_3p+O2bound_3p;   // M (mol O2_3p per L rbc)
real O2free_3m(t,x) = Wrbc*O2_3m; // M (mol O2_3m per L rbc)
real O2bound_3m(t,x) = 4*Hbrbc*SHbO2_3m; // M (mol O2_3m per L rbc)
real O2tot_3m(t,x) = O2free_3m+O2bound_3m;   // M (mol O2_3m per L rbc)
real O2free_4c(t,x) = Wrbc*O2_4c; // M (mol O2_4c per L rbc)
real O2bound_4c(t,x) = 4*Hbrbc*SHbO2_4c; // M (mol O2_4c per L rbc)
real O2tot_4c(t,x) = O2free_4c+O2bound_4c;   // M (mol O2_4c per L rbc)
real O2free_4p(t,x) = Wrbc*O2_4p; // M (mol O2_4p per L rbc)
real O2bound_4p(t,x) = 4*Hbrbc*SHbO2_4p; // M (mol O2_4p per L rbc)
real O2tot_4p(t,x) = O2free_4p+O2bound_4p;   // M (mol O2_4p per L rbc)
real O2free_4m(t,x) = Wrbc*O2_4m; // M (mol O2_4m per L rbc)
real O2bound_4m(t,x) = 4*Hbrbc*SHbO2_4m; // M (mol O2_4m per L rbc)
real O2tot_4m(t,x) = O2free_4m+O2bound_4m;   // M (mol O2_4m per L rbc)
real O2free_5c(t,x) = Wrbc*O2_5c; // M (mol O2_5c per L rbc)
real O2bound_5c(t,x) = 4*Hbrbc*SHbO2_5c; // M (mol O2_5c per L rbc)
real O2tot_5c(t,x) = O2free_5c+O2bound_5c;   // M (mol O2_5c per L rbc)
real O2free_5p(t,x) = Wrbc*O2_5p; // M (mol O2_5p per L rbc)
real O2bound_5p(t,x) = 4*Hbrbc*SHbO2_5p; // M (mol O2_5p per L rbc)
real O2tot_5p(t,x) = O2free_5p+O2bound_5p;   // M (mol O2_5p per L rbc)
real O2free_5m(t,x) = Wrbc*O2_5m; // M (mol O2_5m per L rbc)
real O2bound_5m(t,x) = 4*Hbrbc*SHbO2_5m; // M (mol O2_5m per L rbc)
real O2tot_5m(t,x) = O2free_5m+O2bound_5m;   // M (mol O2_5m per L rbc)
real O2free_6c(t,x) = Wrbc*O2_6c; // M (mol O2_6c per L rbc)
real O2bound_6c(t,x) = 4*Hbrbc*SHbO2_6c; // M (mol O2_6c per L rbc)
real O2tot_6c(t,x) = O2free_6c+O2bound_6c;   // M (mol O2_6c per L rbc)
real O2free_6p(t,x) = Wrbc*O2_6p; // M (mol O2_6p per L rbc)
real O2bound_6p(t,x) = 4*Hbrbc*SHbO2_6p; // M (mol O2_6p per L rbc)
real O2tot_6p(t,x) = O2free_6p+O2bound_6p;   // M (mol O2_6p per L rbc)
real O2free_6m(t,x) = Wrbc*O2_6m; // M (mol O2_6m per L rbc)
real O2bound_6m(t,x) = 4*Hbrbc*SHbO2_6m; // M (mol O2_6m per L rbc)
real O2tot_6m(t,x) = O2free_6m+O2bound_6m;   // M (mol O2_6m per L rbc)
real O2free_7c(t,x) = Wrbc*O2_7c; // M (mol O2_7c per L rbc)
real O2bound_7c(t,x) = 4*Hbrbc*SHbO2_7c; // M (mol O2_7c per L rbc)
real O2tot_7c(t,x) = O2free_7c+O2bound_7c;   // M (mol O2_7c per L rbc)
real O2free_7p(t,x) = Wrbc*O2_7p; // M (mol O2_7p per L rbc)
real O2bound_7p(t,x) = 4*Hbrbc*SHbO2_7p; // M (mol O2_7p per L rbc)
real O2tot_7p(t,x) = O2free_7p+O2bound_7p;   // M (mol O2_7p per L rbc)
real O2free_7m(t,x) = Wrbc*O2_7m; // M (mol O2_7m per L rbc)
real O2bound_7m(t,x) = 4*Hbrbc*SHbO2_7m; // M (mol O2_7m per L rbc)
real O2tot_7m(t,x) = O2free_7m+O2bound_7m;   // M (mol O2_7m per L rbc)
real O2free_8c(t,x) = Wrbc*O2_8c; // M (mol O2_8c per L rbc)
real O2bound_8c(t,x) = 4*Hbrbc*SHbO2_8c; // M (mol O2_8c per L rbc)
real O2tot_8c(t,x) = O2free_8c+O2bound_8c;   // M (mol O2_8c per L rbc)
real O2free_8p(t,x) = Wrbc*O2_8p; // M (mol O2_8p per L rbc)
real O2bound_8p(t,x) = 4*Hbrbc*SHbO2_8p; // M (mol O2_8p per L rbc)
real O2tot_8p(t,x) = O2free_8p+O2bound_8p;   // M (mol O2_8p per L rbc)
real O2free_8m(t,x) = Wrbc*O2_8m; // M (mol O2_8m per L rbc)
real O2bound_8m(t,x) = 4*Hbrbc*SHbO2_8m; // M (mol O2_8m per L rbc)
real O2tot_8m(t,x) = O2free_8m+O2bound_8m;   // M (mol O2_8m per L rbc)
// Compute alphaO2 at given physiological conditions, and set initial guess for pO2.
 //   Input = [pO2_change,40,pHrbc,DPGrbc,Temp,Hbrbc,Hct];
 //   [Output] = SHbO2CO2_EJAP2016(Input); 
 //   alphaO2 = Output{1}(2);  
real pO2old_1(t,x) = pO2new_0; // Initial guess for pO2, otherwise it was previous iteration's value (pO2new_1)
// Newton-Raphson's iterative method for computation of pO2 from TotO2
// while (iter < iter_max && err > err_tol)
// iter:
  pO2c_1 = pO2old_1;
  pO2p_1 = pO2old_1 + 1e-2*pO2old_1;
  pO2m_1 = pO2old_1 - 1e-2*pO2old_1;
//    Inputc = [pO2c1,pCO2,pHrbc,DPGrbc,Temp,Hbrbc,Hct];
//    [Outputc] = SHbO2O2_EJAP2016(Inputc); O2totc1 = Outputc{2}(1);
  dO2tot1c = (O2tot_1p-O2tot_1m)/(pO2p_1-pO2m_1);  // Derivative by central difference
  funcO2_1 = TotO2_Invert_in - O2tot_1c;   // Function to be solved (find pO2 such that funcO2_1 = 0) 
  dfuncO2_1 = -dO2tot1c;        // As TotO2 is an input (constant), dTotO2 = 0
  pO2new_1 = (pO2old_1 - funcO2_1/dfuncO2_1) ;// Newton-Raphson iterative formula
real errO2_1(t,x) = abs(pO2new_1 - pO2old_1)/pO2new_1;   // use this as a check.
real pO2old_2(t,x) = pO2new_1; // Initial guess for pO2, otherwise it was previous iteration's value (pO2new_2)
  pO2c_2 = pO2old_2;
  pO2p_2 = pO2old_2 + 1e-2*pO2old_2;
  pO2m_2 = pO2old_2 - 1e-2*pO2old_2;
  dO2tot2c = (O2tot_2p-O2tot_2m)/(pO2p_2-pO2m_2);  // Derivative by central difference
  funcO2_2 = TotO2_Invert_in - O2tot_2c;   // Function to be solved (find pO2 such that funcO2_2 = 0) 
  dfuncO2_2 = -dO2tot2c;        // As TotO2 is an input (constant), dTotO2 = 0
  pO2new_2 = (pO2old_2 - funcO2_2/dfuncO2_2) ;// Newton-Raphson iterative formula
real errO2_2(t,x) = abs(pO2new_2 - pO2old_2)/pO2new_2;   // use this as a check.
real pO2old_3(t,x) = pO2new_2; // Initial guess for pO2, otherwise it was previous iteration's value (pO2new_3)
  pO2c_3 = pO2old_3;
  pO2p_3 = pO2old_3 + 1e-2*pO2old_3;
  pO2m_3 = pO2old_3 - 1e-2*pO2old_3;
  dO2tot3c = (O2tot_3p-O2tot_3m)/(pO2p_3-pO2m_3);  // Derivative by central difference
  funcO2_3 = TotO2_Invert_in - O2tot_3c;   // Function to be solved (find pO2 such that funcO2_3 = 0) 
  dfuncO2_3 = -dO2tot3c;        // As TotO2 is an input (constant), dTotO2 = 0
  pO2new_3 = (pO2old_3 - funcO2_3/dfuncO2_3) ;// Newton-Raphson iterative formula
real errO2_3(t,x) = abs(pO2new_3 - pO2old_3)/pO2new_3;   // use this as a check.
real pO2old_4(t,x) = pO2new_3; // Initial guess for pO2, otherwise it was previous iteration's value (pO2new_4)
  pO2c_4 = pO2old_4;
  pO2p_4 = pO2old_4 + 1e-2*pO2old_4;
  pO2m_4 = pO2old_4 - 1e-2*pO2old_4;
  dO2tot4c = (O2tot_4p-O2tot_4m)/(pO2p_4-pO2m_4);  // Derivative by central difference
  funcO2_4 = TotO2_Invert_in - O2tot_4c;   // Function to be solved (find pO2 such that funcO2_4 = 0) 
  dfuncO2_4 = -dO2tot4c;        // As TotO2 is an input (constant), dTotO2 = 0
  pO2new_4 = (pO2old_4 - funcO2_4/dfuncO2_4) ;// Newton-Raphson iterative formula
real errO2_4(t,x) = abs(pO2new_4 - pO2old_4)/pO2new_4;   // use this as a check.
real pO2old_5(t,x) = pO2new_4; // Initial guess for pO2, otherwise it was previous iteration's value (pO2new_5)
  pO2c_5 = pO2old_5;
  pO2p_5 = pO2old_5 + 1e-2*pO2old_5;
  pO2m_5 = pO2old_5 - 1e-2*pO2old_5;
  dO2tot5c = (O2tot_5p-O2tot_5m)/(pO2p_5-pO2m_5);  // Derivative by central difference
  funcO2_5 = TotO2_Invert_in - O2tot_5c;   // Function to be solved (find pO2 such that funcO2_5 = 0) 
  dfuncO2_5 = -dO2tot5c;        // As TotO2 is an input (constant), dTotO2 = 0
  pO2new_5 = (pO2old_5 - funcO2_5/dfuncO2_5) ;// Newton-Raphson iterative formula
real errO2_5(t,x) = abs(pO2new_5 - pO2old_5)/pO2new_5;   // use this as a check.
real pO2old_6(t,x) = pO2new_5; // Initial guess for pO2, otherwise it was previous iteration's value (pO2new_6)
  pO2c_6 = pO2old_6;
  pO2p_6 = pO2old_6 + 1e-2*pO2old_6;
  pO2m_6 = pO2old_6 - 1e-2*pO2old_6;
  dO2tot6c = (O2tot_6p-O2tot_6m)/(pO2p_6-pO2m_6);  // Derivative by central difference
  funcO2_6 = TotO2_Invert_in - O2tot_6c;   // Function to be solved (find pO2 such that funcO2_6 = 0) 
  dfuncO2_6 = -dO2tot6c;        // As TotO2 is an input (constant), dTotO2 = 0
  pO2new_6 = (pO2old_6 - funcO2_6/dfuncO2_6) ;// Newton-Raphson iterative formula
real errO2_6(t,x) = abs(pO2new_6 - pO2old_6)/pO2new_6;   // use this as a check.
real pO2old_7(t,x) = pO2new_6; // Initial guess for pO2, otherwise it was previous iteration's value (pO2new_7)
  pO2c_7 = pO2old_7;
  pO2p_7 = pO2old_7 + 1e-2*pO2old_7;
  pO2m_7 = pO2old_7 - 1e-2*pO2old_7;
  dO2tot7c = (O2tot_7p-O2tot_7m)/(pO2p_7-pO2m_7);  // Derivative by central difference
  funcO2_7 = TotO2_Invert_in - O2tot_7c;   // Function to be solved (find pO2 such that funcO2_7 = 0) 
  dfuncO2_7 = -dO2tot7c;        // As TotO2 is an input (constant), dTotO2 = 0
  pO2new_7 = (pO2old_7 - funcO2_7/dfuncO2_7) ;// Newton-Raphson iterative formula
real errO2_7(t,x) = abs(pO2new_7 - pO2old_7)/pO2new_7;   // use this as a check.
real pO2old_8(t,x) = pO2new_7; // Initial guess for pO2, otherwise it was previous iteration's value (pO2new_8)
  pO2c_8 = pO2old_8;
  pO2p_8 = pO2old_8 + 1e-2*pO2old_8;
  pO2m_8 = pO2old_8 - 1e-2*pO2old_8;
  dO2tot8c = (O2tot_8p-O2tot_8m)/(pO2p_8-pO2m_8);  // Derivative by central difference
  funcO2_8 = TotO2_Invert_in - O2tot_8c;   // Function to be solved (find pO2 such that funcO2_8 = 0) 
  dfuncO2_8 = -dO2tot8c;        // As TotO2 is an input (constant), dTotO2 = 0
  pO2new_8 = (pO2old_8 - funcO2_8/dfuncO2_8) ;// Newton-Raphson iterative formula
real errO2_8(t,x) = abs(pO2new_8 - pO2old_8)/pO2new_8;   // use this as a check.
pO2_out = pO2new_8;              // Final Output of inversion
FO2_out = alphaO2Sys*pO2_out; 
SHbO2_out = SHbO2_8c;        // 
//%END O2INVERT
// Set initial values:
Hct_in =0.45;
Temp_in=310;  // Physiological temperature (function of t,x)
TotO2_Invert_in = 0.02827;
pCO2_Invert_in = 40;
pHrbc_in = 7.238;
DPGrbc_in =.00465; 
Hbrbc_in = .0052;
//alphaO2_in = 1.46E-6;
P50_in = 26.8;
real pCO2rbct0= 40 mmHg;
pCO2t0_in= pCO2rbct0;
}
// This MML file generated from TotO2freeO2Invert.mpc using MPC v1.01.
  
