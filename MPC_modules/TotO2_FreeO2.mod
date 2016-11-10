// Calculate free O2 in RBC based on total O2
import nsrunit;
unit conversion on;

math TotO2_FreeO2Invert {

//%START TotO2_FreeO2_Inversion
// *****************************************************************************************
// *****************************************************************************************
// TotO2_FreeO2_Inversion.m

//%START DOMAIN
realDomain t sec; t.min=0; t.max=100; t.delta=0.1;
//%END DOMAIN

//%START SPATIALDOMAIN
realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid;   // BTEX space domain 
//%END SPATIALDOMAIN

// function [pO2,fO2] = TotO2_FreeO2_Inversion(TotO2,pCO2,pHrbc,DPGrbc,Temp,Hbrbc,Hct)
// Mapping:
// Input:
// 1. TotO2 ->TO2PulCap, pCO2 ->PCO2PulCap, pHrbc->pHPulCapRBC
// 2. TotO2 ->TO2SysCap, pCO2 ->PCO2SysCap, pHrbc->pHSysCapRBC
// Output:
// 1. pO2->O2PulCap, fO2->FCO2PulCap
// 2. pO2->O2SysCap, fO2->FCO2SysCap

//%START INPUTS
real TotO2_Invert_in(t,x) M;   // Tot O2 used for inversion to get pO2 
real pCO2in(t,x) mmHg; // pCO2 used for inversion to get pO2
real pHrbc_in(t,x);    // pHrbc used for inversion to get pO2/CO2
real DPGrbc_in(t,x) M; // DPGrbc used for inversion to get pO2/CO2
real Temp_in(t,x) K;   // Temp used for inversion to get pO2/CO2
real Hbrbc_in(t,x) M;  // Hbrbc used for inversion to get pO2/CO2
real Hct_in;           // Vol fraction of RBCs in blood (Hematocrit), needed?
real alphaO2_in(t,x) M/mmHg; // solubility of O2 in water. Calc can be done in SHbO2CO2_EJAP2016.
real P50_in mmHg;   // Init pO2 guess (P50 for O2) for inversion to get pO2.
real pO2new_0(t,x) = P50_in;
real pCO2t0_in(x) mmHg;  // Initial pCO2 at t=0 for all x
//%END INPUTS


//%START InvertOUTPUTS
// TotO2_FreeO2_Inversion: calcs pO2 given TotO2, pCO2, alphaO2, 
real pO2_out(t,x) mmHg;  // free O2 in rbc
real FO2_out(t,x) M;     // free O2 in rbc
real SHbO2_out(t,x);    
//%END InvertOUTPUTS


//%START DEPENDS
// Get Total pO2 value from iterative approx with SHbO2CO2_EJAP2016
real O2totc(t,x) M;  // from SHbO2CO2_EJAP2016 calc
real O2totp(t,x) M;  // from SHbO2CO2_EJAP2016 calc
real O2totm(t,x) M;  // from SHbO2CO2_EJAP2016 calc
//%END DEPENDS

//%START IterativeVars
// IterativeVars
real pO2c(t,x) mmHg;
real pO2p(t,x) mmHg;
real pO2m(t,x) mmHg;

real dO2totc(t,x);
real funcO2(t,x);
real dfuncO2(t,x);
real pO2new(t,x) mmHg;

//%END IterativeVars


//%START TotO2InvertCalcs
// Compute alphaO2 at given physiological conditions, and set initial guess for pO2.
 //   Input = [pO2_change,40,pHrbc,DPGrbc,Temp,Hbrbc,Hct];
 //   [Output] = SHbO2CO2_EJAP2016(Input); 
 //   alphaO2 = Output{1}(2);  
real pO2old(t,x) = pO2new_0; // Initial guess for pO2, otherwise it was previous iteration's value (pO2new)
// Newton-Raphson's iterative method for computation of pO2 from TotO2
// while (iter < iter_max && err > err_tol)

// iter:
  pO2c = pO2old;
  pO2p = pO2old + 1e-2*pO2old;
  pO2m = pO2old - 1e-2*pO2old;

//    Inputc = [pO2c1,pCO2,pHrbc,DPGrbc,Temp,Hbrbc,Hct];
//    [Outputc] = SHbO2O2_EJAP2016(Inputc); O2totc1 = Outputc{2}(1);

  dO2totc = (O2totp-O2totm)/(pO2p-pO2m);  // Derivative by central difference
  funcO2 = TotO2in - O2totc;   // Function to be solved (find pO2 such that funcO2 = 0) 
  dfuncO2 = -dO2totc;        // As TotO2 is an input (constant), dTotO2 = 0
  pO2new = (pO2old - funcO2/dfuncO2) ;// Newton-Raphson iterative formula

real errO2(t,x) = abs(pO2new - pO2old)/pO2new;   // use this as a check.
//%END TotO2InvertCalcs

/* Not used as stand-alone module:
//%START CallSHbO2CO2
// Assign an declare variable inputs to SHbO2CO2 calcs to get pO2 given TO2 in rbc
 Hct = Hct_in;
 pO2_c=pO2c;
 pO2_p=pO2p;
 pO2_m=pO2m;
 pHrbc=pHrbc_in;
 DPGrbc = DPGrbc_in;
 TempSHbO2CO2 =Temp_in;  // Used in Dash2016SHbO2CO2MPC
realState pCO2(t,x);
  when (t=t.min) pCO2 = pCO2t0_in;
  event(t>t.min) {pCO2= pCO2in;} // wait untill pCO2rbc calculated
 
//%END CallSHbO2CO2
*/

//%START FinalResult
pO2_out = pO2new;              // Final Output of inversion
FO2_out = alphaO2_in*pO2_out; 
SHbO2_out = SHbO2;        // 
//%END FinalResult


// *********************************
// Make stand-alone model (Inputs):
TotO2_Invert_in = .028;
pCO2in = 45;
pHrbc_in = 7.24;
DPGrbc_in = .00465;
Temp_in = 310;
Hbrbc_in = .0052;
Hct_in =.45;   
alphaO2_in = ;
pCO2t0_in = 40;
//P50_in = 26.8;
real KHbO2 = 2.7328E5 1/M;
real SHbO2(t,x);
real Wpl = 0.94 dimensionless;    // fractional water space of plasma; unitless
real Wrbc = 0.65;                 // fractional water space of rbc; unitless
real Wbl = .8095 dimensionless;   // fractional water space of blood; unitless
real Hpl = 3.9705E-8 M; // conc of protons in plasma
real Rrbc = 0.69;       // Gibbs-Donnan ratio across the RBC membrane
// Calcs: 
SHbO2 = KHbO2*alphaO2_in*pO2_out/(1+KHbO2*alphaO2_in*pO2_out);
real O2free(t,x) = Wrbc*alphaO2_in*pO2_out;  // M 
real O2bound(t,x) = 4*Hbrbc_in*SHbO2; // M 
O2totc = O2free+O2bound;   // from SHbO2CO2_EJAP2016 calc
O2totp = O2totc;  // from SHbO2CO2_EJAP2016 calc
O2totm = O2totc;  // from SHbO2CO2_EJAP2016 calc

}
