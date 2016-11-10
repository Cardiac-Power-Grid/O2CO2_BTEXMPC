// Calculate free CO2 in rbc based on total CO2
import nsrunit;
unit conversion on;

math TotCO2_FreeCO2Invert {

//%START TotCO2_FreeCO2_Inversion
// *****************************************************************************************
// *****************************************************************************************
// TotCO2_FreeCO2_Inversion.m

//%START DOMAIN
realDomain t sec; t.min=0; t.max=100; t.delta=0.1;
//%END DOMAIN

//%START SPATIALDOMAIN
realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid;   // BTEX space domain 
//%END SPATIALDOMAIN

// function [pO2,fO2] = TotO2_FreeO2_Inversion(TotO2,pCO2,pHrbc,DPGrbc,Temp,Hbrbc,Hct)
// Input:
// 1. TotCO2 ->TCO2, pO2 ->PCO2, pHrbc->pHRBC
// Output:
// 1. pCO2->PO2, fCO2->FCO2

//%START INPUTS
real TotCO2_Invert_in(t,x) M;  // Tot CO2 used for inversion to get pCO2  
real pO2in(t,x) mmHg;  // pO2 used for inversion to get pCO2
real pHrbc_in(t,x);    // pHrbc used for inversion to get pO2/CO2
real DPGrbc_in(t,x) M; // DPGrbc used for inversion to get pO2/CO2
real Temp_in(t,x) K;   // Temp used for inversion to get pO2/CO2
real Hbrbc_in(t,x) M;  // Hbrbc used for inversion to get pO2/CO2
real Hct_in;           // Vol fraction of RBCs in blood (Hematocrit), needed?
real alphaCO2_in(t,x) M/mmHg; // solubility of CO2 in plasma.
real alphaO2_in(t,x) M/mmHg; // solubility of O2 in water. Calc can be done in SHbO2CO2_EJAP2016.
real pCO2initguess =40 mmHg;   // Init guess.
real pCO2new_0(t,x) = pCO2initguess;  
//%END INPUTS


//%START InvertOUTPUTS
// TotCO2_FreeCO2_Inversion: calcs pCO2 given TotCO2, pO2, alphaCO2, 
real pCO2_out(t,x) mmHg;  // free CO2 in rbc
real FCO2_out(t,x) M;     // free CO2 in rbc
real SHbCO2_out(t,x);     // Value from TCO2 inversion calcs
//%END InvertOUTPUTS


//%START DEPENDS
// Get Total pO2 value from iterative approx with SHbO2CO2_EJAP2016
real CO2totc(t,x) M;  // from SHbO2CO2_EJAP2016 calc
real CO2totp(t,x) M;  // from SHbO2CO2_EJAP2016 calc
real CO2totm(t,x) M;  // from SHbO2CO2_EJAP2016 calc
//%END DEPENDS

//%START IterativeVars
// IterativeVars
real pCO2c(t,x) mmHg;
real pCO2p(t,x) mmHg;
real pCO2m(t,x) mmHg;

real dCO2totc(t,x);
real funcCO2(t,x);
real dfuncCO2(t,x);
real pCO2new(t,x) mmHg;

//%END IterativeVars

// ----------- TotCO2InvertCalcs ----------------------------------------
//%START TotCO2InvertCalcs
real pCO2old(t,x) = pCO2new_0; // Initial guess for pCO2, otherwise it was previous iteration's value (pCO2new)
// Newton-Raphson's iterative method for computation of pCO2 from TotCO2
// iter:
  pCO2c = pCO2old;
  pCO2p = pCO2old + 1e-2*pCO2old;
  pCO2m = pCO2old - 1e-2*pCO2old;

//    Inputc = [pCO2c1,pO2,pHrbc,DPGrbc,Temp,Hbrbc,Hct];
//    [Outputc] = SHbO2CO2_EJAP2016(Inputc); CO2totc1 = Outputc{2}(1);

  dCO2totc = (CO2totp-CO2totm)/(pCO2p-pCO2m);  // Derivative by central difference
  funcCO2 = TotCO2in - CO2totc;   // Function to be solved (find pCO2 such that funcO2 = 0) 
  dfuncCO2 = -dCO2totc;        // As TotCO2 is an input (constant), dTotCO2 = 0
  pCO2new = (pCO2old - funcCO2/dfuncCO2) ;// Newton-Raphson iterative formula

real errCO2(t,x) = abs(pCO2new - pCO2old)/pCO2new;   // use this as a check.
//%END TotCO2InvertCalcs

/* Not used as stand-alone module:
//%START CallSHbO2CO2
// Assign inputs to SHbO2CO2 calcs to get pCO2 given TCO2 in rbc
 Hct = Hct_in;
 pCO2_c=pCO2c;
 pCO2_p=pCO2p;
 pCO2_m=pCO2m;
 pHrbc=pHrbc_in;
 DPGrbc = DPGrbc_in;
 pO2= pO2in;
 TempSHbO2CO2 =Temp_in;  // Used in Dash2016SHbO2CO2MPC
//%END CallSHbO2CO2
*/

//%START FinalResult
pCO2_out = pCO2new;              // Final Output of inversion
FCO2_out = alphaCO2_in*pCO2_out; // 
SHbCO2_out = SHbCO2;             // Final SHbCO2 from inversion
//%END FinalResult


// *********************************
// Make stand-alone model (Inputs):
TotCO2_Invert_in = .028;
pO2in = 100;
pHrbc_in = 7.24;
DPGrbc_in = .00465;
Temp_in = 310;
Hbrbc_in = .0052;
Hct_in =.45;   
alphaCO2_in = ;
alphaO2_in =;
//P50_in = 26.8;
real KHbCO2 = 2.7328E5 1/M;
real SHbCO2(t,x);
real Wpl = 0.94 dimensionless;    // fractional water space of plasma; unitless
real Wrbc = 0.65;                 // fractional water space of rbc; unitless
real Wbl = .8095 dimensionless;   // fractional water space of blood; unitless
real Hpl = 3.9705E-8 M; // conc of protons in plasma
real Rrbc = 0.69;       // Gibbs-Donnan ratio across the RBC membrane
// Calcs: 
SHbCO2 = KHbCO2*alphaCO2_in*pCO2_out/(1+KHbCO2*alphaCO2_in*pCO2_out);
real CO2free(t,x) = Wbl*alphaCO2_in*pCO2_out;  // M 
real CO2bound(t,x) = 4// Calculate free CO2 in blood based on total CO2
import nsrunit;
unit conversion on;

}

