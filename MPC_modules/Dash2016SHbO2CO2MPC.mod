// MODEL NAME: Dash2016SHbO2CO2MPC
// SHORT DESCRIPTION: Simulation of oxyhemoglobin (HbO2) and carbaminohemoglobin (HbCO2) 
// dissociation curves and computation of total O2 and CO2 contents in RBCs, 
// Modified from Dash's 2016 Matlab version. Annotated for use with MPC.

import nsrunit;
unit conversion on;

math SHbO2CO2_mod {
// function [Output] = SHbO2CO2_EJAP2016(Input)
// Input physiological state variables for calculations of SHbO2CO2 and blood O2CO2 contents
// pO2 = Input(1); pCO2 = Input(2); pHrbc = Input(3); DPGrbc = Input(4);
// Temp = Input(5); Hbrbc = Input(6); Hct = Input(7);

//%START TIMEDOMAIN
realDomain t sec; t.min=0; t.max=100; t.delta=0.1;
//%END TIMEDOMAIN

//%START SPATIALDOMAIN
realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid; // BTEX space domain 
//%END SPATIALDOMAIN

//%START PO2INPUT
real pO2(t,x) mmHg;    // partial pressure in rbc
//%END PO2INPUT

//%START PCO2INPUT
real pCO2(t,x) mmHg;   // partial pressure in rbc
//%END PCO2INPUT

//%START OTHERINPUTS
real pHrbc(t,x) dimensionless; real DPGrbc(t,x) M ; // pH and DPG conc in RBCs
real Temp(t,x) K;          // Physiological temperature, Used in Dash2016SHbO2CO2MPC
real Hct dimensionless ;   // Volume fraction of RBCs in blood (Hematocrit), unitless
//%END OTHERINPUTS

//%START CONSTANTS
real mol2ml = 22400 mol/ml;          // Conversion factor from mol of gas to ml of gas at STP
real AmHbBl = Hct*(150 g)/(0.45 L);  // Amount of hemoglobin in gm per liter of blood, gm/L
real MWHb = 64100 g/mol;             // Molecular weight of hemoglobin, gm/mol
real HbBl M;
     HbBl = AmHbBl/MWHb;             // Concentration of hemoglobin in blood, mol/L or M
//%END CONSTANTS

//%START PARAMETERS 
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
//%END PARAMETERS

//%START STDCONDITIONS
// Variables those are fixed in the model with values at standard physiological conditions
// (i.e. pO20, pCO20, pHpl0, pHrbc0, DPGrbc0, Temp0)
real pO20 = 100 mmHg;               // standard O2 partial pressure in blood; mmHg
real pCO20 = 40 mmHg;               // standard CO2 partial pressure in blood; mmHg
real pHrbc0 = 7.24 dimensionless;   // standard pH in RBCs; unitless
real pHpl0 = pHrbc0-log(Rrbc);      // standard pH in plsama; unitless
real DPGrbc0 = 4.65e-3;             // standard 2,3-DPG concentration in RBCs; M
real Temp0 = 310 K;                 // standard temperature in blood; deg K
real P500 = 26.8 mmHg;              // standard pO2 for 50% SHbO2; mmHg
//%END STDCONDITIONS

//%START INTERMEDIATE_VARS
// Calculation of intermediate variables in the computations of SHbO2 and SHbCO2
      // Rrbc = Hpl/Hrbc = 10^-(pHpl-pHrbc)
real delpHrbc(t,x) = pHrbc-pHrbc0; 
real delpCO2(t,x) = pCO2-pCO20; 
real delDPGrbc(t,x) = DPGrbc-DPGrbc0;
real Hrbc(t,x) M;
     Hrbc  = 10^(-pHrbc); 
real delTemp(t,x) = Temp-Temp0;
//%END INTERMEDIATE_VARS

//%START O2CO2solubilities
real alphaO20 = 1.46e-6 M/mmHg;     // solubility of O2 in water/plasma at 37 C; M/mmHg
real alphaCO20 = 32.66e-6 M/mmHg;   // solubility of CO2 in water/plasma at 37 C; M/mmHg
real alphaO2(t,x) M/mmHg;        // solubility of O2 in water/plasma at experimental temp
real alphaCO2(t,x) M/mmHg;       // solubility of CO2 in water/plasma at experimental temp
     alphaO2 = alphaO20*(1 - 1e-2*delTemp*(1 K^-1) + 4.234e-4*(1 K^-2)*delTemp^2);  // Corrected solubility of O2
     alphaCO2 = alphaCO20*(1 - 1.86e-2*delTemp*(1 K^-1) + 6.515e-4*(1 K^-2)*delTemp^2);  // Corrected solubility of CO2
//%END O2CO2solubilities

//%START InputConversions
real O2(t,x) M;             // free O2
     O2 = alphaO2*pO2; 
real CO2(t,x) M;            // free CO2
     CO2 = alphaCO2*pCO2;  
//%END InputConversions

//%START P50_adj
// p50 adjustments:
real P501(t,x) = P500 + 1.2*(1 mmHg)*(-21.279*delpHrbc + 8.872*delpHrbc^2 - 1.47*delpHrbc^3); // all standard conditions, except pH
real P502(t,x) = P500 + 1.7*(4.28e-2*delpCO2 + 3.64e-5*(1 mmHg^-1)*delpCO2^2); // all standard conditions, except CO2
real P503(t,x) = P500 + 1.0*(1 mmHg)*(795.633533*(1 M^-1)*delDPGrbc - 19660.8947*(1 M^-2)*delDPGrbc^2); // all standard conditions, except DPG
real P504(t,x) = P500 + 0.98*(1 mmHg)*(1.4945*(1 K^-1)*delTemp + 4.335e-2*(1 K^-2)*delTemp^2 + 7e-4*(1 K^-3)*delTemp^3); // all standard conditions, except T
real P50(t,x) mmHg;             // Final P50 value for 50% O2 binding to Hb
     P50 = P500*(P501/P500)*(P502/P500)*(P503/P500)*(P504/P500); 
real C50(t,x) M;
     C50 = alphaO2*P50; 
//%END P50_adj


//%START O2CO2DISSOC
// *** START O2 CO2 dissoc calcs:
//%START EquilConsts
// equilibirum consts are not directly tied to CO2 and O2 conc:
real BPH1(t,x) dimensionless;
     BPH1 = 1+K2dp/Hrbc;         // Binding polynomial involving K2dp and Hrbc 
real BPH2(t,x) = 1+K3dp/Hrbc;      // Binding polynomial involving K3dp and Hrbc 
real BPH3(t,x) = 1+Hrbc/K5dp;      // Binding polynomial involving K5dp and Hrbc 
real BPH4(t,x) = 1+Hrbc/K6dp;      // Binding polynomial involving K6dp and Hrbc 
//%END EquilConsts

//%START HillCoeff
// Hill coefficient; unitless (redefined as a function pO2)
real alpha = 2.8 dimensionless; real beta = 1.2 dimensionless; real gamma = 29.2 mmHg; // Roughton et al data
real nH(t,x) dimensionless;
     nH = alpha-beta*10^(-pO2/gamma); // pO2 dependent variable nH
//%END HillCoeff

//%START EQUIL_Calcs
//%START EQUILK4p
// Compute the apparent equilibrium constant of Hb with O2 and CO2 (KHbO2 and KHbCO2); 
// O2 and CO2 saturations of Hb (SHbO2 and SHbCO2); and O2 and CO2 contents in blood. 
real K4p(t,x) 1/M;  
     K4p = (1 M^-1)*((O2*(1 M^-1))^(nH-1)*(K2p*BPH1*CO2+BPH3))/( ((1 M^-1)*C50)^nH*(K3p*BPH2*CO2+BPH4));
//%END EQUILK4p

//%START HbO2Saturate
real KHbO2(t,x) 1/M;               // Apparent equilibrium constants of Hb with O2
     KHbO2 = K4p*(K3p*BPH2*CO2+BPH4)/(K2p*BPH1*CO2+BPH3); 
real SHbO2(t,x) dimensionless;     // Saturation of hemoglobin with oxygen
     SHbO2 = KHbO2*O2/(1+KHbO2*O2);
//%END HbO2Saturate

//%START HbCO2Saturate
real KHbCO2(t,x) 1/M;              // Apparent equilibrium constants of Hb with CO2
     KHbCO2 = (K2p*BPH1+K3p*K4p*BPH2*O2)/(BPH3+K4p*BPH4*O2);
real SHbCO2(t,x) dimensionless;
     SHbCO2 = KHbCO2*CO2/(1+KHbCO2*CO2);
//%END HbCO2Saturate
//%END EQUIL_Calcs

//%START O2calcs
real O2free(t,x) = Wrbc*O2; // M (mol O2 per L rbc)
real O2bound(t,x) = 4*Hbrbc*SHbO2; // M (mol O2 per L rbc)
real O2tot(t,x) = O2free+O2bound;   // M (mol O2 per L rbc)
//%END O2calcs

//%START CO2calcs
real CO2free(t,x) = Wrbc*CO2;  // M (mol CO2 per L rbc)
real CO2bicarb(t,x) = (Wrbc*Rrbc)*(K1*CO2/Hrbc); // M (mol CO2 per L rbc)
real CO2bound(t,x) = 4*Hbrbc*SHbCO2; // M (mol CO2 per L rbc)
real CO2tot1(t,x) = CO2free+CO2bound; // M (mol CO2 per L rbc)
real CO2tot2(t,x) = CO2free+CO2bicarb+CO2bound; // M (mol CO2 per L rbc)
//%END CO2calcs

//Output{1} = [alphaO2,alphaCO2,P50,K4p,KHbO2,KHbCO2,SHbO2,SHbCO2];
//Output{2} = [O2tot,CO2tot1,CO2tot2,O2cont,CO2cont1,CO2cont2];
// *** END of O2 CO2 DISSOC calcs
//%END O2CO2DISSOC

//%START O2CO2DISSOC_ALTCALC
// ALTERNATIVE CALCULATIONS FOR SHBO2 AND SHBCO2 BASED ON DIFFERENT HB-BOUND SPECIES
real HbNH2(t,x) = Hbrbc/((K2p*CO2*BPH1+BPH3)+K4p*O2*(K3p*CO2*BPH2+BPH4));
real HbNH3p(t,x) = HbNH2*Hrbc/K5dp;
real O2HbNH2(t,x) = K4p*O2*HbNH2;
real O2HbNH3p(t,x) = O2HbNH2*Hrbc/K6dp;
real HbNHCOOH(t,x) = K2p*CO2*HbNH2;
real HbNHCOOm(t,x) = K2dp*HbNHCOOH/Hrbc;
real O2HbNHCOOH(t,x) = K3p*CO2*O2HbNH2;
real O2HbNHCOOm(t,x) = K3dp*O2HbNHCOOH/Hrbc;
real SHbO2kin(t,x) = (O2HbNH2+O2HbNH3p+O2HbNHCOOH+O2HbNHCOOm)/Hbrbc;
real SHbCO2kin(t,x) = (HbNHCOOH+HbNHCOOm+O2HbNHCOOH+O2HbNHCOOm)/Hbrbc;
//%END O2CO2DISSOC_ALTCALC

// Make stand-alone (Inputs):
pO2 = 100;  pCO2 =40 ;     // free gas in RBCs
pHrbc = 7.24; DPGrbc = 4.65e-3;
Temp = 310;   // Physiological temperature (function of t and x),
Hct = 0.45 ;  // Volume fraction of RBCs in blood (Hematocrit), unitless

}




//---------------------------------------------------------------------------------------
// Simulation of oxyhemoglobin (HbO2) and carbaminohemoglobin (HbCO2) dissociation
// curves and computation of total O2 and CO2 contents in whole blood, revised from
// the original model of Dash and Bassingthwaighte, ABME 38(4):1683-1701, 2010. The
// revision makes the model further simplified, as it bipasses the computations of
// the indices n1, n2, n3 and n4, which are complex expressions. Rather the revision
// necessiates the computations of K4p in terms of P50. Also the calculations of P50
// in terms of pH is enhanced based on a 3rd degree polynomial interpolation. 
//---------------------------------------------------------------------------------------
// Developed by: Ranjan Dash, PhD (Last modified: 2/29/2016 from 3/15/2015 version)
// Department of Physiology and Biotechnology and Bioengineering Center
// Medical College of Wisconsin, Milwaukee, WI-53226
//---------------------------------------------------------------------------------------


//---------------------------------------------------------------------------------------
// BIOCHEMICAL REACTIONS FOR DERIVATION OF THE SHBO2 AND SHBCO2 EQUATIONS----------------
// The equations for O2 and CO2 saturations of hemoglobin (SHbO2 and SHbCO2) are  
// derived by considering the various kinetic reactions involving the binding of
// O2 and CO2 with hemoglobin in RBCs:
//
//            kf1p       K1dp
// 1. CO2+H2O <--> H2CO3 <--> HCO3- + H+;  K1=(kf1p/kb1p)*K1dp
//            kb1p		K1 = 7.43e-7 M; K1dp = 5.5e-4 M
//
//              kf2p          K2dp
// 2. CO2+HbNH2 <--> HbNHCOOH <--> HbNHCOO- + H+;  K2=(kf2p/kb2p)*K2dp
//              kb2p		K2 = 21.5e-6; K2dp = 1.0e-6 M
//
//                kf3p            K3dp
// 3. CO2+O2HbNH2 <--> O2HbNHCOOH <--> O2HbNHCOO- + H+; K3=(kf3p/kb3p)*K3dp
//                kb3p		K3 = 11.3e-6; K3dp = 1.0e-6 M
//
//              kf4p          
// 4. O2+HbNH2 <--> O2HbNH2;  K4p=func([O2];[H+];[CO2];[DPG];T)
//              kb4p
//
//           K5dp
// 5. HbNH3+ <--> HbNH2 + H+; K5dp = 2.4e-8 M
//
//             K6dp
// 6. O2HbNH3+ <--> O2HbNH2 + H+; K6dp = 1.2e-8 M
//---------------------------------------------------------------------------------------





