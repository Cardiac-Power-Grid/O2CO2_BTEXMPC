// MODEL NAME: Bloodpde_MPC
// SHORT DESCRIPTION: Simple model of Blood region that calculates the pO2 and pCO2
// in RBC based on Total O2, CO2, H+, and HCO3. Takes into account temperature
// dependence of O2/CO2 solubility.


//%START BLOOD_MODULE_INFO

// * RBC Variables *:
// TCO2_rbc(t,x) M;     // Total CO2 in rbc
// TO2_rbc(t,x) M;     // Total O2 in rbc
// pCO2_rbc(t,x) mmHg; // free pCO2 in rbc
// CO2_rbc(t,x) M;     // free CO2 in rbc
// O2_rbc(t,x) M;      // free O2 in rbc
// pO2_rbc(t,x) mmHg;  // free pO2 in rbc
// Hp_rbc(t,x) M;      // H+ conc in RBCs
// HCO3m_rbc(t,x) M;   // HCO3- in RBCs

// * RBC Input variables *:
// CO2_pl(t,x) M;     // Free CO2 in plasma
// O2_pl(t,x) M;      // Free O2 in plasma
// HCO3m_pl(t,x) M;   // HCO3m in plasma
// Hp_pl(t,x) M;      // H+ in plasma 

// * RBC Flow input variables (Conc of substrate coming into pipe at x=x.min):
// HCO3m_rbc_Fin(t) M;  // HCO3- coming into the rbc region (x=x.min)
// Hp_rbc_Fin(t) M;     // H+ coming into the rbc region (x=x.min)
// pCO2_rbc_Fin(t) mmHg; // pCO2 coming into the rbc region (x.min), 
                         // TCO2_rbc_Fin calculted from this and initial SHbCO2_rbc0.
// pO2_rbc_Fin(t) mmHg;  // pO2 coming into the rbc region (x.min), 
                         // TO2_rbc_Fin calculted from this and initial SHbO2_rbc0.

// pH_pl_Fin(t) dimensionless;  // pH of plasma coming into region at x.min
// * Plasma Variables *:
// CO2_pl(t,x) M;     // CO2 in plasma
// pCO2_pl(t,x) mmHg; // pCO2 in plasma
// O2_pl(t,x) M;      // O2 in plasma
// pO2_pl(t,x) mmHg;  // pO2 in plasma
// Hp_pl(t,x) M;      // H+ conc in plasma
// pH_pl(t,x) dimensionless;      // pH in plasma
// HCO3m_pl(t,x) M;   // HCO3- in plasma

// * Plasma Input variables *:
// CO2_rbc(t,x) M
// CO2_isf(t,x) M
// O2_rbc(t,x) M
// O2_isf(t,x) M
// Hp_rbc(t,x) M
// Hp_isf(t,x) M

// * Plasma Flow input variables (Conc of substrate coming into pipe at x=x.min):
// HCO3m_rbc_pl(t) M;  // HCO3- coming into the plasma region (x=x.min)
// Hp_rbc_pl(t) M;     // H+ coming into the plasma region (x=x.min)
// pCO2_pl_Fin(t) mmHg; // pCO2 coming into the plasma region (x.min), CO2_pl_Fin calculated from this.
// pO2_pl_Fin(t) mmHg;  // pO2 coming into the plasma region (x.min), O2_pl_Fin calculated from this.
                        
// pH_pl_Fin(t) dimensionless;  // pH of plasma coming into region at x.min
// * Temp Cap (blood) Variables :
// TempSys_cap(t,x) K;    Temperature in capillary

// * Temp Cap (blood) Input variables :
// TempSys_isf_in(t,x) K; Temp in isf

// * Temp Cap (blood) Flow input variables:
// TempSys_capFin(t) K; Temp 'flow' into capillary (blood) region at x=x.min.

//%END BLOOD_MODULE_INFO

import nsrunit;
unit conversion on;
 unit celsius = fundamental;

math Blood_pde_MPC {

 realDomain t sec; t.min=0; t.max=300; t.delta=1;   // time domain
 realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid;  
//%START SYS_PARAMS
// System parameters:   *******************************
real  Rvel = 1.6 dimensionless; // velocity ratio (vrbc/vpl) (unitless)
real  H = 0.45 dimensionless;   // large vessel hematocrit (unitless)
real  Hcap = H/(H+(1-H)*Rvel);  // capillary hematocrit (unitless)
real  Fb = 1.0 ml/(min*g);	// total blood flow  (blood flow per gram tissue)
real  Fpl = (1-H)*Fb;		// flow of plasma, ml/(min*g)
real  Frbc = H*Fb;		// flow of RBCs, ml/(min*g)
real  Vcap = 0.07 ml/g;		// anatomical volume of CAP
real  Vpl = (1-Hcap)*Vcap;	// anatomical volume of plasma
real  Vrbc = Hcap*Vcap;		// anatomical volume of RBC
//%END SYS_PARAMS
 
real TempSys_isf(t,x) = 310;
//%START BLOOD_REGION
real MolV  = 22.414 L/mol;     // Liters of gas per mole at std temp.
real specificHeat = 1 kilocal/(g*K*min);   // Heat required to raise mass one K per min
real ThermCoeff = 0.0001 g/sec;  // Thermal coefficient
real HCO3m_rbc(t,x) M;   // HCO3m in RBCs
real Hp_rbc(t,x) M;      // H+ conc in RBCs
real HCO3m_isf(t,x) M;    // HCO3m in isf
real Hp_isf(t,x) M;       // H+ conc in isf
real CO2_rbc(t,x) M;     // Free CO2 from RBCs
real CO2_isf(t,x) M;      // CO2 in isf
real O2_rbc(t,x) M;     // Free O2 in RBCs
real O2_isf(t,x) M;     // Free O2 in plasma
// ***** *********** *****************************
// ***** Plasma region: *******************************
// *****   *****     ***************************
// Plasma region:
real CFpl = 100;		// catalytic factor in plasma due to CA
real BCpl = 6 mM;		// buffering capacity in plasma
real kp1    = 0.12 1/sec,	// forward rate constant in CO2+H2O reaction
     km1    = 89 1/sec;         // backward rate constant in CO2+H2O reaction
real KH2CO3 = 5.5e-4 M;	        // EQUIL constant in H2CO3 ionization
real K1_calc = (kp1/km1)*KH2CO3; // EQUIL constant in overall CO2+H2O reaction
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
real CO2Dpl = 1e-4 cm^2/sec;   // diffusion coefficient for CO2 in plasma 
real CO2PSrbc ml/(g*min);       // PS for CO2 across RBC Membrane
real CO2PScap ml/(g*min);       // PS for CO2 across Capillary Membrane
real O2Dpl = 1e-4 cm^2/sec;    // diffusion coefficient for O2 in RBCs 
real O2PSrbc ml/(g*min);       // PS for O2 across RBC Membrane
real O2PScap ml/(g*min);       // PS for O2 across Capillary Membrane
real HCO3m_pl_Fin(t) M;   // HCO3m coming into the pl region (x=x.min)
real Hp_pl_Fin(t) M;      // H+ coming into the pl region (x=x.min)
real CO2_pl_Fin(t) M;   // CO2 coming into the plasma region (x.min)
real O2_pl_Fin(t) M;   // O2 coming into the plasma region (x.min)
real Hp_pl(t,x) M;       // H+ conc in plasma
real pH_pl(t,x) dimensionless; // pH in plasma
     pH_pl = -log(Hp_pl/(1 M));
real HCO3m_pl(t,x) M;    // HCO3m in plasma
real CO2_pl(t,x) M;      // CO2 in plasma
real O2_pl(t,x) M;      // O2 in plasma
real HCO3m_plOut(t) M;  // HCO3m from HCO3plasma_MPC
real Hp_plOut(t) M;     // H+ from HCO3plasma_MPC
real CO2_plOut(t) M;    // CO2 coming out of plasma region at x.max
real O2_plOut(t) M;
// Inputs needed for HCO3plasma_MPC
real HpCplt0 M;
real pH_plt0 dimensionless;   // Initial pH in pl
real CO2Cplt0 M;             // Initial free CO2 conc in pl
real HCO3mCplt0 M;           // Initial HCO3m conc in pl
     HpCplt0 = 10^(-pH_plt0)* (1 M);
     HCO3mCplt0 = K1_calc*CO2Cplt0/HpCplt0;
real CO2Cplt0 M;
real O2Cplt0 M;
 when (t=t.min) {	// Plasma INITIAL CONDITIONS
      HCO3m_pl = HCO3mCplt0;
      Hp_pl = HpCplt0;
      CO2_pl = CO2Cplt0;
      O2_pl = O2Cplt0;  // O2 PDE INITIAL CONDITION
 } // end Plasma ICs
 when (x=x.min) {	// LEFT BOUNDARY CONDITIONS
      (-Fpl*L/Vpl)*(HCO3m_pl-HCO3m_pl_Fin)+HCO3mDpl*HCO3m_pl:x =0;
      (-Fpl*L/Vpl)*(Hp_pl-Hp_pl_Fin)+HpDpl*Hp_pl:x =0;
     (-Fpl*L/Vpl)*(CO2_pl-CO2_pl_Fin)+CO2Dpl*CO2_pl:x =0;
     (-Fpl*L/Vpl)*(O2_pl-O2_pl_Fin)+O2Dpl*O2_pl:x =0;
 } // end Plasma left BCs
 when (x=x.max) {	// RIGHT BOUNDARY CONDITIONS
      HCO3mDpl*HCO3m_pl:x = 0;
      HpDpl*Hp_pl:x = 0;
      HCO3m_plOut = HCO3m_pl;
      Hp_plOut = Hp_pl; 
     CO2Dpl*CO2_pl:x = 0;
     CO2_plOut = CO2_pl;
     O2Dpl*O2_pl:x = 0;
     O2_plOut = O2_pl;
 } // end Plasma right BCs
// PDEs:
 CO2_pl:t = - (Fpl/Vpl)*L*(CO2_pl:x) 
                   - (CO2PSrbc/VWpl)*(CO2_pl-CO2_rbc)
                   - (CO2PScap/VWpl)*(CO2_pl-CO2_isf) 
                   + CO2Dpl*(CO2_pl:x:x) 
                   - CFpl*(kp1*CO2_pl-(km1/KH2CO3)*HCO3m_pl*Hp_pl);
 HCO3m_pl:t = - (Fpl/Vpl)*L*(HCO3m_pl:x) + HCO3mDpl*(HCO3m_pl:x:x)
               + (HCO3mPSrbc/VWpl)*(HCO3m_rbc-Rrbc*HCO3m_pl)
               - (HCO3mPScap/VWpl)*(Rcap*HCO3m_pl-HCO3m_isf)
               + CFpl*(kp1*CO2_pl-(km1/KH2CO3)*HCO3m_pl*Hp_pl);
 Hp_pl:t = - (Fpl/Vpl)*L*(Hp_pl:x) + HpDpl*(Hp_pl:x:x) 
            + (HpPSrbc/VWpl)*(Rrbc*Hp_rbc-Hp_pl)
            - (HpPScap/VWpl)*(Hp_pl-Rcap*Hp_isf)
            + (2.303/BCpl)*Hp_pl * CFpl*(kp1*CO2_pl-(km1/KH2CO3)*HCO3m_pl*Hp_pl);
// O2 PDE:
 O2_pl:t = - (Fpl/Vpl)*L*(O2_pl:x) 
           + (O2PSrbc/VWpl)*(O2_rbc-O2_pl)
           - (O2PScap/VWpl)*(O2_pl-O2_isf)
           + O2Dpl*(O2_pl:x:x);             
// ------ End of PLASMA pdes
// ****** ********************************
// ****** RBC region:  ********************************
// ***** **** ****** *****************************
real CFrbc = 13000;		// catalytic factor in RBC due to CA
real BCrbc = 54 mM;		// buffering capacity in RBC
real Frbc ml/(min*g);           // Flow of RBCs.
real Vrbc ml/g;                 // Volume of region
real Wrbc dimensionless;	// fractional water content of RBCs
real VWrbc ml/g;	        // volume of water content in RBCs
real HCO3mDrbc = 1e-4 cm^2/sec; // diffusion coefficient for HCO3- in RBCs
real HCO3mPSrbc ml/(min*g);     // PS for HCO3- across RBC Membrane
real HpDrbc = 1e-4 cm^2/sec;    // diffusion coefficient for H+ in RBCs
real HpPSrbc ml/(min*g);        // PS for H+ across RBC Membrane
real Rrbc dimensionless;               // Gibbs-Donnan ratio [H+]pl/[H+]rbc
real CO2Drbc = 1e-4 cm^2/sec;   // diffusion coefficient for CO2 in RBCs 
real O2Drbc = 1e-4 cm^2/sec;   // diffusion coefficient for O2 in RBCs 
real HCO3m_rbc_Fin(t) M;   // HCO3- coming into the rbc region (x=x.min)
real Hp_rbc_Fin(t) M;      // H+ coming into the rbc region (x=x.min)
real CO2_pl(t,x) M;      // Free CO2 in plasma
real TCO2_rbc_Fin(t) M;   // TCO2 coming into the rbc region (x.min)
real O2_pl(t,x) M;        // Free O2 in plasma
real TO2_rbc_Fin(t) M;   // TO2 coming into the rbc region (x.min)
real TCO2_rbcOut(t) M;
real TO2_rbcOut(t) M;
real TCO2Crbct0 M;
real TO2Crbct0 M;
real HCO3m_rbcOut(t) M;  // HCO3- leaving RBCs at x=x.max
real Hp_rbcOut(t) M;     // H+ leaving RBCs at x=x.max
// ** RBC Variables: **
real TCO2_rbc(t,x) M;    // Total CO2 in rbc
real TO2_rbc(t,x) M;     // Total O2 in rbc
real pCO2_rbc(t,x) mmHg; // free CO2 in rbc
real pO2_rbc(t,x) mmHg;  // free O2 in rbc
real pH_rbc(t,x) dimensionless; // pH in rbc
     pH_rbc = -log(Hp_rbc/(1 M));
// Inputs needed for HCO3PDE_MPC
real HpCrbct0 M;
real pH_rbct0 dimensionless;   // Initial pH in rbc
real CO2Crbct0 M;             // Initial free CO2 conc in rbc
real HCO3mCrbct0 M;           //  Initial HCO3m conc in rbc
     HpCrbct0 = 10^(-pH_rbct0)* (1 M);
     HCO3mCrbct0 = K1_calc*CO2Crbct0/HpCrbct0;
// ------------------------------------------------------
// PDE ICs:
 when (t=t.min) { // PDE RBC ICs
      TO2_rbc = TO2Crbct0;  // TO2 PDE INITIAL CONDITION
      TCO2_rbc = TCO2Crbct0;
      HCO3m_rbc = HCO3mCrbct0;
      Hp_rbc = HpCrbct0;
 } // End RBC ICs
// PDE BCs:
 when (x=x.min) { //  Left RBC BCs
     (-Frbc*L/Vrbc)*(TO2_rbc-TO2_rbc_Fin)+O2Drbc*TO2_rbc:x =0;
     (-Frbc*L/Vrbc)*(TCO2_rbc-TCO2_rbc_Fin)+CO2Drbc*TCO2_rbc:x =0;
      (-Frbc*L/Vrbc)*(HCO3m_rbc-HCO3m_rbc_Fin)+HCO3mDrbc*HCO3m_rbc:x =0;
      (-Frbc*L/Vrbc)*(Hp_rbc-Hp_rbc_Fin)+HpDrbc*Hp_rbc:x =0;
} // end Left RBC BCs
 when (x=x.max) { //  Right RBC BCs
     O2Drbc*TO2_rbc:x = 0;
     TO2_rbcOut = TO2_rbc;
     CO2Drbc*TCO2_rbc:x = 0;
     TCO2_rbcOut = TCO2_rbc;
      HCO3mDrbc*HCO3m_rbc:x = 0;
      HpDrbc*Hp_rbc:x = 0;
      HCO3m_rbcOut = HCO3m_rbc;
      Hp_rbcOut = Hp_rbc; 
} // end Right RBC BCs
 TCO2_rbc:t = - (Frbc/Vrbc)*L*(TCO2_rbc:x) 
                   - (CO2PSrbc/VWrbc)*(CO2_rbc-CO2_pl)
                   + CO2Drbc*(TCO2_rbc:x:x) 
                   - CFrbc*(kp1*CO2_rbc-(km1/KH2CO3)*HCO3m_rbc*Hp_rbc);
// TO2 PDE:
 TO2_rbc:t = - (Frbc/Vrbc)*L*(TO2_rbc:x) 
                   - (O2PSrbc/VWrbc)*(O2_rbc-O2_pl)
                   + O2Drbc*(TO2_rbc:x:x);             
 HCO3m_rbc:t = - (Frbc/Vrbc)*L*(HCO3m_rbc:x) + HCO3mDrbc*(HCO3m_rbc:x:x)
               - (HCO3mPSrbc/VWrbc)*(HCO3m_rbc-Rrbc*HCO3m_pl)
               + CFrbc*(kp1*CO2_rbc-(km1/KH2CO3)*HCO3m_rbc*Hp_rbc);
 Hp_rbc:t = - (Frbc/Vrbc)*L*(Hp_rbc:x) + HpDrbc*(Hp_rbc:x:x) 
            - (HpPSrbc/VWrbc)*(Rrbc*Hp_rbc-Hp_pl)
            + (2.303/BCrbc)*Hp_rbc * CFrbc*(kp1*CO2_rbc-(km1/KH2CO3)*HCO3m_rbc*Hp_rbc);
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
//%START O2CO2SOL_COEFF_PARAMS
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
//%END O2CO2SOL_COEFF_PARAMS
//%START O2CO2SOL_VARCALCS
real O2solSys(t,x)	mL/mL/atm, // ml O2 per ml plasma per atm
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
//%END O2CO2SOL_VARCALCS
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
// ***  CO2 invert to get pCO2:
real TotCO2_Invert_in(t,x) M;  // Tot CO2 used for inversion to get pCO2  
real pO2_CO2_Invert_in(t,x) mmHg;  // pO2 used for inversion to get pCO2
real pCO2initguess =40 mmHg;   // Init guess.
real pCO2new_CO2_Invert_0(t,x) = pCO2initguess;  
// TotCO2_FreeCO2_Inversion: calcs pCO2 given TotCO2, pO2, alphaCO2, 
real pCO2_out(t,x) mmHg;  // free CO2 in rbc
real FCO2_out(t,x) M;     // free CO2 in rbc
real SHbCO2_out(t,x);     // Value from TCO2 inversion calcs
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
 pCO2_CO2_Invert_1c=pCO2c_CO2_Invert_1;
 pCO2_CO2_Invert_1p=pCO2p_CO2_Invert_1;
 pCO2_CO2_Invert_1m=pCO2m_CO2_Invert_1;
 pO2_CO2_Invert_1= pO2_CO2_Invert_in;
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
real delpCO2_CO2_Invert_1c(t,x) = pCO2_CO2_Invert_1c-pCO20; 
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
// Hill coefficient; unitless (redefined as a function pO2_CO2_Invert_1)
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
// -----------------------------------
// Relate module output variables to current variables:
real SHbO2Sys(t,x);  
 SHbO2Sys = SHbO2_out;    // get value from TotO2freeO2Invert
real SHbCO2Sys(t,x);  
 SHbCO2Sys = SHbCO2_out;    // get value from TotCO2freeCO2Invert
// ****** Initial values: *******
real SHbO2_rbc0 = 0.970;   // initial SHbO2 
real SHbCO2_rbc0 = 0.030;   // initial SHbCO2 
real Hb_rbc0 = HbBl/Hct;    // initial Hb conc in RBC 
pH_rbct0 = 7.24;    // pH_rbc at t=t.min
real pO2rbct0 = 90 mmHg;   // pO2 in rbc at t=0 for all x
real pCO2rbct0 = 50 mmHg;  // pCO2 in rbc at t=0 for all x
real alphaCO20 = 2.8472E-5 M/mmHg;
real alphaO20 = 1.1148E-6 M/mmHg;
TCO2Crbct0 =alphaCO20*pCO2rbct0 +4*Hb_rbc0*SHbCO2_rbc0;   // TCO2 in rbc at t=0 for all x
TO2Crbct0 = alphaO20*pO2rbct0 +4*Hb_rbc0*SHbO2_rbc0;   // TO2 in rbc at t=0 for all x
CO2Crbct0 = alphaCO20*pCO2rbct0;      // CO2 in rbc at t=0 to calc HCO3mCrbct0
P50_in = 26;     // mmHg, init guess for pO2 (p50 is used)
real TempSys_capFin(t) K;
real TempSys_cap(t,x) K;      //  Temperature in cap
// Sub-module assignments:
TotCO2_Invert_in = TCO2_rbc;   // Used by TotCO2freeCO2Invert.mpc
TotO2_Invert_in = TO2_rbc;     // Used by TotO2freeO2Invert.mpc
pCO2_rbc = pCO2_out; // from TotCO2freeCO2Invert.mpc
CO2_rbc = FCO2_out;  // from TotCO2freeCO2Invert.mpc
pO2_rbc = pO2_out;   // from TotO2freeO2Invert.mpc
O2_rbc = FO2_out;    // from TotO2freeO2Invert.mpc
pO2_CO2_Invert_in = pO2_rbc;  // pO2 used for CO2 inversion is from O2 inversion. 
pCO2_Invert_in = pCO2_rbc;     // pCO2 used for O2 inversion is from CO2 inversion. 
pHrbc_in = -log(Hp_rbc/(1 M)); // Used in O2/CO2 inversion
DPGrbc_in =.00465;   // DPG conc, used in O2/CO2 inversion
Hbrbc_in = .0052;    // NOT needed, find declaration and remove it.
pCO2t0_in= pCO2rbct0; // pCO2 init conc, used in O2/CO2 inversion
    O2_pl_in = O2_pl;
    Hct_in =0.45;     // Hemotocrit used by model and assigned to Hct
real TempExp K;       // Initial temperature
real TDcap(t,x) s^-1;      // Temp diffusion across cap (capillary) membrane
real TCcap cm^2/sec;       // Thermal conductivity (Diffusion) in cap
     TCcap =  ThermCoeff*Vcap/L;
real TempSys_capOut(t) K;
    pH_plt0 = 7.24;   // Input initial condition
real pO2plt0 = 90 mmHg;
real pCO2plt0 = 40 mmHg;
    O2Cplt0 =  pO2plt0 * alphaO20;
    CO2Cplt0 = pCO2plt0 * alphaCO20; 
real pO2_pl(t,x) mmHg;
    pO2_pl = O2_pl/alphaO2Sys;
real pCO2_pl(t,x) mmHg;
    pCO2_pl = CO2_pl/alphaCO2Sys;
 when (t=t.min) {	// Temp CAP PDE INITIAL CONDITIONS
  TempSys_cap = TempExp;
 }  // END TEMP CAPPILLARY PDE IC
 when (x=x.min) {	// Temp CAP PDE LEFT BOUNDARY CONDITIONS
  (-Fb*L/Vcap)*(TempSys_cap-TempSys_capFin(t)) +TCcap*(TempSys_cap):x = 0; 
 }   //  END Temp Capillary PDE BC
 when (x=x.max) {	// Temp CAP PDE RIGHT BOUNDARY CONDITIONS
  TempSys_cap:x =0;   
  TempSys_capOut = TempSys_cap;
 } // END Temp CAPPILLARY PDE RBC
TempSys_cap:t = -(Fb*L/Vcap)*(TempSys_cap):x 
                +TDcap*(TempSys_isf-TempSys_cap) 
                + TCcap*(TempSys_cap):x:x ; // tissue blood
//%END BLOOD_REGION
// Flow into RBC and plasma region at x=x.min:
//%START RBC_FLOW_INPUTS
// ************** RBC Input at x.min, change as needed: **************************
real pO2_rbc_Fin(t) = 90 mmHg;  // pO2 coming into RBC region at x.min
real pCO2_rbc_Fin(t) = 50 mmHg;  // pCO2 coming into RBC region at x.min
real O2_rbc_Fin(t) = alphaO2Sys(t,x.min)*pO2_rbc_Fin; // free O2 coming into RBC region at x.min needed? 
real CO2_rbc_Fin(t) = alphaCO2Sys(t,x.min)*pCO2_rbc_Fin; // free CO2 coming into RBC region at x.min needed?
TCO2_rbc_Fin = CO2_rbc_Fin+4*Hb_rbc0*SHbCO2_rbc0;  // TCO2 coming into RBC at x.min
TO2_rbc_Fin = O2_rbc_Fin+4*Hb_rbc0*SHbO2_rbc0;     // TO2 coming into RBC at x.min
HCO3m_rbc_Fin = K1*CO2_rbc_Fin/Hp_rbc_Fin;  // HCO3m coming into RBC at x.min
Hp_rbc_Fin = HpCrbct0;   // Hp coming into RBC at x.min
//%END RBC_FLOW_INPUTS
//%START PLASMA_FLOW_INPUTS
// ** Plasma Input at x.min and other regions, change as needed:
real pO2_pl_Fin(t) = 90 mmHg;  // pO2 coming into plasma region at x.min
real pCO2_pl_Fin(t) = 50 mmHg;  // pCO2 coming into plasma region at x.min
real pH_pl_Fin(t) = 7.24;  // pH of plasma coming into region at x.min
    O2_pl_Fin = pO2_pl_Fin*alphaO2Sys(t,x.min);
    CO2_pl_Fin = pCO2_pl_Fin*alphaCO2Sys(t,x.min);
    HCO3m_pl_Fin = K1*CO2_pl_Fin/Hp_pl_Fin;
    Hp_pl_Fin = 10^(-pH_pl_Fin)* (1 M);
//%END PLASMA_FLOW_INPUTS
//%START PS_VALUES
    O2PSrbc = 0;
    O2PScap = 0;
    CO2PSrbc = 0;
    CO2PScap = 0;
    HCO3mPSrbc =0;
    HpPSrbc = 0;
    HCO3mPScap =0;
    HpPScap = 0;
//%END PS_VALUES
//%START ISF_VALUES 
// No iteraction with ISF:
  O2_isf = 0.0;    
  CO2_isf = 0.0;  
  Hp_isf = 0;
  HCO3m_isf = 0;
//%END ISF_VALUES 
    VWrbc = Vrbc*Wrbc;
    VWpl = Vpl*Wpl;           // Water fraction of plasma volume.
    Rcap = 0.63;              // Gibbs-Donnan ratio [H+]pl/[H+]isf
    TDcap =0;                 // Temp diffusion across cap (capillary) membrane
    TempSys_capFin = 310;
    TempExp = 310;
    Temp_in =TempSys_cap;     // Make sure Temp input is using TempSys_cap
} // End of Blood pde model
// This MML file generated from Bloodpde_MPC.mpc using MPC v1.01.
  
