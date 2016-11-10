// MODEL NAME: PlasmaPDE_MPC
// SHORT DESCRIPTION: Simple model of Plasma region that calculates the pO2 and pCO2
// in RBC based on Total O2, CO2, H+, and HCO3. Takes into account temperature
// dependence of O2/CO2 solubility.



//%START PLASMA_MODULE_INFO
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
                        
//%END PLASMA_MODULE_INFO




import nsrunit;
unit conversion on;

math PLASMA_pde_MPC {
realDomain t sec;t.min=0;t.max=10;t.delta=0.1;
realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid; // BTEX space domain 

//%START PLASMA_INPUT
real HCO3m_rbc(t,x) M;   // HCO3m in RBCs
real Hp_rbc(t,x) M;      // H+ conc in RBCs
real HCO3m_isf(t,x) M;    // HCO3m in isf
real Hp_isf(t,x) M;       // H+ conc in isf
real CO2_rbc(t,x) M;     // Free CO2 from RBCs
real CO2_isf(t,x) M;      // CO2 in isf
real O2_rbc(t,x) M;     // Free O2 in RBCs
real O2_isf(t,x) M;     // Free O2 in plasma
//%END PLASMA_INPUT
//%START PLASMA_REGION
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
//%END PLASMA_REGION
real alphaCO2 = 2.8472E-5 M/mmHg;  // Solubility of CO2 in plasma
real alphaO2 = 1.1148E-6 M/mmHg; // Solubility of O2 in plasma
//%START PLASMA_INIT
    pH_plt0 = 7.24;   // Input initial condition
real pO2plt0 = 90 mmHg;
real pCO2plt0 = 40 mmHg;
    O2Cplt0 =  pO2plt0 * alphaO2;
    CO2Cplt0 = pCO2plt0 * alphaCO2; 
//%END PLASMA_INIT
//%START PO2_PCO2
real pO2_pl(t,x) mmHg;
    pO2_pl = O2_pl/alphaO2;
real pCO2_pl(t,x) mmHg;
    pCO2_pl = CO2_pl/alphaCO2;
//%END PO2_PCO2
  
    O2PSrbc = 0;
    O2PScap = 0;
    CO2PSrbc = 0;
    CO2PScap = 0;
    HCO3mPSrbc =0;
    HpPSrbc = 0;
    HCO3mPScap =0;
    HpPScap = 0;
    
//%START PLASMA_FLOW_INPUTS
// ** Plasma Input at x.min and other regions, change as needed:
real pO2_pl_Fin(t) = 90 mmHg;  // pO2 coming into RBC region at x.min
real pCO2_pl_Fin(t) = 50 mmHg;  // pCO2 coming into RBC region at x.min
    O2_pl_Fin = pO2_pl_Fin*alphaO2;
    CO2_pl_Fin = pCO2_pl_Fin*alphaCO2;
    HCO3m_pl_Fin = .01459;
    Hp_pl_Fin = 10^(-pH_plt0)* (1 M);
//%END PLASMA_FLOW_INPUTS
//%START PLASMA_RBC_INPUTS
// Link module var dependencies with common vars:
    CO2_rbc =0;
    O2_rbc =0;
    HCO3m_rbc =0;   
    Hp_rbc = 0; 
//%END PLASMA_RBC_INPUTS
//%START PLASMA_ASSIGN_ISF_INPUTS 
    CO2_isf =0;
    O2_isf = 0;
    HCO3m_isf = 0;   
    Hp_isf = 0;  
//%END PLASMA_ASSIGN_ISF_INPUTS  
//%START SYS_PARAMS
// System parameters:
real  Rvel = 1.6 dimensionless; // velocity ratio (vrbc/vpl) (unitless)
real  H = 0.45 dimensionless;   // large vessel hematocrit (unitless)
real  Hcap = H/(H+(1-H)*Rvel);  // capillary hematocrit (unitless)
real  Fb = 1.0 ml/(min*g),	// total blood flow  (blood flow per gram tissue)
      Fpl = (1-H)*Fb,		// flow of plasma, ml/(min*g)
      Frbc = H*Fb,		// flow of RBCs, ml/(min*g)
      Vcap = 0.07 ml/g,		// anatomical volume of CAP
      Vpl = (1-Hcap)*Vcap,	// anatomical volume of plasma
      Wpl = 0.72;               // fractional water content of plasma
      VWpl = Vpl*Wpl;
      Rcap = 0.63;
      Rrbc = 0.69;
//%END SYS_PARAMS
} // End of Plasma pde model
// This MML file generated from PlasmaPDE_MPC.mpc using MPC v1.01.
