// MODEL NAME: Temp_cap_MPC
// SHORT DESCRIPTION: Modeling Temp cahnge in cap as a function of t and x in cap region.


import nsrunit;   unit conversion on;

//%START TEMPERATURE_CAPILLARY_MOD_INFO
// * Temp Cap (blood) Variables :
// TempSys_cap(t,x) K;    Temperature in capillary

// * Temp Cap (blood) Input variables :
// TempSys_isf_in(t,x) K; Temp in isf

// * Temp Cap (blood) Flow input variables:
// TempSys_capFin(t) K; Temp 'flow' into capillary (blood) region at x=x.min.

//%END TEMPERATURE_CAPILLARY_MOD_INFO

math Temp_cap_MPC {

//%START TIME_DOMAIN
realDomain t sec;t.min=0;t.max=10;t.delta=0.1;
//%END TIME_DOMAIN

//%START SPATIALDOMAIN
realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid; // BTEX space domain 
//%END SPATIALDOMAIN

//%START PARAMS
real Fb ml/(min*g);           // Blood Flow in capillary.
real Vcap ml/g;                 // Volume of region
//%END PARAMS

//%START THERM_STD_PARAMS
real MolV  = 22.414 L/mol;     // Liters of gas per mole at std temp.
real specificHeat = 1 kilocal/(g*K*min);   // Heat required to raise mass one K per min
real ThermCoeff = 0.0001 g/sec;  // Thermal coefficient
//%END THERM_STD_PARAMS

//%START THERM_PARAMETERS
real TDcap(t,x) s^-1;      // Temp diffusion across cap (capillary) membrane
real TCcap cm^2/sec;       // Thermal conductivity (Diffusion) in cap
     TCcap =  ThermCoeff*Vcap/L;
//%END THERM_PARAMETERS

//%START INPUT
real TempSys_isf_in(t,x) K;
//%END INPUT

//%START CAP_FLOW_INPUT
real TempSys_capFin(t) K;
//%END CAP_FLOW_INPUT

//%START HEAT_GAIN

//%END HEAT_GAIN

//%START VARIABLES
real TempSys_cap(t,x) K;      //  Temperature in cap
//%END VARIABLES

//%START VAR_OUTPUT
real TempSys_capOut(t) K;
//%END VAR_OUTPUT

//%START VAR_INITS
real TempExp K;       // Initial temperature
//%END VAR_INITS


 when (t=t.min) {	// Temp CAP PDE INITIAL CONDITIONS
//%START TEMP_CAP_PDE_IC
  TempSys_cap = TempExp;
//%END TEMP_CAP_PDE_IC
 } 

 when (x=x.min) {	// Temp CAP PDE LEFT BOUNDARY CONDITIONS
//%START TEMP_CAP_PDE_LBC
  (-Fb*L/Vcap)*(TempSys_cap-TempSys_capFin(t)) +TCcap*(TempSys_cap):x = 0; 
//%END TEMP_CAP_PDE_LBC
 }   //  END Temp PDE BC

 when (x=x.max) {	// Temp CAP PDE RIGHT BOUNDARY CONDITIONS
//%START TEMP_CAP_PDE_RBC
  TempSys_cap:x =0;   
  TempSys_capOut = TempSys_cap;
//%END TEMP_CAP_PDE_RBC
 }  

//%START TEMP_CAP_PDES
// PDEs:
TempSys_cap:t = -(Fb*L/Vcap)*(TempSys_cap):x 
                +TDcap*(TempSys_isf_in-TempSys_cap) 
                + TCcap*(TempSys_cap):x:x ; // tissue blood
//%END TEMP_CAP_PDES

// To make stand-alone model:
    TempExp = 310;
    TempSys_isf_in = 310;
    TempSys_capFin = 310;
    TDcap =0;
    Vcap = 0.1;      // ml/g
    Fb=1;



}

