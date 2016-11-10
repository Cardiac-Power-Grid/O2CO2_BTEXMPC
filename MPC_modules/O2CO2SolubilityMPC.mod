// MODEL NAME: O2CO2SolubilityMPC 
// SHORT DESCRIPTION: Calculates O2 and CO2 solubilities based on temperature. Time and one 
// spatial dimension dependent. MPC annotated module.

import nsrunit;   unit conversion on;

//%START CELSIUS_DEF
 unit celsius = fundamental;
//%END CELSIUS_DEF

math O2CO2SolubilityMPC {

//%START TIME_DOMAIN
 realDomain t sec; t.min=0; t.max=300; t.delta=1;   // time domain
//%END TIME_DOMAIN

//%START SPATIAL_DOMAIN
 realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid;  
//%END SPATIAL_DOMAIN

//%START CONSTS
real Rconst = 62.36358 L*mmHg/K/mole;  // ideal gas const
//%END CONSTS


//%START INPUTS
real TempSys(t,x) K;
real alphaO2Sys(t,x) M/mmHg;
real alphaCO2Sys(t,x) M/mmHg;
//%END INPUTS

//%START O2CO2SOL_CALCS
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
     TempC_Sys = (TempSys-273.15)* (1 celsius)/(1 K);  // convert kelvin to celsius 

     O2solSys = alphaO20calc*exp(-O2k1*TempC_Sys)+alphaO201calc*exp(-O2k2*TempC_Sys);  
     CO2solSys = alphaCO20calc*exp(-CO2k1*TempC_Sys)+alphaCO201calc*exp(-CO2k2*TempC_Sys);

// convert solubilities to M/mmHg;
     alphaO2Sys =  (1 atm)*O2solSys/Rconst/TempSys;
     alphaCO2Sys = (1 atm)* CO2solSys/Rconst/TempSys;

// **** END of solubility calcs........
//%END O2CO2SOL_VARCALCS
//%END O2CO2SOL_CALCS

// ---------------------------------------------------------
//%START CHECK_SOLUBILITIES
// Check solubility curves:
realDomain Temp celsius; Temp.min=0; Temp.max=60; Temp.delta=1;  // temp domain
real O2solChk(Temp)	mL/mL/atm, // ml O2 per ml plasma per atm
     CO2solChk(Temp)	mL/mL/atm;
O2solChk = alphaO20calc*exp(-O2k1*Temp)+alphaO201calc*exp(-O2k2*Temp);  
CO2solChk = alphaCO20calc*exp(-CO2k1*Temp)+alphaCO201calc*exp(-CO2k2*Temp);
// convert to M/mmHg;
real alphaO2chk(Temp) M/mmHg;
     alphaO2chk = (1 atm)*O2solChk/Rconst/(Temp+273.15 )* (1 celsius)/(1 K);
real alphaCO2chk(Temp) = (1 atm)* CO2solChk/Rconst/(Temp+273.15 )* (1 celsius)/(1 K);
//%END CHECK_SOLUBILITIES

// Needed for stand-alone model:
TempSys = 315;

}

