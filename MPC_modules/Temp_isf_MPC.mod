// MODEL NAME: Temp_isf_MPC
// SHORT DESCRIPTION: Modeling Temp cahnge in isf as a function of t and x in isf region.


import nsrunit;   unit conversion on;

math Temp_isf_MPC {

//%START TIME_DOMAIN
realDomain t sec;t.min=0;t.max=10;t.delta=0.1;
//%END TIME_DOMAIN

//%START SPATIALDOMAIN
realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid; // BTEX space domain 
//%END SPATIALDOMAIN

//%START PARAMS
//real Fisf ml/(min*g);           // Flow of isf. Assume zero.
real Visf ml/g;                 // Volume of region
//%END PARAMS


//%START THERM_STD_PARAMS
real MolV  = 22.414 L/mol;     // Liters of gas per mole at std temp.
real specificHeat = 1 kilocal/(g*K*min);   // Heat required to raise mass one K per min
real ThermCoeff = 0.0001 g/sec;  // Thermal coefficient
//%END THERM_STD_PARAMS

//%START THERM_PARAMETERS
real TDpc(t,x) s^-1;      // Temp diffusion across pc (muscle) membrane
real TDcap(t,x) s^-1;      // Temp diffusion across cap (capillary) membrane
real TCisf cm^2/sec;       // Thermal conductivity (Diffusion) in isf
     TCisf =  ThermCoeff*Visf/L;
//%END THERM_PARAMETERS

//%START INPUTS
real TempSys_cap_in(t,x) K;
real TempSys_pc_in(t,x) K;
//%END INPUTS

//%START HEAT_GAIN

//%END HEAT_GAIN

//%START VARIABLES
real TempSys_isf(t,x) K;      //  Temperature in isf
//%END VARIABLES

//%START VAR_OUTPUT

//%END VAR_OUTPUT

//%START VAR_INITS
real TempExp K;       // Initial temperature
//%END VAR_INITS


 when (t=t.min) {	// Temp isf PDE INITIAL CONDITIONS
//%START TEMP_ISF_PDE_IC
  TempSys_isf = TempExp;
//%END TEMP_ISF_PDE_IC
 } 

 when (x=x.min) {	// Temp isf PDE LEFT BOUNDARY CONDITIONS
//%START TEMP_ISF_PDE_LBC
  TCisf*TempSys_isf:x = 0;  
//%END TEMP_ISF_PDE_LBC
 }   //  END Temp PDE BC

 when (x=x.max) {	// Temp isf PDE RIGHT BOUNDARY CONDITIONS
//%START TEMP_ISF_PDE_RBC
  TempSys_isf:x =0;   
//%END TEMP_ISF_PDE_RBC
 }  

//%START TEMP_ISF_PDES
// PDEs:
TempSys_isf:t =  TDpc*(TempSys_pc_in-TempSys_isf) 
                - TDcap*(TempSys_isf-TempSys_cap_in) 
                + TCisf*(TempSys_isf):x:x ; 
//%END TEMP_ISF_PDES

// To make stand-alone model:
    TempExp = 310;
    TempSys_cap_in = 310;
    TempSys_pc_in = 310;
    TDcap =0;
    TDpc = 0;
    Visf = 0.1;      // ml/g
 


}
