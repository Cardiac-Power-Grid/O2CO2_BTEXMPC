// MODEL NAME: Temp_pc_MPC
// SHORT DESCRIPTION: Modeling Temp cahnge in PC as a function of t and x in pc region.


import nsrunit;   unit conversion on;

math Temp_pc_MPC {

//%START TIME_DOMAIN
realDomain t sec;t.min=0;t.max=10;t.delta=0.1;
//%END TIME_DOMAIN

//%START SPATIALDOMAIN
realDomain x cm;real L=0.1 cm; int Ngrid=51; x.min=0; x.max=L; x.ct=Ngrid; // BTEX space domain 
//%END SPATIALDOMAIN

//%START PARAMS
//real Fpc ml/(min*g);           // Flow of pc. Assume zero.
real Vpc ml/g;                 // Volume of region
//%END PARAMS

//%START THERM_STD_PARAMS
real MolV  = 22.414 L/mol;     // Liters of gas per mole at std temp.
real specificHeat = 1 kilocal/(g*K*min);   // Heat required to raise mass one K per min
real ThermCoeff = 0.0001 g/sec;  // Thermal coefficient
//%END THERM_STD_PARAMS

//%START THERM_PARAMETERS
real TDpc(t,x) s^-1;      // Temp diffusion across pc (muscle) membrane
real TCpc cm^2/sec;       // Thermal conductivity (Diffusion) in pc
     TCpc =  ThermCoeff*Vpc/L;
//%END THERM_PARAMETERS

//%START INPUT
real MRO2pc_in(t,x) mol/min/g;	// O2 consumption rate in PC
real TempSys_isf_in(t,x) K;
//%END INPUT

//%START HEAT_GAIN
real O2energy_cal kilocal/mol/sec;         // Energy released per mole O2 consumed per sec.
real O2energy_calperliter = 1000 kilocal/L/sec; // Energy released per liter of O2 consumed per sec
     O2energy_cal = O2energy_calperliter*MolV;
real HeatGain(t,x) kilocal/min/g/sec;
     HeatGain = O2energy_cal*MRO2pc_in;
real TempGain(t,x) K/sec;
     TempGain = HeatGain/specificHeat;

//%END HEAT_GAIN

//%START VARIABLES
real TempSys_pc(t,x) K;      //  Temperature in pc
//%END VARIABLES

//%START VAR_OUTPUT

//%END VAR_OUTPUT

//%START VAR_INITS
real TempExp K;       // Initial temperature
//%END VAR_INITS


 when (t=t.min) {	// Temp PC PDE INITIAL CONDITIONS
//%START TEMP_PC_PDE_IC
  TempSys_pc = TempExp;
//%END TEMP_PC_PDE_IC
 } 

 when (x=x.min) {	// Temp PC PDE LEFT BOUNDARY CONDITIONS
//%START TEMP_PC_PDE_LBC
  TCpc*TempSys_pc:x = 0;  
//%END TEMP_PC_PDE_LBC
 }   //  END Temp PDE BC

 when (x=x.max) {	// Temp PC PDE RIGHT BOUNDARY CONDITIONS
//%START TEMP_PC_PDE_RBC
  TempSys_pc:x =0;   
//%END TEMP_PC_PDE_RBC
 }  

//%START TEMP_PC_PDES
// PDEs:
TempSys_pc:t = -TDpc*(TempSys_pc-TempSys_isf_in) 
               +(O2energy_cal/specificHeat)*MRO2pc_in  // Temp gain
               + TCpc*(TempSys_pc):x:x ; //pc
//%END TEMP_PC_PDES

// To make stand-alone model:
    TempExp = 310;
    TempSys_isf_in = 310;
    TDpc =0;
    Vpc = 0.01;      // ml/g
    MRO2pc_in = 1;


}
