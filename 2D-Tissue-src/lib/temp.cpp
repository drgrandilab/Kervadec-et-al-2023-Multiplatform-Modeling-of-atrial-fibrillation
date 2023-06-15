// void f( double t, double * y, int p, char * runType, double * ydot , double * currents ){
    
//     //// Model Parameters
//     //// EPI or ENDO?
//     const int epi=1;
//     //// AF
//     const double AF=0;
//     //// ISO
//     const double ISO=0;
//     //// Right ATRIUM
//     const double RA=0;
    
    
//     // Constants
//     const double R = 8314.;       // [J/kmol*K]
//     const double Frdy = 96485.;   // [C/mol]
//     const double Temp = 310.;     // [K]
//     const double FoRT = Frdy / R / Temp;
//     const double Cmem = 1.1e-10;   // [F] membrane capacitance 1.3810e-10;//
//     const double Qpow = ( Temp - 310. ) / 10. ;
    
//     const double pi = 3.141592653589793;
    
//     // Cell geometry
//     const double cellLength = 100.;     // cell length [um]113;//100
//     const double cellRadius = 10.25;   // cell radius [um]12;//10.25
//     const double junctionLength = 160e-3;  //  junc length [um]
//     const double junctionRadius = 15e-3;   //  junc radius [um]
//     const double distSLcyto = 0.45;    //  dist. SL to cytosol [um]
//     const double distJuncSL = 0.5;  //  dist. junc to SL [um]
//     const double DcaJuncSL = 1.64e-6;  //  Dca junc to SL [cm^2/sec]
//     const double DcaSLcyto = 1.22e-6; //  Dca SL to cyto [cm^2/sec]
//     const double DnaJuncSL = 1.09e-5;  //  Dna junc to SL [cm^2/sec]
//     const double DnaSLcyto = 1.79e-5;  //  Dna SL to cyto [cm^2/sec]
//     const double Vcell = pi * cellRadius * cellRadius * cellLength * 1e-15;    //  [L]
//     const double Vmyo = 0.65 * Vcell;
//     const double Vsr = 0.035 * Vcell;
//     const double Vsl = 0.02 * Vcell;
//     const double Vjunc = 1. * 0.0539 * .01 * Vcell;
//     const double SAjunc = 20150. * pi * 2. * junctionLength * junctionRadius;  //  [um^2]
//     const double SAsl = pi * 2. * cellRadius * cellLength;          //  [um^2]
//     // J_ca_juncsl = DcaJuncSL*SAjunc/distSLcyto*1e-10;//  [L/msec] = 1.1074e-13
//     // J_ca_slmyo = DcaSLcyto*SAsl/distJuncSL*1e-10;  //  [L/msec] = 1.5714e-12
//     // J_na_juncsl = DnaJuncSL*SAjunc/distSLcyto*1e-10;//  [L/msec] = 7.36e-13
//     // J_na_slmyo = DnaSLcyto*SAsl/distJuncSL*1e-10;  //  [L/msec] = 2.3056e-11
//     // J_ca_juncsl = DcaJuncSL*SAjunc/distJuncSL*1e-10;//  [L/msec] = 9.9664e-014
//     // J_ca_slmyo = DcaSLcyto*SAsl/distSLcyto*1e-10;  //  [L/msec] = 1.7460e-012
//     // J_na_juncsl = DnaJuncSL*SAjunc/distJuncSL*1e-10;//  [L/msec] = 6.6240e-013
//     // J_na_slmyo = DnaSLcyto*SAsl/distSLcyto*1e-10;  //  [L/msec] = 2.5618e-011
//     //  tau's from c-code, not used here
//     const double J_ca_juncsl = 1. / 1.2134e12; //  [L/msec] = 8.2413e-13
//     const double J_ca_slmyo = 1. / 2.68510e11; //  [L/msec] = 3.2743e-12
//     const double J_na_juncsl = 1. / ( 1.6382e12 / 3. * 100. ); //  [L/msec] = 6.1043e-13
//     const double J_na_slmyo = 1. / ( 1.8308e10 / 3. * 100. );  //  [L/msec] = 5.4621e-11
    
//     //  Fractional currents in compartments
//     const double Fjunc = 0.11;
//     const double Fsl = 1. - Fjunc;
//     const double Fjunc_CaL = 0.9;
//     const double Fsl_CaL = 1. - Fjunc_CaL;
    
//     //  Fixed ion concentrations
//     const double Cli = 15.;   //  Intracellular Cl  [mM]
//     const double Clo = 150.;  //  Extracellular Cl  [mM]
//     const double Ko = 5.4;   //  Extracellular K   [mM]
//     const double Nao = 140.;  //  Extracellular Na  [mM]
//     const double Cao = 1.8;  //  Extracellular Ca  [mM]
//     const double Mgi = 1.;    //  Intracellular Mg  [mM]
    
//     //  Nernst Potentials
//     const double ena_junc = ( 1. / FoRT ) * log( Nao / y[31] );     //  [mV]
//     const double ena_sl = ( 1. / FoRT ) * log( Nao / y[32] );       //  [mV]
//     const double ek = ( 1. / FoRT ) * log( Ko / y[34]);	        //  [mV]
//     const double eca_junc = ( 1. / FoRT / 2. ) * log( Cao / y[35] );   //  [mV]
//     const double eca_sl = ( 1. / FoRT / 2. ) * log( Cao / y[36] );     //  [mV]
//     const double ecl = ( 1. / FoRT ) * log( Cli / Clo );            //  [mV]
    
//     //  Na transport parameters
    
//     // //
//     const double GNa =  para[0] * 23. * ( 1. - 0.1 * AF );  //  [mS/uF]
//     const double GNaB = para[1] * 0.597e-3;    //  [mS/uF]
//     const double IbarNaK = para[2] * 1.26;     //  [uA/uF]
//     const double KmNaip = 11. * ( 1. - 0.25 * ISO );         //  [mM]11
//     const double KmKo = 1.5;         //  [mM]1.5
//     const double Q10NaK = 1.63;
//     const double Q10KmNai = 1.39;
    
//     // //  K current parameters
//     const double pNaK = para[3] *0.01833;
//     const double gkp = para[4] *0.002;
    
//     //  Cl current parameters
//     const double GClCa = para[5] * 0.0548;   //  [mS/uF]
//     const double GClB = para[6] *  9e-3;        //  [mS/uF]
//     const double KdClCa = 100e-3;    //  [mM]
    
//     //  I_Ca parameters
//     const double pNa =  para[7] *( 1. + 0.5 * ISO ) * ( 1. - 0.5 * AF ) * 0.75e-8;       //  [cm/sec]
//     const double pCa =  para[7] *( 1. + 0.5 * ISO ) * ( 1. - 0.5 * AF ) * 2.7e-4;       //  [cm/sec]
//     const double pK = para[7] * ( 1. + 0.5 * ISO ) * ( 1. - 0.5 * AF ) * 1.35e-7;        //  [cm/sec]
//     const double Q10CaL = 1.8;
    
//     // //  Ca transport parameters
//     const double IbarNCX = para[8] *( 1. + 0.4 * AF ) * 3.15;      //  [uA/uF]5.5 before - 9 in rabbit
//     const double KmCai = 3.59e-3;    //  [mM]
//     const double KmCao = 1.3;        //  [mM]
//     const double KmNai = 12.29;      //  [mM]
//     const double KmNao = 87.5;       //  [mM]
//     const double ksat = 0.27;        //  [none]
//     const double nu = 0.35;          //  [none]
//     const double Kdact = 0.384e-3;   //  [mM] 0.256 rabbit
//     const double Q10NCX = 1.57;      //  [none]
//     const double IbarSLCaP =  0.0471; //  IbarSLCaP FEI changed [uA/uF](2.2 umol/L cytosol/sec) jeff 0.093 [uA/uF]
//     const double KmPCa = 0.5e-3;     //  [mM]
//     const double GCaB = 6.0643e-4;    //  [uA/uF] 3
//     const double Q10SLCaP = 2.35;    //  [none]
    
//     //  SR flux parameters
//     const double Q10SRCaP = 2.6;          //  [none]
//     const double Vmax_SRCaP = para[9] * 5.3114e-3;  //  [mM/msec] (286 umol/L cytosol/sec)
//     const double Kmf = ( 2.5 - 1.25 * ISO ) * 0.246e-3;          //  [mM] default
//     const double Kmr = 1.7;               //  [mM]L cytosol
//     const double hillSRCaP = 1.787;       //  [mM]
//     const double ks =  para[10] *25.;                 //  [1/ms]
//     const double koCa = 10. + 20. * AF + 10. * ISO * ( 1. - AF );               //  [mM^-2 1/ms]   // default 10   modified 20
//     const double kom = 0.06;              //  [1/ms]
//     const double kiCa = 0.5;              //  [1/mM/ms]
//     const double kim = 0.005;             //  [1/ms]
//     const double ec50SR = 0.45;           //  [mM]
    
//     //  Buffering parameters
//     //  koff: [1/s] = 1e-3*[1/ms];  kon: [1/uM/s] = [1/mM/ms]
//     const double Bmax_Naj = 7.561;       //  [mM] //  Na buffering
//     const double Bmax_Nasl = 1.65;       //  [mM]
//     const double koff_na = 1e-3;         //  [1/ms]
//     const double kon_na = 0.1e-3;        //  [1/mM/ms]
//     const double Bmax_TnClow = 70e-3;    //  [mM]                      //  TnC low affinity
//     const double koff_tncl = ( 1. + 0.5 * ISO ) * 19.6e-3;    //  [1/ms]
//     const double kon_tncl = 32.7;        //  [1/mM/ms]
//     const double Bmax_TnChigh = 140e-3;  //  [mM]                      //  TnC high affinity
//     const double koff_tnchca = 0.032e-3; //  [1/ms]
//     const double kon_tnchca = 2.37;      //  [1/mM/ms]
//     const double koff_tnchmg = 3.33e-3;  //  [1/ms]
//     const double kon_tnchmg = 3e-3;      //  [1/mM/ms]
//     const double Bmax_CaM = 24e-3;       //  [mM] **? about setting to 0 in c-code**   //  CaM buffering
//     const double koff_cam = 238e-3;      //  [1/ms]
//     const double kon_cam = 34.;           //  [1/mM/ms]
//     const double Bmax_myosin = 140e-3;   //  [mM]                      //  Myosin buffering
//     const double koff_myoca = 0.46e-3;   //  [1/ms]
//     const double kon_myoca = 13.8;       //  [1/mM/ms]
//     const double koff_myomg = 0.057e-3;  //  [1/ms]
//     const double kon_myomg = 0.0157;     //  [1/mM/ms]
//     const double Bmax_SR = 19. * .9e-3;     //  [mM] (Bers text says 47e-3) 19e-3
//     const double koff_sr = 60e-3;        //  [1/ms]
//     const double kon_sr = 100.;           //  [1/mM/ms]
//     const double Bmax_SLlowsl = 37.4e-3 * Vmyo / Vsl;        //  [mM]    //  SL buffering
//     const double Bmax_SLlowj = 4.6e-3 * Vmyo / Vjunc * 0.1;    //  [mM]    // Fei *0.1!!! junction reduction factor
//     const double koff_sll = 1300e-3;     //  [1/ms]
//     const double kon_sll = 100.;          //  [1/mM/ms]
//     const double Bmax_SLhighsl = 13.4e-3 * Vmyo / Vsl;       //  [mM]
//     const double Bmax_SLhighj = 1.65e-3 * Vmyo / Vjunc * 0.1;  //  [mM] // Fei *0.1!!! junction reduction factor
//     const double koff_slh = 30e-3;       //  [1/ms]
//     const double kon_slh = 100.;          //  [1/mM/ms]
//     const double Bmax_Csqn = 140e-3 * Vmyo / Vsr;            //  [mM] //  Bmax_Csqn = 2.6;      //  Csqn buffering
//     const double koff_csqn = 65.;         //  [1/ms]
//     const double kon_csqn = 100.;         //  [1/mM/ms]
    
    
//     // //  Membrane Currents
//     const double mss = 1. / pow( ( 1. + exp( -(56.86 + y[38] ) / 9.03 ) ) , 2 );
//     const double taum = 0.1292 * exp( - pow( ( ( y[38] + 45.79 ) / 15.54 ) , 2 ) ) + 0.06487 * exp( -pow( ( ( y[38] - 4.823 ) / 51.12 ) , 2 ) );
    
//     const double ah = ( ( y[38] >= -40 ) * (0)
//                        + ( y[38] < -40 ) * ( 0.057 * exp( -(  y[38]  + 80 ) / 6.8 ) ) );
//     const double bh = ( ( y[38] >= -40 ) * ( 0.77 / (0.13 * ( 1. + exp( -( y[38] + 10.66 ) / 11.1 ) ) ) )
//                        + (y[38] < -40 ) * ( ( 2.7 * exp( 0.079 * y[38] ) + 3.1 * 1E5 * exp(0.3485 * y[38] ) ) ) );
//     const double tauh = 1. / ( ah + bh );
//     const double hss = 1. / ( pow( ( 1. + exp( ( y[38] + 71.55 ) / 7.43 ) ) , 2 ) );
    
//     const double aj = ( ( y[38]  >= -40 ) * (0)
//                        + ( y[38]  < -40 ) * ( ( ( -2.5428 * 1E4 * exp( 0.2444 * y[38] ) - 6.948 * 1E-6 * exp( -0.04391 * y[38] ) )
//                                                * ( y[38]  + 37.78 ) ) / ( 1. + exp( 0.311 * ( y[38]  + 79.23 ) ) ) ) );
//     const double bj = ( ( y[38]  >= -40 ) * ( ( 0.6 * exp( 0.057 *  y[38] ) ) / ( 1. + exp( -0.1 * ( y[38]  + 32) ) ) )
//                        + ( y[38]  < -40 ) * ( ( 0.02424 * exp( -0.01052 *  y[38]  ) ) / (1. + exp( -0.1378 * ( y[38]  + 40.14) ) ) ) );
//     const double tauj = 1. / (aj + bj);
//     const double jss = 1. / ( pow( ( 1. + exp( ( y[38]  + 71.55 ) / 7.43 ) ) , 2 ) );
    
//     ydot[0]  = (mss -  y[0] ) / taum;
//     ydot[1]  = (hss -  y[1] ) / tauh;
//     ydot[2]  = (jss -  y[2] ) / tauj;
    
//     /*const double */I_Na_junc = Fjunc * GNa * y[0] * y[0] * y[0] * y[1] * y[2] * ( y[38] - ena_junc );
//     /*const double */I_Na_sl = Fsl * GNa * y[0] * y[0] * y[0] * y[1] * y[2] *( y[38] - ena_sl );
//     /*const double */I_Na = I_Na_junc + I_Na_sl;
    
    
//     //  Late I_Na
//     const double GNaL = para[16] *0.0025 * AF;
//     const double aml = 0.32 * ( y[38] + 47.13 ) / ( 1. - exp( -0.1 * ( y[38] + 47.13 ) ) );
//     const double bml = 0.08 * exp( - y[38] / 11. );
//     const double hlinf = 1. / ( 1. + exp( ( y[38] + 91. ) / 6.1 ) );
//     const double tauhl = 600;
//     ydot[59]  = aml * (1. - y[59] ) - bml * y[59] ;
//     ydot[60]  = ( hlinf - y[60] ) / tauhl;
    
//    /* const double */I_NaL_junc = Fjunc * GNaL * y[59]  * y[59]  * y[59] * y[60] * ( y[38] - ena_junc );
//    /* const double */I_NaL_sl = Fsl * GNaL * y[59] * y[59]  * y[59] * y[60] * ( y[38] - ena_sl );
//    /* const double */I_NaL = I_NaL_junc + I_NaL_sl;
    
    
//     /*if( t<9050 ) {
//         ydot[61] =0; }
//     else {
//         ydot[61] =I_NaL; }*/
//     // end
//     //  I_nabk: Na Background Current
//     /*const double*/ I_nabk_junc = Fjunc * GNaB * ( y[38] - ena_junc );
//     /*const double*/ I_nabk_sl = Fsl * GNaB * ( y[38] - ena_sl );
//     /*const double*/ I_nabk = I_nabk_junc + I_nabk_sl;
    
//     //  I_nak: Na/K Pump Current
//     const double sigma = ( exp( Nao / 67.3 ) - 1. ) / 7.;
//     const double fnak = 1. / ( 1. + 0.1245 * exp( -0.1 * y[38] * FoRT ) + 0.0365 * sigma * exp(- y[38] * FoRT ) );
//     /*const double*/ I_nak_junc = 1. * Fjunc * IbarNaK * fnak * Ko / ( 1. + pow( ( KmNaip / y[31] ) , 4 ) ) / ( Ko + KmKo );
//     /*const double*/ I_nak_sl = 1. * Fsl * IbarNaK * fnak * Ko / ( 1. + pow( ( KmNaip / y[32] ) , 4 ) ) / ( Ko + KmKo );
//     const double I_nak = I_nak_junc + I_nak_sl;
    
//     // //  I_kr: Rapidly Activating K Current
//     const double gkr = 0.035 * sqrt( Ko / 5.4 );
//     const double xrss = 1. / ( 1. + exp( -( y[38] + 10. ) / 5. ) ) ;
//     const double tauxr = 550. / ( 1. + exp( ( -22. - y[38] ) / 9. ) ) * 6. / ( 1. + exp( ( y[38] - (-11) ) / 9. ) ) + 230. / ( 1. + exp( ( y[38] - (-40.) ) / 20. ) );
//     ydot[11]  = ( xrss- y[11] ) / tauxr;
//     const double rkr = 1. / ( 1. + exp( ( y[38] + 74. ) / 24. ) );
//     const double I_kr = gkr * y[11] * rkr * ( y[38] - ek );
    
    
//     // //  I_ks: Slowly Activating K Current
//     const double markov_iks = 0;
//     //  pcaks_junc = -log10( y[35] )+3.0;
//     //  pcaks_sl = -log10( y[36] )+3.0;
//     //  gks_junc = 0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_junc)/0.6)));
//     //  gks_sl = 0.07*(0.057 +0.19/(1+ exp((-7.2+pcaks_sl)/0.6)));
    
//     const double eks = (1/FoRT)*log((Ko+pNaK*Nao)/( y[34] +pNaK* y[33] ));
    
//     double  gks_junc, gks_sl, xsss,  tauxs ,I_ks_junc , I_ks_sl , I_ks , alpha, beta, gamma, delta, teta, eta, psi, omega, O2;

//     if( markov_iks == 0 ) {
//          gks_junc =para[11] *  1. * ( 1. + 1. * AF + 2. * ISO ) * 0.0035 * 1. ;
//          gks_sl = para[11] * 1. * ( 1. + 1. * AF + 2. * ISO ) * 0.0035 * 1. ; // FRA
//          xsss = 1. / ( 1. + exp( -( y[38] + 40. * ISO + 3.8 ) / 14.25 ) ); //  fitting Fra
//          tauxs = 990.1 / ( 1. + exp( -( y[38] + 40. * ISO + 2.436 ) / 14.12 ) );
//         ydot[12]  = ( xsss - y[12] ) / tauxs;
//          I_ks_junc = Fjunc * gks_junc * y[12] * y[12] * ( y[38] - eks );
//          I_ks_sl = Fsl * gks_sl * y[12] * y[12] * ( y[38] - eks );
//          I_ks = I_ks_junc + I_ks_sl;
//     } else {
//          gks_junc = para[11] * 1. * 0.0065;
//          gks_sl = para[11] * 1. * 0.0065; // FRA
//          alpha = 3.98e-4 * exp( 3.61e-1 * y[38] * FoRT );
//          beta = 5.74e-5 * exp( -9.23e-2 * y[38] * FoRT );
//          gamma = 3.41e-3 * exp( 8.68e-1 * y[38] * FoRT );
//          delta = 1.2e-3 * exp( -3.3e-1 * y[38] * FoRT );
//          teta = 6.47e-3;
//          eta = 1.25e-2 * exp( -4.81e-1 * y[38] * FoRT );
//          psi = 6.33e-3 * exp( 1.27 * y[38] * FoRT );
//          omega = 4.91e-3 * exp( -6.79e-1 * y[38] * FoRT );
        
//         ydot[41] = -4. * alpha * y[41] + beta * y[42] ;
//         ydot[42] = 4. * alpha * y[41] - ( beta + gamma + 3. * alpha ) * y[42] + 2. * beta * y[43] ;
//         ydot[43] = 3. * alpha * y[42] - ( 2. * beta + 2. * gamma + 2. * alpha ) * y[43] + 3. * beta * y[44] ;
//         ydot[44] = 2. * alpha * y[43] - ( 3. * beta + 3. * gamma + alpha ) * y[44] + 4. * beta * y[45] ;
//         ydot[45] = 1. * alpha * y[43] - ( 4. * beta + 4. * gamma ) * y[45] + delta * y[49] ;
//         ydot[46] = gamma * y[42] - ( delta + 3. * alpha ) * y[46] + beta * y[47] ;
//         ydot[47] = 2. * gamma * y[43] + 3. * alpha * y[46] - ( delta + beta + 2. * alpha + gamma ) * y[47] + 2. * beta * y[48] + 2. * delta * y[50] ;
//         ydot[48] = 3. * gamma * y[44] + 2. * alpha * y[47] - ( delta + 2. * beta + 1. * alpha + 2. * gamma ) * y[48] + 3. * beta * y[49] + 2. * delta * y[51] ;
//         ydot[49] = 4. * gamma * y[45] + 1. * alpha * y[48] - ( delta + 3. * beta + 0 * alpha + 3. * gamma ) * y[49] + 2. * delta * y[52] ;
//         ydot[50] = 1. * gamma * y[47] - ( 2. * delta + 2. * alpha ) * y[50] + beta * y[51] ;
//         ydot[51] = 2. * gamma * y[48] + 2. * alpha * y[50] - ( 2. * delta + beta + 1. * alpha + gamma ) * y[51] + 2. * beta * y[52] + 3. * delta * y[53] ;
//         ydot[52] = 3. * gamma * y[49] + 1. * alpha * y[51] - ( 2. * delta + 2. * beta + 2. * gamma ) * y[52] + 3. * delta * y[54] ;
//         ydot[53] = 1. * gamma * y[51] - ( 3. * delta + 1. * alpha ) * y[53] + beta * y[54] ;
//         ydot[54] = 2. * gamma * y[52] + 1. * alpha * y[53] - ( 3. * delta + 1. * beta + 1. * gamma ) * y[54] + 4. * delta * y[55] ;
//         ydot[55] = 1. * gamma * y[54] - ( 4. * delta + teta ) * y[55] + eta * y[56] ;
//          O2 = 1. - ( y[41] + y[42] + y[43] + y[44] + y[45] + y[46] + y[48] + y[47] + y[49] + y[50] + y[51] + y[52] + y[53] + y[54] + y[55] + y[56] );
//         ydot[56] = 1. * teta * y[55] - ( eta + psi ) * y[56] + omega * O2;
//          I_ks_junc = Fjunc * gks_junc * ( y[56] +O2 ) * ( y[38] - eks );
//          I_ks_sl = Fsl * gks_sl * ( y[56] + O2 ) * ( y[38] - eks );
//          I_ks = I_ks_junc + I_ks_sl;
//     }
//     // end
    
//     // I_kp: Plateau K current
//     const double kp_kp = 1. / ( 1. + exp( 7.488 - y[38] / 5.98 ) );
//     const double I_kp_junc = Fjunc * gkp * kp_kp * ( y[38] - ek );
//     const double I_kp_sl = Fsl * gkp * kp_kp * ( y[38] - ek );
//     /*const double*/ I_kp = I_kp_junc + I_kp_sl;
    
//     // //  I_to: Transient Outward K Current (slow and fast components)
//     //  modified for human myocytes
    
//     const double GtoFast = para[12] * ( 1.0 - 0.7 * AF ) * 0.165 * 1.0; // nS/pF maleckar; // human atrium
    
//     // 11/12/09; changed Itof to that from maleckar/giles/2009; removed I_tos
//     // atrium
//     // equations for activation;
//     const double xtoss = ( (1.) / ( 1. + exp( -( y[38] + 1.0 ) / 11.0 ) ) );
//     const double tauxtof = 3.5 * exp( -( pow( ( y[38] / 30.0 ) , 2.0 ) ) ) + 1.5;
    
//     // equations for inactivation;
//     const double ytoss = ( (1.0) / ( 1. + exp( ( y[38] + 40.5 ) / 11.5 ) ) ) ;
//     const double tauytof = 25.635 * exp( -( pow( ( ( y[38] + 52.45 ) / 15.8827 ) , 2.0 ) ) ) + 24.14;// 14.14
    
//     ydot[9]  = ( xtoss- y[9] ) / tauxtof;
//     ydot[10]  = ( ytoss - y[10] ) / tauytof;
//     const double I_tof = 1.0 * GtoFast * y[9] * y[10] * ( y[38] - ek );
    
//     /*const double*/ I_to = 1. * I_tof;
    
//     // //  I_kur: Ultra rapid delayed rectifier Outward K Current
//     // Equation for IKur; from Maleckar et al. 2009 - EG
//     // atrium
//     // equations for activation;
//     const double Gkur = para[14] * 1. * ( 1.0 - 0.5 * AF ) * ( 1. + 2. * ISO ) * 0.045 * ( 1. + 0.2 * RA ); // nS/pF maleckar 0.045
//     const double xkurss = ( (1.) / ( 1 + exp( ( y[38] + 6. ) / (-8.6) ) ) );
//     const double tauxkur = 9. / ( 1. + exp( ( y[38] + 5. ) / 12.0 ) ) + 0.5;
    
//     // equations for inactivation;
//     const double ykurss = ( (1.) / ( 1. + exp( ( y[38] + 7.5 ) / 10.0 ) ) );
//     const double tauykur = 590. / ( 1. + exp( ( y[38] + 60. ) / 10.0 ) ) + 3050.;
    
//     ydot[57]  = ( xkurss - y[57] ) / tauxkur;
//     ydot[58]  = ( ykurss - y[58] ) / tauykur;
//     /*const double*/ I_kur = 1. * Gkur * y[57] * y[58] * ( y[38] - ek );
    
//     // //  I_ki: Time-Independent K Current
//     const double aki = 1.02 / ( 1. + exp( 0.2385 * ( y[38] - ek - 59.215 ) ) );
//     const double bki = ( 0.49124 * exp( 0.08032 * ( y[38] + 5.476 - ek ) ) + exp( 0.06175 * ( y[38] - ek - 594.31 ) ) ) / ( 1. + exp( -0.5143 * ( y[38] - ek + 4.753 ) ) );
//     const double kiss = aki / ( aki + bki );
    
//     // I_ki =1* 0.35*sqrt(Ko/5.4)*kiss*( y[38] -ek);
//     // SVP 11/11/09
//     // multiplieD IK1 by 0.15 to scale it to single cell isolated atrial cell
//     // resting potential
//     /*const double */I_ki = para[15] *( 1. + 1. * AF ) * 0.0525 * sqrt( Ko / 5.4 ) * kiss * ( y[38] - ek );
    
//     //  I_ClCa: Ca-activated Cl Current, I_Clbk: background Cl Current
//     /*const double */I_ClCa_junc = Fjunc * GClCa / ( 1. + KdClCa / y[35] ) * ( y[38] - ecl );
//     /*const double */I_ClCa_sl = Fsl * GClCa / ( 1. + KdClCa / y[36] ) * ( y[38] - ecl );
//     /*const double */I_ClCa = I_ClCa_junc + I_ClCa_sl;
//     /*const double */I_Clbk = GClB * ( y[38] - ecl );
    
//     const double GClCFTR = 0;// 4.9e-3*ISO;     //  [mS/uF]
//     /*const double */I_ClCFTR = GClCFTR * ( y[38] - ecl );
    
//     // //  I_Ca: L-type Calcium Current
//     const double dss = 1. / ( 1. + exp( -( y[38] + 3 * ISO + 9 ) / 6 ) ); // in Maleckar v1/2=-9 S=6 (mV); Courtemanche v1/2=-9 S=5.8 (mV)
//     const double taud = 1. * dss * ( 1. - exp( -( y[38] + 3 * ISO + 9 ) / 6 ) ) / ( 0.035 * ( y[38] + 3 * ISO + 9 ) );
//     const double fss = 1. / ( 1. + exp( ( y[38] + 3. * ISO + 30. ) / 7. ) ) + 0.2 / ( 1. + exp( ( 50. - y[38] - 3. * ISO ) / 20. ) ); //  in Maleckar v1/2=-27.4 S=7.1 (mV); Courtemanche v1/2=-28 S=6.9 (mV)
//     const double tauf = 1. / ( 0.0197 * exp( - pow( ( 0.0337 * ( y[38] + 3. * ISO + 25. ) ) , 2 ) ) + 0.02 );
//     ydot[3]  = (dss- y[3] )/taud;
//     ydot[4]  = (fss- y[4] )/tauf;
//     ydot[5]  = 1.7 * y[35] * ( 1. - y[5] ) - 1. * 11.9e-3 * y[5] ; //  fCa_junc   koff!!!!!!!!
//     ydot[6]  = 1.7 * y[36] * ( 1. - y[6] ) - 1. * 11.9e-3 * y[6] ; //  fCa_sl
//     // const double fcaCaMSL = 0.1 / ( 1. + ( 0.01 / y[36] ) );
//     // const double fcaCaj = 0.1 / ( 1. + ( 0.01 / y[35] ) );
//     const double fcaCaMSL = 0;
//     const double fcaCaj = 0;
//     const double ibarca_j = pCa * 4. * ( y[38] * Frdy * FoRT ) * ( 0.341 * y[35] * exp( 2 * y[38] * FoRT ) - 0.341 * Cao ) / ( exp( 2 * y[38] * FoRT ) - 1. );
//     const double ibarca_sl = pCa * 4. * ( y[38] * Frdy * FoRT ) * ( 0.341 * y[36] * exp( 2 * y[38] * FoRT ) - 0.341 * Cao ) / ( exp( 2 * y[38] * FoRT ) - 1. );
//     const double ibark = pK * ( y[38] * Frdy * FoRT ) * ( 0.75 * y[34] * exp( y[38] * FoRT ) - 0.75 * Ko ) / ( exp( y[38] * FoRT ) - 1. );
//     const double ibarna_j = pNa * ( y[38] * Frdy * FoRT ) * ( 0.75 * y[31] * exp( y[38] * FoRT ) - 0.75 * Nao ) / ( exp( y[38] * FoRT ) - 1. );
//     const double ibarna_sl = pNa * ( y[38] * Frdy * FoRT ) * ( 0.75 * y[32] * exp( y[38] * FoRT ) - 0.75 * Nao ) / ( exp( y[38] * FoRT ) - 1. );
//     /*const double*/ I_Ca_junc = ( Fjunc_CaL * ibarca_j * y[3] * y[4] *( ( 1. - y[5] ) + fcaCaj ) * pow( Q10CaL , Qpow ) ) * 0.45;
//     /*const double*/ I_Ca_sl = ( Fsl_CaL * ibarca_sl * y[3] * y[4] * ( ( 1. - y[6] ) + fcaCaMSL ) * pow( Q10CaL , Qpow ) ) * 0.45;
//     /*const double*/ I_Ca = I_Ca_junc + I_Ca_sl;
//     /*const double*/ I_CaK = ( ibark * y[3] * y[4] * ( Fjunc_CaL * ( fcaCaj + ( 1. - y[5] ) ) + Fsl_CaL * ( fcaCaMSL + ( 1. - y[6] ) ) ) * pow( Q10CaL , Qpow ) ) * 0.45;
//     /*const double*/ I_CaNa_junc = ( Fjunc_CaL * ibarna_j * y[3] * y[4] * ( ( 1. - y[5] ) + fcaCaj ) * pow( Q10CaL , Qpow ) ) * 0.45;
//     /*const double*/ I_CaNa_sl = ( Fsl_CaL * ibarna_sl * y[3] * y[4] * ( ( 1. - y[6] ) + fcaCaMSL ) * pow( Q10CaL , Qpow ) ) * 0.45;
//     /*const double*/ I_CaNa = I_CaNa_junc + I_CaNa_sl;
//     /*const double*/ I_Catot = I_Ca + I_CaK + I_CaNa;
    
//     //  I_ncx: Na/Ca Exchanger flux
//     const double Ka_junc = 1. / ( 1. + ( Kdact * Kdact / ( y[35] * y[35] ) ) );
//     const double Ka_sl = 1. / ( 1. + ( Kdact * Kdact / ( y[36] * y[36] ) ) );
//     const double s1_junc = exp( nu * y[38] * FoRT ) * y[31] * y[31] * y[31] * Cao;
//     const double s1_sl = exp( nu * y[38] * FoRT ) * y[32] * y[32] * y[32] * Cao;
//     const double s2_junc = exp( ( nu - 1. ) * y[38] * FoRT ) * Nao * Nao * Nao * y[35] ;
//     const double s3_junc = ( KmCai * Nao * Nao * Nao * ( 1. + ( y[31] * y[31] * y[31] / ( KmNai * KmNai * KmNai ) ) )
//                             + KmNao * KmNao * KmNao * y[35] * ( 1. + y[35] / KmCai )
//                             + KmCao * y[31] * y[31] * y[31] + y[31] * y[31] * y[31] * Cao
//                             + Nao * Nao * Nao * y[35] );
//     const double s2_sl = exp( ( nu - 1. ) * y[38] * FoRT ) * Nao * Nao * Nao * y[36] ;
//     // const double s3_sl = KmCai * Nao * Nao * Nao * ( 1. + ( y[32] / KmNai ) ^3) + KmNao ^3 * y[36] * ( 1. + y[36] / KmCai ) + KmCao * y[32] ^3+ y[32] ^3 * Cao + Nao ^3 * y[36] ;
//     const double s3_sl = ( KmCai * Nao * Nao * Nao * ( 1. + ( y[32] * y[32] * y[32] / ( KmNai * KmNai * KmNai ) ) )
//                           + KmNao * KmNao * KmNao * y[36] * ( 1. + y[36] / KmCai )
//                           + KmCao * y[32] * y[32] * y[32]
//                           + y[32] * y[32] * y[32] * Cao
//                           + Nao * Nao * Nao * y[36] );
    
    
//     /*const double */I_ncx_junc = Fjunc * IbarNCX * pow( Q10NCX , Qpow ) * Ka_junc * ( s1_junc - s2_junc ) / s3_junc / ( 1. + ksat * exp( ( nu - 1. ) * y[38] * FoRT ) );
//     /*const double */I_ncx_sl = Fsl * IbarNCX * pow( Q10NCX , Qpow ) * Ka_sl * ( s1_sl - s2_sl ) / s3_sl / ( 1. + ksat * exp( ( nu - 1. ) * y[38] * FoRT ) );
//     /*const double */I_ncx = I_ncx_junc + I_ncx_sl;
    
//     //  I_pca: Sarcolemmal Ca Pump Current
//     /*const double */I_pca_junc = Fjunc * pow( Q10SLCaP , Qpow ) * IbarSLCaP * pow( y[35] , 1.6 ) / ( pow( KmPCa , 1.6 ) + pow( y[35] , 1.6 ) );
//     /*const double */I_pca_sl = Fsl * pow( Q10SLCaP , Qpow ) * IbarSLCaP * pow( y[36] , 1.6 ) / ( pow( KmPCa , 1.6 ) + pow( y[36] , 1.6 ) );
//     /*const double */I_pca = I_pca_junc + I_pca_sl;
    
//     //  I_cabk: Ca Background Current
//     /*const double */I_cabk_junc = Fjunc * GCaB * ( y[38] - eca_junc );
//     /*const double */I_cabk_sl = Fsl * GCaB * ( y[38] - eca_sl);
//     /*const double */I_cabk = I_cabk_junc + I_cabk_sl;
    
//     // //  SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
//     const double MaxSR = 15.;
//     const double MinSR = 1.;
//     const double kCaSR = MaxSR - ( MaxSR - MinSR ) / ( 1. + pow( ( ec50SR / y[30] ) , 2.5 ) );
//     const double koSRCa = (1.) * koCa / kCaSR;//
//     const double kiSRCa = kiCa * kCaSR;
//     const double RI = 1. - y[13] - y[14] - y[15] ;
//     ydot[13]  = ( kim * RI - kiSRCa * y[35] * y[13] ) - ( koSRCa * y[35] * y[35] * y[13] - kom * y[14] );   //  R
//     ydot[14]  = ( koSRCa * y[35] * y[35] * y[13] - kom * y[14] ) - ( kiSRCa * y[35] * y[14] - kim * y[15] );//  O
//     ydot[15]  = ( kiSRCa * y[35] * y[14] - kim * y[15] ) - ( kom * y[15] - koSRCa * y[35] * y[35] * RI);   //  I
//     /*const double*/ J_SRCarel = ks * y[14] * ( y[30] - y[35] );          //  [mM/ms]
    
//     /*const double*/ J_serca = ( 1.0 * pow( Q10SRCaP, Qpow ) * Vmax_SRCaP * ( pow( ( y[37] / Kmf ) , hillSRCaP ) - pow( ( y[30] / Kmr ) , hillSRCaP ) )
//                             / ( 1. + pow( ( y[37] / Kmf ) , hillSRCaP ) + pow( ( y[30] / Kmr ) , hillSRCaP ) ) );
//     /*const double*/ J_SRleak = para[13] *(1.) * ( 1.0 + 0.25 * AF ) * 5.348e-6 * ( y[30] - y[35] );           //    [mM/ms]
    
    
//     // //  Sodium and Calcium Buffering
//     ydot[16]  = kon_na* y[31] *(Bmax_Naj- y[16] )-koff_na* y[16] ;        //  NaBj      [mM/ms]
//     ydot[17]  = kon_na* y[32] *(Bmax_Nasl- y[17] )-koff_na* y[17] ;       //  NaBsl     [mM/ms]
    
//     //  Cytosolic Ca Buffers
//     ydot[18]  = kon_tncl* y[37] *(Bmax_TnClow- y[18] )-koff_tncl* y[18] ;            //  TnCL      [mM/ms]
//     ydot[19]  = kon_tnchca* y[37] *(Bmax_TnChigh- y[19] - y[20] )-koff_tnchca* y[19] ; //  TnCHc     [mM/ms]
//     ydot[20]  = kon_tnchmg*Mgi*(Bmax_TnChigh- y[19] - y[20] )-koff_tnchmg* y[20] ;   //  TnCHm     [mM/ms]
//     ydot[21]  = kon_cam* y[37] *(Bmax_CaM- y[21] )-koff_cam* y[21] ;                 //  CaM       [mM/ms]
//     ydot[22]  = kon_myoca* y[37] *(Bmax_myosin- y[22] - y[23] )-koff_myoca* y[22] ;    //  Myosin_ca [mM/ms]
//     ydot[23]  = kon_myomg*Mgi*(Bmax_myosin- y[22] - y[23] )-koff_myomg* y[23] ;      //  Myosin_mg [mM/ms]
//     ydot[24]  = kon_sr* y[37] *(Bmax_SR- y[24] )-koff_sr* y[24] ;                    //  SRB       [mM/ms]
//     const double J_CaB_cytosol = ydot[18] + ydot[19] + ydot[20] + ydot[21] + ydot[22] + ydot[23] + ydot[24]; // sum(ydot(19:25));
    
//     //  Junctional and SL Ca Buffers
//     ydot[25]  = kon_sll* y[35] *(Bmax_SLlowj- y[25] )-koff_sll* y[25] ;       //  SLLj      [mM/ms]
//     ydot[26]  = kon_sll* y[36] *(Bmax_SLlowsl- y[26] )-koff_sll* y[26] ;      //  SLLsl     [mM/ms]
//     ydot[27]  = kon_slh* y[35] *(Bmax_SLhighj- y[27] )-koff_slh* y[27] ;      //  SLHj      [mM/ms]
//     ydot[28]  = kon_slh* y[36] *(Bmax_SLhighsl- y[28] )-koff_slh* y[28] ;     //  SLHsl     [mM/ms]
//     /*const double */J_CaB_junction =  ydot[25] + ydot[27] ;
//     const double J_CaB_sl =  ydot[26] + ydot[28] ;
    
//     // //  Ion concentrations
//     //  SR Ca Concentrations
//     ydot[29]  = kon_csqn * y[30] * ( Bmax_Csqn - y[29] ) - koff_csqn * y[29] ;       //  Csqn      [mM/ms]
//     ydot[30]  = J_serca - ( J_SRleak * Vmyo / Vsr + J_SRCarel ) - ydot[29] ;         //  Ca_sr     [mM/ms] // Ratio 3 leak current
//     //   ydot[30] =0;
    
//     //  Sodium Concentrations
//     /*const double*/ I_Na_tot_junc = I_Na_junc + I_nabk_junc + 3 * I_ncx_junc + 3 * I_nak_junc + I_CaNa_junc + I_NaL_junc;   //  [uA/uF]
//     /*const double*/ I_Na_tot_sl = I_Na_sl + I_nabk_sl + 3 * I_ncx_sl + 3 * I_nak_sl + I_CaNa_sl + I_NaL_sl;   //  [uA/uF]
//     /*const double*/ I_Na_tot_sl2 = 3 * I_ncx_sl + 3 * I_nak_sl + I_CaNa_sl;   //  [uA/uF]
//     /*const double*/ I_Na_tot_junc2 = 3 * I_ncx_junc + 3 * I_nak_junc + I_CaNa_junc;   //  [uA/uF]
    
//     ydot[31]  = -I_Na_tot_junc * Cmem / ( Vjunc * Frdy ) + J_na_juncsl / Vjunc * ( y[32] - y[31] ) - ydot[16] ;
//     ydot[32]  = ( -I_Na_tot_sl * Cmem / ( Vsl * Frdy ) + J_na_juncsl / Vsl * ( y[31] - y[32] )
//                  + J_na_slmyo / Vsl * ( y[33] - y[32] ) - ydot[17] );
//     // FluxNaSL= ydot[32] ;
//     //   ydot[31]  = 0;
//     //   ydot[32]  = 0;
//     ydot[33]  = J_na_slmyo/Vmyo*( y[32] - y[33] );             //  [mM/msec]
//     //   ydot[33] =0;
    
//     //  Potassium Concentration
//     /*const double */I_K_tot = I_to + I_kr + I_ks + I_ki - 2 * I_nak + I_CaK + I_kp + I_kur;     //  [uA/uF] // SVP: added IKur
//     //   ydot[34]  = 0; // -I_K_tot*Cmem/(Vmyo*Frdy);           //  [mM/msec]
//     ydot[34]  = 0; //  -I_K_tot*Cmem/(Vmyo*Frdy);
    
//     //  Calcium Concentrations
//     /*const double */I_Ca_tot_junc = I_Ca_junc + I_cabk_junc + I_pca_junc - 2 * I_ncx_junc;                   //  [uA/uF]
//     /*const double*/ I_Ca_tot_sl = I_Ca_sl + I_cabk_sl + I_pca_sl - 2 * I_ncx_sl;            //  [uA/uF]
//     ydot[35]  = ( -I_Ca_tot_junc * Cmem / ( Vjunc * 2 * Frdy ) + J_ca_juncsl / Vjunc * ( y[36] - y[35] )
//                  - J_CaB_junction + ( J_SRCarel ) * Vsr / Vjunc + J_SRleak * Vmyo / Vjunc );  //  Ca_j
//     ydot[36]  = ( -I_Ca_tot_sl * Cmem / ( Vsl * 2 * Frdy ) + J_ca_juncsl / Vsl * ( y[35] - y[36] )
//                  + J_ca_slmyo / Vsl * ( y[37] - y[36] ) - J_CaB_sl );   //  Ca_sl
//     //   ydot[35] =0;
//     //   ydot[36] =0;
//     //   ydot[37]  = -J_serca*Vsr/Vmyo-J_CaB_cytosol;// +J_ca_slmyo/Vmyo*( y[36] - y[37] );    //  [mM/msec]
//     ydot[37]  = -J_serca * Vsr / Vmyo - J_CaB_cytosol + J_ca_slmyo / Vmyo * ( y[36] - y[37] );
//     //   ydot[37] =0;
    
//     double rate, factor, V_hold, V_test, V_clamp, R_clamp;
// 	double tdelay[4] = { 10, 1, 5, 80 };
    

    
//     // end
    
//     // //  Membrane Potential
//     // //
//     /*const double*/ I_Na_tot = I_Na_tot_junc + I_Na_tot_sl;          //  [uA/uF]
//     /*const double*/ I_Cl_tot = I_ClCa + I_Clbk + I_ClCFTR;                        //  [uA/uF]
//     /*const double*/ I_Ca_tot = I_Ca_tot_junc+I_Ca_tot_sl;
//     /*const double*/ I_tot = I_Na_tot+I_Cl_tot+I_Ca_tot+I_K_tot;
//     //  ydot[38]  = -(I_Ca_tot+I_K_tot+I_Na_tot-I_app);



//     double i_Stim = S1(20,  12.5,  BCL,  t,  5.0);
//     // dY[38] = -(I_tot - i_Stim);

//     I_app = i_Stim;

//     dV = ydot[38];
//     CaSR = y[30];
//     Caj = y[35];
//     Casl = y[36];
//     Cai = y[37];
//     Nai = y[31];
//     Naj = y[32];
//     Nasl = y[33];


//     ydot[38]  = -(I_tot-I_app);
//     dV =  ydot[38] ;


    
//     //  ----- END EC COUPLING MODEL ---------------
//     //  adjust output depending on the function call
    
//     //    if (nargin == 3)
//     //        output = ydot;
//     //    elseif (nargin == 4) & strcmp(runType,'ydot')
//     //        output = ydot;
//     //    elseif (nargin == 4) & strcmp(runType,'rates')
//     //        output = r;
//     //    elseif (nargin == 4) & strcmp(runType,'currents')
//     //        // currents = [I_Na I_nabk I_nak I_kr I_ks I_kp I_tos I_tof I_ki I_ClCa I_Clbk I_Catot I_ncx I_pca I_cabk J_serca*Vmyo/Vsr];
//     //        // currents = [I_Na I_tof I_tos I_kr I_ks I_ClCa I_Catot J_SRCarel*Vsr/Vmyo J_SRleak RI I_ncx];
//     //        //      currents= [I_Catot J_SRCarel*2*Vsr*Frdy/Cmem];
//     //        //  currents=[I_Na I_Catot I_ncx I_nak I_kr I_ks I_tof I_tos I_ki];
//     //        currents=[];
//     //        output = currents;
//     //    end
    
//     // currents[0] = I_Na + I_NaL;
//     // currents[1] = I_Catot;
//     // currents[2] = I_kr;
//     // currents[3] = I_ks;
//     // currents[4] = I_to;
//     // currents[5] = I_kur;
//     // currents[6] = I_ki;
//     // currents[7] = I_nak;
//     // currents[8] = I_ncx;
//     // currents[9] = ;

// }
