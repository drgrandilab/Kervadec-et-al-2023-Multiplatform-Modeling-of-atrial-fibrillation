/* Grandi model with ISO effect
* with Latest I_kur model here.
* Haibo Ni <haibo.ni0822@gmail.com>

* Wed 08 Jul 2015 10:27:26 BST
*/
/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/

#include <stdio.h>
#include <cmath>
#include <malloc.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
#include "cvode_solver.hpp"
#include "stimulus.h"
#include "APInfo.hpp"

#include "GB_ECC.hpp"
#include <iostream>


/**********************************************************************************************/
/**********************************************************************************************/
/**********************************************************************************************/


#define num_beats 300
#define NEQ 42                /* number of equations  */ //without AF



#define dt 0.1

/* Functions Called by the Solver */
// static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

/*
 *--------------------------------------------------------------------------------
 * Main Program
 *--------------------------------------------------------------------------------
 */



// FILE *current_output;
// #pragma optimize("", off)
int main (int argc, char *argv[])
{
	realtype t, tout;

	int iout;

	clock_t time_begin, time_end;



	FILE *output, *initial_conditions, *current_output;
	int init_counter = 0;

	/* The measurement arrays */
	int NOUT;


	/*	double para[17] = {0};
		for (int i = 0; i < 17; ++i)
		{para[i] = 1.0;}*/





	/**
	* GNa - para[0]
	* INaL - para[16];
	* I_kur - para[14];
	* Gto - para[12]
	* IK1 - para[15]
	* ICaL- para[7]
	* INCX - para[8];
	* SERCA - para[9];

	* RyR - para[10]
	*/

	/*para[0] = atof(argv[1]);
	para[16] = atof(argv[2]);
	para[14] = atof(argv[3]);
	para[12] = atof(argv[4]);
	para[15] = atof(argv[5]);
	para[7] = atof(argv[6]);
	para[8] = atof(argv[7]);
	para[9] = atof(argv[8]);
	para[10] = atof(argv[9]);*/



	// // std::string thread_ID;
	// if (argc == 17) {
	// 	for (int i = 0; i < 17; ++i)
	// 	{

	// 		para[i] = para[i] * atof(argv[i + 1]); // (1:16);
	// 		// std::cout <<  "P_init[" << i << "] = " << p0[i] << endl;
	// 	}
	// }



	/* Start the automation loops here */

	double i_Stim = 0.0;
	double outputTime = 0.0;
	double pcl = 1000.0;

	// y = NULL; // for solver
	// cvode_mem = NULL; // for solver


	// y = N_VNew_Serial(NEQ); // for solver
	// ydot = N_VNew_Serial(NEQ); // for solver


	time_begin = clock();

	GB_ECC mycell;
	cvode_solver cvode(NEQ, 1);
	cvode.set_IC(mycell.y);
	cvode.initialise_mem(fnew);
	cvode.set_user_data(&mycell);
	pcl = atof(argv[1]);
	mycell.kmf_scale = atof(argv[2]);
	mycell.ISO = atof(argv[3]);
	mycell.BCL = pcl;

	if (argc == 18 + 4) {
		for (int i = 0; i < 18; ++i)
		{
			mycell.para[i] = mycell.para[i] * atof(argv[i + 4]); // (1:16);
			// std::cout <<  "P_init[" << i << "] = " << mycell.para[i] << std::endl;  // for debug
		}

	}

	NOUT = int( ((int)(((num_beats + 1) * 1000)) / (int) (pcl) * pcl + pcl - 2) / dt);


	// if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);
	iout = 0;  tout = dt;

	// output = fopen(str, "w+");
	// current_output = fopen("current.dat", "w+");
	APInfor AP("APD_measure.dat", false, true);
	std::ofstream output_file( "HAM_wrap_out.dat");                      // output filename
		int skip = 5;
	mycell.allow_stimulation_flag = true;
	for (iout = 0; iout < NOUT; iout++)
	{
		// CVode(cvode_mem, tout, y, &t, CV_NORMAL); // 1 time step solution.
		cvode.solve_single_step(tout);


		skip = 150;

		if(fabs(mycell.dV) > 0.7)
			skip = 5;
		if (tout >= 0*num_beats * 1000 - 1100/*9000 and tout <= 9650 + 0.05*/)
			if (iout % skip == 0 ) {
				// fprintf(output, "%f %f %f %f %f %f\n", tout, Ith(y, 39), Ith(y, 38), Ith(y, 31), Ith(y, 32), Ith(y, 34));
				// fprintf(output, "%f\t%f\t%f\n", tout, Ith(cvode.y, 39), Ith(cvode.y, 38));
				// printf("%f\t%f\t%f\n", tout, Ith(cvode.y, 39), Ith(cvode.y, 1));
				mycell.print_to_file(tout, output_file);

				// fprintf(current_output, "%f %f %f %f %f %f %f %f %f %f %f\n", tout, I_Ca, I_ki, I_kr, I_ks, I_Na, I_kur, I_ncx, I_to, I_nak, I_nabk);
			}


		AP.MeasureAPD90_INa(tout, mycell.I_app, pcl, dt, mycell.y[38], 0, mycell.y[37]);
		tout += dt;

		if( tout > 280*1000-1-dt/2.0 and tout < 280*1000-1+dt/2.0) {
			mycell.output_inital_condition("Ini.bin");

			// std::cout << "outputing file" << std::endl;
		}


		if(tout > 290*1000-10)
			mycell.allow_stimulation_flag = false;
	} /* end of time loop */
	time_end = clock();

	AP.ReportLastThree();
	output_file.close();

	return (0);
} // end of main.


   