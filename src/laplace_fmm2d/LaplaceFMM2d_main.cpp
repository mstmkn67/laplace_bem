//Boundary Element Method and Fast Multipole Method
// for 2 dimensional Laplace equation solver
//                                     Masato MAKINO 

#include "../common/udf/gourmain.h"
#include "LaplaceFMM2d.h"

void udfHeaderCheck()
{
	string version("1.1"),engine("laplace_fmm2d");
	cout << "**************************************************************" << endl;
	cout <<  "      " <<  engine << "  " << version << "           " << endl;
	cout << "                                        Masato MAKINO         " << endl;
	cout << "**************************************************************" << endl;
	cout <<  endl;
}

void error_massage(){
	cout << "usage: laplace_fmm2d -I input_udf [-O output_udf] " << endl;
}


int gourmain(UDFManager* in,UDFManager* out){
	udfHeaderCheck();
	LaplaceFMM2d* sim=new LaplaceFMM2d(in,out);
	sim->run();
	delete sim;
	return 0;
}
