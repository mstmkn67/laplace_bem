//Boundary Element Method and Fast Multipole Method
// for 3 dimensional Laplace equation solver
//                                     Masato MAKINO 

#include "../common/udf/gourmain.h"
#include "LaplaceFMM3d.h"

void udfHeaderCheck()
{
	string version("1.1"),engine("laplace_fmm3d");
	cout << "**************************************************************" << endl;
	cout <<  "      " <<  engine << "  " << version << "           " << endl;
	cout << "                                        Masato MAKINO         " << endl;
	cout << "**************************************************************" << endl;
	cout <<  endl;
}

void error_massage(){
	cout << "usage: laplace_fmm3d -I input_udf [-O output_udf] " << endl;
}


int gourmain(UDFManager* in,UDFManager* out){
	udfHeaderCheck();
	LaplaceFMM3d* sim=new LaplaceFMM3d(in,out);
	sim->run();
	delete sim;
	return 0;
}
