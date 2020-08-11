//Boundary Element Method
// for 2 dimensional Laplace equation solver
//                                     Masato MAKINO 

#include "../common/udf/gourmain.h"
#include "Laplace2d.h"

void udfHeaderCheck()
{
	string version("1.0"),engine("laplace2d");
	cout << "**************************************************************" << endl;
	cout <<  "      " <<  engine << "  " << version << "           " << endl;
	cout << "                                        Masato MAKINO         " << endl;
	cout << "**************************************************************" << endl;
	cout <<  endl;
}

void error_massage(){
	cout << "usage: laplace2d -I input_udf [-O output_udf] " << endl;
}


int gourmain(UDFManager* in,UDFManager* out){
	udfHeaderCheck();
	Laplace2d* sim=new Laplace2d(in,out);
	sim->run();
	delete sim;
	return 0;
}
