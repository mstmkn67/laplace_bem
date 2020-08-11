//Boundary Element Method
// for 3 dimensional Laplace equation solver
//                                     Masato MAKINO 

#include "../common/udf/gourmain.h"
#include "Laplace3d.h"

void udfHeaderCheck()
{
	string version("1.0"),engine("laplace3d");
	cout << "**************************************************************" << endl;
	cout <<  "      " <<  engine << "  " << version << "           " << endl;
	cout << "                                        Masato MAKINO         " << endl;
	cout << "**************************************************************" << endl;
	cout <<  endl;
}

void error_massage(){
	cout << "usage: laplace3d -I input_udf [-O output_udf] " << endl;
}


int gourmain(UDFManager* in,UDFManager* out){
	udfHeaderCheck();
	Laplace3d* sim=new Laplace3d(in,out);
	sim->run();
	delete sim;
	return 0;
}
