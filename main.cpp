/*
 * main.cpp
 *
 *  Created on: 16.01.2015
 *      Author: philipphannibal
 */


#include <iostream>
#include "myga.h"
#include "Individual.h"
#include "EnvironmentSimulator.h"

using namespace std;



int main() {

/*
	myga* myTestGa = new myga(24,1000,0.05,0.7);
	myTestGa->run();
*/

/**/
	std::cout << "#MAIN###################### START SINGLE ENVIRONMENT SIMULATION ######################\n";
	EnvironmentSimulator* mySim = new EnvironmentSimulator(30,0.01,10,1);//T, dt, L
	double* testgene = new double[4];
	testgene[0] = 0.7096; testgene[1] = -1.102; testgene[2] =-0.6733; testgene[3] = 0.2608;
	//-0.4839---1.831--0---0.6544-- //-1.125-- -0.3455 - 0 - 0
	//-0.384|0.4812|0|0|
	//-0.07776|-0.4787|-0.7312|0.02682| f~=31
	//-0.9311|0.9551|-1.829|1.946| video12
	//0.8564|-1.049|0.001659|-0.1852| IND15 GA8
	//0.7096|-1.102|-0.6733|0.2608|IND6 GA8
	mySim->changeRobotgene(testgene);
	double fitness = 0;
	mySim->setTrace(false);
	mySim->setOnline(true);
	fitness = mySim->run();
	std::cout << "#MAIN:FITNESS: " << fitness << "\n";
/**/


	return 0;
}

