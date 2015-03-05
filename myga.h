/*
 * myga.h
 *
 *  Created on: 16.01.2015
 *      Author: philipphannibal
 */



#ifndef MYGA_H_
#define MYGA_H_

#include "Individual.h"
#include <vector>

class myga {
public:
	//Constructors:

	myga(int /* #individuals */, int /* #generations */, double /* pm */, double /* pc */);
	myga();
	//Destructors:
	virtual ~myga();
	//Functions:
	void run();

	//Variables:
	bool running;

private:
	//Functions:
	void prepare();
	void printState();

	//Variables
	int N;	//nummber of Individuals
	int Gmax; //maximum numbor of generation iterations
	double pm;	//chance of mutation
	double pc; //chance of crossover
	std::vector< Individual* > G;	//generation
	std::vector< Individual* > temp_G;	//temporal copy
	double* temp_genoms;
	double* fitnessOf;	//stores the fitness of the individuals of the current generation
	double faverage;



};

#endif /* MYGA_H_ */
