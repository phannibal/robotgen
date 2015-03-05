/*
 * Individual.cpp
 *
 *  Created on: 16.01.2015
 *      Author: philipphannibal
 */

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <random>
#include <time.h>
//#include <functional>
//#include <algorithm>	//sort
#include "Individual.h"
#include "EnvironmentSimulator.h"


double Individual::pm = 0;
double Individual::pc = 0;
int Individual::N = 0;
int Individual::size = 4;	//size = #weights
EnvironmentSimulator* Individual::mySim = new EnvironmentSimulator(30,0.01,10,0);	//..or(T, dt, L, szenario); are set here
individualPair* Individual::ip;
Individual* Individual::indpair = new Individual[2];

Individual::Individual() {
	//ip = new individualPair();
	//std::cout << "Individual is being constructed..." << "\n";

	//zero-set initiation of the size-bit-genom
	genom = new double[size];
	double rz = 0.0;
	for(int i = 0; i < size; i++) {
		//rz = ((double) rand())/RAND_MAX;
		genom[i] = rz;
	}

}

Individual::~Individual() {
	// TODO Auto-generated destructor stub
}

/*
 * FUNCTIONS for INDIVIDUAL.CLASS
 */

double Individual::fitness(Individual* ri) {	//ri for "reference on individual"
	double f=1;
	//std::cout << "#IND: " << &(ri->genom) << "\n";
	mySim->changeRobotgene( ((ri)->genom) );
	f = mySim->run();
	return f;
}


int* Individual::selection(double* rfitnessOf) {	//r for "reference on"
	int* parentIndexList = new int[Individual::N];
	//sum of all fitnesses:
	double Sigma =0;
	for(int i = 0; i<N; i++) {
		Sigma += rfitnessOf[i];
	}
	//borders {0...1} for fitness proportional selection:
	double* selectionBorder = new double[Individual::N];
	selectionBorder[0] = rfitnessOf[0] / Sigma;
	for(int i = 1; i<N; i++) {
		selectionBorder[i] = rfitnessOf[i] / Sigma + selectionBorder[i-1];
		//std::cout << selectionBorder[i] << ", ";
	}

	double rn; //random number
	bool selected = false;
	for(int i = 0; i<N; i++) {	//searches for N new parents
		selected = false;
		rn = rand();
		for(int j = 0; j<N && !selected; j++) {	//runs over all existing individuals' border widths
			if(rn <= selectionBorder[j] * RAND_MAX) {
				parentIndexList[i] = j;
				selected = true;
			}
		}
	}

	return parentIndexList;
}

/* Individual* Individual::crossover(Individual* rip) {	//rip for "reference on individual pair"
	//single-point-crossover:

	int rn = rand() % size;	//defines crossover cut position
	//std::cout <<"Crossover Called! Cut at rn=" << rn << "\n";
	for(int i = rn; i<size; i++) { //switch the bit values from the rn-th to the last bit
		bool temp_b = rip[0].genom[i];
		rip[0].genom[i] = rip[1].genom[i];
		rip[1].genom[i] = temp_b;
	}

	return rip;
}
*/

void Individual::crossover(Individual* p1, Individual* p2) {
	double* child1Genom = new double[Individual::size];
	double* child2Genom = new double[Individual::size];
	int rn = rand() % size;	//defines crossover cut position
	//std::cout <<"#Crossover Called! Cut at rn=" << rn << "\n";
	for(int i = 0; i<rn; i++) {
		child1Genom[i] = p1->genom[i];
		child2Genom[i] = p2->genom[i];
	}
	for(int i = rn; i<size; i++) { //switch the gene values from the rn-th to the last gene
		child1Genom[i] = p2->genom[i];
		child2Genom[i] = p1->genom[i];
	}
	p1->setGenom(child1Genom);
	p2->setGenom(child2Genom);
	//std::cout << "#MYGA parents' genoms: " << &(p1->genom) << ", "<< &(p2->genom) << "\n";
}

Individual* Individual::mutation(Individual* ri) {	//ri for "reference on individual"
	double rn; //random number
	double rz= 0.0; //another random number, to be added to genom[i]
	for(int i = 0; i < size; i++) {
		rn = rand();
		if(rn <= Individual::pm * RAND_MAX) {
			rz = ((double) rand())/(double)RAND_MAX;//random number between 0...1
			ri->genom[i] += rz-0.5;	//random number between -0.5...0.5
			//std::cout << "#sth mutated!! (w_"<<i<<" by value "<< rz-0.5 << ")\n";
		}
	}

	return ri;
}

void Individual::printGenom() {
	std::cout << "#";
	for(int i = 0; i<size; i++) {
		std::cout << genom[i] << "|";
	}
	std::cout << std::endl;
}

void Individual::setGenom(double* aGenom) {
	for(int i = 0; i<size; i++) {
		genom[i] = aGenom[i];
	}
}


