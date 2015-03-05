/*
 * Individual.h
 *
 *  Created on: 16.01.2015
 *      Author: philipphannibal
 */

#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_

#include <stdlib.h>
//#include <random>
//#include <functional>
#include "EnvironmentSimulator.h"

struct individualPair;
//typedef bool* gen;



class Individual {
public:


	Individual();
	virtual ~Individual();

	//TODO: maybe decide  the following four or five functions on virtual or static:
	static /*virtual*/ double fitness(Individual*);
	static /*virtual*/ int* selection(double* fitnesses);
	static /*virtual*/ Individual* crossover(Individual*);
	static void crossover(Individual*, Individual*);	//actual in use

	static EnvironmentSimulator* mySim;

	static /*virtual*/ Individual* mutation(Individual*);

	virtual void printGenom();
	double* genom;

	static Individual* indpair;
	static int size; 	//number of bits in genom
	static individualPair* ip;	//for push-arounds
	static void setPm(double apm) {
		pm = apm;
	}
	static void setPc(double apc) {
		pc = apc;
	}
	static void setN(int aN) {
		N = aN;
	}
	void setGenom(double*);



private:
	//char* sgenom;	//string of the genom
	//TESTINDIVIDUAL:

	static double pm;
	static double pc;
	static int N;	//Number of created individuals

};

/*
 * STRUCTS
 */

struct individualPair {
	Individual i1;
	Individual i2;
};

#endif /* INDIVIDUAL_H_ */
