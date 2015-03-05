/*
 * myga.cpp
 *
 *  Created on: 16.01.2015
 *      Author: philipphannibal
 */

#include <iostream>
#include "myga.h"
#include <random>
#include <cmath>
#include <time.h>
#include "Individual.h"

//CONSTRUCTORS:

myga::myga(int aN, int aGmax, double apm, double apc) {	//a for adjust
	std::cout << "#myga is being constructed..." << "\n";
	srand ( time(NULL));	//initialize random numbers
	N = aN; Gmax = aGmax; pm = apm; pc = apc;
	running = false;

	fitnessOf = new double[N];
	temp_genoms = new double[N*Individual::size];
	faverage = 0;

	this->prepare();
}

myga::myga() {
	// TODO Auto-generated constructor stub
	srand ( time(NULL));
	N = 1; Gmax = 1; pm = 0; pc = 0;
	running = false;
	fitnessOf = new double[N];
	faverage = 0;
	temp_genoms = new double[N*Individual::size];

	this->prepare();
}

myga::~myga() {
	// TODO Auto-generated destructor stub
}

/* ====================================================================
 * FUNCTIONS for MYGA.CLASS
 * ====================================================================*/

void myga::run() {
	std::cout << "# TABLE HAS FOLLOWING STRUCTURE:\n";
	std::cout << "# Generation, f_average, f0, w00, w01, w02, w03, f1, w10, w11, w12, w13, f2, w20,..., fN, wN0, wN1, wN2, wN3\n";
	std::cout << "#FITNESSPLOT\n";

	for(int g = 0; g<Gmax ; g++) {
		//EVALUATION:
		for(int i = 0; i<N; i++) {
			//std::cout << "#MYGA: i = " << i << " with " << &(G.at(i)) << " and " << &(G.at(i)->genom) << "\n";
			fitnessOf[i] = Individual::fitness(G.at(i));
		}
			//average fitness:
		faverage = 0;
		for(int i = 0; i<N; i++) {
			faverage += fitnessOf[i];
		}
		faverage = faverage/N;
			//outstream:
		//std::cout << "#Generation " << g <<". Following fitnesses of weights evaluated: \n";
		std::cout << g << " " << faverage << " ";
		for(int i = 0; i< N; i++) {
			std::cout <<  fitnessOf[i] << " " << G.at(i)->genom[0] << " " << G.at(i)->genom[1] << " " << G.at(i)->genom[2] << " " << G.at(i)->genom[3] << " ";
		}
		std::cout << "\n";

		//SELECTION:
		int* parentIndexList = Individual::selection(fitnessOf);
		std::cout << "#Selected parents' indices: ";
		for(int i = 0; i<N; i++) {
			std::cout << parentIndexList[i] << ", ";
		}
		std::cout << "\n";

		//apply that selection:
		//(switch the new generations


		for(int i = 0; i < N; i++) {
			for(int b = 0; b<Individual::size; b++) {
				temp_genoms[i*Individual::size+b] = G.at(parentIndexList[i])->genom[b];
			}
		}

		for(int i = 0; i < N; i++) {
			for(int b = 0; b <Individual::size; b++) {
				G.at(i)->genom[b] = temp_genoms[i*Individual::size+b];
			}
		}

		/*
		temp_G.reserve(N);
		for(int i = 0; i<N; i++) {
			temp_G.push_back(new Individual());
			std::cout << "#MYGAvorher: " << &(G.at(parentIndexList[i])->genom) << "\n";
			temp_G.at(i)->setGenom(G.at(parentIndexList[i])->genom);
			std::cout << "#MYGAnachher: " << &(temp_G.at(i)->genom) << "\n";
		}

		G = temp_G;
		for(int i = 0; i < N; i++) {
			std::cout << "#MYGA" << &(G.at(i)->genom) << " and " << &(G.at(i)) << "\n";
		}
		temp_G.clear();
		std::cout << "#MYGA after temp clear\n";
		for(int i = 0; i < N; i++) {
			std::cout << "#MYGA" << &(G.at(i)->genom) << " and " << &(G.at(i)) << "\n";
		}
		*/


		//CROSSOVER:
		//printState();
		//std::cout << "CROSSOVER\n";
		double rn;
		for(int i = 0; 2*i + 1 < N; i++) {
			rn = rand();
			if(rn <= this->pc*RAND_MAX) {
				Individual::crossover(G.at(2*i), G.at(2*i+1));
			}
		}
		//printState();



		//MUTATION:
		//std::cout << "MUTATION\n";
		for(int i = 0; i<N; i++) {
			Individual::mutation(G.at(i));
		}

		//NEW GENERATION:
		//this->printState();

	}

	std::cout << "#Final result:\n";
	this->printState();
}

void myga::prepare() {
	std::cout << "#prepared..." << "\n";
	Individual::setPc(pc);	//pc and pm are static
	Individual::setPm(pm);
	Individual::setN(N);
	G.clear();
	G.reserve(N);	//allocates storage for N Individual*
	double* allZero = new double[Individual::size];
	for(int b = 0; b<Individual::size; b++) allZero[b] = 0.0;
	for(int i = 0; i<N; i++) {
		G.push_back(new Individual());

		G.at(i)->setGenom(allZero);

		G.at(i)->printGenom();
	}
	std::cout.precision( 4 );
	std::cout << "# +++ myga is ready for calculations +++\n";
	std::cout << "# + the parameters are set to:\n# + Number Individuals N = " << N << "\n# + Chance of mutation pm = " << pm << "\n# + Chance of crossover pc = " << pc << "\n# + Calculated Generations Gmax = " << Gmax<< "\n";
}

void myga::printState() {
	Individual::mySim->setTrace(true);
	std::cout << "\n\n#The Generation looks this way now (+trace):\n";
	for(int i = 0; i<N; i++) {
		G.at(i)->printGenom();
		std::cout << "#INDIVIDUAL" << i << "\n";
		Individual::fitness(G.at(i));
	}
	Individual::mySim->setTrace(false);
}
