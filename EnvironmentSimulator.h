/*
 * EnvironmentSimulator.h
 *
 *  Created on: 14.02.2015
 *      Author: philipphannibal
 */

#ifndef ENVIRONMENTSIMULATOR_H_
#define ENVIRONMENTSIMULATOR_H_


#include "RoundRobot.h"
#include "CircleBlock.h"
#include <vector>
#include <fstream>

class EnvironmentSimulator {
public:
	EnvironmentSimulator(int aT/*Simulation time span*/, double aDt/*Simulation step time length*/, int aL/*arena size*/, int scenario/*chose arena*/);
	virtual ~EnvironmentSimulator();

	double run();
	void changeRobotgene(double* robotgene);
	void generateBlocks(int howmany);
	void generateBlock(double X, double Y, double R);
	bool doCollide(RoundRobot*, CircleBlock*);
	bool doesContain(CircleBlock*, double, double);	//block + coordinates x, y
	double angle(double, double, double, double);
	bool isStruckAtBlock(int indexOfBlock);
	bool isStruckAtWall(int indexOfWall);
	void setTrace(bool);
	void setOnline(bool);
	void setVisOptions();

private:
	int Lx;
	int Ly;
	int L2;
	double dt;
	int T;
	int stepMax;
	std::vector< CircleBlock* > blocks;
	RoundRobot* Rob;
	bool* visited;	//counts the visited areals
	int commonSensorReach;
	int scenario;	//0: Koh's problem, 1:
	bool setupready;
	//Temporarily used params for robot movement (on collision)
	double v_temp;
	double phi_temp;
	double ctc;
	bool debug;
	bool trace;
	bool online;
	//collision params:
	double* delta;
	bool* colliding;	//for blocks
	bool* wallcolliding;	//for walls: 0 for SOUTH, 1 for EAST, 2 for NORTH, 3 for WEST
	double* walldelta;
	bool ctemp;
	double cNoise;
	//stream
	std::ofstream myfile;

};

#endif /* ENVIRONMENTSIMULATOR_H_ */
