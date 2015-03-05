/*
 * RoundRobot.cpp
 *
 *  Created on: 14.02.2015
 *      Author: philipphannibal
 */


#include "RoundRobot.h"
#include <math.h>
#include <iostream>

RoundRobot::RoundRobot() {
	//Options:
	debug = false;
	trace = false;
	online = false;
	//Scalars:
	v = 1;
	r = 1;
	sensorReach = 2*r;
	phi = 0;
	w = 0;
	x = 0;
	y = 0;
	crash = false;
	genes = 4;	//for four weights w22, w23, w32, w33 (in this direction)
	theta = M_PI * 0.25;	//one eighth of 2pi

	//Matrices:
	gweights = new double[genes];
	for(int i = 0; i<genes; i++) {
		gweights[i] = 0;
	}
	u = new double[6];
	for(int i = 0; i<6; i++) {
		u[i] = 0;
	}
	a = new double[6];
	for(int i = 0; i<6; i++) {
		a[i] = 0.5;
	}
	W = new double[6*6];
	for(int i = 0; i<6; i++) {//TO
		for(int j = 0; j<6; j++) {//FROM
			W[j+6*i] = 0;
		}
	}
	W[0 + 6*2] = 1; W[2 + 6*4] = 1;
	W[1 + 6*3] = 1; W[3 + 6*5] = 1;

	if(debug) printState();



}

RoundRobot::~RoundRobot() {
	// TODO Auto-generated destructor stub
}

/*
 * Functions:
 */


void RoundRobot::runANNetworkStep() {	//does: a = sigmoid(u); u = SUM( Wij * aj);
	//map sensor u to -1...1 TODO: test this option
	u[0] = (u[0]-1) * 2;
	u[1] = (u[1]-1) * 2;
	//calculate a_t from u_t-1
	for(int i = 0; i<6; i++) {
		a[i] = sigmoid(u[i]);
		u[i] = 0; //for later new summation (*)
	}
	//ATTENTION: at this point are all u_i = 0!
	//calculate W * a_t = u_t
	for(int i = 0; i<6; i++) {	//TO
		for(int j = 0; j<6; j++) {	//FROM
			u[i] += W[j + 6*i] * a[j]; //(*) summation here
		}
	}
	//now u_i are recalculated, but u0 and u1 are zero, because sensors will be updated in next step
	//print this
	if(debug) {
		std::cout << "#";
		for(int i = 0; i < 6; i++) {
			std::cout << i << ":";
			std::cout << u[i] << "_" << a[i] << " -- ";
		}
		std::cout << "\n";
	}

}

//ATTENTION this function uses the ANN results...
void RoundRobot::calcMovementparameters(double dt) {//calc v, phi, w
	w = a[5] - a[4];//...here
	if(debug) std::cout << "#NEW OMEGA: w= " << w ;
	//v = 1;//v is constant
	phi = phi +  w;
	if(debug) std::cout << " AND NEW PHI: phi = " << phi << "\n";
}

double RoundRobot::sigmoid(double u_calc) {
	double sig;
	sig = 1.0 / (1.0 + exp(-u_calc));
	return sig;
}
/*
 * PRINTS
 */

void RoundRobot::printNetwork(int iNetwork) {
	if(online) {
		if(iNetwork < 0) iNetwork = 0;
		if(trace) std::cout << "#";
		std::cout << "print '"<< iNetwork << " RNTW data 6";
		for(int i = 0; i<6; i++) {//TO
			std::cout << "'\nprint '";
			if(trace) std::cout << "#";
			for(int j = 0; j<6; j++) {//FROM
				std::cout << W[j+6*i];
				if(j<5) {std::cout << " ";}

			}
		}
		std::cout << "'\n";
	} else {
		//std::cout << gweights[0] << " " << gweights[1] << " " << gweights[2] << " " << gweights[3];
	}
}

void RoundRobot::printActivity(int iNetwork) {
	if(online) {
		if(iNetwork < 0 ) iNetwork = 0;
		if(trace) std::cout << "#";
		std::cout << "print '"<< iNetwork << " ACTL";
		for(int i = 0; i <6; i++) {
			std::cout << " " << a[i];
		}
		std::cout << "'\npause 0.1\n";
	}
}

void RoundRobot::printState() {
	std::cout << "#STATE OF ROB: \n #Position (x, y, phi): (" << x << ", " << y << ", " << phi << ") \n #Radius r = " << r << "\n #Sensorstate: theta = "
			<< theta << ", reach = " << sensorReach << "\n";
}


/*
 * SETs and GETs
 */
void RoundRobot::setGweights(double* aGweights) {
	//set gweights
	//std::cout << "#RR:SET GWEIGHTS";
	for(int i = 0; i<genes; i++) {
		//std::cout << "..." << i ;
		gweights[i] = aGweights[i];

	}
	//std::cout << "\n#RR:ADJUST MATRIX";
	//adjust matrix
	for(int i = 2; i < 4; i++) {//TO
		for(int j = 2; j<4; j++) {//FROM
			//std::cout << "... " << i << ", " << j;
			W[j+6*i] = gweights[2*(i-2) + (j-2)];
		}
	}
	//std::cout << "\n";
}

void RoundRobot::setX(double aX) { x = aX; }
void RoundRobot::setY(double aY){ y = aY; }
void RoundRobot::setPhi(double aPhi){ phi = aPhi; }
void RoundRobot::setV(double aV){ v = aV; }
void RoundRobot::setW(double aW){ w = aW; }
void RoundRobot::setTheta(double aTheta){ theta = aTheta; }
void RoundRobot::setR(double aR){ r = aR; }
void RoundRobot::setU(int index, double aU) {
	if(index >= 0 && index <=6) {
		u[index] = aU;
	}
}
double RoundRobot::getX() { return x;}
double RoundRobot::getY(){ return y;}
double RoundRobot::getPhi(){ return phi;}
double RoundRobot::getV(){ return v;}
double RoundRobot::getW(){ return w;}
double RoundRobot::getTheta(){ return theta;}
double RoundRobot::getR(){ return r;}
double RoundRobot::getA(int index) {
	if(index >=0 && index <= 6) {
		return a[index];
	} else {
		return 0;
	}
}
double RoundRobot::getU(int index) {
	if(index >=0 && index <= 6) {
		return u[index];
	} else {
		return 0;
	}
}
double RoundRobot::getSensorReach() {
	return sensorReach;
}

void RoundRobot::setOnline(bool aOnline) {
	online = aOnline;
}
