/*
 * RoundRobot.h
 *
 *  Created on: 14.02.2015
 *      Author: philipphannibal
 */

#ifndef ROUNDROBOT_H_
#define ROUNDROBOT_H_

class RoundRobot {
public:
	RoundRobot();
	virtual ~RoundRobot();

	void runANNetworkStep();
	void calcMovementparameters(double dt);
	double sigmoid(double);

	void printNetwork(int iNetwork=0);
	void printState();
	void printActivity(int iNetwork = 0);


	void setGweights(double*);


	void setX(double);
	void setY(double);
	void setPhi(double);
	void setV(double);
	void setW(double);
	void setTheta(double);
	void setR(double);
	void setU(int, double);
	void setOnline(bool);
	double getX();
	double getY();
	double getPhi();
	double getV();
	double getW();
	double getTheta();
	double getR();
	double getA(int);
	double getU(int);
	double getSensorReach();





	int genes;



	bool crash;

private:
	double* gweights; //genomial weights;
	double* u;//0, 1 are sensors; 2,3 are geneweights; 4,5 are output
	double* a;
	double* W;//W[FROM + 6* TO]
	double x;
	double y;
	double v;	//heading speed
	double phi;	//heading direction angle
	double w;	//omega, change of phi
	double theta;
	double r;
	double sensorReach;

	bool debug;
	bool trace;
	bool online;





};

#endif /* ROUNDROBOT_H_ */
