/*
 * EnvironmentSimulator.cpp
 *
 *  Created on: 14.02.2015
 *      Author: philipphannibal
 */

#include "EnvironmentSimulator.h"
#include "RoundRobot.h"
#include "CircleBlock.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>//usleep
#include <fstream>//filestream

EnvironmentSimulator::EnvironmentSimulator(int aT, double aDt, int aLength, int aScenario) {
	// TODO Auto-generated constructor stub
	debug = false;
	trace = false;
	online = false;
	Lx = Ly = aLength;	//creates squared area of edgelength Lx
	L2 = Lx*Ly;
	T = aT;
	dt = aDt;
	v_temp = phi_temp = 0;
	ctc = 0;
	stepMax = (int) floor(T/dt);
	commonSensorReach = Lx;
	Rob = new RoundRobot();

	delta = new double[1];
	colliding = new bool[1];
	wallcolliding = new bool[4]; //0 for SOUTH, 1 for EAST, 2 for NORTH, 3 for WEST, their walldeltas are fix
	walldelta = new double[4];
	walldelta[3] = M_PI;	//180°
	walldelta[0] = 1.5*M_PI;	//270°
	walldelta[1] = 0;	//0°
	walldelta[2] = 0.5*M_PI;	//90°
	visited = new bool[L2];
	for(int area = 0; area<L2; area++) {
		visited[area] = false;
	}
	blocks.clear();
	ctemp = false;
	cNoise = 0.001*Lx;

	scenario = aScenario;
	setupready = false;

	std::cout << "#ENVIRONMENT SIMULATION IS SET UP. PARAMETERS: " << T << ", " << dt << ", " << Lx << "\n";

}

EnvironmentSimulator::~EnvironmentSimulator() {
	// TODO Auto-generated destructor stub
}

/*
 * Functions
 */

double EnvironmentSimulator::run() {	//shall return fitness f of momentaneaous robotgene
	double f =1.0;
	//prepare run:
	//depending on scenario:
	if(scenario == 0 || scenario < 0) {// || scenario > 0) {	//"Koh's problem"
		//generate blocks,
		blocks.clear();
		this->generateBlock(3,7,1);
		this->generateBlock(7,3,1);

		//setup Rob (gene, position, visited)
		Rob->setX(2.0);	Rob->setY(2.0);
		if(trace && !online) std::cout << Rob->getX() << " " << Rob->getY() << "\n";
		Rob->setR(0.4);	Rob->setTheta(M_PI*0.25);
		Rob->setV(1.0); Rob->setW(0.0);
		Rob->setPhi(M_PI*0.25);
		//Rob->printNetwork(0);
	}
	/* this is bugged due to some reason*/
	if(scenario == 1) {	//grid
		//setup blocks and arena
		//std::cout << "#OH NO A BUG\n";
		if(!setupready) {
			if(debug) std::cout << "#SET UP SCENARIO 1\n";
			blocks.clear();
			Lx = 2*Lx;
			Ly = 2*Ly;
			for(int bx = 0; bx <= Lx; bx = bx + 4) {
				for(int by = 0; by <= Ly; by = by + 4) {
					if(abs(bx+by%3-Lx/2) >=2 && abs(by-Lx/2) >=2 ) {
						double rr = ((double)rand())/RAND_MAX;
						this->generateBlock(bx+by%3,by,rr*2);
					}
				}
			}
			L2 = Lx * Ly;
			visited = new bool[L2];
			for(int area = 0; area<L2; area++) {
				visited[area] = false;
			}
			setupready = true;
		}

		//setup Rob
		Rob->setR(0.5); Rob->setX(Lx/2); Rob->setY(Ly/2);
		//if(trace) std::cout << Rob->getX() << " " << Rob->getY() << "\n";
		double rn = 0.125;
		Rob->setPhi(M_PI*2*rn); Rob->setTheta(M_PI*0.25);
		Rob->setV(1.0); Rob->setW(0.0);
	}
	/*end of probably bugged part*/
		//reset collision time counter and area checker
		ctc = 0;
		for(int area = 0; area<L2; area++) {
			visited[area] = false;
		}


		//std::cout << "#START TIME LOOP\n";
	if(online) {
		//usleep(1*1000000);//to ensure vis has started
		setVisOptions();

		//myfile.open ("online.dat", std::ios::out | std::ios::trunc );//app for appending
	}
	//START OF time loop of simulation:
	for(int step = 0; step < stepMax; step++) {
		//wait for real time, if online, and update the online-trace each step:
		if(online) {
			//system("gnuplot ~/workspace/robotgen/Debug/online.load");
			//myfile << Rob->getX() << " " << Rob->getY() << "\n";
			//usleep(dt*1000000);
		}
		//check on sensors
		for(int i = 0; i<2; i++) {	//Rob is limited to two sensors, loops sensors i
			double rho = Rob->getPhi() + (i==0? Rob->getTheta() : -Rob->getTheta());	//sensor direction angle
			double sigma = Rob->getR();
			double sensorStepLength = 0.1; //in m. 10 asks per meter
			bool hit = false;
			double sx, sy;
			for(int s = 0; s<(double)Rob->getSensorReach()/sensorStepLength && !hit; s++) {//loops sensor steps s
				sx = Rob->getX()+cos(rho)*(sigma+s*sensorStepLength);
				sy = Rob->getY()+sin(rho)*(sigma+s*sensorStepLength);
				//check blocks:
				//if(debug) std::cout << "#RIGHT before block loop, blocksize=" << blocks.size() << "\n";
				for(int b=0; b<(int) blocks.size(); b++) {	//loops blocks b
					//if(debug) std::cout << "#INSIDE LOOP: b=" << b << ", " << sx << "," << sy << "\n";
					if(doesContain( blocks.at(b), sx, sy)) {
						hit = true;
						Rob->setU(i,(double)1.0/(sigma+s*sensorStepLength) );	//sensor u value is set here
						if(debug) std::cout << "#SENSOR " << i << " REGISTERED BLOCK. u = " << ((double)1.0)/(sigma+(double)s*sensorStepLength) << "\n";
						//if(debug) std::cout << Rob->getU(i);
					}
				}
				//check walls:
				if(sx <= 0 || sy <= 0 || sx >= Lx || sy >= Ly) {
					if(debug) std::cout << "#SENSOR " << i << " REGISTERED WALL. u = " << ((double)1.0)/(sigma+(double)s*sensorStepLength) << ", phi = " << Rob->getPhi() << "\n";
					hit = true;
					Rob->setU(i,(double)1.0/(sigma+s*sensorStepLength) );
				}

			}
			//if(debug && hit) std::cout << "#SOMETHING WAS hit.\n";
			if(!hit) {
				Rob->setU(i,0);	//the offset of the sensors is chosen here
				if(debug) std::cout << "#SENSOR " << i << " REGISTERED NO BLOCK.\n";
			}

		}

		//calculate Rob's v, phi,
		Rob->runANNetworkStep();
		Rob->calcMovementparameters(dt);//regards w and thus phi
		v_temp = Rob->getV();
		phi_temp = Rob->getPhi();

		if(debug) std::cout << "#temps: " << v_temp << ", " << phi_temp << ";sin=" << sin(phi_temp) << "\n";

		//check on collision. casually adjust v_temp phi_temp
			//check collision on blocks:
		for(int i = 0; i<(int) blocks.size(); i++) {
			colliding[i] = false;
			if( doCollide(Rob, blocks.at(i)) ) {
				//what happens if block i and Rob collide
				colliding[i] = true;
				delta[i] = angle(Rob->getX(), Rob->getY(), blocks.at(i)->X, blocks.at(i)->Y);
				if(Rob->getX() > blocks.at(i)->X) {
					delta[i] += M_PI;	//property of atan() and my geometry: if robot is right of obstacle: add 180°
				}
				if(debug) std::cout << "#COLLISION WITH BLOCK" << i << "\n";
				if(debug) std::cout << "#ANGLE = " << delta[i] << ", phi_temp - delta = " << phi_temp - delta[i] << "\n";
			}
			//all the time: recalc the angle:
			//delta[i] = angle(Rob->getX(), Rob->getY(), blocks.at(i)->X, blocks.at(i)->Y);

		}
			//check collision on walls:, the walls' deltas are fixed in constructor
		ctemp = false;	//temp variable to increase ctc

		for(int i = 0; i<4; i++) {
			wallcolliding[i] = false;
		}
		if(Rob->getX() - Rob->getR() + cNoise < 0 ) {
			wallcolliding[3] = true;//WEST wall
			if(debug) std::cout << "#Collision on WEST wall\n";
			if(!ctemp) {
				ctemp = true;
				ctc += dt;
			}
		}
		if(Rob->getX() + Rob->getR() - cNoise > Lx) {
			wallcolliding[1] = true;//EAST wall
			if(debug) std::cout << "#Collision on EAST wall\n";
			if(!ctemp) {
				ctemp = true;
				ctc += dt;
			}
		}
		if(Rob->getY() - Rob->getR() +cNoise < 0) {
			wallcolliding[0] = true;//SOUTH wall
			if(debug) std::cout << "#Collision on SOUTH wall\n";
			if(!ctemp) {
				ctemp = true;
				ctc += dt;
			}
		}
		if(Rob->getY() + Rob->getR() -cNoise> Ly) {
			wallcolliding[2] = true;//NORTH wall
			if(debug) std::cout << "#Collision on NORTH wall\n";
			if(!ctemp) {
				ctemp = true;
				ctc += dt;
			}
		}

			//adjust:
		for(int i = 0; i<(int) blocks.size(); i++) {	//loop blocks
			if(colliding[i]) {
				if(!ctemp){
					ctc += dt;
					ctemp = true;
				}
				if(!isStruckAtBlock(i)) {
					if(cos(phi_temp-delta[i])>0) {
						v_temp = fabs(v_temp*sin(phi_temp-delta[i]));
						if(phi_temp-delta[i] != 0) phi_temp = delta[i] + 0.5*M_PI * (phi_temp-delta[i])/(fabs(phi_temp-delta[i]));
						if(debug) std::cout << "#DELTAi = " << delta[i] << "\n";
					}
				} else {
					v_temp = 0+cNoise;
					phi_temp = Rob->getPhi();
					if(debug) std::cout << "#NOISE PARAMS SET.\n";
				}
			}
		}
		//TODO check here
		for(int i = 0; i < 4; i++) {	//loop walls
			if(wallcolliding[i]) {
				if(!isStruckAtWall(i) ){
					if(cos(phi_temp-walldelta[i])>0) {
						v_temp = fabs(v_temp*sin(phi_temp-walldelta[i]));
						if(phi_temp-walldelta[i] != 0) phi_temp = walldelta[i] + 0.5*M_PI * (phi_temp-walldelta[i])/(fabs(phi_temp-walldelta[i]));
						if(debug) std::cout << "#WALL" << i <<" V AND PHI ADJUSTED: " << v_temp << ", " << phi_temp << "\n";
					}
				} else {
					v_temp = 0+cNoise;
					phi_temp = Rob->getPhi();
					if(debug) std::cout << "#NOISE PARAMS SET.\n";
				}
			}
		}

		if(debug) std::cout << "#temps(after correction): " << v_temp << ", " << phi_temp << ";sin=" << sin(phi_temp) << "\n";

		//move Rob, print coordinates
		Rob->setX(Rob->getX()+ dt * v_temp * cos(phi_temp));
		Rob->setY(Rob->getY()+ dt * v_temp * sin(phi_temp));
		if(online && step % 10 == 0) {
			std::cout << "f3x(t) = " << Rob->getX() << " + " << Rob->getR() << "* cos(t)\nf3y(t) = " << Rob->getY() << " + " << Rob->getR() << "* sin(t)\n";
			std::cout << "plot '-',f3x(t), f3y(t) lt 1 title 'point agent'\n";//fx(t),fy(t) lt -1 title 'Obstacle', f2x(t),f2y(t) lt -1 title 'Obstacle',
			std::cout << Rob->getX() << " " << Rob->getY() << "\n";
			std::cout << Rob->getX()+Rob->getR()*cos(Rob->getPhi()) << " " << Rob->getY()+Rob->getR()*sin(Rob->getPhi()) << "\ne\n";
			Rob->printActivity();
		}
		if(trace && !online) {
			if(step % 10/*k*/ == 0) std::cout << Rob->getX() << " " << Rob->getY() << "\n";	//print every k-th step
		}

		if(debug) std::cout << Rob->getX() << " " << Rob->getY() << " "<< Rob->getV() <<"\n";


		//check on visited
		visited[(int) Rob->getX() +Ly* (int) Rob->getY()] = true;//visited[index of current field] = true;
		//if(trace) std::cout << "#ON FIELD: "<< (int) Rob->getX() +Ly* (int) Rob->getY() << "\n";
		//visited[abs((int) Rob->getX() +Ly* (int) Rob->getY())%L2] = true;
	}//END OF time loop of simulation
	//if(online) myfile.close();

	//check for visited, calculate fitness
	if(debug || trace) std::cout << "#";
	for(int area = 0; area<L2; area++) {
		if(visited[area]) {
			f++;
			if(debug || trace) std::cout << area <<  ",";
		}
	}

	//only time not colliding counts:
	//f = f * (1-ctc) / (double) T;
	if(debug || trace) std::cout << "#CTC PERCENTAGE IS: ctc/T = " << ctc <<"/"<< T << "=" << ctc/T << "\n";
	f = f* (1.0-ctc/T);
	if(debug || trace) std::cout << "#FITNESS: " << f << "\n\n\n";//2 empty rows for easier plotting with index "INDIVIDUALi"

	return f;
}//END: run();

void EnvironmentSimulator::changeRobotgene(double* aRobotgene) {
	//std::cout << "#CHANGING ROBOTGENE\n";
	//std::cout << "#aRobotgene = " << &aRobotgene << "\n";
	Rob->setGweights(aRobotgene);
	//std::cout << "#ROBOTGENE CHANGED\n";
}

void EnvironmentSimulator::generateBlocks(int amount) {
	blocks.clear();
	blocks.reserve(amount);
	for(int i = 0; i< amount; i++) {
		blocks.push_back(new CircleBlock);
		blocks.at(i)->setParameters(0,0,5);	//X,Y,R
	}
}

void EnvironmentSimulator::generateBlock(double aX, double aY, double aR) {
	blocks.push_back(new CircleBlock);
	blocks.back()->setParameters(aX, aY, aR);
	if(trace) std::cout << "#GENERATED BLOCK AT (x, y): (" << aX << ", " << aY << ")\n";
	colliding = new bool[(int) blocks.size()];
	delta = new double[(int) blocks.size()];
}

bool EnvironmentSimulator::doCollide(RoundRobot* rRob, CircleBlock* rBlock) {	//r for "reference on
	bool collision = false;
	//double delta = atan( (rRob->getY()-rBlock->Y) / (rRob->getX()-rBlock->X) );
	double distance = sqrt( pow(rRob->getY()-rBlock->Y,2) + pow(rRob->getX()-rBlock->X,2) );
	if(distance + cNoise < rRob->getR() + rBlock->R) {
		collision = true;
		//if(debug) std::cout << "#COLLISION WITH BLOCK.\n";//TODO: check, if you want to push back the robot on collision
	}

	return collision;
}

bool EnvironmentSimulator::doesContain(CircleBlock* rBlock, double x, double y) {
	bool contain = false;
	double squaredDistance = pow(x-rBlock->X,2) + pow (y-rBlock->Y,2);
	contain = (squaredDistance <= rBlock->R2);
	//if(debug && contain) std::cout << "#DOESCONTAIN() was called.\n";
	return contain;
}

double EnvironmentSimulator::angle(double x1, double y1, double x2, double y2) { //calculates the angle of line P1(x1,y1) to P2(x2,y2) to x-axis
	return atan((y2-y1)/(x2-x1));
}

bool EnvironmentSimulator::isStruckAtBlock(int index) {	//only called on collision with index. checks if others cause struck
	bool struck = false;
	for(int i = 0; i < index; i++) {
		if(colliding[i]) {
			if( fabs(delta[i]-delta[index]) >= M_PI * 0.49) {//angle of >90° is "spitz" in view of Rob and causes struck, use 49 instead of exact 0.5
				struck = true;
			}
		}
	}
	return struck;
}

bool EnvironmentSimulator::isStruckAtWall(int index) {
	bool struck = false;
	for(int i = 0; i <(int) blocks.size(); i++) {	//loop blocks
		if(colliding[i]) {
			if( fabs(delta[i]-walldelta[index]) >= M_PI * 0.5) {//angle of >90° is "spitz" in view of Rob and causes struck
				struck = true;
			}
		}
	}
	for(int i = 0; i < 4; i++) {
		if(i != index) {
			if(wallcolliding[i]) {
				//if( fabs(walldelta[i]-walldelta[index]) > M_PI * 0.5) {//angle of >90° is "spitz" in view of Rob and causes struck
					struck = true;
				//}
			}
		}
	}
	return struck;
}

/*
 * SETS AND GETS
 */

void EnvironmentSimulator::setTrace(bool aTrace) {
	trace = aTrace;
}

void EnvironmentSimulator::setOnline(bool aOnline) {
	online = aOnline;
	Rob->setOnline(aOnline);
}

void EnvironmentSimulator::setVisOptions() {	//this configures the first gnuplot load file options
	//usleep(1*1000000);
	std::cout << "reset\n";
	std::cout << "set print '-'\n";
	//std::cout << "set multiplot\n";
	std::cout << "set xrange [0:"<< Lx <<"]\n";
	std::cout << "set yrange [0:"<< Ly <<"]\n";
	std::cout << "set size ratio 1\n";
	std::cout << "set xlabel 'x-axis (in [m])'\n";
	std::cout << "set ylabel 'y-axis (in [m])'\n";
	std::cout << "set grid\n";
	std::cout << "set key outside bottom\n";
	std::cout << "set param\n";
	std::cout << "set trange [0:2*pi]\n";
	std::cout << "set tmargin 1\n";
	/*
	std::cout << "fx(t) = 4+cos(t)\n";
	std::cout << "fy(t) = 6+sin(t)\n";
	std::cout << "f2x(t) = 6+cos(t)\n";
	std::cout << "f2y(t) = 4+sin(t)\n";
	*/
	for(int b = 0; b <(int) blocks.size(); b++) {
		std::cout << "set object circle at " << blocks.at(b)->X << "," << blocks.at(b)->Y << " size " << blocks.at(b)->R << "\n";
	}
	//std::cout << "plot fx(t),fy(t) lt -1 title 'Obstacle', f2x(t),f2y(t) lt -1 title 'Obstacle'\n";

	Rob->printNetwork(0);
	std::cout << "print '0 MODE random'" << std::endl;
	std::cout << "print '0 AMOD opacity 30'" << std::endl;
	std::cout << "print '0 NLBL on'" << std::endl;
	std::cout << "print '0 MODE ff 2 2 2'" << std::endl;

	std::cout << "pause -1" << std::endl;

}
