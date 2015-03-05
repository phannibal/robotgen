/*
 * CircleBlock.cpp
 *
 *  Created on: 14.02.2015
 *      Author: philipphannibal
 */

#include "CircleBlock.h"

int CircleBlock::amount = 0;

CircleBlock::CircleBlock() {
	X = Y = 0;
	R = 1;
	R2 = R*R;

	CircleBlock::amount += 1;
}

CircleBlock::~CircleBlock() {
	CircleBlock::amount -= 1;
}

void CircleBlock::setParameters(double aX, double aY, double aR) {
	X=aX;
	Y=aY;
	R=aR;
	R2= R*R;
}

