/*
 * CircleBlock.h
 *
 *  Created on: 14.02.2015
 *      Author: philipphannibal
 */

#ifndef CIRCLEBLOCK_H_
#define CIRCLEBLOCK_H_

class CircleBlock {
public:
	CircleBlock();
	virtual ~CircleBlock();


	void setParameters(double aX, double aY, double aR);


	static int amount;

	double X;
	double Y;
	double R;
	double R2;

private:

};

#endif /* CIRCLEBLOCK_H_ */
