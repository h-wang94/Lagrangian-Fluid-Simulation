/*
 * Sample.h
 *
 *  Created on: Dec 7, 2013
 *      Author: Owner
 */

#ifndef SAMPLE_H_
#define SAMPLE_H_

#include "Vector.h"

struct Sample {
	Sample();
	Sample(float x, float y);
	float x, y;
	float opacity;
	Vector colorVals;

	void addRGBAVals(Sample s);
};



#endif /* SAMPLE_H_ */
