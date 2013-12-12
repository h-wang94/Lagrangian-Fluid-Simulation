/*
 * Sample.cpp
 *
 *  Created on: Dec 7, 2013
 *      Author: Owner
 */

#include "Sample.h"

Sample::Sample() {
	x = 0;
	y = 0;
	colorVals = Vector();
	opacity = 0;
}

Sample::Sample(float a, float b) {
	x = a;
	y = b;
	colorVals = Vector();
	opacity = 0;
}


void Sample::addRGBAVals(Sample sample) {
	colorVals += sample.colorVals * (1 - opacity);
	opacity += (1 - opacity) * sample.opacity;
}
