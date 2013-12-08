/*
 * Film.cpp
 *
 *  Created on: Dec 7, 2013
 *      Author: Owner
 */

#include "Film.h"

Film::Film() {
	imgWidth = 640;
	imgHeight = 480;
}

Film::Film(float x, float y) {
	imgWidth = x;
	imgHeight = y;
}

void Film::commitSample(Sample s) {
	pixelRGBAs.push_back((__int8) s.colorVals.getX());
	pixelRGBAs.push_back((__int8) s.colorVals.getY());
	pixelRGBAs.push_back((__int8) s.colorVals.getZ());
	pixelRGBAs.push_back((__int8) s.opacity);
}

__int8* Film::getBMP() {
	return &(pixelRGBAs[0]);
}
