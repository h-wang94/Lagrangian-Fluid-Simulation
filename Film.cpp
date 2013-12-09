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
	pixelRGBAs.push_back((unsigned char) s.colorVals.getX());
	pixelRGBAs.push_back((unsigned char) s.colorVals.getY());
	pixelRGBAs.push_back((unsigned char) s.colorVals.getZ());
	pixelRGBAs.push_back((unsigned char) s.opacity);
}

unsigned char* Film::getBMP() {
	return &(pixelRGBAs[0]);
}

void Film::renderScene(Camera &c, SpatialGrid &sg, float stepSize) {
	Sampler sampler = Sampler(c);
	Sample s;
	while (sampler.hasNext()) {
		sampler.next();
		s = sampler.getPixelRGBA(sg, stepSize);
		commitSample(s);
	}
}
