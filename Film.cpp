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
	pixelRGBAs.push_back(s.colorVals.getX());
	pixelRGBAs.push_back(s.colorVals.getY());
	pixelRGBAs.push_back(s.colorVals.getZ());
	pixelRGBAs.push_back(s.opacity);
}

float* Film::getBMP() {
	return &(pixelRGBAs[0]);
}

void Film::renderScene(Camera &c, float stepSize) {
	Sampler sampler = Sampler(c);
	Sample s;
	while (sampler.hasNext()) {
		sampler.next();
		if (sampler.getSample().x < 300 && sampler.getSample().y > 100) {
			commitSample(Sample());
			continue;
		}
		s = sampler.getPixelRGBA(stepSize);
		commitSample(s);
	}
}
