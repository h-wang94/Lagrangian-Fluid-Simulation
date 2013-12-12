/*
 * Film.h
 *
 *  Created on: Dec 7, 2013
 *      Author: Owner
 */

#ifndef FILM_H_
#define FILM_H_

#include "Sampler.h"
#include <vector>

class Film {
public:
	Film();
	Film(float imgW, float imgH);
	float* getBMP();
	void commitSample(Sample s);
	void renderScene(Camera &c, float stepSize);

private:
	vector<float> pixelRGBAs;
	float imgWidth, imgHeight;
};



#endif /* FILM_H_ */
