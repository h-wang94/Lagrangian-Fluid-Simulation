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
	unsigned char* getBMP();
	void commitSample(Sample s);

private:
	vector<unsigned char> pixelRGBAs;
	float imgWidth, imgHeight;
};



#endif /* FILM_H_ */
