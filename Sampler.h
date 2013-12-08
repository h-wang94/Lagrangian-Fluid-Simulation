/*
 * Sampler.h
 *
 *  Created on: Dec 7, 2013
 *      Author: Owner
 */

#ifndef SAMPLER_H_
#define SAMPLER_H_

#include "Sample.h"
#include "Camera.h"
#include "Ray.h"
#include "SpatialGrid.h"


class Sampler {
	Sampler();
	Sampler(Camera c);
	bool hasNext();
	Sample next();
	Sample getSample();
	Ray getRayThroughSample();
	Sample getPixelRGBA(SpatialGrid &, float);
private:
	Sample s;
	Camera c;
};



#endif /* SAMPLER_H_ */
