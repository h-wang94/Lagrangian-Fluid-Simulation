/*
 * Sampler.cpp
 *
 *  Created on: Dec 7, 2013
 *      Author: Owner
 */

#include "Sampler.h"

Sampler::Sampler() {
	s = Sample(-1, 0); //On the first iteration of a loop, it will have moved forward to (0, 0)
	c = Camera();
}

Sampler::Sampler(Camera cam) {
	s = Sample(-1, 0);
	c = cam;
}

Sample Sampler::getSample() {
	return s;
}

bool Sampler::hasNext() {
	Sample temp = Sample(s.x, s.y);
	temp.x += 1;
	if (temp.x == c.imgWidth) {
		temp.y += 1;
	}
	if (temp.y == c.imgHeight) {
		return false;
	}
	return true;
}

/* Returns the next sample, provided it exists. Like an iterator. */
Sample Sampler::next() {
	if (!hasNext()) {
		std::cerr << "No next value for Sampler.\n" << std::endl;
		exit(1);
	}
	s.x += 1;
	if (s.x == c.imgWidth) {
		s.x = 0;
		s.y += 1;
	}
	return s;
}

Ray Sampler::getRayThroughSample() {
	Ray temp;
	float px = ((2 * c.R) * (s.x + .5)/c.imgWidth - c.R);
	float py = ((2 * c.T) * (s.y + .5)/c.imgHeight - c.T);
	Vector d = c.xAxis * px + c.yAxis * py + c.zAxis;
	d.normalize();
	temp.setDirection(d);
	temp.setOrigin(c.lookFrom);
	return temp;
}

/* Finds the RBGA value of the current sample for a given spatial grid, taking
 * samples along the ray through the pixel at every stepSize. Returns the now
 * calculated sample.
 *
 * Maybe edit this later so it commits to an array or put this somewhere else or idk. */
Sample Sampler::getPixelRGBA(SpatialGrid &sg, float stepSize) {
	Ray r = getRayThroughSample();
	Particle p, temp;
	vector<Particle> neighbors;
	Sample newSample;
	Point3D testPoint;
	float dist = 0, farClip = 3.0;
	for (int i = 0; dist < farClip && s.opacity < 1.0; i++) { //while the color is still opaque...probably want to set i to something else later btw
		testPoint = r.getPointAtT(stepSize * i);
		p.setPosition(testPoint);
		neighbors = sg.getNeighbors(p);
		if (neighbors.empty() || pSystem.colorFunction(p) < .5) {
			continue; //if there are no neighbors, we're outside the fluid. Probably.
		}
		//find some way to average the opacities and color at this "particle"
		//should these be weighted somehow?
		newSample = Sample();
		for (int j = 0; j < (int) neighbors.size(); j++) {
			temp = neighbors[j];
			newSample.colorVals += temp.getColor();
			newSample.opacity += temp.getOpacity();
		}
		newSample.colorVals = newSample.colorVals / (float) neighbors.size();
		newSample.opacity /= (float) neighbors.size();
		s.addRGBAVals(newSample);
		dist = stepSize * (i + 1);
	}
	return s;
}
