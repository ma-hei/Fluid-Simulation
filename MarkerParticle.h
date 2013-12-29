/*
 * MarkerParticle.h
 *
 *  Created on: 12.08.2013
 *      Author: martenheidemeyer
 */

#ifndef MARKERPARTICLE_H_
#define MARKERPARTICLE_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
class MarkerParticle {
public:
	MarkerParticle();
	virtual ~MarkerParticle();
	void setPosition(float x, float y);
	void setVelocity(float x, float y);
	void setcolor(float a, float b, float c);
	float* getPosition();
	float* getDisplayPosition();
	float* getVelocity();
	float* getcolor();
	float* color;

private:
	float* velocity;
	float* position;

};

#endif /* MARKERPARTICLE_H_ */
