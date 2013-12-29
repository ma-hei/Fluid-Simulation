/*
 * MarkerParticle.cpp
 *
 *  Created on: 12.08.2013
 *      Author: martenheidemeyer
 */

#include "MarkerParticle.h"


MarkerParticle::MarkerParticle() {

	velocity=(float*)malloc(2*sizeof(float));
	position=(float*)malloc(2*sizeof(float));
	color=(float*)malloc(3*sizeof(float));

	velocity[0]=0;
	velocity[1]=0;

	position[0]=0;
	position[1]=0;

}

MarkerParticle::~MarkerParticle() {
	// TODO Auto-generated destructor stub
}

void MarkerParticle::setPosition(float x, float y) {

	this -> position[0]=x;
	this -> position[1]=y;


}

void MarkerParticle::setVelocity(float x, float y) {

	this -> velocity[0] = x;
	this -> velocity[1] = y;
}


void MarkerParticle::setcolor(float a, float b, float c){

	this->color[0]=a;
	this->color[1]=b;
	this->color[2]=c;

}

float* MarkerParticle::getcolor(){
	return this->color;
}

float* MarkerParticle::getPosition() {

	return this->position;
}


float* MarkerParticle::getVelocity() {

	return this->velocity;
}
