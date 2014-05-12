/*
 * Grid2.h
 *
 *  Created on: 07.12.2013
 *      Author: martenheidemeyer
 */

#include "MarkerParticle.h"

#ifndef GRID2_H_
#define GRID2_H_

class Grid2 {
private:
	int cellSize;
	int displaySize;
	float timestep;
	float* vValues;
	float* advecteduValues;
	float* advectedvValues;
	float* uPositionsx;
	float* uPositionsy;
	float* vPositionsx;
	float* vPositionsy;
	float* uValues;
	float** pressure;
public:

	int gridSize;
	Grid2(int gridSize);
	virtual ~Grid2();
	void advectAll();
	void init();
	void printAlluValues();
	void printAllvValues();
	void printAllAdvecteduValues();
	void printPressure();

	void getVelocityAtVerticalCellFace(int i, int k, float* oldVelocity);
	void getVelocityAtHorizontalCellFace(int i, int k, float* oldVelocity);
	void getOldPositionOfParticleAtVerticalCellFace(int i, int k, float* oldPosition);
	void getOldPositionOfParticleAtHorizontalCellFace(int i, int k, float* oldPosition);

	float* getVelocityAtPosition(float x, float y, int verbose);
	void getVelocityAtPosition2(float x, float y, float* velo);
	void getVelocityAtCellPosition(int x, int y, int which, float* velocity);
	void getCellOfPosition(float x, float y, int* cell);
	int getSectorOfPosition(float leftFacex, float lowerFacey, float x, float y);
	void solvePressure();
	void updateVelocities();
	void moveMarkers();
	void setu(int i, int k, float value);
	void setv(int i, int k, float value);
	int getSize();
	MarkerParticle getMarker(int i);

	MarkerParticle* markers;

};

#endif /* GRID2_H_ */
