/*
 * Grid2.cpp
 *
 *  Created on: 07.12.2013
 *      Author: martenheidemeyer
 */

#include "Grid2.h"
#include "pcg_solver.h"
#include <GL/glut.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

Grid2::Grid2(int gridSize) {

	this->gridSize=gridSize;
	timestep=0.8;
	displaySize=60;
	uValues=(float*)malloc((gridSize+1)*gridSize*sizeof(float));
	vValues=(float*)malloc((gridSize+1)*gridSize*sizeof(float));
	uPositionsx=(float*)malloc((gridSize+1)*gridSize*sizeof(float));
	uPositionsy=(float*)malloc((gridSize+1)*gridSize*sizeof(float));
	vPositionsx=(float*)malloc((gridSize+1)*gridSize*sizeof(float));
	vPositionsy=(float*)malloc((gridSize+1)*gridSize*sizeof(float));
	advecteduValues=(float*)malloc((gridSize+1)*gridSize*sizeof(float));
	advectedvValues=(float*)malloc((gridSize+1)*gridSize*sizeof(float));
	pressure = (float**)malloc(gridSize*sizeof(float*));
	markers = (MarkerParticle*)malloc(gridSize*gridSize*4*sizeof(MarkerParticle));
	for (int i=0;i<gridSize;i++)
		pressure[i] = (float*)malloc(gridSize*sizeof(float));

	for (int i=0;i<(gridSize+1)*gridSize;i++){
			uValues[i]=0;
			vValues[i]=0;
			advecteduValues[i]=0;
			advectedvValues[i]=0;
		}
	for (int i=0;i<gridSize;i++)
				for (int k=0;k<gridSize;k++)
						pressure[i][k]=0;
	for (int i=0;i<gridSize*gridSize*4;i++)
			markers[i]=*new MarkerParticle();
	init();

}

void Grid2::printAlluValues(){

	cout<<"u: "<<uValues[39*(gridSize+1)+38]<<"\n";
	cout<<"v: "<<vValues[39*(gridSize+1)+38]<<"\n\n";

}

void Grid2::printAllAdvecteduValues(){

	cout<<"advected u: "<<advecteduValues[39*(gridSize+1)+38]<<"\n";
	cout<<"advected v: "<<advectedvValues[39*(gridSize+1)+38]<<"\n\n";
}

void Grid2::init(){

	float cellWidth = displaySize/gridSize;
	float xCoordinate=0;
	float yCoordinate=0.5*cellWidth;

	for (int i=0;i<gridSize;i++){
		xCoordinate=0;
		for (int k=0;k<gridSize+1;k++){

			uPositionsx[i*(gridSize+1)+k]=xCoordinate;
			uPositionsy[i*(gridSize+1)+k]=yCoordinate;
			vPositionsx[i*(gridSize+1)+k]=yCoordinate;
			vPositionsy[i*(gridSize+1)+k]=xCoordinate;
			xCoordinate+=cellWidth;

		}
		yCoordinate+=cellWidth;
	}

	int counter=0;
	for (int i=0;i<gridSize;i++)
			for (int k=0;k<gridSize;k++){

				markers[counter].setPosition(k+0.5-0.25,i+0.5-0.25);
				markers[counter].setVelocity(0,0);
				counter++;
				markers[counter].setPosition(k+0.5+0.25,i+0.5-0.25);
				markers[counter].setVelocity(0,0);
				counter++;
				markers[counter].setPosition(k+0.5-0.25,i+0.5+0.25);
				markers[counter].setVelocity(0,0);
				counter++;
				markers[counter].setPosition(k+0.5+0.25,i+0.5+0.25);
				markers[counter].setVelocity(0,0);
				counter++;
			}


}

void Grid2::setu(int i, int k, float value){

	uValues[i*(gridSize+1)+k]=value;

}

void Grid2::setv(int i, int k, float value){

	vValues[i*(gridSize+1)+k]=value;

}

void Grid2::advectAll(){

	float* oldPosition=(float*)malloc(2*sizeof(float));
	float* oldVelocity=(float*)malloc(2*sizeof(float));

	for (int i=0;i<gridSize;i++)
		for (int k=0;k<gridSize+1;k++){


			oldPosition=getOldPositionOfParticleAtVerticalCellFace(i, k);
			oldVelocity=getVelocityAtPosition2(oldPosition[0], oldPosition[1]);


			advecteduValues[i*(gridSize+1)+k]=oldVelocity[0];


			oldPosition=getOldPositionOfParticleAtHorizontalCellFace(i, k);
			oldVelocity=getVelocityAtPosition2(oldPosition[0], oldPosition[1]);
			advectedvValues[i*(gridSize+1)+k]=oldVelocity[1];

		}


}

void Grid2::solvePressure(){

	SparseMatrix<float> m(gridSize*gridSize,5);
	vector<float> rhs(gridSize*gridSize);

	int hauptKoeff=4;
	int temp;
	float loesung;
	float indice;


	for (int i=0;i<gridSize;i++)
		for (int k=0;k<gridSize;k++){

			hauptKoeff=4;
			loesung=0;
			indice=i*gridSize+k;

			if (k - 1 < 0)	//links ist solid
				hauptKoeff--;
			else {
				temp = indice - 1;
				m.set_element(indice, temp, -1);
				loesung += advecteduValues[i*(gridSize+1)+k];
			}
			if (i - 1 < 0) //unten ist solid
				hauptKoeff--;
			else {
				temp = (i - 1) * gridSize + k;
				m.set_element(indice, temp, -1);
				loesung += advectedvValues[k*(gridSize+1)+i];
			}
			if (k + 1 > gridSize - 1) //rechts ist solid
				hauptKoeff--;
			else {
				temp = indice + 1;
				m.set_element(indice, temp, -1);
				loesung -= advecteduValues[i*(gridSize+1)+k+1];
			}
			if (i + 1 > gridSize - 1) //oben solid
				hauptKoeff--;
			else {
				temp = (i + 1) * gridSize + k;
				m.set_element(indice, temp, -1);
				loesung -= advectedvValues[k*(gridSize+1)+i+1];
			}
			temp=indice;
			m.set_element(indice,indice,hauptKoeff);
			rhs.at(indice)=loesung;
	}
	PCGSolver<float> solver;
	solver.set_solver_parameters(.0000005,100,0.01,250);
	vector<float> result(gridSize*gridSize);
	float out=0;
	int iter=1;
	solver.solve(m,rhs,result,out,iter);
	int counter = 0;
	int counter2;
	for (int i = 0; i < gridSize; i++) {
		counter2 = i * gridSize;
		for (int k = 0; k < gridSize; k++) {
			counter = counter2 + k;
			pressure[i][k] = (float) result.at(counter);
		}
	}
}

void Grid2::updateVelocities(){

	for (int i=0;i<gridSize;i++)
		for (int k=1;k<gridSize;k++){
			uValues[i*(gridSize+1)+k]=advecteduValues[i*(gridSize+1)+k]-pressure[i][k]+pressure[i][k-1];
		}

	for (int i=0;i<gridSize;i++)
		for (int k=1;k<gridSize;k++){
			vValues[i*(gridSize+1)+k]=advectedvValues[i*(gridSize+1)+k]-pressure[k][i]+pressure[k-1][i];
		}


}

int* Grid2::getCellOfPosition(float x, float y){

	int* cell = (int*) malloc(2*sizeof(int));

	float cellWidth = displaySize/gridSize;
	int faceleft = floor (x/cellWidth);
	int facelower = floor (y/cellWidth);

	if (faceleft == gridSize)
		faceleft--;

	if (facelower == gridSize)
		facelower--;

	cell[0]=faceleft;
	cell[1]=facelower;

	return cell;
}

//give position of left face x and lower face y as params
int Grid2::getSectorOfPosition(float leftFacex, float lowerFacey, float x, float y){

	int sector;
	float cellWidth = displaySize/gridSize;
	float centerx = leftFacex + 0.5*cellWidth;
	float centery = lowerFacey + 0.5*cellWidth;

	if (x<=centerx && y<=centery)
		sector=1;
	else if (x>centerx && y<=centery)
		sector=2;
	else if (x>centerx && y>centery)
		sector=3;
	else if (x<=centerx && y>centery)
		sector=4;

	return sector;
}

void Grid2::printPressure() {

	cout<<pressure[0][0]<<"; "<<pressure[0][1]<<"; "<<pressure[0][2]<<"..."<<pressure[0][37]<<"; "<<pressure[0][38]<<"; "<<pressure[0][39]<<"\n";

}

float* Grid2::getVelocityAtCellPosition(int x, int y, int which){

	float* velo = (float*)malloc(2*sizeof(float));
	velo[0]=-1;
	velo[1]=-1;

	if (which==1){
		if (x==0){
			if (y==0){
			velo[0]=0.5*(uValues[0]);
			velo[1]=0.5*(vValues[0]);
			}
			else{
			velo[0]=0.5*(uValues[(y-1)*(gridSize+1)]+uValues[y*(gridSize+1)]);
			velo[1]=0.5*(vValues[y]);
			}
		}
		else{
			if (y==0){
			velo[0]=0.5*uValues[x];
			velo[1]=0.5*(vValues[x*(gridSize+1)]+vValues[(x-1)*(gridSize+1)]);
			}
			else{
			velo[0]=0.5*(uValues[(y-1)*(gridSize+1)+x]+uValues[y*(gridSize+1)+x]);
			velo[1]=0.5*(vValues[x*(gridSize+1)+y]+vValues[(x-1)*(gridSize+1)+y]);
			}
		}
	}
	else if (which==2){
		if (x==(gridSize-1)){
			if (y==0){
			velo[0]=0.5*uValues[gridSize];
			velo[1]=0.5*vValues[(gridSize-1)*(gridSize+1)];
			}
			else{
			velo[0]=0.5*(uValues[y*(gridSize+1)+gridSize]+uValues[(y-1)*(gridSize+1)+gridSize]);
			velo[1]=0.5*vValues[(gridSize-1)*(gridSize+1)+y];
			}
		}
		else{
			if (y==0){
			velo[0]=0.5*uValues[x+1];
			velo[1]=0.5*(vValues[x*(gridSize+1)]+vValues[(x+1)*(gridSize+1)]);
			}
			else{
			velo[0]=0.5*(uValues[(y-1)*(gridSize+1)+x+1]+uValues[y*(gridSize+1)+x+1]);
			velo[1]=0.5*(vValues[x*(gridSize+1)+y]+vValues[(x+1)*(gridSize+1)+y]);
			}
		}
	}
	else if (which==3){
		if (x==(gridSize-1)){
			if (y==gridSize-1){
				velo[0]=0.5*(uValues[(gridSize-1)*(gridSize+1)+gridSize]);
				velo[1]=0.5*(vValues[(gridSize-1)*(gridSize+1)+gridSize]);
			}
			else{
				velo[0]=0.5*(uValues[y*(gridSize+1)+gridSize]+uValues[(y+1)*(gridSize+1)+gridSize]);
				velo[1]=0.5*(vValues[(gridSize-1)*(gridSize+1)+y+1]);
			}
		}
		else{
			if (y==gridSize-1){
				velo[0]=0.5*(uValues[(gridSize-1)*(gridSize+1)+x+1]);
				velo[1]=0.5*(vValues[x*(gridSize+1)+gridSize]+vValues[(x+1)*(gridSize+1)+gridSize]);
			}
			else{
				velo[0]=0.5*(uValues[y*(gridSize+1)+x+1]+uValues[(y+1)*(gridSize+1)+x+1]);
				velo[1]=0.5*(vValues[x*(gridSize+1)+y+1]+vValues[(x+1)*(gridSize+1)+y+1]);
			}
		}
	}
	else if (which==4){
		if (x==0){
			if (y==gridSize-1){
				velo[0]=0.5*(uValues[(gridSize-1)*(gridSize+1)]);
				velo[1]=0.5*(vValues[gridSize]);
			}
			else{
				velo[0]=0.5*(uValues[y*(gridSize+1)]+uValues[(y+1)*(gridSize+1)]);
				velo[1]=0.5*(vValues[y+1]);
			}
		}
		else{
			if (y==gridSize-1){
				velo[0]=0.5*(uValues[(gridSize-1)*(gridSize+1)+x]);
				velo[1]=0.5*(vValues[(x-1)*(gridSize+1)+gridSize]+vValues[x*(gridSize+1)+gridSize]);
			}
			else{
				velo[0]=0.5*(uValues[y*(gridSize+1)+x+1]+uValues[(y+1)*(gridSize+1)+x+1]);
				velo[1]=0.5*(vValues[x*(gridSize+1)+y+1]+vValues[(x-1)*(gridSize+1)+y+1]);
			}
		}
	}
	else if (which==5){
		velo[0]=0.5*(uValues[y*(gridSize+1)+x]+uValues[y*(gridSize+1)+x+1]);
		velo[1]=0.5*(vValues[x*(gridSize+1)+y]+vValues[x*(gridSize+1)+y+1]);
	}

	return velo;

}

float* Grid2::getVelocityAtPosition(float x, float y, int verbose){

	int* cell = getCellOfPosition(x,y);

	float* velo = (float*)malloc(2*sizeof(float));


	float leftFacex = uPositionsx[cell[1]*(gridSize+1)+cell[0]];
	float lowerFacey = vPositionsy[cell[0]*(gridSize+1)+cell[1]];
	float distLeft = x-leftFacex;
	float distLower = y-lowerFacey;
	float distLeftDivWidth = distLeft;
	float distLowerDivWidth = distLower;
	float factorLeft = 1 - distLeftDivWidth;
	float factorLower = 1 - distLowerDivWidth;
	float factorRight = distLeftDivWidth;
	float factorUpper = distLowerDivWidth;



	velo[0]=factorLeft*uValues[cell[1]*(gridSize+1)+cell[0]]+factorRight*uValues[cell[1]*(gridSize+1)+cell[0]+1];

	velo[1]=factorLower*vValues[cell[0]*(gridSize+1)+cell[1]]+factorUpper*vValues[cell[0]*(gridSize+1)+cell[1]+1];

	return velo;



}

float* Grid2::getVelocityAtPosition2(float x, float y){

	float* velo = (float*)malloc(2*sizeof(float));
	float* veloUpLeft = (float*)malloc(2*sizeof(float));
	float* veloUpRight = (float*)malloc(2*sizeof(float));
	float* veloDownLeft = (float*)malloc(2*sizeof(float));
	float* veloDownRight = (float*)malloc(2*sizeof(float));

	int* cell = getCellOfPosition(x,y);

	int sector = getSectorOfPosition(uPositionsx[cell[1]*(gridSize+1)+cell[0]],vPositionsy[cell[0]*(gridSize+1)+cell[1]],x,y);


	if (sector==1){
		veloUpLeft=getVelocityAtVerticalCellFace(cell[1],cell[0]);
		veloUpRight=getVelocityAtCellPosition(cell[0],cell[1],5);
		veloDownLeft=getVelocityAtCellPosition(cell[0],cell[1],1);
		veloDownRight=getVelocityAtHorizontalCellFace(cell[0],cell[1]);
	}
	else if (sector==2){
		veloUpLeft=getVelocityAtCellPosition(cell[0],cell[1],5);
		veloUpRight=getVelocityAtVerticalCellFace(cell[1],cell[0]+1);
		veloDownLeft=getVelocityAtHorizontalCellFace(cell[0],cell[1]);
		veloDownRight=getVelocityAtCellPosition(cell[0],cell[1],2);
	}
	else if (sector==3){
		veloUpLeft=getVelocityAtHorizontalCellFace(cell[0],cell[1]+1);
		veloUpRight=getVelocityAtCellPosition(cell[0],cell[1],3);
		veloDownLeft=getVelocityAtCellPosition(cell[0],cell[1],5);
		veloDownRight=getVelocityAtVerticalCellFace(cell[1],cell[0]+1);
	}
	else if (sector==4){
		veloUpLeft=getVelocityAtCellPosition(cell[0],cell[1],4);
		veloUpRight=getVelocityAtHorizontalCellFace(cell[0],cell[1]+1);
		veloDownLeft=getVelocityAtVerticalCellFace(cell[1],cell[0]);
		veloDownRight=getVelocityAtCellPosition(cell[0],cell[1],5);
	}

	velo[0]=0.25*(veloUpLeft[0]+veloUpRight[0]+veloDownLeft[0]+veloDownRight[0]);
	velo[1]=0.25*(veloUpLeft[1]+veloUpRight[1]+veloDownLeft[1]+veloDownRight[1]);
	return velo;

}

float* Grid2::getVelocityAtVerticalCellFace(int i, int k){

	float* velocity=(float*)malloc(2*sizeof(float));

	float velx;
	float vely=0;

	velx=uValues[i*(gridSize+1)+k];

	if (k == 0) {
		vely += vValues[k * (gridSize + 1) + i];
		vely += vValues[k * (gridSize + 1) + i + 1];
		vely = vely / 4;
	} else if (k == gridSize) {
		vely += vValues[(k - 1) * (gridSize + 1) + i];
		vely += vValues[(k - 1) * (gridSize + 1) + i + 1];
		vely = vely / 4;
	} else {

		vely += vValues[(k - 1) * (gridSize + 1) + i];
		vely += vValues[(k - 1) * (gridSize + 1) + i + 1];
		vely += vValues[k * (gridSize + 1) + i];
		vely += vValues[k * (gridSize + 1) + i + 1];
		vely = vely / 4;

	}

	velocity[0]=velx;
	velocity[1]=vely;
	return velocity;

}

float* Grid2::getVelocityAtHorizontalCellFace(int i, int k){

	float* velocity=(float*)malloc(2*sizeof(float));


	 float vely=vValues[i*(gridSize+1)+k];
	 float velx=0;

	    if (k==0){
	        velx=uValues[k*(gridSize+1)+i];
	        velx+=uValues[k*(gridSize+1)+i+1];
	        velx=velx/4;
	    }
	    else if (k==gridSize){
	        velx=uValues[(k-1)*(gridSize+1)+i];
	        velx+=uValues[(k-1)*(gridSize+1)+i+1];
	        velx=velx/4;
	    }
	    else{
	        velx+=uValues[(k-1)*(gridSize+1)+i];
	        velx+=uValues[(k-1)*(gridSize+1)+i+1];
	        velx+=uValues[k*(gridSize+1)+i];
	        velx+=uValues[k*(gridSize+1)+i+1];
	        velx=velx/4;

	    }
	    velocity[0]=velx;
	    velocity[1]=vely;
	    return velocity;
}

float* Grid2::getOldPositionOfParticleAtVerticalCellFace(int i,int k){

	float* oldPosition=(float*)malloc(2*sizeof(float));

	float* oldVelocity=getVelocityAtVerticalCellFace(i,k);


	oldPosition[0]=uPositionsx[i*(gridSize+1)+k]-oldVelocity[0]*timestep;
	oldPosition[1]=uPositionsy[i*(gridSize+1)+k]-oldVelocity[1]*timestep;


	return oldPosition;

}

float* Grid2::getOldPositionOfParticleAtHorizontalCellFace(int i, int k){

	float* oldPosition=(float*)malloc(2*sizeof(float));

	float* oldVelocity = getVelocityAtHorizontalCellFace(i, k);

	oldPosition[0]=vPositionsx[i*(gridSize+1)+k] - oldVelocity[0] * timestep;
	oldPosition[1]=vPositionsy[i*(gridSize+1)+k] - oldVelocity[1] * timestep;

	return oldPosition;
}

MarkerParticle Grid2::getMarker(int i){
	return markers[i];
}

int Grid2::getSize(){
	return gridSize;
}

void Grid2::moveMarkers(){

	float* temp;
	float* temp2=(float*)malloc(2*sizeof(float));

	for (int i=0;i<gridSize*gridSize*4;i++){
		temp=getVelocityAtPosition2(markers[i].getPosition()[0], markers[i].getPosition()[1]);
		markers[i].setVelocity(temp[0],temp[1]);
		temp2[0]=markers[i].getPosition()[0]+timestep*temp[0];
		temp2[1]=markers[i].getPosition()[1]+timestep*temp[1];
		markers[i].setPosition(temp2[0],temp2[1]);

	}
}

Grid2::~Grid2() {
	// TODO Auto-generated destructor stub
}

