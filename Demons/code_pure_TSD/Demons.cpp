// Demons.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "MeshObject.h"
#include "GlutManager.h"
#include <time.h>

#include <iostream>
#include <fstream>
using namespace std;

/* surface variables */ 
char* patientNum;
char* caseNum;

BasicMesh* CTSurface;
BasicMesh* RCSurface;
BasicMesh* DFSurface;

link* linkList;
bool* attractor;
bool* attractee;
int* linkTarget;
bool* occluded;

/* registration variables */
double* affinity;
int* bestMatch;
float* matchWeight;
int affinityM = 0;
int affinityN = 0;
float geodesicDistance = GEODECISDIS;

lbfgsfloatval_t* u;
double bendweight;
double regweight;

/*debug*/
double facetBending = 0;
double facetStretching = 0;
double distantLink = 0;
LARGE_INTEGER t1,t2,tc;
double thetaDif1,thetaDif2,thetaDif3;
double facetTrace,facetDet;

/* camera variables */
float camera_Scale = 1;
float camera_Up = 0;
float camera_Right = 0;

float rotate_X = 0;
float rotate_Y = 0;
float rotate_Z = 0;

bool switch1 = false;
bool switch2 = false;
bool switch3 = false;
bool switch4 = false;

float minx,miny,minz;
float maxx,maxy,maxz;
float max_all;

pthread_mutex_t lock;

/* validation */
bool* statisticBool;
bool* isBoundary;
int statisticBoundaryNum = 0;
int statisticVertexNum = 0;

float initialError;
float initialVariance;
float initialBoundError;
float initialBoundVar;

float maxInitialError;

/* initialize memeory */
void initialMemeoryForGlobalVariables(){
	//initial affinity map
	affinity = new double[affinityM * affinityN];
	memset(affinity,0,affinityM * affinityN * sizeof(double));
	//initial best match location
	bestMatch = new int[affinityM];
	memset(bestMatch,0,sizeof(int)*affinityM);
	//initial match confidence
	matchWeight = new float[affinityM];
	memset(matchWeight,0,sizeof(float)*affinityM);
	//initial attractor
	attractor = new bool[affinityM];
	memset(attractor,false,affinityM * sizeof(bool));
	//initial attractee
	attractee = new bool[affinityM];
	memset(attractee,false,affinityM * sizeof(bool));
	//initial linkTarget
	linkTarget = new int[affinityM];
	memset(linkTarget,0,affinityM * sizeof(int));
	//initial occluded region
	occluded = new bool[affinityM];
	memset(occluded,false,affinityM * sizeof(bool));
	//initial statisticalBool
	statisticBool = new bool[affinityN];
	memset(statisticBool,false,sizeof(bool)*affinityN);
	isBoundary = new bool[affinityN];
	memset(isBoundary,false,sizeof(bool)*affinityN);
	//initial DVF
	u = lbfgs_malloc(affinityM * 3);
	memset(u,0,sizeof(lbfgsfloatval_t)*affinityM*3);
}

int main(int argc, char* argv[])
{
	if (argc >= 2){
		patientNum = argv[1];
		caseNum = argv[2];		
	}else{
		patientNum = new char[20]; //manully designate object name
		caseNum = new char[20];

		patientNum = "5705";//manually designate object name
		caseNum = "2";//manually designate index number
	}

	////////////////////////////////////////* read files */////////////////////////////////
	char filename1[50] = "data/";
	strcat(filename1,patientNum);
	strcat(filename1,"/");
	strcat(filename1,caseNum);
	strcat(filename1,"_origin.off");

	char filename2[50] = "data/";
	strcat(filename2,patientNum);
	strcat(filename2,"/");
	strcat(filename2,caseNum);
	strcat(filename2,"_deformed.off");

	cout<<"Begin Constructing First Mesh: "<<filename1<<endl;
	RCSurface = new BasicMesh(filename1,"origin",patientNum,caseNum);
	RCSurface->ComputeMeshProperty(filename1);

	cout<<"Begin Construction Second Mesh: "<<filename2<<endl;
	CTSurface = new BasicMesh(filename2,"deformed",patientNum,caseNum);
	CTSurface->ComputeMeshProperty(filename2);

	affinityM = RCSurface->vertexNum;
	affinityN = CTSurface->vertexNum;

	////////////////////////////////////*compute Geometric Features*/////////////////////////////////////////////
	initialMemeoryForGlobalVariables();

	cout<<"Begin Computing Geometric Features for First Mesh..."<<endl;
	RCSurface->findSignatureAll();
	cout<<"Begin Computing Geometric Features for Second Mesh..."<<endl;
	CTSurface->findSignatureAll();

	////////////////////////////////////*output initial errors*/////////////////////////////////////////////
	computeStatisticBool(RCSurface,CTSurface,initialError,initialBoundError,initialVariance,initialBoundVar, maxInitialError);
	cout<<initialError<<' '<<initialVariance<<' '<<initialBoundError<<' '<<initialBoundVar<<' '<<maxInitialError<<endl;
	cout<<statisticVertexNum<<' '<<statisticBoundaryNum<<endl;

	RCSurface->constructLink();
	//RCSurface->constructOccluded();
	//RCSurface->constructAttractor();
	//RCSurface->constructEdgeList();

	/////////////////////////////* prepare registration *////////////////////////////////////////////////////////////
	RCSurface->findCorrespondenceBothWay(CTSurface,0.9);
	RCSurface->findMatch(CTSurface);
	//RCSurface->findClosestPoint(CTSurface);

	beginGlut(argc,argv);
	return 0;
}

