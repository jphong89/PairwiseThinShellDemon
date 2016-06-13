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

/* registration params */
double YOUNG = 2;
double POISSON = 0.05;
double REGWEIGHT = 40; //stretching weight
double BENDWEIGHT = 400;  //bending weight
double DISTHRESHOLD = 12.0; // we don't want to match things too far away
double EUCLWEIGHT = 1;
int MESHLABOPTION = 1;
double LANDMARKWEIGHT = 50;
bool is_orthotropic = false;

/* the two surfaces we want to register. note!! the variable name doesn't reflect if it's CT or RC */
BasicMesh* CTSurface;
BasicMesh* RCSurface;

link* linkList;
bool* attractor;
bool* attractee;
int* linkTarget;
bool* occluded;

/* registration variables */
double* affinity;
Vector_3* bestMatch;
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
	bestMatch = new Vector_3[affinityM];
	for (int i = 0; i < affinityM; i++) bestMatch[i]=Vector_3(0,0,0);	
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
	char* filename1;
	char* filename2;
	char* landmark1; //landmark1
	char* landmark2; //landmark2

	//format: demons.exe surface1.off surface2.off weight_output.txt landmark1.txt landmark2.txt
	if (argc > 2){
		filename1 = argv[1];
		filename2 = argv[2];
		
		if (argc > 3) {
			REGWEIGHT = atof(argv[3]);
			BENDWEIGHT = atof(argv[4]);
		}

		if (argc > 5) DISTHRESHOLD = atof(argv[5]);
		if (argc > 6) EUCLWEIGHT = atof(argv[6]);
		
		if (argc > 7) {
			YOUNG = atof(argv[7]);
			POISSON = atof(argv[8]);
			if (YOUNG < 0 || POISSON < 0) is_orthotropic = true;
		}

		if (argc > 9) {
			MESHLABOPTION = atoi(argv[9]);
		}

		if (argc > 10) {
			landmark1 = argv[10];
			landmark2 = argv[11];
			LANDMARKWEIGHT = atof(argv[12]);
		}

		patientNum = new char[20]; //manually designate object name
		caseNum = new char[20];	//manually designate case number

		cout<<"STRETCHING WEIGHT: "<<REGWEIGHT<<" -BENDING WEIGHT: "<<BENDWEIGHT<<" -EUCLIDEAN THRESHOLD: "
			<<DISTHRESHOLD<<" -GEOMETRIC FEATURE WEIGHT: "<<EUCLWEIGHT
			<<" -YOUNG: "<<YOUNG<<" -POSSION: "<<POISSON<<" -MESHLAB: "<<MESHLABOPTION<<endl;
	}else{
		patientNum = new char[20]; //manully designate object name
		caseNum = new char[20];

		patientNum = "5708_partial2";//manually designate object name
		caseNum = "2";//manually designate index number

		////////////////////////////////////////* read files */////////////////////////////////

		strcat(filename1,patientNum);
		strcat(filename1,"/");
		strcat(filename1,caseNum);
		strcat(filename1,"_origin.off");

		strcat(filename2,patientNum);
		strcat(filename2,"/");
		strcat(filename2,caseNum);
		strcat(filename2,"_deformed.off");
	}

	////////////////////////////////////*read & translate surfaces */////////////////////////////////////////////
	cout<<"Begin Constructing First Mesh: "<<filename1<<endl;
	RCSurface = new BasicMesh(filename1,"origin",patientNum,caseNum);
	RCSurface->ComputeMeshProperty(filename1);
	if (argc > 10) RCSurface->constructLandmark(landmark1);
	if (is_orthotropic){
		cout<<"computing orthotropic case..."<<endl;
		RCSurface->constructOrthotropic(); 
	
	}
	cout<<"Begin Construction Second Mesh: "<<filename2<<endl;
	CTSurface = new BasicMesh(filename2,"deformed",patientNum,caseNum);
	CTSurface->ComputeMeshProperty(filename2);
	if (argc > 10) CTSurface->constructLandmark(landmark2);

	affinityM = RCSurface->vertexNum;
	affinityN = CTSurface->vertexNum;

	////////////////////////////////////*compute Geometric Features*/////////////////////////////////////////////
	initialMemeoryForGlobalVariables();

	cout<<"Begin Computing Geometric Features for First Mesh..."<<endl;
	RCSurface->findSignatureAll();
	cout<<"Begin Computing Geometric Features for Second Mesh..."<<endl;
	CTSurface->findSignatureAll();

	////////////////////////////////////*output initial errors*/////////////////////////////////////////////
	/* synthetic test validation;   doesn't make sense for real data registration */
	computeStatisticBool(RCSurface,CTSurface,initialError,initialBoundError,initialVariance,initialBoundVar, maxInitialError);
	//cout<<initialError<<' '<<initialVariance<<' '<<initialBoundError<<' '<<initialBoundVar<<' '<<maxInitialError<<endl;
	//cout<<statisticVertexNum<<' '<<statisticBoundaryNum<<endl;

	/* prepare data structures for registration */
	//RCSurface->constructLink();
	RCSurface->constructInitialFrameInfo();
	//RCSurface->constructOccluded();
	//RCSurface->constructAttractor();
	RCSurface->constructEdgeList();
	RCSurface->constructVertexList();

	/////////////////////////////* 3compute correspondences for registration *////////////////////////////////////////////////////////////
	RCSurface->findCorrespondenceBothWay(CTSurface,EUCLWEIGHT);
	RCSurface->findMatch2(CTSurface);
	RCSurface->outputFeatures();
	
	/* output force weight */
	//ofstream fout;
	//fout.open(filename3);
	//for (int i = 0; i<affinityM; i++) fout<<matchWeight[i]<<endl;
	//fout.close();
	//RCSurface->findClosestPoint(CTSurface);

	/////////////////////////////* rendering and optimization *////////////////////////////////////////////////////////////
// 	beginGlut(argc,argv);
	if (pthread_mutex_init(&lock, NULL) != 0)
	{
		cout<<"mutex init failed\n"<<endl;
		return 1;
	}
	int* tmp;
	startOptimization(tmp); //direct begin for script
	return 0;
}

