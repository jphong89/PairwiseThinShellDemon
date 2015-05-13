#include "stdafx.h"
#include "cmaes.h"
#include "GlutManager.h"

#include <iostream>
using namespace std;

void startCMAES(){
	char filename[50];

	regweight = REGWEIGHT;
	bendweight = BENDWEIGHT;

	for (int i = 0; i < 25; i++){
		cout<<"Velocity Iteration: "<<i<<endl;

		if (i > 0){
			cout<<"Reset surface..."<<endl;

			RCSurface->ComputeMeshProperty(filename);

			//cout<<"Recompute geometric features..."<<endl;
			RCSurface->findSignatureAll();
			RCSurface->constructEdgeList();

			memset(u,0,sizeof(lbfgsfloatval_t) * affinityM * 3);
			memset(affinity,0,affinityM * affinityN * sizeof(double));
			memset(adjacent,false,affinityM * affinityM * sizeof(bool));
			memset(bestMatch,0,sizeof(int)*affinityM);
			memset(matchWeight,0,sizeof(double)*affinityM);

			cout<<"Recompute affinity map..."<<endl;
			/* forward reg */
			RCSurface->findCorrespondenceBothWay(CTSurface,EUCLWEIGHT /*+ i * 0.04*/);
			RCSurface->ComputeAdjacency();
			RCSurface->findMatch(CTSurface);

		}

		cout<<"Begin optimization..."<<endl;
		CMAES<double> evo;
		double *arFunvals, *const*pop, *xfinal;

		// Initialize everything
		const int dim = affinityM * 3;
		double* xstart = new double[dim];
		for(int i=0; i<dim; i++) xstart[i] = 0.0;
		double* stddev = new double[dim];
		for(int i=0; i<dim; i++) stddev[i] = 0.01;
		Parameters<double> parameters;

		// TODO Adjust parameters here
		parameters.init(dim, xstart, stddev);
		parameters.lambda = 64;
		parameters.mu = 16;
		arFunvals = evo.init(parameters);
		
		std::cout << evo.sayHello() << std::endl;

		int iterCount = 0;
		// Iterate until stop criterion holds
		while(!evo.testForTermination()){
			iterCount++;
			cout<<"iteration number: "<<iterCount<<endl;
			cout<<"Lambda: "<<evo.get(CMAES<double>::Lambda)<<endl;
			// Generate lambda new search points, sample population
			pop = evo.samplePopulation(); // Do not change content of pop

			// evaluate the new search points using fitfun from above
			double minVal = 300000000;
			for (int i = 0; i < evo.get(CMAES<double>::Lambda); ++i){
				arFunvals[i] = evaluate_cmaes(pop[i], (int) evo.get(CMAES<double>::Dimension));
				if (arFunvals[i] < minVal) minVal = arFunvals[i];
				//cout<<" func val: "<<arFunvals[i];
			}
			
			cout<<"func val: "<<minVal<<endl;
			// update the search distribution used for sampleDistribution()
			evo.updateDistribution(arFunvals);
		}

		std::cout << "Stop:" << std::endl << evo.getStopMessage();
		evo.writeToFile(CMAES<double>::WKResume, "resumeevo1.dat"); // write resumable state of CMA-ES

		// get best estimator for the optimum, xmean
		xfinal = evo.getNew(CMAES<double>::XMean); // "XBestEver" might be used as well


		// do something with final solution and finally release memory
		for (int i = 0; i < dim; i++)
			u[i] = xfinal[i];

		filename[0] = '\0';
		strcat(filename,"temp");
		char iterNum[10];
		itoa(i,iterNum,10);
		strcat(filename,iterNum);
		strcat(filename,".off");
		writeOFF(RCSurface,filename);

		delete[] xfinal;
		delete []xstart;
		delete []stddev;
	}
}

double evaluate_cmaes(double const *u, int N){
	//cout<<"begin optimization"<<endl;
	double data = penalizeData_cmaes(u,N);
	/* stretching */
	double stretching = penalizeStretch_cmaes(u,N);
	/* bending */
	double bending = penalizeBend_cmaes(u,N);

	float linkStretching = 0;
	for (int i = 0; i < RCSurface->linkNum; i++){
		int idx1 = linkList[i].index1;
		int idx2 = linkList[i].index2;

		Vector_3 deformVec(u[idx2*3] - u[idx1*3],u[idx2*3+1] - u[idx1*3+1],u[idx2*3+2] - u[idx1*3+2]);

		linkStretching += LINKWEIGHT * deformVec.squared_length();
	}

	return data + /*stretching + */ bending + linkStretching;
}

double penalizeData_cmaes(double const *u, int N){
	PolyhedralSurf::Vertex_const_iterator vb = RCSurface->P.vertices_begin();
	PolyhedralSurf::Vertex_const_iterator ve = RCSurface->P.vertices_end();
	PolyhedralSurf::Vertex_const_handle vh;

	double fx = 0.0;

	for (int i = 0; vb != ve; vb++,i++){
		vh = vb;

		Point_3 deformed = vh->point() + Vector_3(u[i*3],u[i*3+1],u[i*3+2]);
		Point_3 cor;
		double weight = matchWeight[i];

		if ((matchWeight[i] <= 1) && (matchWeight[i] > 0.01))
			cor = CTSurface->vertexIndex[bestMatch[i]]->point();
		else if (matchWeight[i] > 1)
			cor = RCSurface->vertexIndex[bestMatch[i]]->point()
			+ Vector_3(u[bestMatch[i]*3],u[bestMatch[i]*3+1],u[bestMatch[i]*3+2]);
		else 
			weight = 0;

		Vector_3 dis(deformed,cor);

		fx = fx + weight * dis.squared_length();
	}

	return fx;
}

double penalizeStretch_cmaes(double const *u, int N){
	double stretching = 0;

	if (regweight == REGWEIGHT - 0.5){ // edge-based  stretching energy
		for (int i=0; i < RCSurface->edgeNum/2;i++){


			edge he = RCSurface->edgeList[i];

			double weight =  regweight * 2 * he.stretchWeight;
			int idx1 = he.index1;
			int idx2 = he.index2;

			Vector_3 deformVec(u[idx2*3] - u[idx1*3],u[idx2*3+1] - u[idx1*3+1],u[idx2*3+2] - u[idx1*3+2]);

			stretching += (deformVec.squared_length() * weight);
		}
	}else{ //triangle-based stretching energy
		facet* faceList = RCSurface->faceList;

		double Dtemp_v1_x,Dtemp_v1_y,Dtemp_v1_z;
		double Dtemp_v2_x,Dtemp_v2_y,Dtemp_v2_z;

		Point_3 newV1,newV2,newV3;
		Vector_3 temp_v1,temp_v2;
		double v11,v12,v21,v22,J11,J12,J21,J22,S11,S12,S21,S22,k;
		double Dv11,Dv12,Dv21,Dv22,DJ11,DJ12,DJ21,DJ22,DS11,DS12,DS21,DS22,dk;
		double weight,trace,det,det_inv,Ddet;

		for (int idx=0; idx < RCSurface->faceNum;idx++){
			facet f = faceList[idx];
			weight = f.area * regweight;

			newV1 = f.v1->point() + Vector_3(u[f.index[1]*3],u[f.index[1]*3+1],u[f.index[1]*3+2]);
			newV2 = f.v2->point() + Vector_3(u[f.index[2]*3],u[f.index[2]*3+1],u[f.index[2]*3+2]);
			newV3 = f.v3->point() + Vector_3(u[f.index[3]*3],u[f.index[3]*3+1],u[f.index[3]*3+2]);

			temp_v1 = Vector_3(newV1,newV2);
			temp_v2 = Vector_3(newV1,newV3);

			v11 = 0;
			v12 = sqrt(temp_v1.squared_length());
			if (abs(v12)<0.0000001) v12 = 0.00001;

			k = temp_v1*temp_v2;
			v22 = k / v12;

			double sqv21 = temp_v2.squared_length() - v22*v22;

			if (sqv21 > 0.000000001)
				v21 = - sqrt(sqv21);
			else
				v21 = - 0.00001;

			J11 = 0*f.inverse[0][0] + v21*f.inverse[1][0];
			J12 = 0*f.inverse[0][1] + v21*f.inverse[1][1];
			J21 = v12*f.inverse[0][0] + v22*f.inverse[1][0];
			J22 = v12*f.inverse[0][1] + v22*f.inverse[1][1];

			S11 = J11*J11 + J21*J21; S12 = J11*J12+J21*J22; S21 = S12; S22 = J12*J12 + J22*J22;
			trace = (S11 + S22);
			det = S11*S22 - S12*S21;
			if (det < 0.0000001) det = 0.0000001;
			det_inv = 1 / det;

			stretching += weight * (STRETCH_MIU / 2 * trace + (STRETCH_LAMBDA - STRETCH_MIU*2)/8 * det + (STRETCH_LAMBDA + STRETCH_MIU*2)/8 * det_inv);
		}
	}

	facetStretching = stretching;
	return stretching;
}

double penalizeBend_cmaes(double const *u, int N){
	double bending = 0;
	facet* faceList = RCSurface->faceList;

	Point_3 newV1,newV2,newV3,newNB1,newNB2,newNB3;
	double dv1x, dv1y, dv1z, dv2x, dv2y, dv2z, dv3x, dv3y, dv3z,
		dnb1x,dnb1y,dnb1z,dnb2x,dnb2y,dnb2z,dnb3x,dnb3y,dnb3z;
	double Dtheta1,Dtheta2,Dtheta3,DSOD_00,DSOD_01,DSOD_10,DSOD_11;
	double weight,Dbend;

	/*
	for (int idx=0; idx < RCSurface->faceNum;idx++){
		facet f = faceList[idx];
		if (f.isBorder) continue;
		weight = bendweight;

		newV1 = f.v1->point() + Vector_3(u[f.index[1]*3],u[f.index[1]*3+1],u[f.index[1]*3+2]);
		newV2 = f.v2->point() + Vector_3(u[f.index[2]*3],u[f.index[2]*3+1],u[f.index[2]*3+2]);
		newV3 = f.v3->point() + Vector_3(u[f.index[3]*3],u[f.index[3]*3+1],u[f.index[3]*3+2]);
		newNB1 = f.nb1->point() + Vector_3(u[f.nbIdx1*3],u[f.nbIdx1*3+1],u[f.nbIdx1*3+2]);
		newNB2 = f.nb2->point() + Vector_3(u[f.nbIdx2*3],u[f.nbIdx2*3+1],u[f.nbIdx2*3+2]);
		newNB3 = f.nb3->point() + Vector_3(u[f.nbIdx3*3],u[f.nbIdx3*3+1],u[f.nbIdx3*3+2]);

		double* det1 = computeDeterminant(newV1,newV2,newV3,newNB1);
		double* det2 = computeDeterminant(newV1,newV2,newV3,newNB2);
		double* det3 = computeDeterminant(newV1,newV2,newV3,newNB3);

		double theta1 = - det1[0] / ((f.area*f.sideArea1)/f.l1);
		double theta2 = - det2[0] / ((f.area*f.sideArea2)/f.l2);
		double theta3 = - det3[0] / ((f.area*f.sideArea3)/f.l3);

		double SOD_00 = (theta1-f.theta1) * f.SO1[0][0] + (theta2-f.theta2) * f.SO2[0][0] + (theta3-f.theta3) * f.SO3[0][0];
		double SOD_01 = (theta1-f.theta1) * f.SO1[0][1] + (theta2-f.theta2) * f.SO2[0][1] + (theta3-f.theta3) * f.SO3[0][1];
		double SOD_10 = (theta1-f.theta1) * f.SO1[1][0] + (theta2-f.theta2) * f.SO2[1][0] + (theta3-f.theta3) * f.SO3[1][0];
		double SOD_11 = (theta1-f.theta1) * f.SO1[1][1] + (theta2-f.theta2) * f.SO2[1][1] + (theta3-f.theta3) * f.SO3[1][1];

		bending += weight * (SOD_00*SOD_00 + SOD_01*SOD_01 + SOD_10*SOD_10 + SOD_11*SOD_11); 

		delete[] det1;
		delete[] det2;
		delete[] det3;
	}
	*/
	facetBending = bending;
	return bending;
}