#include "stdafx.h"
#include "GlutManager.h"

#include <iostream>
#include <fstream>
using namespace alglib;

void startOptimization_cg(){
	/* initialize optimization */
	char filename[50];
	char funcValueFileName[50];

	QueryPerformanceCounter(&t1);
	for (int i = 0; i < 25; i++){
		cout<<"Velocity Iteration: "<<i<<endl;

		if (i > 0){
			cout<<"Reset surface..."<<endl;

			RCSurface->ComputeMeshProperty(filename);

			cout<<"Recompute geometric features..."<<endl;
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
		regweight = REGWEIGHT;
		bendweight = BENDWEIGHT;
		int opIterNum = 100000;
		funcValue = new double[opIterNum + 2];
		memset(funcValue,0,sizeof(double)*(opIterNum+2));

		real_1d_array x;
		x.setlength(affinityM*3);
		for (int i = 0; i < affinityM*3; i++) x[i] = 0;

		double epsg = 0.0000000001;
		double epsf = 0;
		double epsx = 0;
		double stpmax = 0.5;
		ae_int_t maxits = opIterNum;
		mincgstate state;
		mincgreport rep;
		
		// first run
		mincgcreate(x, state);
		mincgsetcond(state, epsg, epsf, epsx, maxits);
		mincgsetstpmax(state, stpmax);
		alglib::mincgoptimize(state, evaluate_cg);
		mincgresults(state, x, rep);

		for (int i = 0; i < affinityM*3; i++) u[i] = x[i];
		//output file
		filename[0] = '\0';
		strcat(filename,"temp");
		char iterNum[10];
		itoa(i,iterNum,10);
		strcat(filename,iterNum);
		strcat(filename,".off");
		writeOFF(RCSurface,filename);

		//output funcValue
		funcValueFileName[0] = '\0';
		strcat(funcValueFileName,"temp");
		itoa(i,iterNum,10);
		strcat(funcValueFileName,iterNum);
		strcat(funcValueFileName,".txt");

		ofstream fout;
		fout.open(funcValueFileName);

		for (int i = 0; i < opIterNum+2; i++)
			fout<<funcValue[i]<<endl;

		fout.close();

	
	}

}

void evaluate_cg(const real_1d_array &u, double &func, real_1d_array &g, void *ptr) 
{
	//cout<<"begin optimization"<<endl;
	double data = penalizeData(u,g);
	//cout<<"data: "<<data<<endl;

	/* stretching */
	double stretching = penalizeStretch(u,g);
	//cout<<" stretching: "<<stretching;

	/* bending */
	double bending = penalizeBend(u,g);
	//cout<<" bending: "<<bending<<endl;

	float linkStretching = 0;
	for (int i = 0; i < RCSurface->linkNum; i++){
		int idx1 = linkList[i].index1;
		int idx2 = linkList[i].index2;

		Vector_3 deformVec(u[idx2*3] - u[idx1*3],u[idx2*3+1] - u[idx1*3+1],u[idx2*3+2] - u[idx1*3+2]);

		g[idx1*3] += -2 * LINKWEIGHT * deformVec.x();
		g[idx1*3+1] += -2 * LINKWEIGHT * deformVec.y();
		g[idx1*3+2] += -2 * LINKWEIGHT * deformVec.z();
		g[idx2*3] += 2 * LINKWEIGHT * deformVec.x();
		g[idx2*3+1] += 2 * LINKWEIGHT * deformVec.y();
		g[idx2*3+2] += 2 * LINKWEIGHT * deformVec.z();

		linkStretching += LINKWEIGHT * deformVec.squared_length();
	}

	func = data + stretching +  bending + linkStretching;

}

double penalizeData(const real_1d_array &u, real_1d_array &g){
	PolyhedralSurf::Vertex_const_iterator vb = RCSurface->P.vertices_begin();
	PolyhedralSurf::Vertex_const_iterator ve = RCSurface->P.vertices_end();
	PolyhedralSurf::Vertex_const_handle vh;

	lbfgsfloatval_t fx = 0.0;

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

		g[i*3] = -2 * dis.x() * weight;
		g[i*3+1] = -2 * dis.y() * weight;
		g[i*3+2] = -2 * dis.z() * weight;
		fx = fx + weight * dis.squared_length();
	}

	return fx;
}

double penalizeStretch(const real_1d_array &u, real_1d_array &g){
	double stretching = 0;

	if (regweight == REGWEIGHT - 0.5){ // edge-based  stretching energy
		for (int i=0; i < RCSurface->edgeNum/2;i++){


			edge he = RCSurface->edgeList[i];

			double weight =  regweight * 2 * he.stretchWeight;
			int idx1 = he.index1;
			int idx2 = he.index2;

			Vector_3 deformVec(u[idx2*3] - u[idx1*3],u[idx2*3+1] - u[idx1*3+1],u[idx2*3+2] - u[idx1*3+2]);

			g[idx1*3] += -2 * deformVec.x() * weight;
			g[idx1*3+1] += -2 * deformVec.y() * weight;
			g[idx1*3+2] += -2 * deformVec.z() * weight;
			g[idx2*3] += 2 * deformVec.x() * weight;
			g[idx2*3+1] += 2 * deformVec.y() * weight;
			g[idx2*3+2] += 2 * deformVec.z() * weight;

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

			k = temp_v1*temp_v2;
			v22 = k / v12;
			if (_isnan(v22)){
				cout<<"~~~v22: isnan "<<temp_v1.squared_length()<<endl;
				continue;
			}

			double sqv21 = temp_v2.squared_length() - v22*v22;

			if (sqv21 > 0.000000001)
				v21 = - sqrt(sqv21);
			else
				v21 = - 0.00001;
			//if (_isnan(v21)) v21 = 0;

			J11 = 0*f.inverse[0][0] + v21*f.inverse[1][0];
			J12 = 0*f.inverse[0][1] + v21*f.inverse[1][1];
			J21 = v12*f.inverse[0][0] + v22*f.inverse[1][0];
			J22 = v12*f.inverse[0][1] + v22*f.inverse[1][1];

			S11 = J11*J11 + J21*J21; S12 = J11*J12+J21*J22; S21 = S12; S22 = J12*J12 + J22*J22;
			trace = (S11 + S22);
			det = S11*S22 - S12*S21;
			det_inv = 1 / det;

			if (_isnan(det_inv)) {
				cout<<J11<<' '<<J12<<' '<<J21<<' '<<J22<<endl;
				cout<<f.v1->point()<<endl;
				cout<<f.v2->point()<<endl;
				cout<<u[f.index[1]*3]<<' '<<u[f.index[1]*3+1]<<' '<<u[f.index[1]*3+2]<<endl;
				cout<<u[f.index[2]*3]<<' '<<u[f.index[2]*3+1]<<' '<<u[f.index[2]*3+2]<<endl;
				cout<<f.index[1]<<' '<<f.index[2]<<endl;
				cout<<idx<<endl;
				cin>>k;
			}

			stretching += weight * (STRETCH_MIU / 2 * trace + (STRETCH_LAMBDA - STRETCH_MIU*2)/8 * det + (STRETCH_LAMBDA + STRETCH_MIU*2)/8 * det_inv);

			if ((f.index[1] == 2640)&& (f.index[2] == 2638) && (f.index[3] == 2570)){
				facetStretching = weight * (STRETCH_MIU / 2 * trace + (STRETCH_LAMBDA - STRETCH_MIU*2)/8 * det + (STRETCH_LAMBDA + STRETCH_MIU*2)/8 * det_inv);
			}
			//stretching += weight * ((S11 - 1)*(S11 - 1) + (S22 - 1)*(S22 - 1) + S12*S12 + S21*S21);

			/* COMPUTE GRADIENT */
			for (int i = 1; i < 4; i++)
				for (int j = 1; j < 4; j++){
					Dtemp_v1_x = computeStretchGradient(1,1,i,j);Dtemp_v1_y = computeStretchGradient(1,2,i,j);Dtemp_v1_z = computeStretchGradient(1,3,i,j);
					Dtemp_v2_x = computeStretchGradient(2,1,i,j);Dtemp_v2_y = computeStretchGradient(2,2,i,j);Dtemp_v2_z = computeStretchGradient(2,3,i,j);

					dk = (Dtemp_v1_x*temp_v2.x()  +temp_v1.x()*Dtemp_v2_x) + (Dtemp_v1_y*temp_v2.y()  +temp_v1.y()*Dtemp_v2_y) + (Dtemp_v1_z*temp_v2.z()  +temp_v1.z()*Dtemp_v2_z);
					Dv12 = 0.5 * 1 / v12 * (2*temp_v1.x()*Dtemp_v1_x + 2*temp_v1.y()*Dtemp_v1_y + 2*temp_v1.z()*Dtemp_v1_z);
					Dv22 = (dk * v12 - Dv12 * k) / (v12*v12);
					Dv21 = 0.5 * 1 / v21 * (2*temp_v2.x()*Dtemp_v2_x + 2*temp_v2.y()*Dtemp_v2_y + 2*temp_v2.z()*Dtemp_v2_z - 2*v22*Dv22);

					DJ11 = Dv21*f.inverse[1][0];
					DJ12 = Dv21*f.inverse[1][1];
					DJ21 = Dv12*f.inverse[0][0] + Dv22*f.inverse[1][0];
					DJ22 = Dv12*f.inverse[0][1] + Dv22*f.inverse[1][1];

					DS11 = 2*J11*DJ11 + 2*J21*DJ21;
					DS22 = 2*J12*DJ12 + 2*J22*DJ22;
					DS12 = (DJ11*J12+J11*DJ12) + (DJ21*J22+J21*DJ22);
					DS21 = DS12;

					Ddet = (DS11*S22+S11*DS22) - (DS12*S21+S12*DS21);

					g[f.index[i]*3+j-1] += weight * STRETCH_MIU / 2 * (DS11 + DS22);
					g[f.index[i]*3+j-1] += weight * (STRETCH_LAMBDA - STRETCH_MIU*2)/8 * Ddet;
					g[f.index[i]*3+j-1] += weight * (STRETCH_LAMBDA + STRETCH_MIU*2)/8 * (-1/(det*det)*Ddet);

					//g[f.index[i]*3+j-1] += weight * (2*(S11 - 1)*DS11 + 2*(S22 - 1)*DS22 + 2*S12*DS12 + 2*S21*DS21);
				}


		}
	}

	facetStretching = stretching;
	return stretching;
}

double penalizeBend(const real_1d_array &u, real_1d_array &g){
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

		if ((f.index[1] == 2640)&& (f.index[2] == 2638) && (f.index[3] == 2570)){
			facetBending = weight * (SOD_00*SOD_00 + SOD_01*SOD_01 + SOD_10*SOD_10 + SOD_11*SOD_11);
			thetaDif1 = theta1 - f.theta1;
			thetaDif2 = theta2 - f.theta2;
			thetaDif3 = theta3 - f.theta3;
		}
		//COMPUTE GRADIENT

		for (int i = 1; i < 4; i++)
			for (int j = 1; j < 4; j++){
				Dtheta1 = - computeDetGradient(i,j,newV1,newV2,newV3,newNB1,det1) / ((f.area*f.sideArea1)/f.l1);
				Dtheta2 = - computeDetGradient(i,j,newV1,newV2,newV3,newNB2,det2) / ((f.area*f.sideArea2)/f.l2);
				Dtheta3 = - computeDetGradient(i,j,newV1,newV2,newV3,newNB3,det3) / ((f.area*f.sideArea3)/f.l3);

				DSOD_00 = Dtheta1 * f.SO1[0][0] + Dtheta2 * f.SO2[0][0] + Dtheta3 * f.SO3[0][0];
				DSOD_01 = Dtheta1 * f.SO1[0][1] + Dtheta2 * f.SO2[0][1] + Dtheta3 * f.SO3[0][1];
				DSOD_10 = Dtheta1 * f.SO1[1][0] + Dtheta2 * f.SO2[1][0] + Dtheta3 * f.SO3[1][0];
				DSOD_11 = Dtheta1 * f.SO1[1][1] + Dtheta2 * f.SO2[1][1] + Dtheta3 * f.SO3[1][1];

				Dbend = weight * (2*SOD_00*DSOD_00 + 2*SOD_01*DSOD_01 + 2*SOD_10*DSOD_10 + 2*SOD_11*DSOD_11);

				g[f.index[i]*3+j-1] += Dbend;
			}

			for (int j = 1; j < 4; j++){
				Dtheta1 = - computeDetGradient(4,j,newV1,newV2,newV3,newNB1,det1) / ((f.area*f.sideArea1)/f.l1);
				Dtheta2 = 0;
				Dtheta3 = 0;

				DSOD_00 = Dtheta1 * f.SO1[0][0] + Dtheta2 * f.SO2[0][0] + Dtheta3 * f.SO3[0][0];
				DSOD_01 = Dtheta1 * f.SO1[0][1] + Dtheta2 * f.SO2[0][1] + Dtheta3 * f.SO3[0][1];
				DSOD_10 = Dtheta1 * f.SO1[1][0] + Dtheta2 * f.SO2[1][0] + Dtheta3 * f.SO3[1][0];
				DSOD_11 = Dtheta1 * f.SO1[1][1] + Dtheta2 * f.SO2[1][1] + Dtheta3 * f.SO3[1][1];

				Dbend = weight * (2*SOD_00*DSOD_00 + 2*SOD_01*DSOD_01 + 2*SOD_10*DSOD_10 + 2*SOD_11*DSOD_11);

				g[f.nbIdx1*3+j-1] += Dbend;
			}
			for (int j = 1; j < 4; j++){
				Dtheta1 = 0;
				Dtheta2 = - computeDetGradient(4,j,newV1,newV2,newV3,newNB2,det2) / ((f.area*f.sideArea2)/f.l2);
				Dtheta3 = 0;

				DSOD_00 = Dtheta1 * f.SO1[0][0] + Dtheta2 * f.SO2[0][0] + Dtheta3 * f.SO3[0][0];
				DSOD_01 = Dtheta1 * f.SO1[0][1] + Dtheta2 * f.SO2[0][1] + Dtheta3 * f.SO3[0][1];
				DSOD_10 = Dtheta1 * f.SO1[1][0] + Dtheta2 * f.SO2[1][0] + Dtheta3 * f.SO3[1][0];
				DSOD_11 = Dtheta1 * f.SO1[1][1] + Dtheta2 * f.SO2[1][1] + Dtheta3 * f.SO3[1][1];

				Dbend = weight * (2*SOD_00*DSOD_00 + 2*SOD_01*DSOD_01 + 2*SOD_10*DSOD_10 + 2*SOD_11*DSOD_11);

				g[f.nbIdx2*3+j-1] += Dbend;
			}
			for (int j = 1; j < 4; j++){
				Dtheta1 = 0;
				Dtheta2 = 0;
				Dtheta3 = - computeDetGradient(4,j,newV1,newV2,newV3,newNB3,det3) / ((f.area*f.sideArea3)/f.l3);

				DSOD_00 = Dtheta1 * f.SO1[0][0] + Dtheta2 * f.SO2[0][0] + Dtheta3 * f.SO3[0][0];
				DSOD_01 = Dtheta1 * f.SO1[0][1] + Dtheta2 * f.SO2[0][1] + Dtheta3 * f.SO3[0][1];
				DSOD_10 = Dtheta1 * f.SO1[1][0] + Dtheta2 * f.SO2[1][0] + Dtheta3 * f.SO3[1][0];
				DSOD_11 = Dtheta1 * f.SO1[1][1] + Dtheta2 * f.SO2[1][1] + Dtheta3 * f.SO3[1][1];

				Dbend = weight * (2*SOD_00*DSOD_00 + 2*SOD_01*DSOD_01 + 2*SOD_10*DSOD_10 + 2*SOD_11*DSOD_11);

				g[f.nbIdx3*3+j-1] += Dbend;
			}

			delete[] det1;
			delete[] det2;
			delete[] det3;
	}
	*/

	facetBending = bending;
	return bending;
}

