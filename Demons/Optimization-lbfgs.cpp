#include "stdafx.h"
#include "GlutManager.h"
#include <omp.h>
#include <iostream>
#include <fstream>
using namespace std;

static int progress(
        void *instance,
        const lbfgsfloatval_t *tu,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
{
    pthread_mutex_lock(&lock);
    if ((k % 5 == 0) || (k <=1)){
        memcpy(u,tu,sizeof(lbfgsfloatval_t)*affinityM*3);

        cout<<"Iteration: "<<k<<" -------------------------------"<<endl;
        cout<<"fx = "<<fx<<endl;
        cout<<"xnorm = "<<xnorm<<", gnorm = "<<gnorm<<", step = "<<step<<endl;
        cout<<endl;
        cout<<"Bending: "<<facetBending<<" Stretching: "<<facetStretching<<" Link: "<<distantLink<<endl;
        cout<<thetaDif1<<' '<<thetaDif2<<' '<<thetaDif3<<endl;

        //QueryPerformanceFrequency(&tc);
        //QueryPerformanceCounter(&t2);
        //cout<<"Time Elapse:"<<(t2.QuadPart - t1.QuadPart)*1.0/tc.QuadPart * 1000<<endl;
        //cout<<endl;
        ///t1 = t2;
    }
    pthread_mutex_unlock(&lock);
    return 0;
}

void initialDVF(double* initialU){
    PolyhedralSurf::Vertex_const_iterator vb = RCSurface->P.vertices_begin();
    PolyhedralSurf::Vertex_const_iterator ve = RCSurface->P.vertices_end();
    PolyhedralSurf::Vertex_const_handle vh;

    for (int i = 0; vb != ve; vb++,i++){
        vh = vb;

        initialU[i*3]   = CTSurface->vertexIndex[i]->point().x() - vh->point().x();
        initialU[i*3+1] = CTSurface->vertexIndex[i]->point().y()  - vh->point().y();
        initialU[i*3+2] = CTSurface->vertexIndex[i]->point().z()  - vh->point().z();
    }
}

double evaluateError(double bestError){
    /* output errer */
    float finalErr,finalVar,finalBoundErr,finalBoundVar,maxErr;
    computeStatisticBool(RCSurface,CTSurface,finalErr,finalBoundErr,finalVar,finalBoundVar,maxErr);

    if (finalBoundErr < bestError){
        bestError = finalBoundErr;

        ofstream fout;
        fout.open("error.txt");

        fout<<initialError<<' '<<initialVariance<<endl;
        fout<<initialBoundError<<' '<<initialBoundVar<<endl;


        fout<<finalErr<<' '<<finalVar<<endl;
        fout<<finalBoundErr<<' '<<finalBoundVar<<endl;

        fout.close();
    }

    return bestError;
    /* output errer */
}

void *startOptimization(void * arg){
    /* initialize optimization */
    char filename[50];

    int ret;
    lbfgsfloatval_t *tempU = new lbfgsfloatval_t[affinityM*3];
    memset(tempU,0,sizeof(lbfgsfloatval_t) * affinityM * 3);

    // 	omp_set_dynamic(0);
    // 	omp_set_num_threads(2);

    double bestError = 100000;
    regweight = REGWEIGHT;
    bendweight = BENDWEIGHT;

    for (int i = 0; i < 10; i++){
        cout<<"Velocity Iteration: "<<i<<endl;

        if (i > 0){
            cout<<"Reset surface..."<<endl;

            RCSurface->ComputeMeshProperty(filename);
            bestError = evaluateError(bestError);

            cout<<"Recompute geometric features..."<<endl;
            RCSurface->findSignatureAll();

            memset(tempU,0,sizeof(lbfgsfloatval_t) * affinityM * 3);
            memset(affinity,0,affinityM * affinityN * sizeof(double));
            for (int i = 0; i < affinityM; i++) bestMatch[i]=Vector_3(0,0,0);	
            memset(matchWeight,0,sizeof(float)*affinityM);

            cout<<"Recompute affinity map..."<<endl;

            /* prepare registration */
            RCSurface->constructEdgeList();
            //RCSurface->constructLink();
            RCSurface->constructVertexList();
            RCSurface->findCorrespondenceBothWay(CTSurface,EUCLWEIGHT);
            cout<<"Recompute affinity map..."<<endl;
            RCSurface->findMatch2(CTSurface);
            //RCSurface->findClosestPoint(CTSurface);
        }

        cout<<"Begin optimization..."<<endl;    /* actual registration (optimization) */

        int opIterNum = 20000;

        lbfgs_parameter_t param;
        ret = 0;

        clock_t start, finish;
        double duration;

        start = clock();
        //QueryPerformanceCounter(&t1);

        lbfgs_parameter_init(&param);
        param.epsilon = 1e-5f;
        int CurrentIter = 0;
        lbfgsfloatval_t fx = 0;

        {
            param.max_iterations = 2000;
            double tempWeight = bendweight;
            bendweight = 0;
            ret = lbfgs(affinityM * 3, tempU, &fx, evaluate, progress, NULL, &param);

            fx = 0;
            param.max_iterations = 20000;
            bendweight = tempWeight;
            ret = lbfgs(affinityM * 3, tempU, &fx, evaluate, progress, NULL, &param);
        }


        finish = clock();
        duration = (double)(finish - start) / CLOCKS_PER_SEC;
        cout<<"function value: "<<fx<<"  Elapsed Time: "<<duration<<" s"<<endl;
        memcpy(u,tempU,sizeof(lbfgsfloatval_t)*affinityM*3);
        // What is this function???
        //if (is_orthotropic) transportFrame();

        //output file
        filename[0] = '\0';
        strcat(filename,"temp");
        char iterNum[10];
        //itoa(i,iterNum,10);
        sprintf(iterNum, "%d", i);;
        strcat(filename,iterNum);
        strcat(filename,".off");
        writeOFF(RCSurface,filename);

        cout<<"L-BFGS optimization terminated with status code: "<<ret<<"fx: "<<fx<<endl;

        if ((MESHLABOPTION > 0) && (i % MESHLABOPTION == 0) && (i > 0)) {
            char cmdLine[100];
            sprintf(cmdLine,"meshlabserver -i %s -o %s -s meshlabscript_demons.mlx\n",filename,filename);
            system(cmdLine);
        }

        regweight *= 0.95;
        bendweight *= 0.95;
        EUCLWEIGHT *= 0.8;
    }

    //final iteration evaluation
    RCSurface->ComputeMeshProperty(filename);
    bestError = evaluateError(bestError);

    return NULL;
}

static lbfgsfloatval_t evaluate(
        void *instance,
        const lbfgsfloatval_t *u,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
{
    pthread_mutex_lock(&lock);
    //cout<<"begin optimization"<<endl;
    double data = penalizeData(u,g);

    /* stretching */
    //double stretching = penalizeStretchQuadratic(u,g);
    double stretching = penalizeStretch(u,g);
    //double stretching = 0;
    /* bending */
    double bending = penalizeBendQuadratic(u,g);
    //double bending = penalizeBendAngleLP(u,g);
    //double bending = 0;

    //double landmark = penalizaeLandmark(u,g);

    //double linkStretching = penalizeLink(u,g);
    // 	cout<<data<<endl;
    // 	cout<<stretching<<endl;
    // 	cout<<bending<<endl;

    pthread_mutex_unlock(&lock);
    return data + stretching +  bending /*+landmark + linkStretching*/;
}

lbfgsfloatval_t penalizeData(const lbfgsfloatval_t *u, lbfgsfloatval_t *g){
    PolyhedralSurf::Vertex_const_iterator vb = RCSurface->P.vertices_begin();
    PolyhedralSurf::Vertex_const_iterator ve = RCSurface->P.vertices_end();
    PolyhedralSurf::Vertex_const_handle vh;

    lbfgsfloatval_t fx = 0.0;

    int i = 0;
    //#pragma omp parallel for private(i)
    for (i = 0; vb != ve; vb++,i++){
        vh = vb;

        if (occluded[i]){
            g[i*3] = 0;
            g[i*3+1] = 0;
            g[i*3+2] = 0;
            continue; // occluded region
        }

        Point_3 deformed = vh->point() + Vector_3(u[i*3],u[i*3+1],u[i*3+2]);
        Point_3 cor = Point_3(bestMatch[i].x(),bestMatch[i].y(),bestMatch[i].z());
        double weight = matchWeight[i];

        Vector_3 dis(deformed,cor);

        g[i*3] = -2 * dis.x() * weight;
        g[i*3+1] = -2 * dis.y() * weight;
        g[i*3+2] = -2 * dis.z() * weight;
        fx = fx + weight * dis.squared_length();
    }

    return fx;
}

lbfgsfloatval_t penalizeStretchQuadratic(const lbfgsfloatval_t *u, lbfgsfloatval_t *g){
    double stretching = 0;

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
    facetStretching = stretching;
    return stretching;
}
lbfgsfloatval_t penalizeStretch(const lbfgsfloatval_t *u, lbfgsfloatval_t *g){
    double stretching = 0;

    { //triangle-based stretching energy
        facet* faceList = RCSurface->faceList;

        int idx;
        //#pragma omp parallel for private(idx)
        for ( idx=0; idx < RCSurface->faceNum;idx++){
            facet f = faceList[idx];
            double weight = f.area * regweight;

            Point_3 newV1 = f.v1->point() + Vector_3(u[f.index[1]*3],u[f.index[1]*3+1],u[f.index[1]*3+2]);
            Point_3 newV2 = f.v2->point() + Vector_3(u[f.index[2]*3],u[f.index[2]*3+1],u[f.index[2]*3+2]);
            Point_3 newV3 = f.v3->point() + Vector_3(u[f.index[3]*3],u[f.index[3]*3+1],u[f.index[3]*3+2]);

            Vector_3 temp_v1 = Vector_3(newV1,newV2);
            Vector_3 temp_v2 = Vector_3(newV1,newV3);

            double v11 = 0;
            double v12 = sqrt(temp_v1.squared_length()+EPS);

            double k = temp_v1*temp_v2;
            double v22 = k / v12;

            double sqv21 = temp_v2.squared_length() - v22*v22;

            double v21 = - sqrt(sqv21+EPS);

            double J11 = 0*f.inverse[0][0] + v21*f.inverse[1][0];
            double J12 = 0*f.inverse[0][1] + v21*f.inverse[1][1];
            double J21 = v12*f.inverse[0][0] + v22*f.inverse[1][0];
            double J22 = v12*f.inverse[0][1] + v22*f.inverse[1][1];

            double S11 = J11*J11 + J21*J21 - 1; double S12 = J11*J12+J21*J22; double S21 = S12; double S22 = J12*J12 + J22*J22 - 1;

            double trace,trace_2;
            if (is_orthotropic){
                stretching += weight * (RCSurface->C[idx][0] * S11*S11 + 2 * RCSurface->C[idx][1] * S11*S22 + RCSurface->C[idx][2] * S22*S22 + RCSurface->C[idx][3] * S12*S12);
            }else{
                trace = (S11 + S22);
                trace_2 = S11*S11 + 2*S12*S21 + S22*S22;
                stretching += weight *YOUNG/(1-POISSON*POISSON)*((1-POISSON)*trace_2+POISSON*trace*trace);
            }

            /* COMPUTE GRADIENT */
            for (int i = 1; i < 4; i++)
                for (int j = 1; j < 4; j++){
                    double Dtemp_v1_x = computeStretchGradient(1,1,i,j);double Dtemp_v1_y = computeStretchGradient(1,2,i,j);double Dtemp_v1_z = computeStretchGradient(1,3,i,j);
                    double Dtemp_v2_x = computeStretchGradient(2,1,i,j);double Dtemp_v2_y = computeStretchGradient(2,2,i,j);double Dtemp_v2_z = computeStretchGradient(2,3,i,j);

                    double dk = (Dtemp_v1_x*temp_v2.x()  +temp_v1.x()*Dtemp_v2_x) + (Dtemp_v1_y*temp_v2.y()  +temp_v1.y()*Dtemp_v2_y) + (Dtemp_v1_z*temp_v2.z()  +temp_v1.z()*Dtemp_v2_z);
                    double Dv12 = 0.5 * 1 / v12 * (2*temp_v1.x()*Dtemp_v1_x + 2*temp_v1.y()*Dtemp_v1_y + 2*temp_v1.z()*Dtemp_v1_z);
                    double Dv22 = (dk * v12 - Dv12 * k) / (v12*v12);
                    double Dv21 = 0.5 * 1 / v21 * (2*temp_v2.x()*Dtemp_v2_x + 2*temp_v2.y()*Dtemp_v2_y + 2*temp_v2.z()*Dtemp_v2_z - 2*v22*Dv22);

                    double DJ11 = Dv21*f.inverse[1][0];
                    double DJ12 = Dv21*f.inverse[1][1];
                    double DJ21 = Dv12*f.inverse[0][0] + Dv22*f.inverse[1][0];
                    double DJ22 = Dv12*f.inverse[0][1] + Dv22*f.inverse[1][1];

                    double DS11 = 2*J11*DJ11 + 2*J21*DJ21;
                    double DS22 = 2*J12*DJ12 + 2*J22*DJ22;
                    double DS12 = (DJ11*J12+J11*DJ12) + (DJ21*J22+J21*DJ22);
                    double DS21 = DS12;

                    if (is_orthotropic){
                        g[f.index[i]*3+j-1] += weight * (2*RCSurface->C[idx][0] * S11*DS11 + 
                                2*RCSurface->C[idx][1] * DS11*S22 + 
                                2*RCSurface->C[idx][1] * S11*DS22 +  
                                2*RCSurface->C[idx][2] * S22*DS22 + 
                                2*RCSurface->C[idx][3] * S12*DS12);
                    }else{
                        g[f.index[i]*3+j-1] += weight * YOUNG/(1-POISSON*POISSON)*
                            ((1-POISSON)*(2*S11*DS11+2*S12*DS21+2*DS12*S21+2*S22*DS22)+
                             POISSON*2*trace*(DS11+DS22));
                    }
                }


        }	
    }

    facetStretching = stretching;
    return stretching;
}

lbfgsfloatval_t penalizeBendAngle(const lbfgsfloatval_t *u, lbfgsfloatval_t *g){
    double bending = 0;
    facet* faceList = RCSurface->faceList;

    Point_3 newV1,newV2,newV3,newNB1,newNB2,newNB3;
    double theta1,theta2,theta3;
    threeTuple db1,db2,db3,ds1,ds2,ds3,dnm,dn1,dn2,dn3;
    double dx1,dx2,dx3,Dtheta1,Dtheta2,Dtheta3,DSOD_00,DSOD_01,DSOD_10,DSOD_11;
    double weight,Dbend;

    for (int idx=0; idx < RCSurface->faceNum;idx++){
        facet f = faceList[idx];
        if (f.isBorder) continue;
        weight = f.area / RCSurface->averageArea * bendweight;

        newV1 = f.v1->point() + Vector_3(u[f.index[1]*3],u[f.index[1]*3+1],u[f.index[1]*3+2]);
        newV2 = f.v2->point() + Vector_3(u[f.index[2]*3],u[f.index[2]*3+1],u[f.index[2]*3+2]);
        newV3 = f.v3->point() + Vector_3(u[f.index[3]*3],u[f.index[3]*3+1],u[f.index[3]*3+2]);
        newNB1 = f.nb1->point() + Vector_3(u[f.nbIdx1*3],u[f.nbIdx1*3+1],u[f.nbIdx1*3+2]);
        newNB2 = f.nb2->point() + Vector_3(u[f.nbIdx2*3],u[f.nbIdx2*3+1],u[f.nbIdx2*3+2]);
        newNB3 = f.nb3->point() + Vector_3(u[f.nbIdx3*3],u[f.nbIdx3*3+1],u[f.nbIdx3*3+2]);

        /* -------------- compute angle ------------------------- */
        /* ------------------------------------------------------ */
        Vector_3 b1(newV1,newV2);
        Vector_3 b2(newV2,newV3);
        Vector_3 b3(newV3,newV1);
        Vector_3 s1(newV1,newNB1);
        Vector_3 s2(newV2,newNB2);
        Vector_3 s3(newV3,newNB3);

        Vector_3 nm = cross_product(b1,b2);
        Vector_3 n1 = cross_product(s1,b1);
        Vector_3 n2 = cross_product(s2,b2);
        Vector_3 n3 = cross_product(s3,b3);

        double sign1 = s1*nm;
        double sign2 = s2*nm;
        double sign3 = s3*nm;

        double x1 = n1*nm / (sqrt(n1.squared_length()+EPS)*sqrt(nm.squared_length()+EPS));
        double x2 = n2*nm / (sqrt(n2.squared_length()+EPS)*sqrt(nm.squared_length()+EPS));
        double x3 = n3*nm / (sqrt(n3.squared_length()+EPS)*sqrt(nm.squared_length()+EPS));

        theta1 = acos(x1);
        theta2 = acos(x2);
        theta3 = acos(x3);

        if (sign1 > 0) theta1 = -theta1;
        if (sign2 > 0) theta2 = -theta2;
        if (sign3 > 0) theta3 = -theta3;

        double SOD_00 = (theta1-f.theta1) * f.SO1[0][0] + (theta2-f.theta2) * f.SO2[0][0] + (theta3-f.theta3) * f.SO3[0][0];
        double SOD_01 = (theta1-f.theta1) * f.SO1[0][1] + (theta2-f.theta2) * f.SO2[0][1] + (theta3-f.theta3) * f.SO3[0][1];
        double SOD_10 = (theta1-f.theta1) * f.SO1[1][0] + (theta2-f.theta2) * f.SO2[1][0] + (theta3-f.theta3) * f.SO3[1][0];
        double SOD_11 = (theta1-f.theta1) * f.SO1[1][1] + (theta2-f.theta2) * f.SO2[1][1] + (theta3-f.theta3) * f.SO3[1][1];

        if (is_orthotropic)
            bending += weight * (RCSurface->C[idx][0] * SOD_00*SOD_00 + 
                    2 * RCSurface->C[idx][1] * SOD_00*SOD_11 + 
                    RCSurface->C[idx][2] * SOD_11*SOD_11 + 
                    RCSurface->C[idx][3] * SOD_01*SOD_10); 
        else
            bending += weight * (SOD_00*SOD_00 + SOD_01*SOD_01 + SOD_10*SOD_10 + SOD_11*SOD_11); 

        /* -------------- compute gradient ------------------------- */
        /* --------------------------------------------------------- */
        // 		Dtheta1 = computeDtheta(newV1,newV2,newV3,newNB1,1);
        // 		Dtheta2 = computeDtheta(newV2,newV3,newV1,newNB2,3);
        // 		Dtheta3 = computeDtheta(newV3,newV1,newV2,newNB3,2);
        // 		upgradeGradient();

    }
}

lbfgsfloatval_t penalizeBendAngleLP(const lbfgsfloatval_t *u, lbfgsfloatval_t *g){
    double bending = 0;
    facet* faceList = RCSurface->faceList;

    Point_3 newV1,newV2,newV3,newNB1,newNB2,newNB3;
    double theta1,theta2,theta3;
    threeTuple db1,db2,db3,ds1,ds2,ds3,dnm,dn1,dn2,dn3;
    double dx1,dx2,dx3,Dtheta1,Dtheta2,Dtheta3,DSOD_00,DSOD_01,DSOD_10,DSOD_11;
    double weight,Dbend,trace,trace_2;

    for (int idx=0; idx < RCSurface->faceNum;idx++){
        facet f = faceList[idx];
        if (f.isBorder) continue;
        weight = f.area / RCSurface->averageArea * bendweight;

        newV1 = f.v1->point() + Vector_3(u[f.index[1]*3],u[f.index[1]*3+1],u[f.index[1]*3+2]);
        newV2 = f.v2->point() + Vector_3(u[f.index[2]*3],u[f.index[2]*3+1],u[f.index[2]*3+2]);
        newV3 = f.v3->point() + Vector_3(u[f.index[3]*3],u[f.index[3]*3+1],u[f.index[3]*3+2]);
        newNB1 = f.nb1->point() + Vector_3(u[f.nbIdx1*3],u[f.nbIdx1*3+1],u[f.nbIdx1*3+2]);
        newNB2 = f.nb2->point() + Vector_3(u[f.nbIdx2*3],u[f.nbIdx2*3+1],u[f.nbIdx2*3+2]);
        newNB3 = f.nb3->point() + Vector_3(u[f.nbIdx3*3],u[f.nbIdx3*3+1],u[f.nbIdx3*3+2]);

        /* -------------- compute angle ------------------------- */
        /* ------------------------------------------------------ */
        Vector_3 b1(newV1,newV2);
        Vector_3 b2(newV2,newV3);
        Vector_3 b3(newV3,newV1);
        Vector_3 s1(newV1,newNB1);
        Vector_3 s2(newV2,newNB2);
        Vector_3 s3(newV3,newNB3);

        Vector_3 nm = cross_product(b1,b2);
        Vector_3 n1 = cross_product(s1,b1);
        Vector_3 n2 = cross_product(s2,b2);
        Vector_3 n3 = cross_product(s3,b3);

        double sign1 = s1*nm;
        double sign2 = s2*nm;
        double sign3 = s3*nm;

        double x1 = n1*nm / (sqrt(n1.squared_length()+EPS)*sqrt(nm.squared_length()+EPS));
        double x2 = n2*nm / (sqrt(n2.squared_length()+EPS)*sqrt(nm.squared_length()+EPS));
        double x3 = n3*nm / (sqrt(n3.squared_length()+EPS)*sqrt(nm.squared_length()+EPS));

        if (sign1 > 0) x1 = (2 - x1);
        if (sign2 > 0) x2 = (2 - x2);
        if (sign3 > 0) x3 = (2 - x3);

        if (x1 <= 1) theta1 = (-0.69813170079773212 * x1 * x1 - 0.87266462599716477) * x1 + 1.5707963267948966;
        else	     theta1 = (-0.69813170079773212 * (x1-2) * (x1-2) - 0.87266462599716477) * (x1-2) - 1.5707963267948966;
        if (x2 <= 1) theta2 = (-0.69813170079773212 * x2 * x2 - 0.87266462599716477) * x2 + 1.5707963267948966;
        else	     theta2 = (-0.69813170079773212 * (x2-2) * (x2-2) - 0.87266462599716477) * (x2-2) - 1.5707963267948966;
        if (x3 <= 1) theta3 = (-0.69813170079773212 * x3 * x3 - 0.87266462599716477) * x3 + 1.5707963267948966;
        else	     theta3 = (-0.69813170079773212 * (x3-2) * (x3-2) - 0.87266462599716477) * (x3-2) - 1.5707963267948966;

        /* */
        double SOD_00 = (theta1-f.theta1) * f.SO1[0][0] + (theta2-f.theta2) * f.SO2[0][0] + (theta3-f.theta3) * f.SO3[0][0];
        double SOD_01 = (theta1-f.theta1) * f.SO1[0][1] + (theta2-f.theta2) * f.SO2[0][1] + (theta3-f.theta3) * f.SO3[0][1];
        double SOD_10 = (theta1-f.theta1) * f.SO1[1][0] + (theta2-f.theta2) * f.SO2[1][0] + (theta3-f.theta3) * f.SO3[1][0];
        double SOD_11 = (theta1-f.theta1) * f.SO1[1][1] + (theta2-f.theta2) * f.SO2[1][1] + (theta3-f.theta3) * f.SO3[1][1];

        if (is_orthotropic)
            bending += weight * (RCSurface->C[idx][0] * SOD_00*SOD_00 + 
                    2 * RCSurface->C[idx][1] * SOD_00*SOD_11 + 
                    RCSurface->C[idx][2] * SOD_11*SOD_11 + 
                    RCSurface->C[idx][3] * SOD_01*SOD_10); 
        else{
            trace = (SOD_00 + SOD_11);
            trace_2 = SOD_00*SOD_00 + 2*SOD_01*SOD_10 + SOD_11*SOD_11;
            bending += weight *YOUNG/(1-POISSON*POISSON)*((1-POISSON)*trace_2+POISSON*trace*trace);
            //bending += weight * (SOD_00*SOD_00 + SOD_01*SOD_01 + SOD_10*SOD_10 + SOD_11*SOD_11); 
        }
        /* -------------- compute gradient ------------------------- */
        /* --------------------------------------------------------- */

        for (int i = 1; i < 7; i++)
            for (int j = 1; j < 4; j++){
                /* */
                db1 = computeGradient(1,i,j);db2 = computeGradient(2,i,j);db3 = computeGradient(3,i,j);
                ds1 = computeGradient(4,i,j);ds2 = computeGradient(5,i,j);ds3 = computeGradient(6,i,j);

                dnm = computeCrossProductGradient(b1,b2,db1,db2);
                dn1 = computeCrossProductGradient(s1,b1,ds1,db1);
                dn2 = computeCrossProductGradient(s2,b2,ds2,db2);
                dn3 = computeCrossProductGradient(s3,b3,ds3,db3);

                dx1 = computeDotProductGradient(n1,nm,dn1,dnm);
                dx2 = computeDotProductGradient(n2,nm,dn2,dnm);
                dx3 = computeDotProductGradient(n3,nm,dn3,dnm);

                if (sign1 > 0) dx1 = -dx1;
                if (sign2 > 0) dx2 = -dx2;
                if (sign3 > 0) dx3 = -dx3;

                if (x1 <= 1) Dtheta1 = (-3 * 0.69813170079773212 * x1 * x1 - 0.87266462599716477)  * dx1;
                else		 Dtheta1 = (-3 * 0.69813170079773212 * (x1-2) * (x1-2) - 0.87266462599716477)  * dx1;
                if (x2 <= 1) Dtheta2 = (-3 * 0.69813170079773212 * x2 * x2 - 0.87266462599716477)  * dx2;
                else		 Dtheta2 = (-3 * 0.69813170079773212 * (x2-2) * (x2-2) - 0.87266462599716477)  * dx2;
                if (x3 <= 1) Dtheta3 = (-3 * 0.69813170079773212 * x3 * x3 - 0.87266462599716477)  * dx3;
                else		 Dtheta3 = (-3 * 0.69813170079773212 * (x3-2) * (x3-2) - 0.87266462599716477)  * dx3;

                DSOD_00 = Dtheta1 * f.SO1[0][0] + Dtheta2 * f.SO2[0][0] + Dtheta3 * f.SO3[0][0];
                DSOD_01 = Dtheta1 * f.SO1[0][1] + Dtheta2 * f.SO2[0][1] + Dtheta3 * f.SO3[0][1];
                DSOD_10 = Dtheta1 * f.SO1[1][0] + Dtheta2 * f.SO2[1][0] + Dtheta3 * f.SO3[1][0];
                DSOD_11 = Dtheta1 * f.SO1[1][1] + Dtheta2 * f.SO2[1][1] + Dtheta3 * f.SO3[1][1];

                if (is_orthotropic){
                    Dbend = weight * (2 * RCSurface->C[idx][0] * SOD_00*DSOD_00 + 
                            2 * RCSurface->C[idx][1] * SOD_00*DSOD_11 + 
                            2 * RCSurface->C[idx][1] * DSOD_00*SOD_11 + 
                            2 * RCSurface->C[idx][2] * SOD_11*DSOD_11 + 
                            2 * RCSurface->C[idx][3] * SOD_01*DSOD_10); 
                }else{
                    Dbend = weight * YOUNG/(1-POISSON*POISSON)*
                        ((1-POISSON)*(2*SOD_00*DSOD_00+2*SOD_01*DSOD_10+2*DSOD_01*SOD_10+2*SOD_11*DSOD_11)+
                         POISSON*2*trace*(DSOD_00+DSOD_11));
                }
                //Dbend = weight * (2*SOD_00*DSOD_00 + 2*SOD_01*DSOD_01 + 2*SOD_10*DSOD_10 + 2*SOD_11*DSOD_11);

                if (i <= 3)
                    g[f.index[i]*3+j-1] += Dbend;
                else if (i == 4)
                    g[f.nbIdx1*3+j-1] += Dbend;
                else if (i == 5)
                    g[f.nbIdx2*3+j-1] += Dbend;
                else if (i == 6)
                    g[f.nbIdx3*3+j-1] += Dbend;
            }

    }

    facetBending = bending;
    return bending;
}

lbfgsfloatval_t penalizeBendDet(const lbfgsfloatval_t *u, lbfgsfloatval_t *g){
    double bending = 0;
    facet* faceList = RCSurface->faceList;

    Point_3 newV1,newV2,newV3,newNB1,newNB2,newNB3;
    double dv1x, dv1y, dv1z, dv2x, dv2y, dv2z, dv3x, dv3y, dv3z,
           dnb1x,dnb1y,dnb1z,dnb2x,dnb2y,dnb2z,dnb3x,dnb3y,dnb3z;
    double Dtheta1,Dtheta2,Dtheta3,DSOD_00,DSOD_01,DSOD_10,DSOD_11;
    double weight,Dbend;

    for (int idx=0; idx < RCSurface->faceNum;idx++){
        facet f = faceList[idx];
        if (f.isBorder) continue;
        weight = f.area / RCSurface->averageArea * bendweight;

        newV1 = f.v1->point() + Vector_3(u[f.index[1]*3],u[f.index[1]*3+1],u[f.index[1]*3+2]);
        newV2 = f.v2->point() + Vector_3(u[f.index[2]*3],u[f.index[2]*3+1],u[f.index[2]*3+2]);
        newV3 = f.v3->point() + Vector_3(u[f.index[3]*3],u[f.index[3]*3+1],u[f.index[3]*3+2]);
        newNB1 = f.nb1->point() + Vector_3(u[f.nbIdx1*3],u[f.nbIdx1*3+1],u[f.nbIdx1*3+2]);
        newNB2 = f.nb2->point() + Vector_3(u[f.nbIdx2*3],u[f.nbIdx2*3+1],u[f.nbIdx2*3+2]);
        newNB3 = f.nb3->point() + Vector_3(u[f.nbIdx3*3],u[f.nbIdx3*3+1],u[f.nbIdx3*3+2]);

        /* */
        double* det1 = computeDeterminant(newV1,newV2,newV3,newNB1);
        double* det2 = computeDeterminant(newV1,newV2,newV3,newNB2);
        double* det3 = computeDeterminant(newV1,newV2,newV3,newNB3);
        double theta1 = - det1[0] / ((f.area*f.sideArea1)/f.l1);
        double theta2 = - det2[0] / ((f.area*f.sideArea2)/f.l2);
        double theta3 = - det3[0] / ((f.area*f.sideArea3)/f.l3);
        /* */
        /* angle based 
           double theta1 = computeAngle(newV1,newV2,newV3,newNB1);
           double theta2 = computeAngle(newV2,newV3,newV1,newNB2);
           double theta3 = computeAngle(newV3,newV1,newV2,newNB3);
           */

        double SOD_00 = (theta1-f.theta1) * f.SO1[0][0] + (theta2-f.theta2) * f.SO2[0][0] + (theta3-f.theta3) * f.SO3[0][0];
        double SOD_01 = (theta1-f.theta1) * f.SO1[0][1] + (theta2-f.theta2) * f.SO2[0][1] + (theta3-f.theta3) * f.SO3[0][1];
        double SOD_10 = (theta1-f.theta1) * f.SO1[1][0] + (theta2-f.theta2) * f.SO2[1][0] + (theta3-f.theta3) * f.SO3[1][0];
        double SOD_11 = (theta1-f.theta1) * f.SO1[1][1] + (theta2-f.theta2) * f.SO2[1][1] + (theta3-f.theta3) * f.SO3[1][1];

        if (is_orthotropic)
            bending += weight * (RCSurface->C[idx][0] * SOD_00*SOD_00 + 
                    2 * RCSurface->C[idx][1] * SOD_00*SOD_11 + 
                    RCSurface->C[idx][2] * SOD_11*SOD_11 + 
                    RCSurface->C[idx][3] * SOD_01*SOD_10); 
        else
            bending += weight * (SOD_00*SOD_00 + SOD_01*SOD_01 + SOD_10*SOD_10 + SOD_11*SOD_11); 

        //COMPUTE GRADIENT

        for (int i = 1; i < 4; i++)
            for (int j = 1; j < 4; j++){
                /* */
                Dtheta1 = - computeDetGradient(i,j,newV1,newV2,newV3,newNB1,det1) / ((f.area*f.sideArea1)/f.l1);
                Dtheta2 = - computeDetGradient(i,j,newV1,newV2,newV3,newNB2,det2) / ((f.area*f.sideArea2)/f.l2);
                Dtheta3 = - computeDetGradient(i,j,newV1,newV2,newV3,newNB3,det3) / ((f.area*f.sideArea3)/f.l3);

                DSOD_00 = Dtheta1 * f.SO1[0][0] + Dtheta2 * f.SO2[0][0] + Dtheta3 * f.SO3[0][0];
                DSOD_01 = Dtheta1 * f.SO1[0][1] + Dtheta2 * f.SO2[0][1] + Dtheta3 * f.SO3[0][1];
                DSOD_10 = Dtheta1 * f.SO1[1][0] + Dtheta2 * f.SO2[1][0] + Dtheta3 * f.SO3[1][0];
                DSOD_11 = Dtheta1 * f.SO1[1][1] + Dtheta2 * f.SO2[1][1] + Dtheta3 * f.SO3[1][1];

                if (is_orthotropic){
                    Dbend = weight * (2 * RCSurface->C[idx][0] * SOD_00*DSOD_00 + 
                            2 * RCSurface->C[idx][1] * SOD_00*DSOD_11 + 
                            2 * RCSurface->C[idx][1] * DSOD_00*SOD_11 + 
                            2 * RCSurface->C[idx][2] * SOD_11*DSOD_11 + 
                            2 * RCSurface->C[idx][3] * SOD_01*DSOD_10); 
                }else
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

            if (is_orthotropic){
                Dbend = weight * (2 * RCSurface->C[idx][0] * SOD_00*DSOD_00 + 
                        2 * RCSurface->C[idx][1] * SOD_00*DSOD_11 + 
                        2 * RCSurface->C[idx][1] * DSOD_00*SOD_11 + 
                        2 * RCSurface->C[idx][2] * SOD_11*DSOD_11 + 
                        2 * RCSurface->C[idx][3] * SOD_01*DSOD_10); 
            }else
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

            if (is_orthotropic){
                Dbend = weight * (2 * RCSurface->C[idx][0] * SOD_00*DSOD_00 + 
                        2 * RCSurface->C[idx][1] * SOD_00*DSOD_11 + 
                        2 * RCSurface->C[idx][1] * DSOD_00*SOD_11 + 
                        2 * RCSurface->C[idx][2] * SOD_11*DSOD_11 + 
                        2 * RCSurface->C[idx][3] * SOD_01*DSOD_10); 
            }else
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

            if (is_orthotropic){
                Dbend = weight * (2 * RCSurface->C[idx][0] * SOD_00*DSOD_00 + 
                        2 * RCSurface->C[idx][1] * SOD_00*DSOD_11 + 
                        2 * RCSurface->C[idx][1] * DSOD_00*SOD_11 + 
                        2 * RCSurface->C[idx][2] * SOD_11*DSOD_11 + 
                        2 * RCSurface->C[idx][3] * SOD_01*DSOD_10); 
            }else
                Dbend = weight * (2*SOD_00*DSOD_00 + 2*SOD_01*DSOD_01 + 2*SOD_10*DSOD_10 + 2*SOD_11*DSOD_11);

            g[f.nbIdx3*3+j-1] += Dbend;
        }

        delete[] det1;
        delete[] det2;
        delete[] det3;
    }

    facetBending = bending;
    return bending;
}

lbfgsfloatval_t penalizeBendQuadratic(const lbfgsfloatval_t *tu, lbfgsfloatval_t *g){
    PolyhedralSurf::Vertex_const_iterator vb = RCSurface->P.vertices_begin();
    PolyhedralSurf::Vertex_const_iterator ve = RCSurface->P.vertices_end();
    PolyhedralSurf::Vertex_const_handle vh,nb;
    vertex* vl = RCSurface->vertexList;

    double fv = 0;
    int idx1,idx2;
    double laplacian,l1,l2,l3;
    for (int i = 0; vb != ve; vb++,i++){
        vh = vb;

        l1 = 0;l2 = 0;l3 = 0;

        for (int j = 0; j < vl[i].neighbourNum; j++){
            idx1 = vl[i].idx;
            idx2 = vl[i].nbIdx[j];

            l1 += vl[i].weight[j]*(tu[idx2*3] -   tu[idx1*3]);
            l2 += vl[i].weight[j]*(tu[idx2*3+1] - tu[idx1*3+1]);
            l3 += vl[i].weight[j]*(tu[idx2*3+2] - tu[idx1*3+2]);
        }

        double weight = bendweight / (vl[i].area*vl[i].area);
        fv += weight * (l1*l1 + l2*l2 + l3*l3);

        for (int j = 0; j < vl[i].neighbourNum; j++){
            idx1 = vl[i].idx;
            idx2 = vl[i].nbIdx[j];

            g[idx2*3]   += weight * 2 * l1 * vl[i].weight[j];
            g[idx2*3+1] += weight * 2 * l2 * vl[i].weight[j];
            g[idx2*3+2] += weight * 2 * l3 * vl[i].weight[j];

            g[idx1*3]   -= weight * 2 * l1 * vl[i].weight[j];
            g[idx1*3+1] -= weight * 2 * l2 * vl[i].weight[j];
            g[idx1*3+2] -= weight * 2 * l3 * vl[i].weight[j];
        }
    }

    facetBending = fv;
    return fv;
}

lbfgsfloatval_t penalizeLink(const lbfgsfloatval_t *u, lbfgsfloatval_t *g){
    // 	float linkEnergy = 0;
    // 	float linkWeight = 50;
    // 	
    // 	for (int i = 0; i < affinityM; i++){
    // 		if (!attractor[i]) continue;
    // 
    // 		int idx1 = i;
    // 		int idx2 = linkTarget[i];
    // 
    // 		Point_3 P1 = RCSurface->vertexIndex[idx1]->point() + Vector_3(u[idx1*3],u[idx1*3+1],u[idx1*3+2]);
    // 		Point_3 P2 = RCSurface->vertexIndex[idx2]->point() + Vector_3(u[idx2*3],u[idx2*3+1],u[idx2*3+2]);
    // 		Vector_3 dis(P1,P2);
    // 
    // 		float linkLength = sqrt(dis.squared_length());
    // 		if (linkLength < 0.00001) return 0;
    // 		float a = 4;
    // 		if (linkLength < a){
    // 			linkEnergy += linkWeight * (a*a/6) * (1-pow((1-linkLength*linkLength/a/a),3));
    // 
    // 			float dB = linkWeight * linkLength * pow((1-linkLength*linkLength/a/a),2);
    // 			g[idx1*3] += dB * (- 0.5* 1 / linkLength * dis.x());
    // 			g[idx1*3+1] += dB * (- 0.5* 1 / linkLength * dis.y());
    // 			g[idx1*3+2] += dB * (- 0.5* 1 / linkLength * dis.z());
    // 
    // 			g[idx2*3] += dB * (0.5* 1 / linkLength * dis.x());
    // 			g[idx2*3+1] += dB * (0.5* 1 / linkLength * dis.y());
    // 			g[idx2*3+2] += dB * (0.5* 1 / linkLength * dis.z());
    // 		}
    // 		else{ 
    // 			linkEnergy += linkWeight * (a*a/6);
    // 		}
    // 	}
    // 	
    // 	distantLink = linkEnergy;
    float linkEnergy = 0;
    float linkWeight = 50;

    for (int i = 0; i < RCSurface->linkNum; i++){
        int idx1 = linkList[i].index1;
        int idx2 = linkList[i].index2;

        Point_3 P1 = RCSurface->vertexIndex[idx1]->point() + Vector_3(u[idx1*3],u[idx1*3+1],u[idx1*3+2]);
        Point_3 P2 = RCSurface->vertexIndex[idx2]->point() + Vector_3(u[idx2*3],u[idx2*3+1],u[idx2*3+2]);
        Vector_3 dis(P1,P2);

        double inc = sqrt(dis.squared_length()) - linkList[i].length;

        linkEnergy += linkWeight * inc * inc;

        double dB = linkWeight * 2 * inc * (1/sqrt(dis.squared_length()));
        g[idx1*3] +=   dB * 2*dis.x()* -1;
        g[idx1*3+1] += dB * 2*dis.y()* -1;
        g[idx1*3+2] += dB * 2*dis.z()* -1;

        g[idx2*3] +=   dB * 2*dis.x();
        g[idx2*3+1] += dB * 2*dis.y();
        g[idx2*3+2] += dB * 2*dis.z();
    }

    return linkEnergy;
}

lbfgsfloatval_t penalizaeLandmark(const lbfgsfloatval_t *u, lbfgsfloatval_t *g){
    PolyhedralSurf::Vertex_const_handle vh;
    double weight = LANDMARKWEIGHT;

    lbfgsfloatval_t fx = 0.0;
    for (int i = 0; i<RCSurface->landmarkNum; i++){
        int idx = RCSurface->landmark[i];
        vh = RCSurface->vertexIndex[idx];
        Point_3 deformed = vh->point() + Vector_3(u[idx*3],u[idx*3+1],u[idx*3+2]);
        Vector_3 dis(deformed,CTSurface->vertexIndex[CTSurface->landmark[i]]->point());

        fx = fx + weight*dis.squared_length();
        g[idx*3] += -2 * dis.x() * weight;
        g[idx*3+1] += -2 * dis.y() * weight;
        g[idx*3+2] += -2 * dis.z() * weight;

    }

    return fx;
}

inline double computeStretchGradient(int n,int m,int x,int y){
    if (x == 1){
        if (m != y) return 0; else return -1;
    }else if (x == 2){
        if (n == 2) return 0;
        if (m != y) return 0; else return 1;
    }if (x == 3){
        if (n == 1) return 0;
        if (m != y) return 0; else return 1;
    }
}

inline threeTuple computeGradient(int v,int i,int j){
    if (v <= 3){ // computing derivatives w.r.t. b1,b2,b3
        if (i>3) return threeTuple(0,0,0);

        if (v == i){
            if (j == 1) return threeTuple(-1,0,0);
            if (j == 2) return threeTuple(0,-1,0);
            if (j == 3) return threeTuple(0,0,-1);
        }

        if ((v+1 == i) || (v-2 == i)){
            if (j == 1) return threeTuple(1,0,0);
            if (j == 2) return threeTuple(0,1,0);
            if (j == 3) return threeTuple(0,0,1);
        }

        return threeTuple(0,0,0);
    }

    if (v-3 == i){
        if (j == 1) return threeTuple(-1,0,0);
        if (j == 2) return threeTuple(0,-1,0);
        if (j == 3) return threeTuple(0,0,-1);
    }

    if (v == i){
        if (j == 1) return threeTuple(1,0,0);
        if (j == 2) return threeTuple(0,1,0);
        if (j == 3) return threeTuple(0,0,1);
    }

    return threeTuple(0,0,0);
}

inline threeTuple computeCrossProductGradient(Vector_3 b1,Vector_3 b2,threeTuple db1,threeTuple db2){
    double dx = (db1.y*b2.z()+b1.y()*db2.z)-(db1.z*b2.y()+b1.z()*db2.y);
    double dy = (db1.z*b2.x()+b1.z()*db2.x)-(db1.x*b2.z()+b1.x()*db2.z);
    double dz = (db1.x*b2.y()+b1.x()*db2.y)-(db1.y*b2.x()+b1.y()*db2.x);

    return threeTuple(dx,dy,dz);
}

inline double computeDotProductGradient(Vector_3 n1,Vector_3 n2,threeTuple dn1,threeTuple dn2){
    double g = n1*n2;
    double dg = dn1.x*n2.x() + n1.x()*dn2.x +
        dn1.y*n2.y() + n1.y()*dn2.y +
        dn1.z*n2.z() + n1.z()*dn2.z;

    double N1 = n1.squared_length() + EPS;
    double N2 = n2.squared_length() + EPS;
    double N1N2 = N1*N2;

    double dN1 = 2*n1.x()*dn1.x + 2*n1.y()*dn1.y + 2*n1.z()*dn1.z;
    double dN2 = 2*n2.x()*dn2.x + 2*n2.y()*dn2.y + 2*n2.z()*dn2.z;

    double dN1N2 = dN1*N2 + N1*dN2;

    double dSN1N2 = 0.5*dN1N2/sqrt(N1N2);

    return (sqrt(N1N2)*dg - dSN1N2*g)/(N1N2);
}

double computeDetGradient(int i,int j,Point_3 v1,Point_3 v2,Point_3 v3,Point_3 v4,double* det){
    double D1,D2,D3,D4,Ddet;

    if (i == 1){
        switch (j)
        {
            case 1:
                Ddet = -det[2];
                break;
            case 2:
                Ddet = det[3];
                break;
            case 3:
                Ddet = -det[4];
        }
    }

    if (i == 2){
        switch (j)
        {
            case 1:
                D1 = v3.y() * v4.z() -v3.z() * v4.y(); 
                D2 = 0;
                D3 = v3.z() - v4.z();
                D4 = v3.y() - v4.y();
                Ddet = D1 - v1.x()*D2 + v1.y()*D3 - v1.z()*D4;
                break;
            case 2:
                D1 = v3.z() * v4.x() - v3.x() * v4.z();
                D2 = v3.z() - v4.z();
                D3 = 0;
                D4 = v4.x() - v3.x();
                Ddet = D1 - v1.x()*D2 + v1.y()*D3 - v1.z()*D4;
                break;

            case 3:
                D1 = v3.x() * v4.y() - v3.y() * v4.x();
                D2 = v4.y() - v3.y();
                D3 = v4.x() - v3.x();
                D4 = 0;
                Ddet = D1 - v1.x()*D2 + v1.y()*D3 - v1.z()*D4;
        }
    }

    if (i == 3){
        switch (j)
        {
            case 1:
                D1 = v2.z() * v4.y() - v2.y() * v4.z();
                D2 = 0;
                D3 = v4.z() - v2.z();
                D4 = v4.y() - v2.y();
                Ddet = D1 - v1.x()*D2 + v1.y()*D3 - v1.z()*D4;
                break;

            case 2:
                D1 = v2.x() * v4.z() - v2.z() * v4.x();
                D2 = v4.z() - v2.z();
                D3 = 0;
                D4 = v2.x() - v4.x();
                Ddet = D1 - v1.x()*D2 + v1.y()*D3 - v1.z()*D4;
                break;

            case 3:
                D1 = v2.y() * v4.x() - v2.x() * v4.y();
                D2 = v2.y() - v4.y();
                D3 = v2.x() - v4.x();
                D4 = 0;
                Ddet = D1 - v1.x()*D2 + v1.y()*D3 - v1.z()*D4;
        }
    }

    if (i == 4){
        switch (j)
        {
            case 1:
                D1 = v2.y() * v3.z() - v2.z() * v3.y();
                D2 = 0;
                D3 = v2.z() - v3.z();
                D4 = v2.y() - v3.y();
                Ddet = D1 - v1.x()*D2 + v1.y()*D3 - v1.z()*D4;
                break;

            case 2:
                D1 = v2.z() * v3.x() - v2.x() * v3.z();
                D2 = v2.z() - v3.z();
                D3 = 0;
                D4 = v3.x() - v2.x();
                Ddet = D1 - v1.x()*D2 + v1.y()*D3 - v1.z()*D4;
                break;

            case 3:
                D1 = v2.x() * v3.y() - v2.y() * v3.x();
                D2 = v3.y() - v2.y();
                D3 = v3.x() - v2.x();
                D4 = 0;
                Ddet = D1 - v1.x()*D2 + v1.y()*D3 - v1.z()*D4;
        }
    }

    return Ddet;
}
