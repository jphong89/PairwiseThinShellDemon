#include "stdafx.h"

#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#endif

#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

#include "MeshObject.h"

signature* BasicMesh::findSignature(int dirOpt,PolyhedralSurf::Vertex_const_handle vh,double rotation){
	float distance = geodesicDistance;

	Vector_3 d1 = vertex2d1_map[vh];
	Vector_3 d2 = vertex2d2_map[vh];
	Vector_3 normal = cross_product(d1,d2);

	Vector_3 dir;
	if (dirOpt == 0)
		dir = d1;
	else if (dirOpt == 1)
		dir = -d1;
	else if (dirOpt == 2)
		dir = d2;
	else if (dirOpt == 3)
		dir = -d2;
	else if (dirOpt == 4)
		dir = (d1 + d2)/2;
	else if (dirOpt == 5)
		dir = (-d1 + d2)/2;
	else if (dirOpt == 6)
		dir = (-d1 - d2)/2;
	else
		dir = (d1 - d2)/2;

	Vector_3 odir = dir;
	dir = rotate(dir,normal,rotation);

	int edgeNum = vh->degree();

	// find the cooresponding facet/halfedge for a certain walking direction
	PolyhedralSurf::Halfedge_around_vertex_const_circulator hfe = vh->vertex_begin();
	PolyhedralSurf::Vertex_const_handle firstVertex =  hfe->opposite()->vertex();
	PolyhedralSurf::Vertex_const_handle lastVertex,vertex1,vertex2;
	Vector_3 firstProjected;

	Vector_3 projected,lastProjected;
	bool found = false;
	float theta1,theta2,theta3;

	for (int j = 0; j < edgeNum; j++){
		if (hfe->is_border_edge()) 
			if (strcmp(tag.c_str(),"deformed") == 0) isBoundary[indexMap.find(vh)->second] = true;

		PolyhedralSurf::Vertex_const_handle nextV = hfe->opposite()->vertex();

		projected = projectVertex(d1,d2,Vector_3(vh->point(),nextV->point()));

		if (j == 0){
			firstProjected = projected;

			lastProjected = projected;
			lastVertex = nextV;
			hfe++;
			continue;
		}

		theta1 = acos(lastProjected*dir/sqrt(lastProjected.squared_length()*dir.squared_length()));
		theta2 = acos(projected*dir/sqrt(projected.squared_length()*dir.squared_length()));
		theta3 = acos(lastProjected*projected/sqrt(lastProjected.squared_length()*projected.squared_length()));

		if (abs((theta1 + theta2 - theta3))<0.001){
			vertex1 = lastVertex;
			vertex2 = nextV;
			found  =true;
			break;
		}

		lastProjected = projected;
		lastVertex = nextV;
		hfe++; 
	}

	if (!found){
		theta1 = acos(lastProjected*dir/sqrt(lastProjected.squared_length()*dir.squared_length()));
		theta2 = acos(firstProjected*dir/sqrt(firstProjected.squared_length()*dir.squared_length()));
		theta3 = acos(lastProjected*firstProjected/sqrt(lastProjected.squared_length()*firstProjected.squared_length()));

		if ((theta1 + theta2 - theta3)<0.001){
			vertex1 = lastVertex;
			vertex2 = firstVertex;
		}
		else{ // v1 is at an edge
// 			marchList[dirOpt] = NULL;
// 			signature* sig = new signature();
// 			return sig;
			std::map<Vertex_const_handle, int>::iterator iter;
			iter = indexMap.find(vh);

			return NULL;
		}
	}

	// now find the halfedge of vertex1,vertex2
	PolyhedralSurf::Halfedge_around_vertex_const_circulator tempHfe = vertex1->vertex_begin();
	PolyhedralSurf::Halfedge_const_handle nextHalfEdge;
	PolyhedralSurf::Facet_const_handle nextFacet;

	int v1Degree = vertex1->degree();
	int degreeCount = 1;
	while (degreeCount<=v1Degree){
		if (tempHfe->opposite()->vertex() == vertex2){
			nextHalfEdge = tempHfe->opposite();
			nextFacet = nextHalfEdge->facet();
			break;
		}
		tempHfe++;
		degreeCount++;
	}

	if (degreeCount > v1Degree){ // v1 is at an edge
// 		marchList[dirOpt] = NULL;
// 		signature* sig = new signature();
// 		return sig;
		std::map<Vertex_const_handle, int>::iterator iter;
		iter = indexMap.find(vh);

		return NULL;
	}

	//build a marching Path
	marchList[dirOpt] = new threeTuple(vh->point());

	float sx = (vertex1->point().x() * theta2 + vertex2->point().x() * theta1) / (theta1 + theta2);
	float sy = (vertex1->point().y() * theta2 + vertex2->point().y() * theta1) / (theta1 + theta2);
	float sz = (vertex1->point().z() * theta2 + vertex2->point().z() * theta1) / (theta1 + theta2);

	dir = Vector_3(vh->point(),Point_3(sx,sy,sz));

	PolyhedralSurf::Vertex_const_handle endPoint;

	marchList[dirOpt]->next = new threeTuple();
	marchList[dirOpt]->next->next = NULL;
	endPoint = marchOnSurface(nextHalfEdge,nextFacet,dir,theta1,theta2,distance,dirOpt,marchList[dirOpt]->next);

	if (endPoint == NULL) return NULL;

	float k1 = vertex2k1_map[endPoint];
	float k2 = vertex2k2_map[endPoint];

	signature* sig = new signature();
	sig->k1 = k1;
	sig->k2 = k2;
	sig->K = k1*k2;
	sig->H = (k1 + k2)/2;
	sig->C = sqrt((k1 * k1 + k2 * k2)/2);
	sig->S = atan((k1 + k2)/(k1 - k2));

	sig->normal = cross_product(vertex2d1_map[endPoint],vertex2d2_map[endPoint]);
	sig->normal = sig->normal / sqrt(sig->normal.squared_length());

	sig->d1 = vertex2d1_map[endPoint];
	sig->d2 = vertex2d2_map[endPoint];

	sig->point = Vector_3(endPoint->point().x(),endPoint->point().y(),endPoint->point().z());

	std::map<Vertex_const_handle, int>::iterator iter = indexMap.find(vh);

	return sig;
}

PolyhedralSurf::Vertex_const_handle BasicMesh::marchOnSurface(PolyhedralSurf::Halfedge_const_handle hfe,PolyhedralSurf::Facet_const_handle fc,Vector_3 lastD,float theta1,float theta2,float distance,int dirOpt,threeTuple* marchList){
	PolyhedralSurf::Vertex_const_handle v1,v2,v3;

	Vector_3 d1,d2;

	float coveredDistance = 0;
	float sigma = 2;
	PolyhedralSurf::Halfedge_const_handle h1,h2,h3;
	PolyhedralSurf::Vertex_const_handle tv1,tv2,tv3;
	PolyhedralSurf::Vertex_const_handle vertex1,vertex2;

	int step = 0;
	bool valid = false;

	while (coveredDistance < distance){
		v1 = hfe->opposite()->vertex();
		v2 = hfe->vertex();
		vertex1 = v1;
		vertex2 = v2;

		if (abs(theta1) < 0.0001) theta1 = 0.01;
		if (abs(theta2) < 0.0001) theta2 = 0.01;

		float sx = (v1->point().x() * theta2 + v2->point().x() * theta1) / (theta1 + theta2);
		float sy = (v1->point().y() * theta2 + v2->point().y() * theta1) / (theta1 + theta2);
		float sz = (v1->point().z() * theta2 + v2->point().z() * theta1) / (theta1 + theta2);

		marchList->x = sx;
		marchList->y = sy;
		marchList->z = sz;		

		if (hfe->is_border_edge()){
			if (step >= 0)
				valid = true;
			break;
		}

		h1 = fc->halfedge();
		tv1 = h1->vertex();
		h2 = h1->next();
		tv2 = h2->vertex();
		h3 = h2->next();
		tv3 = h3->vertex();

		if ((tv1 != v1)&&(tv1 != v2)) v3 = tv1;
		if ((tv2 != v1)&&(tv2 != v2)) v3 = tv2;
		if ((tv3 != v1)&&(tv3 != v2)) v3 = tv3;

		Point_3 s(sx,sy,sz); // starting point
		Vector_3 nv1(s,v1->point());
		Vector_3 nv2(s,v2->point());
		Vector_3 nv3(s,v3->point());
		Vector_3 projected = projectVertex(nv1,nv3,lastD);

		lastD = projected;

		theta1 = acos(projected*nv1/sqrt(projected.squared_length()*nv1.squared_length()));
		theta2 = acos(projected*nv3/sqrt(projected.squared_length()*nv3.squared_length()));
		float theta3 = acos(nv1*nv3/sqrt(nv1.squared_length()*nv3.squared_length()));


		if ((theta1 + theta2 - theta3)<0.0001){
			vertex1 = v1;
			vertex2 = v3;

			if (abs(theta1+theta2+theta3 - 3.141592654 * 2) < 0.0001)
				lastD = -lastD;

			theta1 = sin(theta1) * sqrt(nv1.squared_length());
			theta2 = sin(theta2) * sqrt(nv3.squared_length());
		}else{		
			theta1 = acos(projected*nv3/sqrt(projected.squared_length()*nv3.squared_length()));
			theta2 = acos(projected*nv2/sqrt(projected.squared_length()*nv2.squared_length()));
			theta3 = acos(nv2*nv3/sqrt(nv2.squared_length()*nv3.squared_length()));

			if (abs(theta1+theta2+theta3 - 3.141592654 * 2) < 0.0001)
				lastD = -lastD;

			vertex1 = v3;
			vertex2 = v2;
			theta1 = sin(theta1) * sqrt(nv3.squared_length());
			theta2 = sin(theta2) * sqrt(nv2.squared_length());
		}

		h1 = fc->halfedge();
		while (true){
			if((h1->vertex()==vertex1)&&(h1->opposite()->vertex()==vertex2)){
				hfe = h1->opposite();
				fc = hfe->facet();
				
				break;
			}

			h1 = h1->next();
		}

		sx = (vertex1->point().x() * theta2 + vertex2->point().x() * theta1) / (theta1 + theta2);
		sy = (vertex1->point().y() * theta2 + vertex2->point().y() * theta1) / (theta1 + theta2);
		sz = (vertex1->point().z() * theta2 + vertex2->point().z() * theta1) / (theta1 + theta2);
		Point_3 e(sx,sy,sz); // ending point
		coveredDistance = coveredDistance + sqrt(Vector_3(s,e).squared_length());

		step++;

		if (step > 100){
			valid = true;
			break;
		}

		if (coveredDistance <distance){
			marchList->next = new threeTuple();
			marchList->next->next = NULL;
			marchList = marchList->next;
		}else{
			valid = true;
			break;
		}
	}

	if (valid){
		if (theta2 > theta1)
			return vertex1;
		else
			return vertex2;
	}else
		return NULL;
}

Vertex_const_handle BasicMesh::findVertexHandle(threeTuple a){
	PolyhedralSurf::Vertex_const_iterator vb = P.vertices_begin();
	PolyhedralSurf::Vertex_const_iterator ve = P.vertices_end();

	float testX = a.x;
	float testY = a.y;
	float testZ = a.z;
	float minDis = 10000;

	PolyhedralSurf::Vertex_const_handle vh = vb;
	for (int i = 0; vb != ve; vb++,i++){
		if (abs(vb->point().x() - testX)+
			abs(vb->point().y() - testY)+
			abs(vb->point().z() - testZ) <minDis){
				vh = vb;
				minDis = abs(vb->point().x() - testX)+
					abs(vb->point().y() - testY)+
					abs(vb->point().z() - testZ);
		}

	}

	centre = vh;

	return vh;
}

/* construct geometric signature for one single vertex */
signature* BasicMesh::constructSignature(Vertex_const_handle vh,double rotation){
	signature* sig = new signature[dirNum + 1];

	/* go along 8 directions and collect features at those end points */
	
	/*
	for (int i = 0; i < dirNum; i++){
		signature* temp = findSignature(i,vh,rotation);

		if (temp == NULL){
			delete[] sig;
			deleteMarchList();
			return NULL;
		}

		sig[i+1] = *temp;
	}
	*/

	float k1 = vertex2k1_map[vh];
	float k2 = vertex2k2_map[vh];

	sig[0].k1 = k1;
	sig[0].k2 = k2;
	sig[0].K = k1*k2;
	sig[0].H = (k1 + k2)/2;
	sig[0].C = sqrt((k1 * k1 + k2 * k2)/2);
	sig[0].S = atan((k1 + k2)/(k1 - k2));
	
	sig[0].normal = cross_product(vertex2d1_map[vh],vertex2d2_map[vh]);
	sig[0].normal = sig[0].normal / sqrt(sig[0].normal.squared_length());
	sig[0].d1 = vertex2d1_map[vh];
	sig[0].d2 = vertex2d2_map[vh];

	sig[0].point = Vector_3(vh->point().x(),vh->point().y(),vh->point().z());

	/*
	for (int i = 1; i <= DIR_NUM; i++){
		sig[i].deltaS = (sig[i].S - sig[0].S)/geodesicDistance;
		sig[i].deltaN = abs(acos(sig[i].normal * sig[0].normal / sqrt(sig[i].normal.squared_length()*sig[0].normal.squared_length())))/geodesicDistance;
	}

	for (int i = 1; i <= DIR_NUM; i++){
		sig[i].rotation = constructQuaternion(sig,sig+i);
	}

	sig[0].deltaN12 = abs(acos(sig[1].normal * sig[2].normal / sqrt(sig[1].normal.squared_length()*sig[2].normal.squared_length())))/geodesicDistance;
	sig[0].deltaN34 = abs(acos(sig[3].normal * sig[4].normal / sqrt(sig[3].normal.squared_length()*sig[4].normal.squared_length())))/geodesicDistance;

	sig[0].dis12 = sqrt((sig[1].point-sig[2].point).squared_length())/geodesicDistance;
	sig[0].dis34 = sqrt((sig[3].point-sig[4].point).squared_length())/geodesicDistance;
	sig[0].dis56 = sqrt((sig[5].point-sig[6].point).squared_length())/geodesicDistance;
	sig[0].dis78 = sqrt((sig[7].point-sig[8].point).squared_length())/geodesicDistance;
	sig[0].dis58 = sqrt((sig[5].point-sig[8].point).squared_length())/geodesicDistance;
	sig[0].dis67 = sqrt((sig[7].point-sig[6].point).squared_length())/geodesicDistance;

	std::map<Vertex_const_handle, int>::iterator iter = indexMap.find(vh);

	deleteMarchList();
	*/
	return sig;
}

/* compute geometric signatures for all vertices */
void BasicMesh::findSignatureAll(){
	PolyhedralSurf::Vertex_const_iterator vb = P.vertices_begin();
	PolyhedralSurf::Vertex_const_iterator ve = P.vertices_end();
	signatureMap.clear();

	matched = vb;
	corresponding = vb;

	PolyhedralSurf::Vertex_const_handle vh;
	for (int i = 0; vb != ve; vb++,i++){
		vh = vb;

		signature* sig;

		sig = constructSignature(vh,0);

		if (sig != NULL)
			signatureMap.insert(pair<Vertex_const_handle, signature*>(vh,sig));
	}

	deleteMarchList();
}

void BasicMesh::outputFeatures(){
	ofstream fout;
	fout.open(if_name+"features.txt");

	fout<<vertexNum<<endl;

	PolyhedralSurf::Vertex_const_iterator vb = P.vertices_begin();
	PolyhedralSurf::Vertex_const_iterator ve = P.vertices_end();
	PolyhedralSurf::Vertex_const_handle vh;
	std::map<Vertex_const_handle, signature*>::iterator iterSig;
	signature* sig;

	for (int i = 0; vb != ve; vb++,i++){
// 		vh = vb;
// 
// 		fout<<i+1<<' ';
// 
// 		Vector_3 d1 = vertex2d1_map[vh];
// 		Vector_3 d2 = vertex2d2_map[vh];
// 		Vector_3 normal = cross_product(d1,d2);
// 
// 		float k1 = vertex2k1_map[vh];
// 		float k2 = vertex2k2_map[vh];
// 		float K = k1*k2;
// 		float H = (k1 + k2)/2;
// 
// 		fout<<normal<<' '<<K<<' '<<H<<endl;

		fout<<matchWeight[i]<<endl;
	}

	fout.close();
}

void BasicMesh::deleteMarchList(){
	for (int i = 0; i < 8; i++){
		threeTuple* head = marchList[i];

		threeTuple* temp = head;
		while (head != NULL){
			if (head->next != NULL)
				temp = head->next;
			else 
				temp = NULL;

			delete head;
			
			head = temp;
		}

		marchList[i] = NULL;
	}
}

float compareSignature(signature* sig1,signature* sig2){
	/*
	Eigen::VectorXd x1(FEATDIM);
	Eigen::VectorXd x2(FEATDIM);

	x1.setZero();
	x2.setZero();

	int textureRatio = 50;

	x1(0) = sig1[0].C;
	x1(1) = sig1[0].deltaN12;
	x1(2) = sig1[0].deltaN34;

	for (int i = 1; i <= DIR_NUM; i++){
		x1((i-1)*10+5+1) = sig1[i].C;
		x1((i-1)*10+5+2) = sig1[i].deltaS;
		x1((i-1)*10+5+3) = sig1[i].deltaN;
		x1((i-1)*10+5+4) = sig1[i].rotation.w();
		x1((i-1)*10+5+5) = sig1[i].rotation.x();
		x1((i-1)*10+5+6) = sig1[i].rotation.y();
		x1((i-1)*10+5+7) = sig1[i].rotation.z();
	}

	x2(0) = sig2[0].C;
	x2(1) = sig2[0].deltaN12;
	x2(2) = sig2[0].deltaN34;

	for (int i = 1; i <= DIR_NUM; i++){
		x2((i-1)*10+5+1) = sig2[i].C;
		x2((i-1)*10+5+2) = sig2[i].deltaS;
		x2((i-1)*10+5+3) = sig2[i].deltaN;
		x2((i-1)*10+5+4) = sig2[i].rotation.w();
		x2((i-1)*10+5+5) = sig2[i].rotation.x();
		x2((i-1)*10+5+6) = sig2[i].rotation.y();
		x2((i-1)*10+5+7) = sig2[i].rotation.z();
	}

	float ans = 0;
	for (int i = 0; i < FEATDIM; i++)
		ans += (x1(i)-x2(i))*(x1(i)-x2(i));
	//cout<<"dif: "<<ans<<endl;
	*/
	
	float dis1 = (sig1[0].normal - sig2[0].normal).squared_length();
	float dis2 = (sig1[0].C - sig2[0].C) * (sig1[0].C - sig2[0].C);
	float dis3 = (sig1[0].S - sig2[0].S) * (sig1[0].S - sig2[0].S);


	return dis1 + dis2*5 + dis3*5;
}

void BasicMesh::findCorrespondenceBothWay(BasicMesh* secondMesh,double featWeight){
	PolyhedralSurf::Vertex_const_iterator vb = secondMesh->P.vertices_begin();
	PolyhedralSurf::Vertex_const_iterator ve = secondMesh->P.vertices_end();
	PolyhedralSurf::Vertex_const_handle vh1,vh2;
	float dis,dif,normalDif;

	float* featDis = new float[affinityM * affinityN];

	for (int i = 0; i < affinityM; i++)
		for (int j = 0; j < affinityN; j++){
			VAL(affinity,i,j,affinityN) = -1;
			VAL(featDis,i,j,affinityN) = -1;
		}

	vb = P.vertices_begin();
	ve = P.vertices_end();

	std::map<Vertex_const_handle, signature*>::iterator iterSig;

	for (int i = 0; vb != ve; vb++,i++){
		//cout<<"Working on vertex: "<<i<<endl;

		vh1 = vb;

		iterSig = signatureMap.find(vh1);
		if (iterSig == signatureMap.end()) continue;
		signature* sig1 = iterSig->second;
		if (sig1 == NULL) continue;
		
		for (int j = 0; j < secondMesh->vertexNum; j++){
			vh2 = secondMesh->vertexIndex[j];

			iterSig = secondMesh->signatureMap.find(vh2);
			if (iterSig == secondMesh->signatureMap.end()) continue;
			signature* sig2 = iterSig->second;
			if (sig2 == NULL) continue;

			dif = compareSignature(sig1, sig2);

 			dis = computeEuclideanDis(vh1->point(),vh2->point());
			if (dis > DISTHRESHOLD) continue;

			float theta = acos(sig1->normal*sig2->normal/sqrt(sig1->normal.squared_length()*sig2->normal.squared_length()));
			if (theta  > PI/2 * 0.8) continue;

			dif = dif * featWeight + dis;

			VAL(affinity,i,j,affinityN) = dif;
			VAL(featDis,i,j,affinityN) = dif;
		}
	}

	//optional affinity
	float minA = 1000;
	float maxA = 0;

	float temp;
	for (int i = 0; i < affinityM; i++){
		for (int j = 0; j < affinityN; j++){
			temp = VAL(affinity,i,j,affinityN);
			if (temp >= 0){
				if (temp > maxA) maxA = temp;
				if (temp < minA) minA = temp;
			}
		}
	}

	for (int i = 0; i < affinityM; i++){
		for (int j = 0; j < affinityN; j++){
			temp = VAL(affinity,i,j,affinityN);
			if (temp >= 0)
				VAL(affinity,i,j,affinityN) = 1-(temp - minA) / (maxA - minA);
			else
				VAL(affinity,i,j,affinityN) = 0;
			
		}
	}

	delete [] featDis;
}

// void BasicMesh::findMatch(BasicMesh* secondMesh){
// 	//use all correspondences
// 	for (int i = 0; i < affinityM; i++){
// 		if (matchWeight[i] > 1) continue;
// 		
// 		float maxA = -1;
// 
// 		double sum = 0;
// 		for (int j = 0; j < affinityN; j++){
// 			if (VAL(affinity,i,j,affinityN) > maxA){
// 				maxA = VAL(affinity,i,j,affinityN);
// 				bestMatch[i] = j;
// 				matchWeight[i] = VAL(affinity,i,j,affinityN);
// 			}
// 			sum += VAL(affinity,i,j,affinityN);
// 		}
// 		
// /*		matchWeight[i] = sum;*/
// 		if (sum > 2.5)
// 			matchWeight[i] = 1 / (1 + exp(-0.5*(matchWeight[i]*100-73)));
// 		else
// 			matchWeight[i] = 0;
// 		
// 		/*syn 5: 0.3 100 35 */
// 		/* 5708 -frame 240: 27,0.5,100,35 */
// 
// 		/* 5708-frame 10,450,760: 35,100,45
// 		   5708-frame 1400: 10,100,40*/
// 	}
// 	
// 	for (int i = 0; i < affinityM; i++){
// 		if (!attractor[i]) continue;
// 
// 		matchWeight[i] = 0;
// 	}
// }

int compareNode (const void * a, const void * b)
{
	if ( ((qnode*)a)->value <  ((qnode*)b)->value ) return -1;
	if ( ((qnode*)a)->value ==  ((qnode*)b)->value ) return 0;
	if ( ((qnode*)a)->value >  ((qnode*)b)->value ) return 1;
}

void BasicMesh::findMatch2(BasicMesh* secondMesh){
	//use all correspondences
	qnode* queue = new qnode[affinityM*affinityN];
	for (int i = 0; i< affinityM; i++)
		for (int j = 0; j < affinityN; j++){
			queue[j+i*affinityN].value = VAL(affinity,i,j,affinityN);
			queue[j+i*affinityN].i = i;
			queue[j+i*affinityN].j = j;
		}
	
	qsort(queue,affinityM*affinityN,sizeof(qnode),compareNode);
	
	bool* ibool = new bool[affinityM];
	bool* jbool = new bool[affinityN];
	memset(ibool,false,sizeof(bool)*affinityM);
	memset(jbool,false,sizeof(bool)*affinityN);

	double maxWeight = 1;
	/* first add these landmarks */
	for (int k = 0; k < landmarkNum; k++){
		ibool[landmark[k]] = true;
		jbool[secondMesh->landmark[k]] = true;

		matchWeight[landmark[k]] = maxWeight;
		Vertex_const_handle vh = secondMesh->vertexIndex[secondMesh->landmark[k]];
		bestMatch[landmark[k]] = Vector_3(vh->point().x(),vh->point().y(),vh->point().z()); 
	}

	/* then use geometric feature to find some other correspondences */
	for (int k = affinityM * affinityN-1; k >= 0; k--)
		if (!ibool[queue[k].i] && !jbool[queue[k].j]){
			ibool[queue[k].i] = true;
			jbool[queue[k].j] = true;
			if (queue[k].value < 0.9) break;

			matchWeight[queue[k].i] = (queue[k].value * maxWeight - 0.9) / 0.1;
			//matchWeight[queue[k].i] = 0.8;
			PolyhedralSurf::Vertex_const_handle vh =  secondMesh->vertexIndex[queue[k].j];
			PolyhedralSurf::Vertex_const_handle source = vertexIndex[queue[k].i];
			PolyhedralSurf::Halfedge_around_vertex_const_circulator tempHfe = vh->vertex_begin();
			Vector_3 minTarget;
			double minDis = 1000;

			for (int i = 0; i< vh->degree(); i++){
				if (tempHfe->is_border()) tempHfe++;
				
				PolyhedralSurf::Halfedge_around_facet_const_circulator shf = tempHfe->facet()->facet_begin();

				PolyhedralSurf::Vertex_const_handle v1 = shf->vertex();
				shf++;
				PolyhedralSurf::Vertex_const_handle v2 = shf->vertex();
				shf++;
				PolyhedralSurf::Vertex_const_handle v3 = shf->vertex();

				Vector_3 tri[3] = {Vector_3(v1->point().x(),v1->point().y(),v1->point().z()),
								   Vector_3(v2->point().x(),v2->point().y(),v2->point().z()),
								   Vector_3(v3->point().x(),v3->point().y(),v3->point().z())};

				Vector_3 sourcePos = Vector_3(source->point().x(),source->point().y(),source->point().z());

				Vector_3 target = closesPointOnTriangle(tri,sourcePos);
				Vector_3 dis = target - sourcePos;
				if (dis.squared_length() < minDis){
					minTarget = target;
					minDis = dis.squared_length();
				}

				tempHfe++;
			}

			bestMatch[queue[k].i] = minTarget;
			
		}
	
	delete[] queue;
	delete[] ibool;
	delete[] jbool;
}

// void BasicMesh::findClosestPoint(BasicMesh* secondMesh){
// 	PolyhedralSurf::Vertex_const_iterator vb = P.vertices_begin();
// 	PolyhedralSurf::Vertex_const_iterator ve = P.vertices_end();
// 	PolyhedralSurf::Vertex_const_handle vh1,vh2;
// 
// 	std::map<Vertex_const_handle, signature*>::iterator iterSig;
// 
// 	for (int i = 0; vb != ve; vb++,i++){
// 		//cout<<"Working on vertex: "<<i<<endl;
// 		double minDis = 100000;
// 		int minJ = 0;
// 
// 		vh1 = vb;
// 
// 		iterSig = signatureMap.find(vh1);
// 		if (iterSig == signatureMap.end()) continue;
// 
// 		signature* sig1 = iterSig->second;
// 		if (sig1 == NULL) continue;
// 
// 		for (int j = 0; j < secondMesh->vertexNum; j++){
// 			vh2 = secondMesh->vertexIndex[j];
// 
// 			iterSig = secondMesh->signatureMap.find(vh2);
// 
// 			if (iterSig == secondMesh->signatureMap.end()) continue;
// 
// 			signature* sig2 = iterSig->second;
// 			if (sig2 == NULL) continue;
// 
// 			float theta = acos(sig1->normal*sig2->normal/sqrt(sig1->normal.squared_length()*sig2->normal.squared_length()));
// 			if (theta > PI/2) continue;
// 
// 			double dis = computeEuclideanDis(vertexIndex[i]->point(),secondMesh->vertexIndex[j]->point());
// 
// 			if (dis < minDis){
// 				minDis = dis;
// 				minJ = j;
// 			}
// 		}
// 
// 		bestMatch[i] = minJ;
// 		
// 		////////////////
// 		int minI;
// 		minDis = 100000;
// 
// 		vh1 = secondMesh->vertexIndex[minJ];
// 		iterSig = secondMesh->signatureMap.find(vh1);
// 		sig1 = iterSig->second;
// 		if (sig1 == NULL) continue;
// 
// 		for (int j = 0; j < vertexNum; j++){
// 			vh2 = vertexIndex[j];
// 
// 			iterSig = signatureMap.find(vh2);
// 
// 			if (iterSig == signatureMap.end()) continue;
// 
// 			signature* sig2 = iterSig->second;
// 			if (sig2 == NULL) continue;
// 
// 			float theta = acos(sig1->normal*sig2->normal/sqrt(sig1->normal.squared_length()*sig2->normal.squared_length()));
// 			if (theta > PI/2) continue;
// 
// 			double dis = computeEuclideanDis(vertexIndex[j]->point(),secondMesh->vertexIndex[minJ]->point());
// 
// 			if (dis < minDis){
// 				minDis = dis;
// 				minI = j;
// 			}
// 		}
// 
// 		double dis = computeEuclideanDis(vertexIndex[minI]->point(),vertexIndex[i]->point());
// 
// // 		if (minDis < 3)
// // 			matchWeight[i] = 0.4;
// // 		else 
// // 			matchWeight[i] = 0;
// 		matchWeight[i] = exp(- 0.5 * dis*dis);
// 	}
// 
// 
// }