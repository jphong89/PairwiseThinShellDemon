#include "stdafx.h"
#include "GlutManager.h"
#include <iostream>
using namespace std;

int beginGlut(int argc, char* argv[]){
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH ); 
	glutInitWindowSize(400, 400); 

	glutCreateWindow("OpenGL"); 
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	init();
	glutKeyboardFunc(ProcessKeyboard);
	glutSpecialFunc(ProcessSpecialKeyboead);
	glutMouseFunc(Mouse);

	cout<<"Begin Rendering..."<<endl;
	glutMainLoop();
	return 0;
}

void ProcessKeyboard(unsigned char key,int x,int y)
{
	if (key == '=')
		camera_Scale = camera_Scale + 0.1;
	else if (key == '-')
		camera_Scale = camera_Scale - 0.1;
	else if (key == 'w')
		rotate_X = rotate_X - 2;
	else if (key == 's')
		rotate_X = rotate_X + 2;
	else if (key == 'a')
		rotate_Y = rotate_Y - 2;
	else if (key == 'd')
		rotate_Y = rotate_Y + 2;
	else if (key == 'q')
		rotate_Z = rotate_Z - 2;
	else if (key == 'e')
		rotate_Z = rotate_Z + 2;
	else if (key == 'b'){
		/* initialize optimization */
		char filename[50];

		lbfgsfloatval_t* tempU = lbfgs_malloc(affinityM * 3);
		int ret;

		for (int i = 0; i < 1; i++){
			cout<<"Velocity Iteration: "<<i<<endl;

			if (i > 0){
				cout<<"Reset surface..."<<endl;
				RCSurface->ComputeMeshProperty(filename);

				cout<<"Recompute geometric features..."<<endl;
				RCSurface->findSignatureAll();
				

				memset(u,0,sizeof(lbfgsfloatval_t) * affinityM * 3);
				
				//cout<<"Recompute affinity map..."<<endl;
				memset(affinity,0,affinityM * affinityN * sizeof(double));
				RCSurface->findCorrespondenceBothWay(CTSurface);
			
				/* surface color */
				paintColor(RCSurface,threeTuple(0.8,0.8,0.8),threeTuple(0.8,0.8,0.8),threeTuple(0.8,0.8,0.8));
				paintColor(CTSurface,threeTuple(0.2,0.2,0.9),threeTuple(0.2,0.2,0.9),threeTuple(0.2,0.2,0.9));
			}

			cout<<"Begin optimization..."<<endl;

			lbfgsfloatval_t lastFx = 100000001;
			lbfgsfloatval_t fx = 100000000;
			

			int alterCount = 0;
			memset(tempU,0,sizeof(lbfgsfloatval_t)*affinityM*3);
			
			while (lastFx > fx){
				memcpy(u,tempU,sizeof(lbfgsfloatval_t)*affinityM*3);
				lastFx = fx;
				
				if (alterCount == 1) break;
				alterCount++;
				cout<<"--------------------------zig iteration: "<<alterCount<<endl;

				lbfgs_parameter_t param;
				bendweight = 0;
				ret = 0;

				lbfgs_parameter_init(&param);
				//param.max_iterations = 5;
				ret = lbfgs(affinityM * 3, tempU, &fx, evaluate, progress, NULL, &param);

				////
				continue;
				cout<<"-------------------------zag iteration: "<<alterCount<<endl;
				lbfgs_parameter_t param1;
				bendweight = BENDWEIGHT;
				ret = 0;

				lbfgs_parameter_init(&param1);
				ret = lbfgs(affinityM * 3, tempU, &fx, evaluate, progress, NULL, &param1);

				
			}

			
			filename[0] = '\0';
			strcat(filename,"temp");
			char iterNum[10];
			itoa(i,iterNum,10);
			strcat(filename,iterNum);
			strcat(filename,".off");
			writeOFF(RCSurface,filename);

			cout<<"L-BFGS optimization terminated with status code: "<<ret<<"fx: "<<lastFx<<endl;
	
		}
	}
	else if (key == '1')
		switch1 = !switch1;
	else if (key == '2')
		switch2 = !switch2;

	display();
}

void ProcessSpecialKeyboead(int key, int x, int y) 
{
	if (key == GLUT_KEY_LEFT)
		camera_Right = camera_Right - 0.1;
	else if (key == GLUT_KEY_RIGHT)
		camera_Right = camera_Right + 0.1;
	else if (key == GLUT_KEY_UP)
		camera_Up = camera_Up + 0.1;
	else if (key == GLUT_KEY_DOWN)
		camera_Up = camera_Up - 0.1;

	display();
}

void Mouse(int button, int state, int x, int y) //处理鼠标点击
{
	GLint viewport[4];   
	GLdouble modelview[16];   
	GLdouble projection[16];   
	GLfloat winX, winY, winZ;   
	GLdouble posX, posY, posZ;   

	glPushMatrix(); 
	glLoadIdentity();    
  
	glTranslatef(camera_Right,camera_Up,-6.0f);

	glRotatef(rotate_X,1.0f,0.0f,0.0f);	
	glRotatef(rotate_Y,0.0f,1.0f,0.0f);	
	glRotatef(rotate_Z,0.0f,0.0f,1.0f);	

	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);   
	glGetDoublev(GL_PROJECTION_MATRIX, projection);   
	glPopMatrix();   

	winX = x;   
	winY = viewport[3] - (float)y;   
	glReadPixels((int)winX, (int)winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);   
	gluUnProject(winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ); 

	threeTuple a;
	a.x = posX / camera_Scale * max_all + (maxx - minx) / 2 + minx;
	a.y = posY / camera_Scale * max_all + (maxy - miny) / 2 + miny;
	a.z = posZ / camera_Scale * max_all + (maxz - minz) / 2 + minz;
	
	cout<<"seleted coordinates:    "<<a.x<<' '<<a.y<<' '<<a.z<<endl;

	PolyhedralSurf::Vertex_const_handle vh = RCSurface->findVertexHandle(a,CTSurface);

	int center = RCSurface->indexMap.find(vh)->second;
	//cout<<center<<endl;
	RCSurface->deleteMarchList();
	RCSurface->vertexBool.clear();
	for (int j = 0; j < RCSurface->dirNum; j++)
		RCSurface->findSignature(j,RCSurface->centre,0);


	CTSurface->meshColor.clear();
	PolyhedralSurf::Vertex_const_iterator vb = CTSurface->P.vertices_begin();
	PolyhedralSurf::Vertex_const_iterator ve = CTSurface->P.vertices_end();

	for (int i = 0; vb != ve; vb++,i++){
		vh = vb;
		CTSurface->meshColor.insert(pair<Vertex_const_handle,float>(vh,VAL(affinity,center,i,affinityN)));
	}
	for (int i = 0; i < CTSurface->faceNum; i++){
		float color = CTSurface->meshColor.find(CTSurface->faceList[i].v1)->second;
		CTSurface->faceList[i].color1 = getJetColor(color);

		color = CTSurface->meshColor.find(CTSurface->faceList[i].v2)->second;
		CTSurface->faceList[i].color2 = getJetColor(color);

		color = CTSurface->meshColor.find(CTSurface->faceList[i].v3)->second;
		CTSurface->faceList[i].color3 = getJetColor(color);
	}

	
	display();
}

void init(){
	glEnable(GL_DEPTH_TEST);			//开启深度测试
	
	glEnable(GL_LIGHT0);								// Enable Light One
	glEnable(GL_LIGHTING); 

	GLfloat light_ambient [] = { 0.3,0.3,0.3, 1.0 };
	GLfloat light_diffuse [] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_specular [] = { 0.2, 0.2, 0.2, 1.0 };
	GLfloat light_position[] = {0.0f, 0.0f, 0.0f, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT,   light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,   light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR,  light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION,  light_position);
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0f);

	glEnable(GL_COLOR_MATERIAL);

	GLfloat ambient [] = { 0.029412, 0.223529, 0.327451, 1.0 };
	GLfloat diffuse [] = { 0.180392, 0.268627, 0.813725, 1.0 };
	GLfloat specular[] = { 0.492157, 0.441176, 0.807843, 1.0 };
	GLfloat shininess[] = { 50.8974f };
	GLfloat emission[]  = {0.3f,0.3f,0.3f,1.0f};

	glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT,   ambient);
	glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE,   diffuse);
	glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR,  specular);
	glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, shininess);
	glMaterialfv( GL_FRONT_AND_BACK, GL_EMISSION,  emission);

	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	cout<<"Initializing..."<<endl;
	paintColor(RCSurface,threeTuple(0.8,0.8,0.8),threeTuple(0.8,0.8,0.8),threeTuple(0.8,0.8,0.8));
}

void display(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear The Screen And The Depth Buffer
	glLoadIdentity();
	//m_Camera.DisplayScene();
	glTranslatef(camera_Right,camera_Up,-6.0f);


	glRotatef(rotate_X,1.0f,0.0f,0.0f);	
	glRotatef(rotate_Y,0.0f,1.0f,0.0f);	
	glRotatef(rotate_Z,0.0f,0.0f,1.0f);	

	if (switch1) drawSurface(RCSurface,1);
	else drawSurface(RCSurface,0);

	if (switch2) drawSurface(CTSurface,0);

	drawExtra();

	glFlush();
	glutSwapBuffers();
}

void reshape(int w,int h)
{
	glViewport( 0, 0, w, h );
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	gluPerspective(75.0f, (float)w/h, 1.0f, 1000.0f);
	glMatrixMode( GL_MODELVIEW );
}

void drawSurface(BasicMesh* mesh,int deform){
	threeTuple color;

	Vector_3 pd1,pd2;

	float colorRatio = 1;

	for (int i = 0; i < mesh->linkNum; i++){
		Point_3 p1 = mesh->linkList[i].v1->point();

		threeTuple a;
		a.x = (p1.x()-  (maxx - minx) / 2 - minx) / max_all * camera_Scale;
		a.y = (p1.y() - (maxy - miny) / 2 - miny) / max_all * camera_Scale;
		a.z = (p1.z() - (maxz - minz) / 2 - minz) / max_all * camera_Scale;
	}

	for (int i = 0; i < mesh->faceNum; i++){
		facet f = mesh->faceList[i];

		f.p1.x += deform*u[f.index1*3];
		f.p1.y += deform*u[f.index1*3+1];
		f.p1.z += deform*u[f.index1*3+2];

		f.p2.x += deform*u[f.index2*3];
		f.p2.y += deform*u[f.index2*3+1];
		f.p2.z += deform*u[f.index2*3+2];

		f.p3.x += deform*u[f.index3*3];
		f.p3.y += deform*u[f.index3*3+1];
		f.p3.z += deform*u[f.index3*3+2];

		threeTuple a;
		a.x = (f.p1.x-  (maxx - minx) / 2 - minx) / max_all * camera_Scale;
		a.y = (f.p1.y - (maxy - miny) / 2 - miny) / max_all * camera_Scale;
		a.z = (f.p1.z - (maxz - minz) / 2 - minz) / max_all * camera_Scale;

		threeTuple b;
		b.x = (f.p2.x - (maxx - minx) / 2 - minx) / max_all * camera_Scale;
		b.y = (f.p2.y - (maxy - miny) / 2 - miny) / max_all * camera_Scale;
		b.z = (f.p2.z - (maxz - minz) / 2 - minz) / max_all * camera_Scale;

		threeTuple c;
		c.x = (f.p3.x - (maxx - minx) / 2 - minx) / max_all * camera_Scale;
		c.y = (f.p3.y - (maxy - miny) / 2 - miny) / max_all * camera_Scale;
		c.z = (f.p3.z - (maxz - minz) / 2 - minz) / max_all * camera_Scale;

		glLineWidth(1.0);
		glBegin(GL_TRIANGLES);
		glNormal3f(f.normal.x,f.normal.y,f.normal.z);
		
		glColor3f(f.color1.x,f.color1.y,f.color1.z);
		glVertex3f(a.x,a.y,a.z);
		
		glColor3f(f.color2.x,f.color2.y,f.color2.z);
		glVertex3f(b.x,b.y,b.z);

		glColor3f(f.color3.x,f.color3.y,f.color3.z);
		glVertex3f(c.x,c.y,c.z);
		glEnd();

		glLineWidth(3.0);

	}

}

void drawExtra(){
	glLineWidth(5.0);

	threeTuple color;
	for (int i = 0; i < RCSurface->dirNum; i++){

		color = getJetColor(0.125*(i+1)); 
		glColor3f(color.x,color.y,color.z);

		threeTuple* marchList = RCSurface->marchList[i];

		glBegin(GL_LINE_STRIP);
		while (marchList != NULL){

			threeTuple a;
			a.x = (marchList->x - (maxx - minx) / 2 - minx) / max_all * camera_Scale;
			a.y = (marchList->y - (maxy - miny) / 2 - miny) / max_all * camera_Scale;
			a.z = (marchList->z - (maxz - minz) / 2 - minz) / max_all * camera_Scale;
			glVertex3f(a.x,a.y,a.z);

			marchList = marchList->next;
		}
		glEnd();
	}
}

void drawSphere(GLfloat xx, GLfloat yy, GLfloat zz, GLfloat radius, GLfloat M, GLfloat N)
{
	float step_z = PI/M;
	float step_xy = 2*PI/N;
	float x[4],y[4],z[4];
	float angle_z = 0.0;
	float angle_xy = 0.0;
	int i=0, j=0;
	glBegin(GL_QUADS);
	for(i=0; i<M; i++)
	{
		angle_z = i * step_z;

		for(j=0; j<N; j++)
		{
			angle_xy = j * step_xy;
			x[0] = radius * sin(angle_z) * cos(angle_xy);
			y[0] = radius * sin(angle_z) * sin(angle_xy);
			z[0] = radius * cos(angle_z);
			x[1] = radius * sin(angle_z + step_z) * cos(angle_xy);
			y[1] = radius * sin(angle_z + step_z) * sin(angle_xy);
			z[1] = radius * cos(angle_z + step_z);
			x[2] = radius*sin(angle_z + step_z)*cos(angle_xy + step_xy);
			y[2] = radius*sin(angle_z + step_z)*sin(angle_xy + step_xy);
			z[2] = radius*cos(angle_z + step_z);
			x[3] = radius * sin(angle_z) * cos(angle_xy + step_xy);
			y[3] = radius * sin(angle_z) * sin(angle_xy + step_xy);
			z[3] = radius * cos(angle_z);
			for(int k=0; k<4; k++)
			{
				glVertex3f(xx+x[k], yy+y[k],zz+z[k]);
			}
		}
	}
	glEnd();
}

void paintColor(BasicMesh* mesh, threeTuple color1,threeTuple color2,threeTuple color3){
	for (int i = 0; i < mesh->faceNum; i++){
		mesh->faceList[i].color1 = color1;
		mesh->faceList[i].color2 = color2;
		mesh->faceList[i].color3 = color3;
	}
}


static int progress(
	void *instance,
	const lbfgsfloatval_t *u,
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

	cout<<"Iteration: "<<k<<endl;
	cout<<"fx = "<<fx<<endl;
	cout<<"xnorm = "<<xnorm<<", gnorm = "<<gnorm<<", step = "<<step<<endl;
	display();
	Sleep(1000);

	return 0;
}

static lbfgsfloatval_t evaluate(
	void *instance,
	const lbfgsfloatval_t *u,
	lbfgsfloatval_t *g,
	const int n,
	const lbfgsfloatval_t step
	)
{
	//cout<<"begin optimization"<<endl;
	lbfgsfloatval_t fx = 0.0;

	PolyhedralSurf::Vertex_const_iterator vb = RCSurface->P.vertices_begin();
	PolyhedralSurf::Vertex_const_iterator ve = RCSurface->P.vertices_end();
	PolyhedralSurf::Vertex_const_handle vh;

	for (int i = 0; vb != ve; vb++,i++){
		vh = vb;
	
		Point_3 deformed = vh->point() + Vector_3(u[i*3],u[i*3+1],u[i*3+2]);
		Point_3 cor = CTSurface->vertexIndex[bestMatch[i]]->point();

		Vector_3 dis(deformed,cor);
		double weight = VAL(affinity,i,bestMatch[i],affinityN);

		g[i*3] = -2 * dis.x() * weight;
		g[i*3+1] = -2 * dis.y() * weight;
		g[i*3+2] = -2 * dis.z() * weight;
		fx = fx + weight * dis.squared_length();

// 		if (sqrt(dis.squared_length()) <= HUBERSIGMA){
// 			fx = fx + weight * dis.squared_length();
// 			g[i*3] = -2 * dis.x() * weight;
// 			g[i*3+1] = -2 * dis.y() * weight;
// 			g[i*3+2] = -2 * dis.z() * weight;
// 		}else{
// 			fx = fx + weight * 2 * HUBERSIGMA * (sqrt(dis.squared_length()) - HUBERSIGMA / 2);
// 			
// 			g[i*3] = - 2 * HUBERSIGMA * weight * 0.5 * 1 / sqrt(dis.squared_length()) * 2 * dis.x();
// 			g[i*3+1] = - 2 * HUBERSIGMA * weight * 0.5 * 1 / sqrt(dis.squared_length()) * 2 * dis.y();
// 			g[i*3+2] = - 2 * HUBERSIGMA * weight * 0.5 * 1 / sqrt(dis.squared_length()) * 2 * dis.z();
// 		}
		
	}

	/* stretching */

	double stretching = 0;
	for (int i=0; i < RCSurface->edgeNum/2;i++){


		edge he = RCSurface->edgeList[i];

		int idx1 = he.index1;
		int idx2 = he.index2;

		Vector_3 deformVec(u[idx2*3] - u[idx1*3],u[idx2*3+1] - u[idx1*3+1],u[idx2*3+2] - u[idx1*3+2]);

		g[idx1*3] += -2 * REGWEIGHT * deformVec.x();
		g[idx1*3+1] += -2 * REGWEIGHT * deformVec.y();
		g[idx1*3+2] += -2 * REGWEIGHT * deformVec.z();
		g[idx2*3] += 2 * REGWEIGHT * deformVec.x();
		g[idx2*3+1] += 2 * REGWEIGHT * deformVec.y();
		g[idx2*3+2] += 2 * REGWEIGHT * deformVec.z();

		stretching += deformVec.squared_length();
	}

	/* bending */
	double bending = 0;

	for (int i=0; i < RCSurface->edgeNum/2;i++){
		edge he = RCSurface->edgeList[i];
		
		Point_3 newV1 = he.v1->point() + Vector_3(u[he.index1*3],u[he.index1*3+1],u[he.index1*3+2]);
		Point_3 newV2 = he.v2->point() + Vector_3(u[he.index2*3],u[he.index2*3+1],u[he.index2*3+2]);
		Point_3 newVl = he.vl->point() + Vector_3(u[he.indexl*3],u[he.indexl*3+1],u[he.indexl*3+2]);
		Point_3 newVr = he.vr->point() + Vector_3(u[he.indexr*3],u[he.indexr*3+1],u[he.indexr*3+2]);

		Vector_3 v1(newV1,newVr);
		Vector_3 v2(newV1,newV2);
		Vector_3 v3(newV1,newVl);

		Vector_3 n1 = cross_product(v1,v2);
		Vector_3 n2 = cross_product(v2,v3);

		double sm1 = n1.squared_length();
		double sq1 = sqrt(sm1);
		double sm2 = n2.squared_length();
		double sq2 = sqrt(sm2);

		double p = n1*n2;
		double gg = sq1*sq2;

		double kos = p/gg;
		
		if (kos > 1) kos = 1;
		if (kos < -1) kos = -1;

		double theta = acos(kos);
		
		double w = (theta - he.angle)*(theta - he.angle);

		bending += w;

		/* compute bending gradient */
		double dn1x,dn1y,dn1z,dn2x,dn2y,dn2z,dsm1,dsm2,dsq1,dsq2,dp,dg,dk,dtheta,dw;
		double dv1x,dv1y,dv1z,dv2x,dv2y,dv2z,dv3x,dv3y,dv3z;

		for (int j = 1; j <= 4; j++)
			for (int k = 1; k <= 3; k++){
				dv1x = computeGradient(1,1,j,k);dv1y = computeGradient(1,2,j,k); dv1z = computeGradient(1,3,j,k);
				dv2x = computeGradient(2,1,j,k);dv2y = computeGradient(2,2,j,k); dv2z = computeGradient(2,3,j,k);
				dv3x = computeGradient(3,1,j,k);dv3y = computeGradient(3,2,j,k); dv3z = computeGradient(3,3,j,k);

				dn1x = (dv1y*v2.z()+v1.y()*dv2z)-(dv1z*v2.y()+v1.z()*dv2y);
				dn1y = (dv1z*v2.x()+v1.z()*dv2x)-(dv1x*v2.z()+v1.x()*dv2z);
				dn1z = (dv1x*v2.y()+v1.x()*dv2y)-(dv1y*v2.x()+v1.y()*dv2x);
				dn2x = (dv2y*v3.z()+v2.y()*dv3z)-(dv2z*v3.y()+v2.z()*dv3y);
				dn2y = (dv2z*v3.x()+v2.z()*dv3x)-(dv2x*v3.z()+v2.x()*dv3z);
				dn2z = (dv2x*v3.y()+v2.x()*dv3y)-(dv2y*v3.x()+v2.y()*dv3x);

				dsm1 = 2*n1.x()*dn1x + 2*n1.y()*dn1y + 2*n1.z()*dn1z;
				dsm2 = 2*n2.x()*dn2x + 2*n2.y()*dn2y + 2*n2.z()*dn2z;

				dsq1 = 0.5 * 1 / sqrt(sm1) * dsm1;
				dsq2 = 0.5 * 1 / sqrt(sm2) * dsm2;

				dg = dsq1*sq2+sq1*dsq2;
				dp = (dn1x*n2.x()+n1.x()*dn2x) + (dn1y*n2.y()+n1.y()*dn2y) + (dn1z*n2.z()+n1.z()*dn2z);

				dk = (dp*gg-p*dg)/(gg*gg);

				if (abs(kos) == 1)
					dtheta = 0;
				else
					dtheta = -1 / sqrt(1-kos*kos)*dk;

				dw = 2*(theta - he.angle)*dtheta;

				g[he.index[j]*3+k-1] += dw * bendweight;
			}
		
	}

	float linkStretching = 0;
	for (int i = 0; i < RCSurface->linkNum; i++){
		int idx1 = RCSurface->linkList[i].index1;
		int idx2 = RCSurface->linkList[i].index2;

		Vector_3 deformVec(u[idx2*3] - u[idx1*3],u[idx2*3+1] - u[idx1*3+1],u[idx2*3+2] - u[idx1*3+2]);

		g[idx1*3] += -2 * LINKWEIGHT * deformVec.x();
		g[idx1*3+1] += -2 * LINKWEIGHT * deformVec.y();
		g[idx1*3+2] += -2 * LINKWEIGHT * deformVec.z();
		g[idx2*3] += 2 * LINKWEIGHT * deformVec.x();
		g[idx2*3+1] += 2 * LINKWEIGHT * deformVec.y();
		g[idx2*3+2] += 2 * LINKWEIGHT * deformVec.z();

		linkStretching += deformVec.squared_length();
	}

	//cout<<"fx:= "<<fx<<" streching:= "<<stretching<<" bending:= "<<bending<<endl;

	return fx + REGWEIGHT * stretching +  bendweight * bending + LINKWEIGHT * linkStretching;
}

double computeGradient(int n,int nd,int e,int ed){
	if (nd != ed) return 0;

	if (n == 1){
		if (e == 1) return -1;
		else if (e == 4) return 1;
		else return 0;
	}
	else if (n == 2){
		if (e == 1) return -1;
		else if (e == 2) return 1;
		else return 0;
	}
	else if (n == 3){
		if (e == 1) return -1;
		else if (e == 3) return 1;
		else return 0;
	}

	return 0;
}