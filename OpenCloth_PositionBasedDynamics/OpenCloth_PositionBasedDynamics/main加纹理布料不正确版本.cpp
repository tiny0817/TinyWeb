#include <iostream>
#include <GL/glew.h>
#include <GL/wglew.h>
#include <GL/freeglut.h>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp> 
#include <glm/gtc/type_ptr.hpp>
#include <gl/glaux.h>
#include <gl/glut.h>
#include <time.h>
using namespace std;

//#undef USE_TRIANGLE_BENDING_CONSTRAINT
//#undef USE_DIHEDRAL_ANGLE_BENDING_CONSTRAINT
//#define USE_TRIANGLE_BENDING_CONSTRAINT  1
//#define USE_DIHEDRAL_ANGLE_BENDING_CONSTRAINT 
#define USE_NEW_BENDING_CONSTRAINT

#pragma comment(lib, "glew32.lib")
#pragma comment(lib, "glaux.lib")

using namespace std;  
const int width = 1024, height = 1024;

#define PI 3.1415926536f
#define EPSILON  0.0000001f
 
int          numX = 20, numY=20;  //numX是宽  numY是长度
const size_t total_points = (numX+1)*(numY+1);
int          size = 10;
float        hsize = size/2.0f;

char   info[MAX_PATH]={0};

float  timeStep = 1.0f/60.0f; 
float  currentTime = 0;
double accumulator = timeStep;
int    selected_index = -1;
float  global_dampening = 0.93f; //0.98

struct DistanceConstraint {	int p1, p2;	float rest_length, k; float k_prime; };

#ifdef USE_TRIANGLE_BENDING_CONSTRAINT
struct BendingConstraint {	int p1, p2, p3;	float rest_length,  w,  k; float k_prime;};
#else
#ifdef USE_NEW_BENDING_CONSTRAINT
struct BendingConstraint { int p1, p2, p3; float k; float k_prime; float phi0;};
#else 
struct BendingConstraint {	int p1, p2, p3, p4;	float rest_length1, rest_length2, w1, w2,  k; float k_prime;};
#endif
#endif


vector<GLushort> indices;
vector<DistanceConstraint> d_constraints;
vector<BendingConstraint> b_constraints;

//particle system
vector<glm::vec3> X; //position
vector<glm::vec3> tmp_X; //predicted position
vector<glm::vec3> V; //velocity
vector<glm::vec3> F;
vector<float>     W; //inverse particle mass 
vector<glm::vec3> Ri; //Ri = Xi-Xcm 

int       oldX=0, oldY=0;
float     rX=15, rY=0;
int       state =1 ;
float     dist=-23;
const int GRID_SIZE=10;

const size_t solver_iterations = 4; //number of solver iterations per step. PBD  

float     kBend = 0.0018f; //0.5
float     kStretch = 0.3f;//0.25f; 
float     kDamp = 0.00125f;
glm::vec3 gravity=glm::vec3(0.0f,-0.00981f,0.0f);  
float     mass = 1.0f/( total_points);
 
GLint    viewport[4];
GLdouble MV[16];
GLdouble P[16];

LARGE_INTEGER frequency;        // ticks per second
LARGE_INTEGER t1, t2;           // ticks
double frameTimeQP=0;
float  frameTime =0 ;

glm::vec3 Up=glm::vec3(0,1,0), Right, viewDir;
float     startTime =0, fps=0;
int       totalFrames=0;
bool      IsStop = false;

void StepPhysics(float dt);

inline int getIndex(int i, int j) {
	return j*(numX+1) + i;
}

void AddDistanceConstraint(int a, int b, float k) {
	DistanceConstraint c;
	c.p1=a;
	c.p2=b;
	c.k =k;
	c.k_prime = 1.0f-pow((1.0f-c.k), 1.0f/solver_iterations); 
	
	if(c.k_prime>1.0) 
		c.k_prime = 1.0;
	 
	glm::vec3 deltaP = X[c.p1]-X[c.p2]; 
	c.rest_length = glm::length(deltaP);  

	d_constraints.push_back(c);
}
#ifdef USE_TRIANGLE_BENDING_CONSTRAINT
void AddBendingConstraint(int pa, int pb, int pc, float k) {
	BendingConstraint c;
	c.p1=pa;
	c.p2=pb;
	c.p3=pc; 
	
	c.w = W[pa] + W[pb] + 2*W[pc];  
	glm::vec3 center = 0.3333f * (X[pa] + X[pb] + X[pc]); 
 	c.rest_length = glm::length(X[pc]-center);
	c.k = k;
	c.k_prime = 1.0f-pow((1.0f-c.k), 1.0f/solver_iterations);  
	if(c.k_prime>1.0) 
		c.k_prime = 1.0;
	b_constraints.push_back(c);
}
#else
#ifdef USE_NEW_BENDING_CONSTRAINT 
void AddBendingConstraint(int pa, int pb, int pc, float k ) {
	BendingConstraint c;
	float d;
	glm::vec3 n1=glm::vec3(0), n2=glm::vec3(0);
	c.p1=pa;
	c.p2=pb;
	c.p3=pc;

	glm::vec3 np1 =X[c.p2]-X[c.p1];
	glm::vec3 np2 =X[c.p3]-X[c.p1];
	n1 = glm::normalize(np1);
	n2 = glm::normalize(np2);

	d = glm::dot(n1, n2);
	c.phi0 = acos(d);

	c.k = k;
	c.k_prime = 1.0f-pow((1.0f-c.k), 1.0f/solver_iterations);  
	if(c.k_prime>1.0) 
		c.k_prime = 1.0;
	b_constraints.push_back(c);
	
}

#else
void AddBendingConstraint(int pa, int pb, int pc,int pd, float k) {
	BendingConstraint c;
	c.p1=pa;
	c.p2=pb;
	c.p3=pc; 
	c.p4=pd; 
	c.w1 = W[pa] + W[pb] + 2*W[pc];
	c.w2 = W[pa] + W[pb] + 2*W[pd]; 
	glm::vec3 center1 = 0.3333f * (X[pa] + X[pb] + X[pc]);
	glm::vec3 center2 = 0.3333f * (X[pa] + X[pb] + X[pd]);
	c.rest_length1 = glm::length(X[pc]-center1);
	c.rest_length2 = glm::length(X[pd]-center2);
	c.k = k;

	c.k_prime = 1.0f-pow((1.0f-c.k), 1.0f/solver_iterations);  //1.0f-pow((1.0f-c.k), 1.0f/ns);
	if(c.k_prime>1.0) 
		c.k_prime = 1.0;
	b_constraints.push_back(c);
}

#endif
#endif

/*捕捉键盘事件响应处理函数*/
void processNormalKeys(unsigned char key,int x,int y) 
{
	if(key == 32)
	{
		IsStop = !IsStop;
	}
}

void OnMouseDown(int button, int s, int x, int y)
{
	if (s == GLUT_DOWN) 
	{
		oldX = x; 
		oldY = y; 
		int window_y = (height - y);
		float norm_y = float(window_y)/float(height/2.0);
		int window_x = x ;
		float norm_x = float(window_x)/float(width/2.0);
		float winZ=0;
		glReadPixels( x, height-y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
		if(winZ==1)
			winZ=0; 
		double objX=0, objY=0, objZ=0;
		gluUnProject(window_x,window_y, winZ,  MV,  P, viewport, &objX, &objY, &objZ);
		glm::vec3 pt(objX,objY, objZ); 
		size_t i=0;
		for(i=0;i<total_points;i++) {			 
			if( glm::distance(X[i],pt)<0.1) {
				selected_index = i;
				printf("Intersected at %d\n",i);
				break;
			}
		}
	}	

	if(button == GLUT_MIDDLE_BUTTON)
		state = 0;
	else
		state = 1;
	if(s==GLUT_UP) {
		selected_index= -1;
		glutSetCursor(GLUT_CURSOR_INHERIT);
	}
}

void OnMouseMove(int x, int y)
{
	if(selected_index == -1) {
		if (state == 0)
			dist *= (1 + (y - oldY)/60.0f); 
		else
		{
			rY += (x - oldX)/5.0f; 
			rX += (y - oldY)/5.0f; 
		} 
	} else {
		float delta = 1500/abs(dist);
		float valX = (x - oldX)/delta; 
		float valY = (oldY - y)/delta; 
		if(abs(valX)>abs(valY))
			glutSetCursor(GLUT_CURSOR_LEFT_RIGHT);
		else 
			glutSetCursor(GLUT_CURSOR_UP_DOWN); 
		V[selected_index] = glm::vec3(0);
		X[selected_index].x += Right[0]*valX ;
		float newValue = X[selected_index].y+Up[1]*valY;
		if(newValue>0)
			X[selected_index].y = newValue;
		X[selected_index].z += Right[2]*valX + Up[2]*valY;		
	}
	oldX = x; 
	oldY = y; 
	glutPostRedisplay(); 
}

void DrawGrid()
{
	glBegin(GL_LINES);
	glColor3f(0.5f, 0.5f, 0.5f);
	for(int i=-GRID_SIZE;i<=GRID_SIZE;i++)
	{
		glVertex3f((float)i,0,(float)-GRID_SIZE);
		glVertex3f((float)i,0,(float)GRID_SIZE);

		glVertex3f((float)-GRID_SIZE,0,(float)i);
		glVertex3f((float)GRID_SIZE,0,(float)i);
	}
	glEnd();
}


void InitGL() { 
	startTime = (float)glutGet(GLUT_ELAPSED_TIME);
	currentTime = startTime;

	// get ticks per second
    QueryPerformanceFrequency(&frequency);
    // start timer
    QueryPerformanceCounter(&t1);

	glEnable(GL_DEPTH_TEST);

	size_t i=0, j=0, count=0;
	int l1=0, l2=0;
	int v = numY+1;
	int u = numX+1;

	indices.resize( numX*numY*2*3);
	X.resize(total_points);
	tmp_X.resize(total_points);
	V.resize(total_points);
	F.resize(total_points); 
	Ri.resize(total_points); 
	W.resize(total_points);

	for(int j=0;j<=numY;j++) {		 
		for(int i=0;i<=numX;i++) {	 
			X[count++] = glm::vec3( ((float(i)/(u-1)) *2-1)* hsize, size+1, ((float(j)/(v-1) )* size)); 
			//X[count++] = glm::vec3( ((float(i)/(u-1)) *2-1)* hsize,-(1-((float(j)/(v-1)) )* size),0); 
		}
	}

	W.resize(total_points); //W中每个元素对应一个元素点 
	for(i=0;i<total_points;i++) {	
		W[i] = 1.0f/mass; //质量权重
	}
	
	W[0] = 0.0;
	W[numX] = 0.0;
	/*W[420] = 0.0;
	W[440] = 0.0;*/

	memcpy(&tmp_X[0].x, &X[0].x, sizeof(glm::vec3)*total_points);  
	//fill in velocities	 
	memset(&(V[0].x),0,total_points*sizeof(glm::vec3)); 
	//fill in indices

	GLushort* id=&indices[0];
	for (int i = 0; i < numY; i++) {        
		for (int j = 0; j < numX; j++) {            
			int i0 = i * (numX+1) + j; 
			int i1 = i0 + 1;            
			int i2 = i0 + (numX+1);     
			int i3 = i2 + 1;           
			if ((j+i)%2) {            //index  index+(numX+1) index+1 index+1 index+(numX+1) index+(numX+1)+1            
				*id++ = i0; *id++ = i2; *id++ = i1;                
				*id++ = i1; *id++ = i2; *id++ = i3;            
			} else {                  //index  index+(numX+1) index+(numX+1)+1 index index+(numX+1)+1 index+1
				*id++ = i0; *id++ = i2; *id++ = i3;                
				*id++ = i0; *id++ = i3; *id++ = i1;            
			}        
		}    
	}
	 
//	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); 
	//glPolygonMode(GL_BACK, GL_LINE);
	glPointSize(5);
    //wglSwapIntervalEXT(0);

	//check the damping values
	if(kStretch>1)
		kStretch=1;  
	if(kStretch<0)
		kStretch=0;
	if(kBend>1)
		kBend=1;   
	if(kBend<0)
		kBend=0;
	if(kDamp>1)
		kDamp=1;  
	if(kDamp<0)
		kDamp=0;
	if(global_dampening>1)
		global_dampening = 1;   

	
	// Horizontal 
	for (l1 = 0; l1 < v; l1++)	// v=(numY)+1
		for (l2 = 0; l2 < (u - 1); l2++) {  //u=(numX)+1
			AddDistanceConstraint((l1 * u) + l2,(l1 * u) + l2 + 1, kStretch);  
		}
	// Vertical 
	for (l1 = 0; l1 < u; l1++)	
		for (l2 = 0; l2 < (v - 1); l2++) {
			AddDistanceConstraint((l2 * u) + l1,((l2 + 1) * u) + l1, kStretch); 
		}	
	// Shearing distance constraint
	for (l1 = 0; l1 < (v - 1); l1++)	
		for (l2 = 0; l2 < (u - 1); l2++) {
			AddDistanceConstraint((l1 * u) + l2,((l1 + 1) * u) + l2 + 1, kStretch);  
			AddDistanceConstraint(((l1 + 1) * u) + l2,(l1 * u) + l2 + 1, kStretch);   
		}
	
	//#ifdef USE_TRIANGLE_BENDING_CONSTRAINT

	//add vertical constraints
	for(int i=0;i<=numX;i++) {
		for(int j=0;j<numY-1 ;j++) {
			//AddBendingConstraint(getIndex(i,j), getIndex(i,(j+1)), getIndex(i,j+2), kBend);
			AddBendingConstraint(getIndex(i,(j+1)),getIndex(i,j),  getIndex(i,j+2), kBend);
		}
	}
	//add horizontal constraints
	for(int i=0;i<numX-1;i++) {
		for(int j=0;j<=numY;j++) {	   
			//AddBendingConstraint(getIndex(i,j), getIndex(i+1,j), getIndex(i+2,j), kBend);
			AddBendingConstraint(getIndex(i+1,j),getIndex(i,j),  getIndex(i+2,j), kBend);
		}
	}
	
}

void OnReshape(int nw, int nh) {
	glViewport(0,0,nw, nh);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, (GLfloat)nw / (GLfloat)nh, 1.f, 100.0f);
	
	glGetIntegerv(GL_VIEWPORT, viewport); 
	glGetDoublev(GL_PROJECTION_MATRIX, P);

	glMatrixMode(GL_MODELVIEW);
}

void OnRender() 
{		
	size_t i=0;
	float newTime = (float) glutGet(GLUT_ELAPSED_TIME);
	frameTime = newTime-currentTime;
	currentTime = newTime;
    QueryPerformanceCounter(&t2);
	
	frameTimeQP = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;
	t1=t2;
	accumulator += frameTimeQP;

	++totalFrames;
	if((newTime-startTime)>1000)
	{		
		float elapsedTime = (newTime-startTime);
		fps = (totalFrames/ elapsedTime)*1000 ;
		startTime = newTime;
		totalFrames=0;
	}

    sprintf_s(info, "FPS: %3.2f, Frame time (GLUT): %3.4f msecs, Frame time (QP): %3.3f FrameNum:%d", fps, frameTime, frameTimeQP,totalFrames);
	//sprintf_s(info, "Frame time (GLUT): %3.4f msecs, Frame time (QP): %3.3f ,FrameNum:%d,Time:%3.4f", frameTime, frameTimeQP,totalFrames,accumulator);
	glutSetWindowTitle(info);

	glClear(GL_COLOR_BUFFER_BIT| GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	glTranslatef(0,0,dist);
	glRotatef(rX,1,0,0);
	glRotatef(rY,0,1,0);

	glGetDoublev(GL_MODELVIEW_MATRIX, MV);
	viewDir.x = (float)-MV[2];
	viewDir.y = (float)-MV[6];
	viewDir.z = (float)-MV[10];
	Right = glm::cross(viewDir, Up);

	//draw grid 
	DrawGrid();

//	glColor3f(1,1,1);
	
	for(i=0;i<indices.size();i+=3) {
		glm::vec3 p1 = X[indices[i]];  
		glm::vec3 p2 = X[indices[i+1]];  
		glm::vec3 p3 = X[indices[i+2]];

		

		GLfloat color[3];
		if( i % 2 == 0 )
		{
			color[0]=0.6f;color[1]=0.2f;color[2]=0.2f;
		}else{
		//	color[0]=0.6f;color[1]=0.2f;color[2]=0.2f;
			color[0]=1.0f;color[1]=1.0f;color[2]=1.0f;
		}
		glColor3fv(color);

		glBegin(GL_TRIANGLES);
		
		 glVertex3f(p1.x,p1.y,p1.z);  
		 glVertex3f(p2.x,p2.y,p2.z);
	     glVertex3f(p3.x,p3.y,p3.z);
		 glEnd();

	//	 glFlush();  
	}
	




#if 0
	//draw points
	glBegin(GL_POINTS);
	for(i=0;i<total_points;i++) {
		glm::vec3 p = X[i];
		int is = (i==selected_index);
		glColor3f((float)!is,(float)is,(float)is);
		glVertex3f(p.x,p.y,p.z);
	}

	glEnd();
#endif	
	glutSwapBuffers();
}

void OnShutdown() {	
	d_constraints.clear();
	b_constraints.clear();
	indices.clear();
	X.clear();
	F.clear();
	V.clear();
	W.clear();
	tmp_X.clear();
	Ri.clear();
}

void ComputeForces( ) {
	size_t i=0;
	for(i=0;i<total_points;i++) {
		F[i] = glm::vec3(0);  
		if(W[i]>0)		 
			F[i] += gravity ; //gravity=glm::vec3(0.0f,-0.00981f,0.0f);
	}	
} 

void IntegrateExplicitWithDamping(float deltaTime) {
	float deltaTimeMass = deltaTime;
	size_t i=0;
 
	glm::vec3 Xcm = glm::vec3(0);
	glm::vec3 Vcm = glm::vec3(0);
	float sumM = 0;
	
	for(i=0;i<total_points;i++) {
		V[i] *= global_dampening; //为了减小速度，个人理解为假设空气中存在阻尼力会对速度有影响	
		V[i] = V[i] + (F[i]*deltaTime)*W[i]; 	 					
		Xcm += (X[i]*mass); 
		Vcm += (V[i]*mass);  
		sumM += mass;  
	} 
	Xcm /= sumM; //sumM = 1
	Vcm /= sumM; //sumM = 1 

	glm::mat3 I = glm::mat3(1);
	glm::vec3 L = glm::vec3(0);
	glm::vec3 w = glm::vec3(0);

	for(i=0;i<total_points;i++) { 
		Ri[i] = (X[i] - Xcm);	
		L += glm::cross(Ri[i],mass*V[i]); 
		glm::mat3 tmp = glm::mat3(0,-Ri[i].z,  Ri[i].y, 
							 Ri[i].z,       0,-Ri[i].x,
							 -Ri[i].y,Ri[i].x,        0);
		I +=(tmp*glm::transpose(tmp))*mass;  
	} 
	w = glm::inverse(I)*L;  
	
	for(i=0;i<total_points;i++) {
		glm::vec3 delVi = Vcm + glm::cross(w,Ri[i])-V[i];		
		V[i] += kDamp*delVi;
	}

	//calculate predicted position
	for(i=0;i<total_points;i++) {
		if(W[i] <= 0.0) { 
			tmp_X[i] = X[i]; //fixed points 
		} else {
			tmp_X[i] = X[i] + (V[i]*deltaTime);	 
		}
	} 
}
 
void Integrate(float deltaTime) {	
	float inv_dt = 1.0f/deltaTime;
	size_t i=0; 

	for(i=0;i<total_points;i++) {	
		V[i] = (tmp_X[i] - X[i])*inv_dt; 	
		X[i] = tmp_X[i]; 
}
}

void UpdateDistanceConstraint(int i){

	DistanceConstraint c = d_constraints[i];  
	glm::vec3 dir = tmp_X[c.p1] - tmp_X[c.p2]; 

	float len = glm::length(dir); 
	if(len <= EPSILON) 
		return;
	
	float w1 = W[c.p1];  
	float w2 = W[c.p2];   
	float invMass = w1+ w2; 
	if(invMass <= EPSILON) 
		return;

	glm::vec3 dP = (1.0f/invMass) * (len-c.rest_length ) * (dir/len)* c.k_prime;
	if(w1 > 0.0)
		tmp_X[c.p1] -= dP*w1;	

	if(w2 > 0.0)
		tmp_X[c.p2] += dP*w2;
}

void UpdateBendingConstraint(int index) {
	size_t i=0;
	BendingConstraint c = b_constraints[index]; 

#ifdef USE_TRIANGLE_BENDING_CONSTRAINT
	float global_k = global_dampening*0.01f; //0.01f
	glm::vec3 center = 0.3333f * (tmp_X[c.p1] + tmp_X[c.p2] + tmp_X[c.p3]);
	glm::vec3 dir_center = tmp_X[c.p3]-center;
	float dist_center = glm::length(dir_center);

	float diff = 1.0f - ((global_k + c.rest_length) / dist_center);
	glm::vec3 dir_force = dir_center * diff;
	glm::vec3 fa = c.k_prime * ((2.0f*W[c.p1])/c.w) * dir_force;
	glm::vec3 fb = c.k_prime * ((2.0f*W[c.p2])/c.w) * dir_force;
	glm::vec3 fc = -c.k_prime * ((4.0f*W[c.p3])/c.w) * dir_force;

	if(W[c.p1] > 0.0)  {
		tmp_X[c.p1] += fa;
	}
	if(W[c.p2] > 0.0) {
		tmp_X[c.p2] += fb;
	}
	if(W[c.p3] > 0.0) {
		tmp_X[c.p3] += fc;
	}
#else
#ifdef USE_NEW_BENDING_CONSTRAINT
	float d = 0, phi=0, i_d=0 ;
	glm::vec3 n1=glm::vec3(0), n2=glm::vec3(0);

	glm::vec3 p1 = tmp_X[c.p1];
	glm::vec3 p2 = tmp_X[c.p2]-p1;
	glm::vec3 p3 = tmp_X[c.p3]-p1;

	float lenp2 = glm::length(p2);
	if(lenp2 == 0.0) { return; } 

	float lenp3 = glm::length(p3);
	if(lenp3 == 0.0) { return; } 

	n1 = glm::normalize(p2);
	n2 = glm::normalize(p3);

	d = glm::dot(n1,n2);
	phi = acos(d);

	if(d<-1.0) 
		d = -1.0; 
	else if(d>1.0) 
		d=1.0;

	if(d == 1.0){ 
		phi = 0.0; 
		if(phi == c.phi0)
		return;  
	}
	if(d == -1.0){ 
		phi = PI;  
		if(phi == c.phi0)
			if(c.p1!=0 && c.p1!=numX)
			tmp_X[c.p2] += n1/100.0f;
			tmp_X[c.p3] += n2/100.0f;
		return; 
	}

	i_d =sqrt(1-(d*d))*(phi-c.phi0) ;

	glm::vec3 q1 = (n2-n1*d)/lenp2;
	glm::vec3 q2 = (n1-n2*d)/lenp3;
	glm::vec3 q3 =-q1-q2;

	float q1_len2 = glm::dot(q1,q1);  
	float q2_len2 = glm::dot(q2,q2);
	float q3_len2 = glm::dot(q3,q3);

	float sum = W[c.p1]*(q1_len2) +W[c.p2]*(q2_len2) +W[c.p3]*(q3_len2);

	glm::vec3 dP1 = -( (W[c.p1] * i_d) /sum)*q1;
	glm::vec3 dP2 = -( (W[c.p2] * i_d) /sum)*q2;
	glm::vec3 dP3 = -( (W[c.p3] * i_d) /sum)*q3;

	if(W[c.p1] > 0.0) {
		tmp_X[c.p1] += dP1*c.k;
	}
	if(W[c.p2] > 0.0) {
		tmp_X[c.p2] += dP2*c.k;
	}
	if(W[c.p3] > 0.0) {
		tmp_X[c.p3] += dP3*c.k;
	}
#else 
	//Using the dihedral angle approach of the position based dynamics		
	float d = 0, phi=0,i_d=0;
	glm::vec3 n1=glm::vec3(0), n2=glm::vec3(0);

	glm::vec3 p1 = tmp_X[c.p1];
	glm::vec3 p2 = tmp_X[c.p2]-p1;
	glm::vec3 p3 = tmp_X[c.p3]-p1;
	glm::vec3 p4 = tmp_X[c.p4]-p1;

	glm::vec3 p2p3 = glm::cross(p2,p3);		
	glm::vec3 p2p4 = glm::cross(p2,p4);		

	float lenp2p3 = glm::length(p2p3);

	if(lenp2p3 == 0.0) { return; } //need to handle this case.

	float lenp2p4 = glm::length(p2p4);

	if(lenp2p4 == 0.0) { return; } //need to handle this case.

	n1 = glm::normalize(p2p3);
	n2 = glm::normalize(p2p4); 

	d	= glm::dot(n1,n2);
	phi = acos(d);

	//try to catch invalid values that will return NaN.
	// sqrt(1 - (1.0001*1.0001)) = NaN 
	// sqrt(1 - (-1.0001*-1.0001)) = NaN 
	if(d<-1.0) 
		d = -1.0; 
	else if(d>1.0) 
		d=1.0; //d = clamp(d,-1.0,1.0);

	//in both case sqrt(1-d*d) will be zero and nothing will be done.
	//0?case, the triangles are facing in the opposite direction, folded together.
	if(d == -1.0){ 
		phi = PI;  //acos(-1.0) == PI
		if(phi == phi0[index]) 
			return; //nothing to do 

		//in this case one just need to push 
		//vertices 1 and 2 in n1 and n2 directions, 
		//so the constrain will do the work in second iterations.
		if(c.p1!=0 && c.p1!=numX)
			tmp_X[c.p3] += n1/100.0f;

		if(c.p2!=0 && c.p2!=numX)
			tmp_X[c.p4] += n2/100.0f;

		return;
	}
	if(d == 1.0){ //180?case, the triangles are planar
		phi = 0.0;  //acos(1.0) == 0.0
		if(phi == phi0[index]) 
			return; //nothing to do 
	}

	i_d = sqrt(1-(d*d))*(phi-phi0[index]) ;

	glm::vec3 p2n1 = glm::cross(p2,n1);
	glm::vec3 p2n2 = glm::cross(p2,n2);
	glm::vec3 p3n2 = glm::cross(p3,n2);
	glm::vec3 p4n1 = glm::cross(p4,n1);
	glm::vec3 n1p2 = -p2n1;
	glm::vec3 n2p2 = -p2n2;
	glm::vec3 n1p3 = glm::cross(n1,p3);
	glm::vec3 n2p4 = glm::cross(n2,p4);

	glm::vec3 q3 =  (p2n2 + n1p2*d)/ lenp2p3;
	glm::vec3 q4 =  (p2n1 + n2p2*d)/ lenp2p4;
	glm::vec3 q2 =  (-(p3n2 + n1p3*d)/ lenp2p3) - ((p4n1 + n2p4*d)/lenp2p4);

	glm::vec3 q1 = -q2-q3-q4;

	float q1_len2 = glm::dot(q1,q1);// glm::length(q1)*glm::length(q1);
	float q2_len2 = glm::dot(q2,q2);// glm::length(q2)*glm::length(q1);
	float q3_len2 = glm::dot(q3,q3);// glm::length(q3)*glm::length(q1);
	float q4_len2 = glm::dot(q4,q4);// glm::length(q4)*glm::length(q1); 

	float sum = W[c.p1]*(q1_len2) +
		W[c.p2]*(q2_len2) +
		W[c.p3]*(q3_len2) +
		W[c.p4]*(q4_len2);	

	glm::vec3 dP1 = -( (W[c.p1] * i_d) /sum)*q1;
	glm::vec3 dP2 = -( (W[c.p2] * i_d) /sum)*q2;
	glm::vec3 dP3 = -( (W[c.p3] * i_d) /sum)*q3;
	glm::vec3 dP4 = -( (W[c.p4] * i_d) /sum)*q4;

	if(W[c.p1] > 0.0) {
		tmp_X[c.p1] += dP1*c.k;
	}
	if(W[c.p2] > 0.0) {
		tmp_X[c.p2] += dP2*c.k;
	}
	if(W[c.p3] > 0.0) {
		tmp_X[c.p3] += dP3*c.k;
	}	
	if(W[c.p4] > 0.0) {
		tmp_X[c.p4] += dP4*c.k;
	}  	
#endif
#endif
}

void GroundCollision() 
{
	for(size_t i=0;i<total_points;i++) {	
		if(tmp_X[i].y<0) 
			tmp_X[i].y=0;
	}
}

void UpdateInternalConstraints(float deltaTime) {
	size_t i=0;
 
	for (size_t si=0;si<solver_iterations;++si) {

		for(i=0;i<d_constraints.size();i++) {
			UpdateDistanceConstraint(i); 
		}

		clock_t start = clock();  //获取当前时间，单位是毫秒
	    for(i=0;i<b_constraints.size();i++) {
			UpdateBendingConstraint(i);
		}
		clock_t finsh = clock(); //获取当前时间，单位是毫秒

		cout<<(double)finsh-start<<"(ms),total="<<b_constraints.size()
			<<" , avg="<<((double)finsh-start)/b_constraints.size()<<endl;

      
	GroundCollision();  
	}
}

void OnIdle() {	
	//printf(" ### OnIdle %f ### \n",accumulator);
	//Fixed time stepping + rendering at different fps	
	if(IsStop == false)
	{
		if ( accumulator >= timeStep ) //时间步长，使得不是每一次调用idle函数都调用，保持在60的帧率
		{	 
			StepPhysics(timeStep );	 
			accumulator -= timeStep;
		}
		glutPostRedisplay();
	}
	Sleep(5); //TODO
}

void StepPhysics(float dt ) {
	ComputeForces();  
	IntegrateExplicitWithDamping(dt);  
	UpdateInternalConstraints(dt);
    Integrate(dt);
}

void main(int argc, char** argv) {
	
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
//	glEnable(GL_DEPTH_TEST);
	glutInitWindowSize(width, height);
	glutCreateWindow("GLUT Cloth Demo [Position based Dynamics]");

	glutDisplayFunc(OnRender);
	glutReshapeFunc(OnReshape);
	glutIdleFunc(OnIdle);

	glutMouseFunc(OnMouseDown);
	glutMotionFunc(OnMouseMove);	
	glutKeyboardFunc(processNormalKeys); 
	glutCloseFunc(OnShutdown);

	glewInit();
	InitGL();

	glutMainLoop();		
}
 
