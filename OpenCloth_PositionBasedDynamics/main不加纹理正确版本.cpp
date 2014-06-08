#include <GL/glew.h>
#include <GL/wglew.h>
#include <GL/freeglut.h>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp> //for matrices
#include <glm/gtc/type_ptr.hpp>
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include "glm.h"
using namespace std;

//#undef USE_TRIANGLE_BENDING_CONSTRAINT
#define USE_TRIANGLE_BENDING_CONSTRAINT

#pragma comment(lib, "glew32.lib")
using namespace std;  
const int width = 1024, height = 1024;

#define PI 3.1415926536f
#define EPSILON  0.0000001f
int size = 4;
float hsize = size/2.0f;

char info[MAX_PATH]={0};

float timeStep = 1.0f/60.0f; //1.0/60.0f;
float currentTime = 0;
double accumulator = timeStep;
int selected_index = -1;
float global_dampening = 0.93f; //DevO: 24.07.2011  //global velocity dampening !!!

struct DistanceConstraint {	int p1, p2;	float rest_length, k; float k_prime; };
#ifdef USE_TRIANGLE_BENDING_CONSTRAINT
struct BendingConstraint {	int p1, p2, p3;	float rest_length,  w,  k; float k_prime;};
#else
struct BendingConstraint { int p1, p2, p3; float k; float k_prime; float phi0;};
#endif

int oldX=0, oldY=0;
float rX=15, rY=0;
int state =1 ;
float dist=-3.5;
const int GRID_SIZE=10;

const size_t solver_iterations = 4; //number of solver iterations per step. PBD  
int isPause = 0;

float kBend = 1.0f; 
float kStretch = 0.99f; 
float kDamp = 0.125f; //0.125f
glm::vec3 gravity=glm::vec3(0.0f,-0.00981f,0.0f);  

GLint viewport[4];
GLdouble MV[16];
GLdouble P[16];

LARGE_INTEGER frequency;        // ticks per second
LARGE_INTEGER t1, t2;           // ticks
double frameTimeQP=0;
float frameTime =0 ;


glm::vec3 Up=glm::vec3(0,1,0), Right, viewDir;
float startTime =0, fps=0;
int totalFrames=0;

GLMmodel* pmodel = NULL;
glm::mat4 bunny_pos,inverse_bunny_pos;

//particle system
vector<glm::vec3> X_bunny; //position
vector<glm::vec3> tmp_X_bunny; //predicted position
vector<glm::vec3> V_bunny; //velocity
vector<glm::vec3> F_bunny;
vector<float> W_bunny; //inverse particle mass 
vector<glm::vec3> Ri_bunny; //Ri = Xi-Xcm 
int total_points_bunny;
float mass_bunny;

vector<DistanceConstraint> d_constraints_bunny; 

typedef struct _triangle_node{
	int ver1_index;
	int ver2_index;
	int ver3_index;
}triangle_node;

vector<triangle_node> triangle_bunny;
vector<BendingConstraint> b_constraints_bunny;
vector<vector<int>> ring_neighborhood_bunny;
vector<vector<int>> vertex_map_to_triangle;
vector<glm::vec3> clothBalloons_p;

float total_volumn = 0.0f;
float kprime_clothBalloons;
float k_clothBalloons = 1.0f;

void StepPhysics(float dt);

void AddClothBalloons()
{
	total_volumn = 0.0f;
	for(int i=0; i<triangle_bunny.size(); i++)
	{
		int ver1 = triangle_bunny[i].ver1_index;
		int ver2 = triangle_bunny[i].ver2_index;
		int ver3 = triangle_bunny[i].ver3_index;
		if(ver1 >= X_bunny.size() || ver2 >= X_bunny.size() || ver3 >= X_bunny.size())
		{
			continue;
		}

		glm::vec3 temp = glm::cross(X_bunny[ver1],X_bunny[ver2]);
		float result = glm::dot(temp, X_bunny[ver3]);
		total_volumn += result;
	}
//	total_volumn = total_volumn*10;

	kprime_clothBalloons = 1.0f-pow((1.0f-k_clothBalloons), 1.0f/solver_iterations); 
	if(kprime_clothBalloons > 1.0)
	{ 
		kprime_clothBalloons = 1.0f;
	}
	kprime_clothBalloons = 1.0f;
}
float tmp_total_volumn =0.0f;

void UpdateClothBalloons()
{
	float kpressure = 1.0f;
	tmp_total_volumn = 0.0f;

	for(int i=0; i<triangle_bunny.size(); i++)
	{
		int ver1 = triangle_bunny[i].ver1_index;
		int ver2 = triangle_bunny[i].ver2_index;
		int ver3 = triangle_bunny[i].ver3_index;
		if(ver1 >= tmp_X_bunny.size() || ver2 >= tmp_X_bunny.size() || ver3 >= tmp_X_bunny.size() )
		{
			//		printf("data error\r\n");
			continue;

		}
		glm::vec3 temp = glm::cross(tmp_X_bunny[ver1], tmp_X_bunny[ver2] );
		float result = glm::dot(temp, tmp_X_bunny[ver3]);
		tmp_total_volumn += result;
	}

	float c_upper = tmp_total_volumn - kpressure*total_volumn;
	//vertex_map_to_triangle;

	clothBalloons_p.clear();
	clothBalloons_p.resize(X_bunny.size());

	float c_down = 0.0f;

	for(int i = 0; i < tmp_X_bunny.size();  i++)
	{
		vector<int> involving_tris = vertex_map_to_triangle[i];
		glm::vec3 vec_pi; 
		for(int j = 0; j < involving_tris.size(); j++)
		{
			int tris_index = involving_tris[j];
			int ver1 = triangle_bunny[tris_index].ver1_index;
			int ver2 = triangle_bunny[tris_index].ver2_index;
			int ver3 = triangle_bunny[tris_index].ver3_index;

			glm::vec3 temp;
			if(ver1 == i)
			{
				temp = glm::cross(tmp_X_bunny[ver2], tmp_X_bunny[ver3]);
			}else if(ver2 == i)
			{
				temp = glm::cross(tmp_X_bunny[ver3], tmp_X_bunny[ver1]);
			}else if(ver3 == i)
			{
				temp = glm::cross(tmp_X_bunny[ver1], tmp_X_bunny[ver2]);
			}else{
				printf("logic error!\r\n");
			}
			vec_pi += temp;		
		}
		clothBalloons_p[i] = vec_pi;
		c_down += W_bunny[i]*glm::length(vec_pi)*glm::length(vec_pi);
	}

	for(int i = 0; i < X_bunny.size(); i++)
	{
		tmp_X_bunny[i] -= kprime_clothBalloons*W_bunny[i]*(c_upper/c_down)*clothBalloons_p[i];
	}
}

void AddDistanceConstraint(int a, int b, float k) {
	DistanceConstraint c;
	c.p1=a;
	c.p2=b;
	c.k =k;
	c.k_prime = 1.0f-pow((1.0f-c.k), 1.0f/solver_iterations);  //1.0f-pow((1.0f-c.k), 1.0f/ns);

	if(c.k_prime>1.0) 
		c.k_prime = 1.0;

	glm::vec3 deltaP = X_bunny[c.p1]-X_bunny[c.p2]; 
	c.rest_length = glm::length(deltaP);  
	d_constraints_bunny.push_back(c);
}

#ifdef USE_TRIANGLE_BENDING_CONSTRAINT
void AddBendingConstraint(int pa, int pb, int pc, float k) {
	BendingConstraint c;
	c.p1=pa;
	c.p2=pb;
	c.p3=pc; 

	c.w = W_bunny[pa] + W_bunny[pb] + 2*W_bunny[pc];  
	glm::vec3 center = 0.3333f * (X_bunny[pa] + X_bunny[pb] + X_bunny[pc]); 
	c.rest_length = glm::length(X_bunny[pc]-center); 
	c.k = k;
	c.k_prime = 1.0f-pow((1.0f-c.k), 1.0f/solver_iterations); 
	if(c.k_prime>1.0) 
		c.k_prime = 1.0;
	b_constraints_bunny.push_back(c);
}
#else
void AddBendingConstraint(int pa, int pb, int pc, float k ) {
	BendingConstraint c;
	float d;
	glm::vec3 n1=glm::vec3(0), n2=glm::vec3(0);
	c.p1=pa;
	c.p2=pb;
	c.p3=pc;

	glm::vec3 np1 =X_bunny[c.p2]-X_bunny[c.p1];
	glm::vec3 np2 =X_bunny[c.p3]-X_bunny[c.p1];
	n1 = glm::normalize(np1);
	n2 = glm::normalize(np2);

	d = glm::dot(n1, n2);
	c.phi0 = acos(d);

	c.k = k;
	c.k_prime = 1.0f-pow((1.0f-c.k), 1.0f/solver_iterations);  
	if(c.k_prime>1.0) 
		c.k_prime = 1.0;
	b_constraints_bunny.push_back(c);

}
#endif

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
		for(i=0;i<total_points_bunny;i++) {			 
			if( glm::distance(X_bunny[i],pt)<0.1) {
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
		{
			dist *= (1 + (y - oldY)/60.0f); 
		}
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

		V_bunny[selected_index] = glm::vec3(0);
		X_bunny[selected_index].x += Right[0]*valX ;
		float newValue =X_bunny[selected_index].y+Up[1]*valY;
		if(newValue>0)
			X_bunny[selected_index].y = newValue;
		X_bunny[selected_index].z += Right[2]*valX + Up[2]*valY;		
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

//inline int getIndex(int i, int j) {
//	return j*(numX+1) + i;
//}

void loadmodel(void)
{
	if (!pmodel) {
		pmodel = glmReadOBJ("data/soccerball.obj");  //bunny.obj
		if (!pmodel) exit(0);
		glmUnitize(pmodel);
		glmFacetNormals(pmodel);
		glmVertexNormals(pmodel, 90.0);
	}
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
	float ypos = 7.0f;

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	//glPolygonMode(GL_BACK, GL_LINE);
	glPointSize(5);

	if(kStretch>1)
		kStretch=1; //全局的拉伸系数,这里使得值得范围在0-1之间 
	if(kStretch<0)
		kStretch=0;

	if(kBend>1)
		kBend=1;   //全局的弯曲系数
	if(kBend<0)
		kBend=0;

	if(kDamp>1)
		kDamp=1;  //全局的阻尼系数
	if(kDamp<0)
		kDamp=0;
	if(global_dampening>1)
		global_dampening = 1;   //全局的阻尼系数

	//加载兔子的数据以及初始化兔子的位置
	loadmodel();
	bunny_pos =  glm::translate(glm::mat4(1),glm::vec3(0,0,0));
	bunny_pos = glm::rotate(bunny_pos, 0.0f ,glm::vec3(1,0,0));
	//bunny_pos = glm::scale(bunny_pos, glm::vec3(fRadius,fRadius,fRadius/2));
	inverse_bunny_pos = glm::inverse(bunny_pos);

	assert(pmodel);
	assert(pmodel->vertices);
	X_bunny.resize(pmodel->numvertices);
	GLfloat *vertices_pos = pmodel->vertices;

	

	for(int i = 0; i < pmodel->numvertices; i++)
	{
		X_bunny[i] = glm::vec3(*vertices_pos, *(vertices_pos+1)+1, *(vertices_pos+2));
		X_bunny[i] = 10.0f*X_bunny[i];
		vertices_pos += 3;

	}
	//缺失的那个点
	X_bunny.push_back(glm::vec3(0.6,0.2,0.2)); 

	total_points_bunny = X_bunny.size();

	GLMtriangle *tmp_triangle_pos = pmodel->triangles;
	triangle_bunny.resize(pmodel->numtriangles); //bunny -40
	for( int i = 0; i < pmodel->numtriangles; i++,tmp_triangle_pos++)  //bunny 40
	{
		triangle_node element;
		element.ver1_index = tmp_triangle_pos->vindices[0];
		element.ver2_index = tmp_triangle_pos->vindices[1];
		element.ver3_index = tmp_triangle_pos->vindices[2];
		
		triangle_bunny[i] = element;
	}

	//这个里面出现重复的情况
	vertex_map_to_triangle.resize(X_bunny.size());

	for(int i=0; i< triangle_bunny.size(); i++)
	{
		int ver1 = triangle_bunny[i].ver1_index;
		int ver2 = triangle_bunny[i].ver2_index;
		int ver3 = triangle_bunny[i].ver3_index;

		if(ver1 >= X_bunny.size() || ver2 >= X_bunny.size() ||ver3 >= X_bunny.size() )
		{
			printf("model is error\r\n");
			continue;
		}
		vertex_map_to_triangle[ver1].push_back(i);
		vertex_map_to_triangle[ver2].push_back(i);
		vertex_map_to_triangle[ver3].push_back(i);
	}

	mass_bunny = 1.0/(total_points_bunny-1);

	tmp_X_bunny.resize(total_points_bunny);
	V_bunny.resize(total_points_bunny);
	F_bunny.resize(total_points_bunny); 
	W_bunny.resize(total_points_bunny);
	Ri_bunny.resize(total_points_bunny);

	for(int i=1; i < X_bunny.size(); i++)
	{
		W_bunny[i] = 1/mass_bunny;
	}	

	memcpy(&tmp_X_bunny[0].x, &X_bunny[0].x, sizeof(glm::vec3)*total_points_bunny);   
	memset(&(V_bunny[0].x),0,total_points_bunny*sizeof(glm::vec3));

	// 遍历全部三角形，增加约束
	for (int i = 0; i < triangle_bunny.size(); i++) 
	{  //u=(numX)+1
		int ver_index1 = triangle_bunny[i].ver1_index;
		int ver_index2 = triangle_bunny[i].ver2_index;
		int ver_index3 = triangle_bunny[i].ver3_index;

		if(ver_index1 >= X_bunny.size() || ver_index2 >= X_bunny.size() ||ver_index3 >= X_bunny.size() )
		{
			continue;
		}
		AddDistanceConstraint(ver_index1,ver_index2, kStretch); 
		AddDistanceConstraint(ver_index1,ver_index3, kStretch); 
		AddDistanceConstraint(ver_index2,ver_index3, kStretch); 
	}


	//	ring_neighborhood_bunny.resize(X_bunny.size());
	vector<int> sort_order;
	map<string,int> check_bend_constraint;


	
	//这个是最耗时的部分，需要进一步的进行改进流程 
	for(int i = 0; i< X_bunny.size(); i++)
	{
		int index = i;

		vector<int> ring_neiborhood;
		
	//	printf("%d------------------------\r\n", i);

		for(int j=0; j<triangle_bunny.size(); j++)
		{
			int ver1_index = triangle_bunny[j].ver1_index;
			int ver2_index = triangle_bunny[j].ver2_index;
			int ver3_index = triangle_bunny[j].ver3_index;
			if(ver1_index == i || ver2_index == i || ver3_index == i)
			{
				if(ver1_index == i)
				{
					ring_neiborhood.push_back(ver2_index);
					ring_neiborhood.push_back(ver3_index);
				}else if(ver2_index == i)
				{
					ring_neiborhood.push_back(ver1_index);
					ring_neiborhood.push_back(ver3_index);
				}else if(ver3_index == i)
				{
					ring_neiborhood.push_back(ver1_index);
					ring_neiborhood.push_back(ver2_index);
				}
			}
		}
		//删除重复的结点信息
		std::sort(ring_neiborhood.begin(), ring_neiborhood.end());
		vector<int>::iterator pos;
		pos = std::unique(ring_neiborhood.begin(), ring_neiborhood.end());
		ring_neiborhood.erase(pos, ring_neiborhood.end());


		ring_neighborhood_bunny.push_back(ring_neiborhood);
	}  //计算相邻的节点


	////将计算出来的结点之间的关系存放在文件中，下次的时候就直接读文件就可以了
	//FILE *fp_bunny_relation = fopen("bunny_relation.txt", "wb");
	//if(fp_bunny_relation == NULL)
	//{
	//	exit(0);
	//}

	//fprintf(fp_bunny_relation, "%d\r\n ", ring_neighborhood_bunny.size());
	//for(int i=0 ; i < ring_neighborhood_bunny.size(); i++)
	//{
	//	vector<int> ring_neiborhood = ring_neighborhood_bunny[i];
	//	for(int j=0; j<ring_neiborhood.size(); j++)
	//	{
	//		fprintf(fp_bunny_relation, "%d ", ring_neiborhood[j]);
	//	}
	//	fprintf(fp_bunny_relation, "\r\n");
	//}

	//fclose(fp_bunny_relation);


	////read the relationship from txt
	//FILE *fp_bunny_relation = fopen("bunny_relation.txt", "rb");
	//if(fp_bunny_relation == NULL)
	//{
	//	exit(0);
	//}

	//
	//fprintf(fp_bunny_relation, "%d\r\n ", ring_neighborhood_bunny.size());
	//for(int i=0 ; i < ring_neighborhood_bunny.size(); i++)
	//{
	//	vector<int> ring_neiborhood = ring_neighborhood_bunny[i];
	//	for(int j=0; j<ring_neiborhood.size(); j++)
	//	{
	//		fprintf(fp_bunny_relation, "%d ", ring_neiborhood[j]);
	//	}
	//	fprintf(fp_bunny_relation, "\r\n");
	//}

	//fclose(fp_bunny_relation);


	for(int i =0; i < ring_neighborhood_bunny.size(); i++)
	{
		vector<int> ring_neighborhood = ring_neighborhood_bunny[i];
	//	int vertex = i;

		for(int j = 0; j < ring_neighborhood.size(); j++)
		{
			float cos_best = 0.0f;
			int v_best = ring_neighborhood[j];
			int v_i = ring_neighborhood[j];
			for(int k=0; k<ring_neighborhood.size(); k++)
			{
				int v_j = ring_neighborhood[k];
			//	if(v_i >= X_bunny.size() || v_j >= X_bunny.size() || i >= X_bunny.size())
			//	{
			//		continue;
			//	}
				glm::vec3 vi_v = X_bunny[v_i]-X_bunny[i];
				glm::vec3 vj_v = X_bunny[v_j]-X_bunny[i];
				float temp = glm::dot(vi_v,vj_v);
				float vi_v_length = glm::length(vi_v);
				float vj_v_length = glm::length(vj_v);

				float cos_temp = temp/(vi_v_length*vj_v_length);

				if(cos_temp <cos_best)
				{
					//printf("valid data");
					cos_best = cos_temp;
					v_best = ring_neighborhood[k];
				}

			}
			if(v_best != ring_neighborhood[j])
			{
				//invaild bend constraint
				//增加避免重复添加弯曲约束的代码实现
				sort_order.clear();
				sort_order.push_back(v_i);
				sort_order.push_back(v_best);
				sort_order.push_back(i);
				sort(sort_order.begin(), sort_order.end());

				vector<int>::iterator sort_iter;
				char temp[128];
				char result[128];
				memset(result , 0x00, 128);

				for(sort_iter = sort_order.begin() ; sort_iter != sort_order.end() ; sort_iter++)
				{

					_snprintf(temp, 128, "%d_", *sort_iter);
					strcat(result, temp);
				}

				string str = result;
				std::map<string, int>::iterator it;

				it = check_bend_constraint.find(str);
				if(it ==  check_bend_constraint.end()){
					check_bend_constraint[str] = 1;
					AddBendingConstraint(v_i, v_best, i, kBend);
				}else{

				}					
			}
		}
	}


	//计算封闭球的体积
	AddClothBalloons();
	printf("success!\r\n");
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

void OnRender() {		
	size_t i=0;
	float newTime = (float) glutGet(GLUT_ELAPSED_TIME);
	frameTime = newTime-currentTime;
	currentTime = newTime;

	QueryPerformanceCounter(&t2);
	// compute and print the elapsed time in millisec
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

	sprintf_s(info, "FPS: %3.2f, Frame time (GLUT): %3.4f msecs, Frame time (QP): %3.3f", fps, frameTime, frameTimeQP);
	glutSetWindowTitle(info);

	glClear(GL_COLOR_BUFFER_BIT| GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	//set viewing transformation
	glTranslatef(0,0,dist);
	glRotatef(rX,1,0,0);
	glRotatef(rY,0,1,0);

	glGetDoublev(GL_MODELVIEW_MATRIX, MV);
	viewDir.x = (float)-MV[2];
	viewDir.y = (float)-MV[6];
	viewDir.z = (float)-MV[10];
	Right = glm::cross(viewDir, Up);

	//draw grid 绘制地面上的方格线
	DrawGrid();

	//绘制兔子
	glColor3f(0,1,0);
	glPushMatrix();
	glMultMatrixf(glm::value_ptr(bunny_pos));  //????


	glBegin(GL_TRIANGLES);
	for(int index=0; index < triangle_bunny.size(); index++)
	{

		int ver1_index = triangle_bunny[index].ver1_index;
		int ver2_index = triangle_bunny[index].ver2_index;
		int ver3_index = triangle_bunny[index].ver3_index;

		if(ver1_index >= X_bunny.size() 
			|| ver2_index >= X_bunny.size()
			||ver3_index >= X_bunny.size() )
		{
			continue;
		}
		glVertex3f(X_bunny[ver1_index].x,X_bunny[ver1_index].y, X_bunny[ver1_index].z);
		glVertex3f(X_bunny[ver2_index].x,X_bunny[ver2_index].y, X_bunny[ver2_index].z);
		glVertex3f(X_bunny[ver3_index].x,X_bunny[ver3_index].y, X_bunny[ver3_index].z);

	}
	glEnd();
	glPopMatrix();
	glutSwapBuffers();
}

void OnShutdown() {	
	d_constraints_bunny.clear();
	b_constraints_bunny.clear();
	X_bunny.clear();
	F_bunny.clear();
	V_bunny.clear();
	W_bunny.clear();
	tmp_X_bunny.clear();
	Ri_bunny.clear();
}

void ComputeForces( ) {
	size_t i=0;
	for(i=1;i<total_points_bunny;i++) {
		F_bunny[i] = glm::vec3(0); 

		//add gravity force
		if(W_bunny[i]>0)		 
			F_bunny[i] += gravity ;  
	}
} 


void IntegrateExplicitWithDamping(float deltaTime) {
	float deltaTimeMass = deltaTime;
	size_t i=0;

	glm::vec3 Xcm = glm::vec3(0);
	glm::vec3 Vcm = glm::vec3(0);
	float sumM = 0;

	//Notice 这里从1开始，是因为bunny的数据集中第一个点的坐标特别的诡异
	for(i=1;i<total_points_bunny;i++) {

		V_bunny[i] *= global_dampening; 
		V_bunny[i] = V_bunny[i] + (F_bunny[i]*deltaTime)*W_bunny[i];  	 					
		Xcm += (X_bunny[i]*mass_bunny); 
		Vcm += (V_bunny[i]*mass_bunny);  
		sumM += mass_bunny;  
	} 
	Xcm /= sumM; //sumM = 1
	Vcm /= sumM; //sumM = 1 

	glm::mat3 I = glm::mat3(1);
	glm::mat3 Test = glm::mat3(1);
	glm::vec3 L = glm::vec3(0);
	glm::vec3 w = glm::vec3(0);//angular velocity 角速度

	//Notice 这里从1开始，是因为bunny的数据集中第一个点的坐标特别的诡异
	for(i=1;i<total_points_bunny;i++) { 
		Ri_bunny[i] = (X_bunny[i] - Xcm);	

		L += glm::cross(Ri_bunny[i],mass_bunny*V_bunny[i]); 
		glm::mat3 tmp = glm::mat3(0,-Ri_bunny[i].z,  Ri_bunny[i].y, 
			Ri_bunny[i].z,       0,-Ri_bunny[i].x,
			-Ri_bunny[i].y,Ri_bunny[i].x,        0);
		Test = (tmp*glm::transpose(tmp))*mass_bunny;  
		I +=(tmp*glm::transpose(tmp))*mass_bunny;  
	} 

	w = glm::inverse(I)*L;  

	//apply center of mass damping
	for(i=1;i<total_points_bunny;i++) {
		glm::vec3 delVi = Vcm + glm::cross(w,Ri_bunny[i])-V_bunny[i];		
		V_bunny[i] += kDamp*delVi;
	}

	//calculate predicted position
	for(i=1;i<total_points_bunny;i++) {
		if(W_bunny[i] <= 0.0) { 
			tmp_X_bunny[i] = X_bunny[i]; //fixed points 在这里表示没有
		} else {
			tmp_X_bunny[i] = X_bunny[i] + (V_bunny[i]*deltaTime);				 
		}
	} 
}

void Integrate(float deltaTime) {	
	float inv_dt = 1.0f/deltaTime;
	size_t i=0; 

	for(i=1;i<total_points_bunny;i++) {	
		V_bunny[i] = (tmp_X_bunny[i] - X_bunny[i])*inv_dt; 
		X_bunny[i] = tmp_X_bunny[i]; 

	}
}

void UpdateDistanceConstraint(int i) {
	DistanceConstraint c = d_constraints_bunny[i]; 
	glm::vec3 dir = tmp_X_bunny[c.p1] - tmp_X_bunny[c.p2];  

	float len = glm::length(dir); 
	if(len <= EPSILON) 
	{
		printf("1:exception");
		return;
	}
	float w1 = W_bunny[c.p1];  
	float w2 = W_bunny[c.p2];   
	float invMass = w1+ w2; 
	if(invMass <= EPSILON) 
	{
		printf("2:exception");
		return;
	}
	glm::vec3 dP = (1.0f/invMass) * (len-c.rest_length ) * (dir/len)* c.k_prime;
	if(len-c.rest_length != 0)
	{
		//	printf("dP is not zero\r\n");
	}
	if(w1 > 0.0)
		tmp_X_bunny[c.p1] -= dP*w1;
	if(w2 > 0.0)
		tmp_X_bunny[c.p2] += dP*w2;	
}

void UpdateBendingConstraint(int index) {
	size_t i=0;
	BendingConstraint c = b_constraints_bunny[index]; 

#ifdef USE_TRIANGLE_BENDING_CONSTRAINT
	float global_k =0.0f; 
	glm::vec3 center = 0.3333f * (tmp_X_bunny[c.p1] + tmp_X_bunny[c.p2] + tmp_X_bunny[c.p3]);
	glm::vec3 dir_center = tmp_X_bunny[c.p3]-center;
	float dist_center = glm::length(dir_center);

	float diff = 1.0f - ((global_k + c.rest_length) / dist_center);
	glm::vec3 dir_force = dir_center * diff;
	glm::vec3 fa = c.k_prime * ((2.0f*W_bunny[c.p1])/c.w) * dir_force;
	glm::vec3 fb = c.k_prime * ((2.0f*W_bunny[c.p2])/c.w) * dir_force;
	glm::vec3 fc = -c.k_prime * ((4.0f*W_bunny[c.p3])/c.w) * dir_force;

	if(W_bunny[c.p1] > 0.0)  {
		tmp_X_bunny[c.p1] += fa;
	}
	if(W_bunny[c.p2] > 0.0) {
		tmp_X_bunny[c.p2] += fb;
	}
	if(W_bunny[c.p3] > 0.0) {
		tmp_X_bunny[c.p3] += fc;
	}
#else
	float d = 0, phi=0, i_d=0 ;
	glm::vec3 n1=glm::vec3(0), n2=glm::vec3(0);

	glm::vec3 p1 = tmp_X_bunny[c.p1];
	glm::vec3 p2 = tmp_X_bunny[c.p2]-p1;
	glm::vec3 p3 = tmp_X_bunny[c.p3]-p1;

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
			return; 
	}

	i_d =sqrt(1-(d*d))*(phi-c.phi0) ;

	glm::vec3 q1 = (n2-n1*d)/lenp2;
	glm::vec3 q2 = (n1-n2*d)/lenp3;
	glm::vec3 q3 =-q1-q2;

	float q1_len2 = glm::dot(q1,q1);  
	float q2_len2 = glm::dot(q2,q2);
	float q3_len2 = glm::dot(q3,q3);

	float sum = W_bunny[c.p1]*(q1_len2) +W_bunny[c.p2]*(q2_len2) +W_bunny[c.p3]*(q3_len2);

	glm::vec3 dP1 = -( (W_bunny[c.p1] * i_d) /sum)*q1;
	glm::vec3 dP2 = -( (W_bunny[c.p2] * i_d) /sum)*q2;
	glm::vec3 dP3 = -( (W_bunny[c.p3] * i_d) /sum)*q3;

	if(W_bunny[c.p1] > 0.0) {
		tmp_X_bunny[c.p1] += dP1*c.k;
	}
	if(W_bunny[c.p2] > 0.0) {
		tmp_X_bunny[c.p2] += dP2*c.k;
	}
	if(W_bunny[c.p3] > 0.0) {
		tmp_X_bunny[c.p3] += dP3*c.k;
	}		
#endif
}
//----------------------------------------------------------------------------------------------------
void GroundCollision() //DevO: 24.07.2011
{
	for(size_t i=0;i<total_points_bunny;i++) {	
		if(tmp_X_bunny[i].y<0) //collision with ground
		{
			//	printf("GroundCollision\r\n");
			tmp_X_bunny[i].y=0;
		}
	}
}
//----------------------------------------------------------------------------------------------------
void UpdateInternalConstraints(float deltaTime) {
	size_t i=0;
	for (size_t si=0;si<solver_iterations;++si) {
		for(i=0;i<d_constraints_bunny.size();i++) {
			UpdateDistanceConstraint(i); 
		} 
		for(i=0;i<b_constraints_bunny.size();i++) {
			UpdateBendingConstraint(i);
		}
		UpdateClothBalloons();
	}
	
	
	GroundCollision();  
}

void OnIdle() {	
	if ( accumulator >= timeStep ) //时间步长，使得不是每一次调用idle函数都调用，保持在60的帧率
	{	 
		if(isPause == 0)
		{
			StepPhysics(timeStep );	 //设置物理约束	
			accumulator -= timeStep;

		}
	}

	glutPostRedisplay(); 
	Sleep(5); //TODO
}

void StepPhysics(float dt ) {

	ComputeForces();  
	IntegrateExplicitWithDamping(dt);  
	UpdateInternalConstraints(dt);	 
	Integrate(dt);
}

void processNormalKeys(unsigned char key, int x, int y)
{

	switch(key){
	case 27://ESC键

		break;
	case 32: //空格键,按下空格键的时候来回的切换
		isPause = !isPause;	

		break;
	}
}


void main(int argc, char** argv) {

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
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

