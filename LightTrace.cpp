#include<iostream>
#include<stdlib.h>
#include <graphics.h>
#include<math.h>
using namespace std;
#define inf 999.0//可视最大距离


#define width_c 800//画布像素宽度
#define high_c 600///画布像素高度


#define PIE 3.1415926//Π
COLORREF BG_COLOR=RGB(7,25,46);//背景颜色
double Z_C = 1;//画布的Z轴,(视角范围，数值越小范围越大）
double width_v = 1.6;//画布现实宽度
double high_v = 1.2;//画布现实高度
int R, G, B;//RGB颜色参数
//点
struct DOT {
	double x;
	double y;
	double z;
};
//球体
struct SPHERE {
	DOT center;//球心
	double radial;//半径
	COLORREF color;//颜色
	int specular;//光滑程度，-1完全粗糙
	float reflective;//反射强度0-1
};
//向量
struct VECTOR {
	double x;
	double y;
	double z;
};
//方向光
struct DIRECTLIGHT {
	double intensify;
	VECTOR direction;
};
//全局光照
struct AMBIENTLIGHT {
	double intensify;
};
//点光源
struct POINTLIGHT {
	double intensify;
	DOT position;
};
//画布点的集合
DOT canvas[high_c][width_c];
//球的集合和数量
SPHERE spheres[1000];
int s_i;//反射强度
int num_sph;
//光源的集合和数量
int pointnum;
int dirnum;
DIRECTLIGHT dirlights[10];
POINTLIGHT pointlights[10];
AMBIENTLIGHT ambientlight;
//交点数组
double t[2];
//函数的声明
double computinglight(DOT P, VECTOR N, VECTOR V, int s);
DOT canvastoviewpoint(double x, double y);
void intersectsphere(DOT O, DOT D, SPHERE sphere);
COLORREF traceray(DOT O, DOT D, double t_min, double t_max,int recursion_depth);
double dot(VECTOR a, VECTOR b);
double length(VECTOR a);
double calcolor(BYTE c, double i);
inline double transx(double a) {//坐标转换
	return a * (width_v / width_c);
}
inline double transy(double a) {
	return a * (high_v / high_c);
}
inline VECTOR reflectray(VECTOR R, VECTOR N) {//计算反射光线方向向量，R是入射光线的反向
	VECTOR Ref = { N.x * 2 * dot(N,R) - R.x,N.y * 2 * dot(N,R) - R.y,N.z * 2 * dot(N,R) - R.z };
	return Ref;
}
inline DOT canvastoviewpoint(double x, double y) {//将像素坐标xy转化为世界坐标
	DOT OV = { transx(x),transy(y),Z_C };
	return OV;
}
inline double calcolor(BYTE c, double i) {//光照参数范围限定
	double temp = c * i;
	if (temp > 255)
		temp = 255;
	if (temp < 0)
		temp = 0;
	return temp;
}
inline double dot(VECTOR a, VECTOR b) {//点乘
	return (a.x * b.x + a.y * b.y + a.z * b.z); /*
	(sqrt(a.x * a.x + a.y * a.y + a.z * a.z) * sqrt(b.x * b.x + b.y * b.y + b.z * b.z));
	*/
}
int main()
{
	DOT O = { 0,-1,0 };	//设置相机位置
	//DOT C = {0,0,0 };//相机位置
	VECTOR R = {0,0,0 };//沿xyz轴旋转弧度,相机朝向
	//反射递归次数
	int recursion_depth=2;
	//创建光源
	s_i = 2;
	dirnum = 0;
	pointnum = 1;
	//方向光
	dirlights[0] = { 0.5,{-1,-1,-1} };
	//点光源
	pointlights[0] = { 0.8,{5,5,4} };
	//全局光照
	ambientlight = { 0.2 };
	//创建球

	num_sph = 4;
	int count = 1;
	spheres[0] = { {3,5,5 }, 0.5, RGB(255, 1, 1),100,0.3};
	spheres[1] = { {0,-5002,0 }, 5000, RGB(255, 255, 255),-1,0.2};
	spheres[2] = { {-1,0,4 }, 0.5, RGB(1, 255, 1),100,1 };
	spheres[3] = { {-5002,0,0 }, 5000, RGB(200, 200, 200),50,0.6};
	


	int x, y;	//创建坐标变量
	//创建画布并转化为标准直角坐标系
	initgraph(width_c, high_c);
	setorigin(width_c / 2, high_c / 2);
	setaspectratio(1, -1);
	int i = 0;
	int j = 0;
	//定义画布的点
	for (y = high_c / 2;y > -high_c / 2;y--) {
		j = 0;
		for (x = -width_c / 2;x < width_c / 2;x++) {
			canvas[i][j].x = x;
			canvas[i][j].y = y;
			canvas[i][j].z = Z_C;
			j++;
		}
		i++;
	}
	DOT D;//画布上的点
	COLORREF color;//颜色
	struct pixel {
		int x;
		int y;
		COLORREF color;
	};
	double ST;
	double i1=1;
	int R1, G1, B1;
		for (int i = 0;i < high_c;i++) {
			for (int j = 0;j < width_c;j++) {
				D = canvastoviewpoint(canvas[i][j].x, canvas[i][j].y);
				D = { O.x + (D.x * cos(R.y) * cos(R.z) + sin(R.y) * cos(R.z) * D.y - sin(R.z) * D.z),
					O.y + (D.x * (sin(R.x) * sin(R.y) * sin(R.z) - sin(R.y) * cos(R.x)) + D.y * (sin(R.x) * sin(R.y) * sin(R.z) + cos(R.x) * cos(R.y)) + D.z * sin(R.x) * cos(R.z)),
						O.z + (D.x * (sin(R.z) * cos(R.x) * cos(R.y) + sin(R.x) * sin(R.y)) + D.y * (sin(R.z) * cos(R.x) * sin(R.y) - sin(R.x) * cos(R.y)) + D.z * cos(R.x) * cos(R.z)) };
				color = traceray(O, D, 1, inf, recursion_depth);
				//color = RGB(R1, G1, B1);
				pixel p = { canvas[i][j].x, canvas[i][j].y, color };
				putpixel(canvas[i][j].x, canvas[i][j].y, color);
			}
	}
		system("pause");
		closegraph();
		return 0;
}

void closestintersection(SPHERE** sphere_closest, double* t_closest, DOT O, DOT D, double t_min, double t_max) {
	//
	*t_closest = inf;
	int i, j;
	*sphere_closest = NULL;
	//这个循环遍历所有的球体
	for (i = 0;i < num_sph;i++) {
		//计算交点的参数t，存入数组
		intersectsphere(O, D, spheres[i]);
		//这个循环用于找小的t
		for (j = 0;j < 2;j++) {
			if (t[j] > t_min && t[j] < t_max && t[j] < *t_closest) {
				//记录下来前面的交点
				*t_closest = t[j];
				//记录下来前面的球体
				*sphere_closest = &spheres[i];
			}
		}
	}
}

COLORREF traceray(DOT O, DOT D, double t_min, double t_max,int recursion_depth) {//计算O点和D点射线上是否有物体，并返回物体的颜色
	//int i,j;
	COLORREF local_color;
	double t_closest=inf;
	SPHERE* sphere_closest;
	//没有球体观察到的是背景
	closestintersection(&sphere_closest, &t_closest, O, D, t_min, t_max);

	if (sphere_closest == NULL)
		return BG_COLOR;
	//局部颜色计算
	VECTOR DO = { O.x - D.x,O.y - D.y,O.z - D.z };//相机到画布一点D的向量的相反向量
	DOT P = { O.x + t_closest * -DO.x,O.y + t_closest * -DO.y,O.z + t_closest * -DO.z };	//计算交点P
	VECTOR N = { sphere_closest->center.x - P.x,
		sphere_closest->center.y - P.y  ,
		sphere_closest->center.z - P.z };	//计算P处的单位法向量N
	N = { N.x / length(N),N.y / length(N) ,N.z / length(N) };
	//提取RGB在进行颜色计算，防止超出RGB颜色定义域
	R = GetRValue(sphere_closest->color);
	G=GetGValue(sphere_closest->color);
	B = GetBValue(sphere_closest->color);
	//返回前面的球体的颜色
	local_color= RGB(calcolor(R, computinglight(P, N, DO, sphere_closest->specular)),
		calcolor(G, computinglight(P, N, DO, sphere_closest->specular)),
		calcolor(B, computinglight(P, N, DO, sphere_closest->specular)));
	/*return RGB(calcolor(R, computinglight(P, N,DO,sphere_closest->specular)), 
		calcolor(G, computinglight(P, N, DO, sphere_closest->specular)), 
		calcolor(B, computinglight(P, N, DO, sphere_closest->specular)));*/
	//return sphere_closest->color*computinglight(P,N);
	double r = sphere_closest->reflective;
	if (recursion_depth <= 0 || r <= 0) {
		return local_color;
	}
	VECTOR Ref;//反射光线
	Ref=reflectray(DO, N);
	COLORREF Ref_color;
	DOT R1 = { P.x + Ref.x, P.y + Ref.y, P.z + Ref.z };//反射光线上的一点
	Ref_color = traceray(P, R1, 0.0001, inf, recursion_depth - 1);
	int r1 = GetRValue(Ref_color);
	R = GetRValue(local_color);
	G = GetGValue(local_color);
	B = GetBValue(local_color);
	return RGB(R * (1 - r) + GetRValue(Ref_color) * r, G * (1 - r) + GetGValue(Ref_color) * r, B * (1 - r) + GetBValue(Ref_color) * r);
}

void intersectsphere(DOT O, DOT D, SPHERE sphere) {//计算OD射线上是否有物体，返回tOD中的t
	double r = sphere.radial;//代表圆的半径
	//根据公式的需要求的两个向量
	VECTOR CO = { O.x - sphere.center.x,O.y - sphere.center.y,O.z - sphere.center.z };
	VECTOR OD = { D.x - O.x,D.y - O.y,D.z - O.z };
	//计算系数
	double a = dot(OD, OD);
	double b = 2 * dot(CO, OD);
	double c = dot(CO, CO) - r * r;

	double discriminant = b * b - 4 * a * c;//解的判别式
	if (discriminant < 0) {//没有交点
		t[0] = inf;
		t[1] = inf;
	}
	else {//有交点
		t[0] = (-b + sqrt(discriminant)) / (2 * a);
		t[1] = (-b - sqrt(discriminant)) / (2 * a);
	}
}

double length(VECTOR a) {//求向量长度
	return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

double computinglight(DOT P, VECTOR N,VECTOR V,int s) {//计算点P的光照强度并返回
	double i = 0;//光照强度
	double t_max;
	DOT L1;
	VECTOR L;//代表光的方向向量
	double n_dot_l;//方向向量和法向量的点乘
	VECTOR R;//镜面反射光线
	float r_dot_v;//相机与画布上的点与反射光线的夹角的余弦
	i += ambientlight.intensify;//全局光照直接加
	SPHERE* sphere_shadow;
	double t_shadow;
	int j;
	for (j = 0;j < pointnum;j++) {//这个循环遍历点光源
		t_max = 1;
		L = { P.x - pointlights[j].position.x , P.y - pointlights[j].position.y, P.z - pointlights[j].position.z };
		L1 = { pointlights[j].position.x, pointlights[j].position.y, pointlights[j].position.z };
		closestintersection(&sphere_shadow, &t_shadow, P, L1, 0.0001, t_max);
		//阴影检测
		if (sphere_shadow != NULL) {
			continue;
		}
		n_dot_l = dot(N, L);
		if (n_dot_l > 0) {//正面才加上光照强度
			i += pointlights[j].intensify * n_dot_l / (length(N) * length(L));
		}

		//镜面反射
		if (s != -1) {
			L = { -L.x,-L.y, -L.z };
			R = { N.x * 2 * dot(N,L) - L.x,N.y * 2 * dot(N,L) - L.y,N.z * 2 * dot(N,L) - L.z };
			r_dot_v = dot(R, V);
			if (r_dot_v > 0) {
				i += pointlights[j].intensify * pow(r_dot_v / (length(R) * length(V)), s)*s_i;
			}
		}
	}
	for (j = 0;j < dirnum;j++) {//这个循环遍历方向光
		t_max = inf;
		L = { dirlights[j].direction.x,dirlights[j].direction.y,dirlights[j].direction.z };
		L1 = { P.x-L.x,P.y-L.y,P.z-L.z };
		closestintersection(&sphere_shadow, &t_shadow, P, L1, 0.001, t_max);
		if (sphere_shadow != NULL) {
			continue;
		}
		n_dot_l = dot(N, L);
		if (n_dot_l > 0) {//正面才加上光照强度
			i += dirlights[j].intensify * n_dot_l / (length(N) * length(L));
		}
		//镜面反射
		if (s != -1) {
			L = { -L.x,-L.y, -L.z };
			R = { N.x * 2 * dot(N,L) - L.x,N.y * 2 * dot(N,L) - L.y,N.z * 2 * dot(N,L) - L.z };
			r_dot_v = dot(R, V);
			if (r_dot_v > 0) {
				i += dirlights[j].intensify * pow(r_dot_v / (length(R) * length(V)), s)*s_i;
			}
		}
	}
	return i;
}
