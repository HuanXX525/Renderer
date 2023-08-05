#include<iostream>
#include<stdlib.h>
#include <graphics.h>
#include<math.h>
using namespace std;
#define inf 999.0//����������


#define width_c 800//�������ؿ��
#define high_c 600///�������ظ߶�


#define PIE 3.1415926//��
COLORREF BG_COLOR=RGB(7,25,46);//������ɫ
double Z_C = 1;//������Z��,(�ӽǷ�Χ����ֵԽС��ΧԽ��
double width_v = 1.6;//������ʵ���
double high_v = 1.2;//������ʵ�߶�
int R, G, B;//RGB��ɫ����
//��
struct DOT {
	double x;
	double y;
	double z;
};
//����
struct SPHERE {
	DOT center;//����
	double radial;//�뾶
	COLORREF color;//��ɫ
	int specular;//�⻬�̶ȣ�-1��ȫ�ֲ�
	float reflective;//����ǿ��0-1
};
//����
struct VECTOR {
	double x;
	double y;
	double z;
};
//�����
struct DIRECTLIGHT {
	double intensify;
	VECTOR direction;
};
//ȫ�ֹ���
struct AMBIENTLIGHT {
	double intensify;
};
//���Դ
struct POINTLIGHT {
	double intensify;
	DOT position;
};
//������ļ���
DOT canvas[high_c][width_c];
//��ļ��Ϻ�����
SPHERE spheres[1000];
int s_i;//����ǿ��
int num_sph;
//��Դ�ļ��Ϻ�����
int pointnum;
int dirnum;
DIRECTLIGHT dirlights[10];
POINTLIGHT pointlights[10];
AMBIENTLIGHT ambientlight;
//��������
double t[2];
//����������
double computinglight(DOT P, VECTOR N, VECTOR V, int s);
DOT canvastoviewpoint(double x, double y);
void intersectsphere(DOT O, DOT D, SPHERE sphere);
COLORREF traceray(DOT O, DOT D, double t_min, double t_max,int recursion_depth);
double dot(VECTOR a, VECTOR b);
double length(VECTOR a);
double calcolor(BYTE c, double i);
inline double transx(double a) {//����ת��
	return a * (width_v / width_c);
}
inline double transy(double a) {
	return a * (high_v / high_c);
}
inline VECTOR reflectray(VECTOR R, VECTOR N) {//���㷴����߷���������R��������ߵķ���
	VECTOR Ref = { N.x * 2 * dot(N,R) - R.x,N.y * 2 * dot(N,R) - R.y,N.z * 2 * dot(N,R) - R.z };
	return Ref;
}
inline DOT canvastoviewpoint(double x, double y) {//����������xyת��Ϊ��������
	DOT OV = { transx(x),transy(y),Z_C };
	return OV;
}
inline double calcolor(BYTE c, double i) {//���ղ�����Χ�޶�
	double temp = c * i;
	if (temp > 255)
		temp = 255;
	if (temp < 0)
		temp = 0;
	return temp;
}
inline double dot(VECTOR a, VECTOR b) {//���
	return (a.x * b.x + a.y * b.y + a.z * b.z); /*
	(sqrt(a.x * a.x + a.y * a.y + a.z * a.z) * sqrt(b.x * b.x + b.y * b.y + b.z * b.z));
	*/
}
int main()
{
	DOT O = { 0,-1,0 };	//�������λ��
	//DOT C = {0,0,0 };//���λ��
	VECTOR R = {0,0,0 };//��xyz����ת����,�������
	//����ݹ����
	int recursion_depth=2;
	//������Դ
	s_i = 2;
	dirnum = 0;
	pointnum = 1;
	//�����
	dirlights[0] = { 0.5,{-1,-1,-1} };
	//���Դ
	pointlights[0] = { 0.8,{5,5,4} };
	//ȫ�ֹ���
	ambientlight = { 0.2 };
	//������

	num_sph = 4;
	int count = 1;
	spheres[0] = { {3,5,5 }, 0.5, RGB(255, 1, 1),100,0.3};
	spheres[1] = { {0,-5002,0 }, 5000, RGB(255, 255, 255),-1,0.2};
	spheres[2] = { {-1,0,4 }, 0.5, RGB(1, 255, 1),100,1 };
	spheres[3] = { {-5002,0,0 }, 5000, RGB(200, 200, 200),50,0.6};
	


	int x, y;	//�����������
	//����������ת��Ϊ��׼ֱ������ϵ
	initgraph(width_c, high_c);
	setorigin(width_c / 2, high_c / 2);
	setaspectratio(1, -1);
	int i = 0;
	int j = 0;
	//���廭���ĵ�
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
	DOT D;//�����ϵĵ�
	COLORREF color;//��ɫ
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
	//���ѭ���������е�����
	for (i = 0;i < num_sph;i++) {
		//���㽻��Ĳ���t����������
		intersectsphere(O, D, spheres[i]);
		//���ѭ��������С��t
		for (j = 0;j < 2;j++) {
			if (t[j] > t_min && t[j] < t_max && t[j] < *t_closest) {
				//��¼����ǰ��Ľ���
				*t_closest = t[j];
				//��¼����ǰ�������
				*sphere_closest = &spheres[i];
			}
		}
	}
}

COLORREF traceray(DOT O, DOT D, double t_min, double t_max,int recursion_depth) {//����O���D���������Ƿ������壬�������������ɫ
	//int i,j;
	COLORREF local_color;
	double t_closest=inf;
	SPHERE* sphere_closest;
	//û������۲쵽���Ǳ���
	closestintersection(&sphere_closest, &t_closest, O, D, t_min, t_max);

	if (sphere_closest == NULL)
		return BG_COLOR;
	//�ֲ���ɫ����
	VECTOR DO = { O.x - D.x,O.y - D.y,O.z - D.z };//���������һ��D���������෴����
	DOT P = { O.x + t_closest * -DO.x,O.y + t_closest * -DO.y,O.z + t_closest * -DO.z };	//���㽻��P
	VECTOR N = { sphere_closest->center.x - P.x,
		sphere_closest->center.y - P.y  ,
		sphere_closest->center.z - P.z };	//����P���ĵ�λ������N
	N = { N.x / length(N),N.y / length(N) ,N.z / length(N) };
	//��ȡRGB�ڽ�����ɫ���㣬��ֹ����RGB��ɫ������
	R = GetRValue(sphere_closest->color);
	G=GetGValue(sphere_closest->color);
	B = GetBValue(sphere_closest->color);
	//����ǰ����������ɫ
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
	VECTOR Ref;//�������
	Ref=reflectray(DO, N);
	COLORREF Ref_color;
	DOT R1 = { P.x + Ref.x, P.y + Ref.y, P.z + Ref.z };//��������ϵ�һ��
	Ref_color = traceray(P, R1, 0.0001, inf, recursion_depth - 1);
	int r1 = GetRValue(Ref_color);
	R = GetRValue(local_color);
	G = GetGValue(local_color);
	B = GetBValue(local_color);
	return RGB(R * (1 - r) + GetRValue(Ref_color) * r, G * (1 - r) + GetGValue(Ref_color) * r, B * (1 - r) + GetBValue(Ref_color) * r);
}

void intersectsphere(DOT O, DOT D, SPHERE sphere) {//����OD�������Ƿ������壬����tOD�е�t
	double r = sphere.radial;//����Բ�İ뾶
	//���ݹ�ʽ����Ҫ�����������
	VECTOR CO = { O.x - sphere.center.x,O.y - sphere.center.y,O.z - sphere.center.z };
	VECTOR OD = { D.x - O.x,D.y - O.y,D.z - O.z };
	//����ϵ��
	double a = dot(OD, OD);
	double b = 2 * dot(CO, OD);
	double c = dot(CO, CO) - r * r;

	double discriminant = b * b - 4 * a * c;//����б�ʽ
	if (discriminant < 0) {//û�н���
		t[0] = inf;
		t[1] = inf;
	}
	else {//�н���
		t[0] = (-b + sqrt(discriminant)) / (2 * a);
		t[1] = (-b - sqrt(discriminant)) / (2 * a);
	}
}

double length(VECTOR a) {//����������
	return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

double computinglight(DOT P, VECTOR N,VECTOR V,int s) {//�����P�Ĺ���ǿ�Ȳ�����
	double i = 0;//����ǿ��
	double t_max;
	DOT L1;
	VECTOR L;//�����ķ�������
	double n_dot_l;//���������ͷ������ĵ��
	VECTOR R;//���淴�����
	float r_dot_v;//����뻭���ϵĵ��뷴����ߵļнǵ�����
	i += ambientlight.intensify;//ȫ�ֹ���ֱ�Ӽ�
	SPHERE* sphere_shadow;
	double t_shadow;
	int j;
	for (j = 0;j < pointnum;j++) {//���ѭ���������Դ
		t_max = 1;
		L = { P.x - pointlights[j].position.x , P.y - pointlights[j].position.y, P.z - pointlights[j].position.z };
		L1 = { pointlights[j].position.x, pointlights[j].position.y, pointlights[j].position.z };
		closestintersection(&sphere_shadow, &t_shadow, P, L1, 0.0001, t_max);
		//��Ӱ���
		if (sphere_shadow != NULL) {
			continue;
		}
		n_dot_l = dot(N, L);
		if (n_dot_l > 0) {//����ż��Ϲ���ǿ��
			i += pointlights[j].intensify * n_dot_l / (length(N) * length(L));
		}

		//���淴��
		if (s != -1) {
			L = { -L.x,-L.y, -L.z };
			R = { N.x * 2 * dot(N,L) - L.x,N.y * 2 * dot(N,L) - L.y,N.z * 2 * dot(N,L) - L.z };
			r_dot_v = dot(R, V);
			if (r_dot_v > 0) {
				i += pointlights[j].intensify * pow(r_dot_v / (length(R) * length(V)), s)*s_i;
			}
		}
	}
	for (j = 0;j < dirnum;j++) {//���ѭ�����������
		t_max = inf;
		L = { dirlights[j].direction.x,dirlights[j].direction.y,dirlights[j].direction.z };
		L1 = { P.x-L.x,P.y-L.y,P.z-L.z };
		closestintersection(&sphere_shadow, &t_shadow, P, L1, 0.001, t_max);
		if (sphere_shadow != NULL) {
			continue;
		}
		n_dot_l = dot(N, L);
		if (n_dot_l > 0) {//����ż��Ϲ���ǿ��
			i += dirlights[j].intensify * n_dot_l / (length(N) * length(L));
		}
		//���淴��
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
