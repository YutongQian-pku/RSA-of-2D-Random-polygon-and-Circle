#ifndef _VOXEL_H_
#define _VOXEL_H_

#include"point.h"

using namespace std;

class Point;

class voxel {
public:
	voxel(double character_L, double box_l);//初始化整体属性
	voxel();
	~voxel();
	Point center;
	Point v1, v2, v3, v4;
	bool is_occupy = false;//默认为未占用，false
	int voxel_n;
	long int voxel_N;
	double length;
};

int get_voxel_id(Point loc, voxel v);
void initial_vox(voxel* v);
void scr_vox(double box_l, voxel* vox, int num, int cycle);

#endif