#include <cstdio>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include "voxel.h"
#include "particle.h"
using namespace std;

voxel::voxel() {}
voxel::voxel(double character_L, double box_l) {
	double a = character_L / 8.0;
	voxel_n = std::max(1, int(box_l / a));
	voxel_N = voxel_n * voxel_n;
	length = box_l / voxel_n;

}
voxel::~voxel() {}

int get_voxel_id(Point loc, voxel v) {
	int ix, iy, id;
	int n = v.voxel_n;

	ix = (int)(loc.x / v.length);
	iy = (int)(loc.y / v.length);
	if (ix > n - 1) ix = n - 1;
	if (iy > n - 1) iy = n - 1;
	if (ix < 0) ix = 0;
	if (iy < 0) iy = 0;

	id = ix + n * iy;
	return id;
}

void initial_vox(voxel* v) {
	for (int i = 0; i < v[0].voxel_N; i++) {
		int iy = int(double(i) / double(v[0].voxel_n));
		int ix = (i - iy * v[0].voxel_n);
		v[i].center.x = (double(ix) + 0.5) * v[i].length;
		v[i].center.y = (double(iy) + 0.5) * v[i].length;
	}
}

void scr_vox(double box_l, voxel* vox, int num, int cycle) {
	char cha[40];
	sprintf(cha, "packing_voxel_%d.scr", cycle);
	ofstream scr(cha);
	scr << fixed << setprecision(8);

	scr << "-osnap off" << endl;
	//scr << "erase all " << endl;
	scr << "vscurrent 2" << endl;
	scr << "color 2" << endl;
	for (int i = 0; i < num; ++i) {
		if (vox[i].is_occupy == false) {
			Point v1, v4;//voxel   ĸ    㣬   ĸ        particle  ʱ  voxel  particle  
			v1.x = vox[i].center.x - vox[i].length * 0.5;
			v1.y = vox[i].center.y - vox[i].length * 0.5;
			v4.x = vox[i].center.x + vox[i].length * 0.5;
			v4.y = vox[i].center.y + vox[i].length * 0.5;
			scr << "rectang " << v1.x << "," << v1.y << " " << v4.x << "," << v4.y << endl;
		}
	}
	scr << "rectang "
		<< "0"
		<< ","
		<< "0"
		<< " " << box_l << "," << box_l << endl;
	scr << "zoom e" << endl;
	scr << "vpoint 0,0,1" << endl;
	scr.close();
}

