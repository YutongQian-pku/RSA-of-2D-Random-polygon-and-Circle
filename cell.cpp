#include "cell.h"
#include "particle.h"
#include "voxel.h"
#include <algorithm>

using namespace std;
//**************cell��غ���**************************
Cell::Cell() { num = 0; }
Cell::Cell(double ccharacter_L, double L) {
	num = 0;
	cell_L = L;
	cellcharacter_L = ccharacter_L;
	cell_n = std::max(1, (int)(cell_L / cellcharacter_L));
	cell_l = cell_L / (double)(cell_n);
	cell_N = cell_n * cell_n;
}
Cell::~Cell() {}

void initial_grid(Cell* cell) {
	int N = cell[0].cell_N;
	int n = cell[0].cell_n;

	int x, y;
	int ix, iy, id, iid;

	for (iy = 0; iy < n; iy++) {
		for (ix = 0; ix < n; ix++) {

			iid = ix + n * iy;

			y = iy - 1; if (y < 0) { y = n - 1; cell[iid].boundary_tag = 4; cell[iid].tag_num++; }

			x = ix - 1; if (x < 0) { x = n - 1; cell[iid].boundary_tag = 2; cell[iid].tag_num++; }
			id = x + n * y;
			cell[iid].cell_id[0] = id;
			//***************************************
			x = ix;
			id = x + n * y;
			cell[iid].cell_id[1] = id;
			//***************************************
			x = ix + 1; if (x > n - 1) { x = 0; cell[iid].boundary_tag = 1; cell[iid].tag_num++; }
			id = x + n * y;
			cell[iid].cell_id[2] = id;
			//***************************************

			y = iy;

			x = ix - 1; if (x < 0) { x = n - 1; cell[iid].boundary_tag = 2; cell[iid].tag_num++; }
			id = x + n * y;
			cell[iid].cell_id[3] = id;
			//***************************************
			//x = ix;
			//id = x + n * y;
			//cell[iid].cell_id[4] = id;
			//***************************************
			x = ix + 1; if (x > n - 1) { x = 0; cell[iid].boundary_tag = 1; cell[iid].tag_num++; }
			id = x + n * y;
			cell[iid].cell_id[4] = id;
			//***************************************

			y = iy + 1; if (y > n - 1) { y = 0; cell[iid].boundary_tag = 3; cell[iid].tag_num++; }

			x = ix - 1; if (x < 0) { x = n - 1; cell[iid].boundary_tag = 2; cell[iid].tag_num++; }
			id = x + n * y;
			cell[iid].cell_id[5] = id;
			//***************************************
			x = ix;
			id = x + n * y;
			cell[iid].cell_id[6] = id;
			//***************************************
			x = ix + 1; if (x > n - 1) { x = 0; cell[iid].boundary_tag = 1; cell[iid].tag_num++; }
			id = x + n * y;
			cell[iid].cell_id[7] = id;
			//***************************************
		}
	}
}

int get_cell_id(Cell &cell, double x, double y) {
	int ix, iy, id;
	double l = cell.cell_l;
	int n = cell.cell_n;

	ix = (int)(x / l);
	iy = (int)(y / l);
	if (ix > n - 1) ix = n - 1;
	if (iy > n - 1) iy = n - 1;
	if (ix < 0) ix = 0;
	if (iy < 0) iy = 0;

	id = ix + n * iy;
	return id;
}

bool check_over(Cell* cell, Particle* p, voxel* vox) {
	if (vox[get_voxel_id(p->center, vox[0])].is_occupy == true)return true;//flagΪtrueʱ��Ϊ�ཻ���Σ����ܷ���
	int id;
	double box_l = cell[0].cell_L;
	
	//����ԭ����ʵ�����������д����У��ô�����ͼ��cell[0]��������Ǵ���ָ�룬��͵����˸��Ʋ�������cell������particleҲ�������ˣ�Υ����unique��ԭ��
	id = get_cell_id(cell[0], p->center.x, p->center.y);

	for (const auto& particlePtr : cell[id].particles) {
		auto p1 = particlePtr->clone();
		if (p->intersects(p1.get())){
			return true;
		}
	}

	//���ֱ߽�λ�� //0���ڲ���1��x+��2��x-��3��y+��4��y-
	int x_min[3] = { 0,3,5 };
	int x_max[3] = { 2,4,7 };
	int y_min[3] = { 0,1,2 };
	int y_max[3] = { 5,6,7 };
	int a = cell[id].boundary_tag;
	int b = cell[id].tag_num;
	if (a == 0) {
		for (int k = 0; k < 8; k++) {
			for (const auto& particlePtr : cell[cell[id].cell_id[k]].particles) {
				auto p1 = particlePtr->clone();
				if (p->intersects(p1.get())) {
					return true;
				}
			}
		}
	}
	else if (a == 1 && b == 3) {
		for (int i = 0; i < 3; i++) {
			for (const auto& particlePtr : cell[cell[id].cell_id[x_max[i]]].particles) {
				auto p1 = particlePtr->clone();
				p1->translate(1.0f * box_l, 0);				
				if (p->intersects(p1.get())) return true;
			}
		}
		for (int k = 0; k < 8; k++) {
			for (const auto& particlePtr : cell[cell[id].cell_id[k]].particles) {
				auto p1 = particlePtr->clone();
				if (p->intersects(p1.get())) return true;
			}
		}
	}
	else if (a == 2 && b == 3) {
		for (int i = 0; i < 3; i++) {
			for (const auto& particlePtr : cell[cell[id].cell_id[x_min[i]]].particles) {
				auto p1 = particlePtr->clone();
				p1->translate(-1.0 * box_l, 0);
				if (p->intersects(p1.get()))return true;
			}
		}
		for (int k = 0; k < 8; k++) {
			for (const auto& particlePtr : cell[cell[id].cell_id[k]].particles) {
				auto p1 = particlePtr->clone();
				if (p->intersects(p1.get())) return true;
			}
		}
	}
	else if (a == 3 && b == 1) {
		for (int i = 0; i < 3; i++) {
			for (const auto& particlePtr : cell[cell[id].cell_id[y_max[i]]].particles) {
				auto p1 = particlePtr->clone();
				p1->translate(0, 1.0 * box_l);
				if (p->intersects(p1.get())) return true;
			}
		}
		for (int k = 0; k < 8; k++) {
			for (const auto& particlePtr : cell[cell[id].cell_id[k]].particles) {
				auto p1 = particlePtr->clone();
				if (p->intersects(p1.get())) return true;
			}
		}
	}
	else if (a == 4 && b == 1) {
		for (int i = 0; i < 3; i++) {
			for (const auto& particlePtr : cell[cell[id].cell_id[y_min[i]]].particles) {
				auto p1 = particlePtr->clone();
				p1->translate(0, -1.0 * box_l);
				if (p->intersects(p1.get())) return true;
			}
		}
		for (int k = 0; k < 8; k++) {
			for (const auto& particlePtr : cell[cell[id].cell_id[k]].particles) {
				auto p1 = particlePtr->clone();
				if (p->intersects(p1.get())) return true;
			}
		}
	}
	else {
		for (int k = 0; k < 8; ++k) {
			for (int m = -1; m < 2; ++m) {
				for (int n = -1; n < 2; ++n) {
					for (const auto& particlePtr : cell[cell[id].cell_id[k]].particles) {
						auto p1 = particlePtr->clone();
						p1->translate(m * box_l, n * box_l);
						if (p->intersects(p1.get())) return true;
					}
				}
			}
		}
	}
	return false;//���ཻ
}