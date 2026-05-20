#ifndef _CELL_H_
#define _CELL_H_

#include <vector>
#include <memory>
using namespace std;
#include"particle.h"

class Cell {
public:
	Cell(double ccharacter_L, double box_l);
	//��һ���������������ȣ��ڶ���������packing�ı߳�
	Cell();
	~Cell();

	//cell parameter
	double cell_L;//cell���ܱ߳�����packing�ռ�ı߳�
	double cell_l;//cell�ı߳�
	int cell_n;//һ������cell�ĸ���
	int cell_N;//�ܵ�cell����
	double cellcharacter_L;

	//for search
	int num; //cell�ڵĿ�����
	vector<unique_ptr<Particle>> particles;

	int cell_id[8];

	int boundary_tag = 0;//0���ڲ���1��x+��2��x-��3��y+��4��y-
	int tag_num = 0;//����tag��ǩ�Ĵ���

	Cell& operator=(Cell&& other)noexcept {
		if (this != &other) {
			particles = std::move(other.particles);
			cell_L = other.cell_L;
			cellcharacter_L = other.cellcharacter_L;
			for (int i = 0; i < 8; i++) {
				cell_id[i] = other.cell_id[i];
			}
			cell_l = other.cell_l;
			cell_N = other.cell_N;
			cell_n = other.cell_n;
			num = other.num;
			tag_num = other.tag_num;
			boundary_tag = other.boundary_tag;
		}
		return *this;
	}
};

void initial_grid(Cell* cell);
int get_cell_id(Cell &cell, double x, double y);
bool check_over(Cell* cell, Particle* originalParticle, voxel* vox);


#endif