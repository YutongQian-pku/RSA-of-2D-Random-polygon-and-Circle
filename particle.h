#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include <cmath>
#include <memory>
#include <vector>
#include "voxel.h"
#include "point.h"

#define PI 3.1415926535898

class voxel;

class Segment {
public:
    Point start;
    Point end;
    Point normal;

    Segment(const Point& s, const Point& e) : start(s), end(e) {
        Point edge = end - start;
        normal.x = -edge.y;
        normal.y = edge.x;
        double length = std::sqrt(normal.x * normal.x + normal.y * normal.y);
        if (length > 0.0) {
            normal.x /= length;
            normal.y /= length;
        }
    }
};

class Particle {
public:
    Point center;
    double character_L = 0.0;
    double R = 0.0;       // bounding/circum radius about center
    double inner_r = 0.0; // inscribed radius about center

    virtual void translate(double x, double y) = 0;
    virtual bool intersects(Particle* p) const = 0;
    virtual std::unique_ptr<Particle> clone() const = 0;
    virtual bool is_onboundary(double box_l) const = 0;
    virtual bool point_is_in(Point v) const = 0;

    virtual void voxel_is_in(voxel* v) const {
        Point v1, v2, v3, v4;
        v1.x = v->center.x - v->length * 0.5;
        v1.y = v->center.y - v->length * 0.5;
        v2.x = v->center.x + v->length * 0.5;
        v2.y = v->center.y - v->length * 0.5;
        v3.x = v->center.x - v->length * 0.5;
        v3.y = v->center.y + v->length * 0.5;
        v4.x = v->center.x + v->length * 0.5;
        v4.y = v->center.y + v->length * 0.5;
        int a = point_is_in(v1) + point_is_in(v2) + point_is_in(v3) + point_is_in(v4);
        if (a == 4) v->is_occupy = true;
    }

    virtual void cal_occupy_single(voxel* v, double /*box_l*/) const {
        double up = center.y + R;
        double down = center.y - R;
        double left = center.x - R;
        double right = center.x + R;

        Point p1(left, down), p2(right, down), p3(left, up), p4(right, up);
        int id1 = get_voxel_id(p1, v[0]);
        int id2 = get_voxel_id(p2, v[0]);
        int id3 = get_voxel_id(p3, v[0]);
        int row = (id3 - id1) / v[0].voxel_n + 1;
        int col = id2 - id1 + 1;
        if (row < 1) row = 1;
        if (col < 1) col = 1;

        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                int iid = id1 + j + v[0].voxel_n * i;
                if (iid >= 0 && iid < v[0].voxel_N) voxel_is_in(&v[iid]);
            }
        }
    }

    virtual void cal_occupy(voxel* v, double box_l) {
        if (is_onboundary(box_l)) {
            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    translate(i * box_l, j * box_l);
                    cal_occupy_single(&v[0], box_l);
                    translate(-1.0 * i * box_l, -1.0 * j * box_l);
                }
            }
        } else {
            cal_occupy_single(&v[0], box_l);
        }
    }

    virtual std::unique_ptr<Particle> enlarge(double add) const = 0;
    virtual ~Particle() {}
};

class Circle : public Particle {
public:
    Circle();
    Circle(double R, Point ccenter);
    ~Circle();

    void translate(double x, double y) override;
    bool intersects(Particle* p) const override;
    std::unique_ptr<Particle> clone() const override {
        return std::make_unique<Circle>(*this);
    }
    bool is_onboundary(double box_l) const override;
    bool point_is_in(Point v) const override;
    std::unique_ptr<Particle> enlarge(double add) const override;
};

class Polygon : public Particle {
public:
    int shape_id = 0;
    double theta = 0.0;
    double buffer_radius = 0.0;        // used only for voxel pruning: exact disk buffer around polygon
    std::vector<Point> local_vertices; // centroid-relative, absolute-size template vertices

    Polygon();
    Polygon(const std::vector<Point>& localVertices,
            double angle,
            Point ccenter,
            int shapeId = 0,
            double bufferRadius = 0.0);
    ~Polygon();

    void translate(double x, double y) override;
    bool intersects(Particle* p) const override;
    std::unique_ptr<Particle> clone() const override {
        return std::make_unique<Polygon>(*this);
    }
    bool is_onboundary(double box_l) const override;
    bool point_is_in(Point v) const override;
    std::unique_ptr<Particle> enlarge(double add) const override;
    std::vector<Point> getVertices() const;
};

double distanceFromPointToSegment(const Point& point, const Point& start, const Point& end);
std::vector<Segment> getEdges(const std::vector<Point>& vertices);
void projectVertices(const std::vector<Point>& vertices, const Point& axis, double& minProj, double& maxProj);
bool intersectProjections(double minProj1, double maxProj1, double minProj2, double maxProj2);
Point get_randloc(double L);
void scr(double box_l, std::vector<std::unique_ptr<Particle>>& particles, int num, int cycle);

#endif
