#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include "particle.h"

using namespace std;

namespace {
constexpr double EPS = 1e-10;

Point rotate_point(const Point& p, double theta) {
    double c = std::cos(theta);
    double s = std::sin(theta);
    return Point(c * p.x - s * p.y, s * p.x + c * p.y);
}

double polygon_base_radius(const vector<Point>& local_vertices) {
    double r = 0.0;
    for (const auto& v : local_vertices) r = max(r, PointMag(v));
    return r;
}

double polygon_base_inner_radius(const vector<Point>& local_vertices) {
    if (local_vertices.size() < 3) return 0.0;
    Point origin(0.0, 0.0);
    double r = numeric_limits<double>::max();
    for (size_t i = 0; i < local_vertices.size(); ++i) {
        size_t j = (i + 1) % local_vertices.size();
        r = min(r, distanceFromPointToSegment(origin, local_vertices[i], local_vertices[j]));
    }
    return r == numeric_limits<double>::max() ? 0.0 : r;
}

bool pointInPolygonRay(const vector<Point>& vertices, Point v) {
    int count = 0;
    const size_t n = vertices.size();
    for (size_t i = 0; i < n; ++i) {
        size_t j = (i + 1) % n;

        double cross = (vertices[j].x - vertices[i].x) * (v.y - vertices[i].y)
                     - (vertices[j].y - vertices[i].y) * (v.x - vertices[i].x);
        double dot1 = (v.x - vertices[i].x) * (v.x - vertices[j].x)
                    + (v.y - vertices[i].y) * (v.y - vertices[j].y);
        if (std::fabs(cross) < EPS && dot1 <= EPS) return true;

        if ((vertices[i].y <= v.y && vertices[j].y > v.y) ||
            (vertices[j].y <= v.y && vertices[i].y > v.y)) {
            double xIntersect = vertices[i].x +
                (v.y - vertices[i].y) * (vertices[j].x - vertices[i].x) /
                (vertices[j].y - vertices[i].y);
            if (xIntersect > v.x) count++;
        }
    }
    return (count % 2 == 1);
}

bool circlePolygonIntersects(const Point& circle_center,
                             double circle_radius,
                             const Polygon& polygon) {
    const auto vertices = polygon.getVertices();
    const double effective_radius = circle_radius + polygon.buffer_radius;

    // Any polygon vertex inside the circle.
    for (const auto& vertex : vertices) {
        if (PointMag(vertex - circle_center) < effective_radius - EPS) return true;
    }

    // Circle intersects one of polygon edges, including vertex neighborhoods.
    for (size_t i = 0; i < vertices.size(); ++i) {
        size_t j = (i + 1) % vertices.size();
        double dist = distanceFromPointToSegment(circle_center, vertices[i], vertices[j]);
        if (dist < effective_radius - EPS) return true;
    }

    // Circle center inside polygon means containment/intersection.
    if (polygon.point_is_in(circle_center)) return true;

    return false;
}
}

Circle::Circle() {}
Circle::Circle(double cR, Point ccenter) {
    center = ccenter;
    R = cR;
    inner_r = cR;
    character_L = 2.2 * R;
}
Circle::~Circle() {}

Polygon::Polygon() {}
Polygon::Polygon(const vector<Point>& localVertices,
                 double angle,
                 Point ccenter,
                 int shapeId,
                 double bufferRadius)
    : shape_id(shapeId), theta(angle), buffer_radius(bufferRadius), local_vertices(localVertices) {
    center = ccenter;
    const double base_r = polygon_base_radius(local_vertices);
    const double base_inner = polygon_base_inner_radius(local_vertices);
    R = base_r + buffer_radius;
    inner_r = base_inner + buffer_radius;
    character_L = 2.2 * R;
}
Polygon::~Polygon() {}

void Circle::translate(double x, double y) {
    center.x += x;
    center.y += y;
}

bool Circle::intersects(Particle* p) const {
    if (PointMag(center - p->center) >= (R + p->R) - EPS) return false;

    const Circle* c = dynamic_cast<const Circle*>(p);
    if (c) {
        double dist = PointMag(center - c->center);
        return dist < (R + c->R) - EPS;
    }

    const Polygon* poly = dynamic_cast<const Polygon*>(p);
    if (poly) {
        return circlePolygonIntersects(center, R, *poly);
    }
    return false;
}

bool Circle::is_onboundary(double box_l) const {
    double up = center.y + R;
    double down = center.y - R;
    double left = center.x - R;
    double right = center.x + R;
    return (up >= box_l || down <= 0 || left <= 0 || right >= box_l);
}

bool Circle::point_is_in(Point v) const {
    return PointMag(v - center) <= R + EPS;
}

unique_ptr<Particle> Circle::enlarge(double add) const {
    return make_unique<Circle>(R + add, center);
}

void Polygon::translate(double x, double y) {
    center.x += x;
    center.y += y;
}

bool Polygon::intersects(Particle* p) const {
    if (PointMag(center - p->center) >= (R + p->R) - EPS) return false;

    const Circle* c = dynamic_cast<const Circle*>(p);
    if (c) {
        return circlePolygonIntersects(c->center, c->R, *this);
    }

    const Polygon* poly = dynamic_cast<const Polygon*>(p);
    if (poly) {
        const auto vertices1 = getVertices();
        const auto vertices2 = poly->getVertices();
        const double expand1 = buffer_radius;
        const double expand2 = poly->buffer_radius;

        for (const auto& edge1 : getEdges(vertices1)) {
            double minProj1, maxProj1, minProj2, maxProj2;
            projectVertices(vertices1, edge1.normal, minProj1, maxProj1);
            projectVertices(vertices2, edge1.normal, minProj2, maxProj2);
            minProj1 -= expand1; maxProj1 += expand1;
            minProj2 -= expand2; maxProj2 += expand2;
            if (!intersectProjections(minProj1, maxProj1, minProj2, maxProj2)) return false;
        }
        for (const auto& edge2 : getEdges(vertices2)) {
            double minProj1, maxProj1, minProj2, maxProj2;
            projectVertices(vertices1, edge2.normal, minProj1, maxProj1);
            projectVertices(vertices2, edge2.normal, minProj2, maxProj2);
            minProj1 -= expand1; maxProj1 += expand1;
            minProj2 -= expand2; maxProj2 += expand2;
            if (!intersectProjections(minProj1, maxProj1, minProj2, maxProj2)) return false;
        }
        return true;
    }
    return false;
}

bool Polygon::is_onboundary(double box_l) const {
    double up = center.y + R;
    double down = center.y - R;
    double left = center.x - R;
    double right = center.x + R;
    return (up >= box_l || down <= 0 || left <= 0 || right >= box_l);
}

bool Polygon::point_is_in(Point v) const {
    const auto vertices = getVertices();
    if (vertices.size() < 3) return false;

    if (pointInPolygonRay(vertices, v)) return true;

    // Exact membership in polygon dilated by a disk of radius buffer_radius.
    // This is used by voxel pruning; it avoids over-rejecting valid RSA trials.
    if (buffer_radius > EPS) {
        for (size_t i = 0; i < vertices.size(); ++i) {
            size_t j = (i + 1) % vertices.size();
            if (distanceFromPointToSegment(v, vertices[i], vertices[j]) <= buffer_radius + EPS) return true;
        }
    }
    return false;
}

unique_ptr<Particle> Polygon::enlarge(double add) const {
    return make_unique<Polygon>(local_vertices, theta, center, shape_id, buffer_radius + add);
}

vector<Point> Polygon::getVertices() const {
    vector<Point> vertices;
    vertices.reserve(local_vertices.size());
    for (const auto& local : local_vertices) {
        Point r = rotate_point(local, theta);
        vertices.emplace_back(center.x + r.x, center.y + r.y);
    }
    return vertices;
}

double distanceFromPointToSegment(const Point& point, const Point& start, const Point& end) {
    Point seg = end - start;
    double segLen2 = seg * seg;
    if (segLen2 <= EPS) return PointMag(point - start);
    double t = ((point - start) * seg) / segLen2;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
    Point projection = start + seg * t;
    return PointMag(point - projection);
}

vector<Segment> getEdges(const vector<Point>& vertices) {
    vector<Segment> edges;
    edges.reserve(vertices.size());
    for (size_t i = 0; i < vertices.size(); ++i) {
        size_t j = (i + 1) % vertices.size();
        edges.emplace_back(vertices[i], vertices[j]);
    }
    return edges;
}

void projectVertices(const vector<Point>& vertices, const Point& axis, double& minProj, double& maxProj) {
    minProj = vertices[0] * axis;
    maxProj = minProj;
    for (const auto& vertex : vertices) {
        double proj = vertex * axis;
        if (proj < minProj) minProj = proj;
        if (proj > maxProj) maxProj = proj;
    }
}

bool intersectProjections(double minProj1, double maxProj1, double minProj2, double maxProj2) {
    return std::max(minProj1, minProj2) < std::min(maxProj1, maxProj2) - EPS;
}

Point get_randloc(double L) {
    Point result;
    result.x = L * double(rand()) / double(RAND_MAX);
    result.y = L * double(rand()) / double(RAND_MAX);
    return result;
}

void scr(double box_l, vector<unique_ptr<Particle>>& particles, int num, int cycle) {
    char cha[40];
    sprintf(cha, "packing_%d.scr", cycle);
    ofstream scr(cha);
    scr << fixed << setprecision(8);

    scr << "-osnap off" << endl;
    scr << "erase all " << endl;
    scr << "vscurrent 2" << endl;

    for (int i = 0; i < num; ++i) {
        const Circle* c = dynamic_cast<const Circle*>(particles[i].get());
        if (c) {
            scr << "color " << 1 << endl;
            scr << "circle " << c->center.x << "," << c->center.y << " " << c->R << endl;
            continue;
        }

        const Polygon* poly = dynamic_cast<const Polygon*>(particles[i].get());
        if (poly) {
            scr << "color " << 2 << endl;
            const auto vertices = poly->getVertices();
            scr << "pline ";
            for (const auto& vertex : vertices) scr << vertex.x << "," << vertex.y << " ";
            scr << "c" << endl;
        }
    }

    scr << "color 255" << endl;
    scr << "rectang 0,0 " << box_l << "," << box_l << endl;
    scr << "zoom e" << endl;
    scr << "vpoint 0,0,1" << endl;
    scr.close();
}
