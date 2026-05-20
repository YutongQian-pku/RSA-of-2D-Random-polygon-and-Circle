#ifndef _PARTICLE_TYPE_H_
#define _PARTICLE_TYPE_H_

#include <string>
#include <vector>
#include "point.h"

// Particle shape enum
enum class ParticleShape {
    CIRCLE,
    POLYGON
};

// Particle type information.  For polygons, vertices are the absolute-size
// template coordinates translated to the polygon centroid, i.e. local vertices.
struct ParticleType {
    ParticleShape shape = ParticleShape::CIRCLE;
    int shape_id = 0;                 // polygon file id, e.g. shapes/1.txt -> 1
    double radius = 0.0;              // circle radius, or polygon circumradius about centroid
    double prob = 0.0;                // selection probability
    std::vector<Point> vertices;      // local polygon vertices, clockwise, without repeated last point
    double area = 0.0;                // actual particle area
    double inner_r = 0.0;             // circle radius fully contained around center/centroid
    double character_L = 0.0;         // 2.2 * radius
    std::string label;                // display/output label
};

#endif
