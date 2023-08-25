#ifndef RAY_H
#define RAY_H

#include "vec3.h"

class ray {
    public:
        ray() {}
        ray(const point3& eye_pos, const vec3& eye_dir)
            : Eye_position(eye_pos), Eye_direction(eye_dir)
        {}

        point3 show_eye_position() const  { return Eye_position; }
        vec3 show_eye_dir() const { return Eye_direction; }

        point3 at(double t) const {
            return Eye_position + t*Eye_direction;
        }

    public:
        point3 Eye_position;
        vec3 Eye_direction;
};

#endif