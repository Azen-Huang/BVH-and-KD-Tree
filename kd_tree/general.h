#ifndef GENERAL_H
#define GENERAL_H


#include <cmath>
#include <limits>
#include <memory>

#include <cstdlib>




// Constants

const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;


// Usings
using std::shared_ptr;
using std::make_shared;
using std::sqrt;



// Utility Functions
inline double random_double() {
    // Returns a random real in [0,1).
    return rand() / (RAND_MAX + 1.0);
}

inline double random_double(double min, double max) {
    // Returns a random real in [min,max).
    return min + (max-min)*random_double();
}

inline double clamp(double x, double min, double max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

inline double degrees_to_radians(double degrees) {
    return degrees * pi / 180.0;
}

inline vec3 random_in_unit_disk() {
    while (true) {
        auto p = vec3(random_double(-1,1), random_double(-1,1), 0);
        if (p.length_squared() >= 1) continue;
        return p;
    }
}

// Common Headers

#include "ray.h"
#include "vec3.h"






#endif