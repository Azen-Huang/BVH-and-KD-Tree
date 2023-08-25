#ifndef HITTABLE_H
#define HITTABLE_H


//#include "rtweekend.h"
#include "general.h"


class material;

struct hit_record {
    point3 hit_position;    //hit position
    vec3 normal;            //hit position's normal
    std::shared_ptr<material> mat_ptr; //hit object's material
    double t;    //
    bool front_face; //the direction of normal: point outward or inward




    inline void set_face_normal(const ray& r, const vec3& outward_normal) {
        front_face = dot(r.show_eye_dir(), outward_normal) < 0; //判斷法線指向外為1, 內為0
        normal = front_face ? outward_normal :-outward_normal;  // 對法線加+-號
    }
};

class hittable {
    public:
        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
        
        virtual char get_kind()=0;
        virtual point3 get_center()=0;
        virtual void get_outer_bound(point3* p)=0;
};

#endif