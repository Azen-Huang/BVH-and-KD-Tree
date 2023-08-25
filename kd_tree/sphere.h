#ifndef SPHERE_H
#define SPHERE_H


#include "general.h"
#include "hittable.h"

//----------------------------------------------------sphere----------------------------------------------------------------------
class sphere: public hittable {
    public:
        sphere(){}
        sphere(point3 cen, double r, std::shared_ptr<material> m)
            : center(cen), radius(r), mat_ptr(m) {};
        
        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;

        virtual char get_kind(){return 's';};
        virtual point3 get_center(){return center;};
        virtual void get_outer_bound(point3* p){
            p[0] = center+radius*point3(1, 0, 0);
            p[1] = center+radius*point3(-1, 0, 0);
            p[2] = center+radius*point3(0, 1, 0);
            p[3] = center+radius*point3(0, -1, 0);
            p[4] = center+radius*point3(0, 0, 1);
            p[5] = center+radius*point3(0, 0, -1);
        };

    public:
        point3 center;
        double radius;
        std::shared_ptr<material> mat_ptr;
    
};


bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    vec3 oc = r.show_eye_position() - center;
    auto a = r.show_eye_dir().length_squared();
    auto half_b = dot(oc, r.show_eye_dir());
    auto c = oc.length_squared() - radius*radius;

    auto discriminant = half_b*half_b - a*c;
    if (discriminant < 0) return false;
    auto sqrtd = sqrt(discriminant);

    // Find the nearest root that lies in the acceptable range.
    auto root = (-half_b - sqrtd) / a;
    if (root < t_min || t_max < root) {
        root = (-half_b + sqrtd) / a;
        if (root < t_min || t_max < root)
            return false;
    }

    rec.t = root;
    rec.hit_position = r.at(rec.t);
    vec3 outward_normal = (rec.hit_position - center) / radius;
    rec.set_face_normal(r, outward_normal);
    rec.mat_ptr = mat_ptr;

    return true;
}









#endif