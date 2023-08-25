#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "general.h"
#include "hittable.h"

class triangle : public hittable{
    public:
        triangle() {}

        triangle(point3 p1, point3 p2, point3 p3, std::shared_ptr<material> m)
            : point_1(p1), point_2(p2), point_3(p3), mat_ptr(m) {};

        virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;

        virtual char get_kind(){return 't';};
        virtual point3 get_center(){return (point_1+point_2+point_3)/3;};
        virtual void get_outer_bound(point3* p){
            p[0] = point_1;
            p[1] = point_2;
            p[2] = point_3;
        };

        
    public:
        point3 point_1;
        point3 point_2;
        point3 point_3;
        std::shared_ptr<material> mat_ptr;

};


bool triangle::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    bool isIn = false;
    
    //build edge
    vec3 edge12 = point_2 - point_1;
    vec3 edge13 = point_3 - point_1;
    
    vec3 S = r.show_eye_position() - point_1;
    vec3 S1 = cross(r.show_eye_dir(), edge13);
    vec3 S2 = cross(S, edge12);

    double coeff = 1.0f / dot(S1, edge12);
    double t = coeff * dot(S2, edge13);
    double b1 = coeff * dot(S1, S);
    double b2 = coeff * dot(S2, r.show_eye_dir());

    if (t>=0 && b1>=0 && b2 >=0 && (1-b1-b2)>=0){
        isIn = true;
        // tnear = t;
        // u = b1;
        // v = b2;
        rec.t = t;
        rec.hit_position = r.at(rec.t);
        vec3 outward_normal = cross(edge12, edge13);
        rec.set_face_normal(r, outward_normal);
        rec.mat_ptr = mat_ptr;
    }
    
    if (t < t_min || t_max < t) {
        return false;
    }

    return isIn;
    
    
    
    
    
    /*
    vec3 pvec = cross(r.show_eye_dir(), edge13); //pvec == r.direction() x edge13
    double det = dot(edge12, pvec);

    if (det < 1e-8f){
        return false;
    }

    auto inv_det = 1.0f/det;
    // 先計算T
    auto tvec = r.show_eye_position() - point_1; 
    // 求u解
    double u = dot(tvec, pvec) * inv_det;
    if (u < 0.0 || u > 1.0){
        return false;
    }

    //準備Q == T‧E1
    auto qvec = cross(tvec, edge12);
    // 求v解
    double v = dot(r.show_eye_dir(), qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0f){
        return false;
    }
    // 求t解
    float t = dot(edge13, qvec) * inv_det;
    if(t <= 0){
        return false;
    }
    */
    return true;
}



#endif