#ifndef MATERIAL_H
#define MATERIAL_H


#include "vec3.h"


class material{
    public:
        material(){};
        material(color object_color, double ka, double kd, double ks, double exp, double reflection):
            Object_color(object_color), Ka(ka), Kd(kd), Ks(ks), Exp(exp), Reflection(reflection) {}; 
        
        void check_material_data();

    public:
        color Object_color;
        double Ka;
        double Kd;
        double Ks;
        double Exp;
        double Reflection;
};


void material::check_material_data(){
    printf("material data: color(r:%lf, g:%lf, b:%lf), ka:%lf, kd:%lf, ks:%lf, exp:%lf, reflect:%lf\n", Object_color.x(), Object_color.y(), Object_color.z(), Ka, Kd, Ks, Exp, Reflection);
}


















#endif