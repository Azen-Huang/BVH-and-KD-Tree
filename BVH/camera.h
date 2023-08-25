#ifndef CAMERA_H
#define CAMERA_H

#include "vec3.h"
const double pi = 3.1415926535897932385;
class camera {
    public:
        camera(){}
        //camera() : camera(point3(0,0,-1), point3(0,0,0), vec3(0,1,0), 40, 1, 0, 10) {}

        camera(
            point3 lookfrom,
            point3 lookat,
            vec3   eye_dir_up,
            double view_angle, // vertical field-of-view in degrees
            double aspect_ratio,
            double aperture,
            double focus_dist//,
            //double _time0 = 0,
            //double _time1 = 0
        ) {

            double viewport_height = focus_dist*tan(view_angle*pi/180)*2;
            double viewport_width = viewport_height*aspect_ratio;
            //printf("Windows Property: viewport(W:%lf, H:%lf), focal length: %lf, resolution:(W:%d, H:%d)\n", viewport_width, viewport_height, focus_dist, resolution[0], resolution[1]);
            horizontal = viewport_width * unit_vector(cross(lookat-lookfrom, eye_dir_up));
            vertical = viewport_height * unit_vector(-eye_dir_up);
            
            // auto h = tan(theta/2);
            // auto viewport_height = 2.0 * h;
            // auto viewport_width = aspect_ratio * viewport_height;

            // w = unit_vector(lookfrom - lookat);
            // u = unit_vector(cross(vup, w));
            // v = cross(w, u);
            u = unit_vector(horizontal);
            v = unit_vector(eye_dir_up);

            origin = lookfrom;
            top_left_corner = origin + focus_dist*unit_vector(lookat-lookfrom) + viewport_height/2 * unit_vector(eye_dir_up) - viewport_width/2 * unit_vector(horizontal);
            // horizontal = focus_dist * viewport_width * u;
            // vertical = focus_dist * viewport_height * v;
            // top_left_corner = origin - horizontal/2 + vertical/2 + focus_dist*w;

            lens_radius = aperture / 2;
            //time0 = _time0;
            //time1 = _time1;
        }

        ray get_ray(double s, double t) const {
            vec3 rd = lens_radius * random_in_unit_disk();
            vec3 offset = u * rd.x() + v * rd.y();
            return ray(
                origin , //+ offset,
                top_left_corner + s*horizontal + t*vertical - origin //- offset
            );
        }

    private:
        point3 origin;
        point3 top_left_corner;
        vec3 horizontal;
        vec3 vertical;
        vec3 u, v;//, w;
        double lens_radius;
        //double time0, time1;  // shutter open/close times
};

#endif
