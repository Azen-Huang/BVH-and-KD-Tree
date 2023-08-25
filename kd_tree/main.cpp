#include <iostream>
#include <fstream>
#include <memory>
#include <omp.h>

//#include "vec3.h"
#include "material.h"
#include "sphere.h"
#include "triangle.h"
#include "outputPPM.h"
//#include "ray.h"
#include "hittable_list.h"
#include "general.h"
#include "camera.h"

//for obj read
# include <string.h>
# define FILENAME "dragon.obj"

struct V3{  //暫存頂點或法線
    float x;
    float y;
    float z;
};
struct V2{  // 暫存vt順序
    float x;
    float y;
};

#define FOCAL_LENGTH 60
# define ia_default 0.1

color ray_color(const ray& Ray, const hittable& world, int depth, const vec3 &light_position){
    hit_record Ray_hit_record;

    if (depth < 0)
        return vec3(0, 0, 0);

    if(world.hit(Ray, 0.000001, infinity, Ray_hit_record)){
        color light_intensity(1, 1, 1);

        //=========================shadow==============================================
        ray shadow_ray(light_position, (Ray_hit_record.hit_position - light_position));
        hit_record shadow_hit_record;
        if(world.hit(shadow_ray, 0.001, 1-0.000001, shadow_hit_record)){
            light_intensity = color(0, 0, 0);
        };
        //==========================phong lighting=====================================
        vec3 N = unit_vector(Ray_hit_record.normal);  //printf("Normal: (%lf, %lf, %lf}\n", N.x(), N.y(), N.z());   
        vec3 V = unit_vector(-Ray.show_eye_dir());    //reversed eye ray from hit position   
        vec3 L = unit_vector(light_position - Ray_hit_record.hit_position);  //Light vector from hit position
        vec3 H = unit_vector((V + L) / 2);            //H for quickly calculating

        //ambient = Ka*ia
        color ia(ia_default, ia_default, ia_default); 
        //diffuse = Kd*id;  Id = Ii * N·L
        color id = Ray_hit_record.mat_ptr->Object_color*light_intensity*std::max(0.0, dot(N, L));
        //specular = Ks*is;  Is = Ii(N · H)^n
        color is = light_intensity * pow(std::max(0.0, dot(N, H)), Ray_hit_record.mat_ptr->Exp);

        color phong_liting = Ray_hit_record.mat_ptr->Ka*ia + Ray_hit_record.mat_ptr->Kd*id + Ray_hit_record.mat_ptr->Ks*is;
        //===========================reflection==========================================
        ray reflect_Ray = ray(Ray_hit_record.hit_position, reflect(-V, N));
        color reflect_color = ray_color(reflect_Ray, world, depth - 1, light_position);

        color output_color  = phong_liting + Ray_hit_record.mat_ptr->Reflection * reflect_color;
        for(int i=0;i<3;i++){
            if(output_color.e[i]>1) output_color.e[i] = 1;
            if(output_color.e[i]<0) output_color.e[i] = 0;
        }

        return output_color;

    }
    return color(0, 0, 0);
}




int main(){

    //-----------------------------------------------READ DATA-------------------------------------------------------------------------------
    std::ifstream input("hw3_input.txt");

    std::string read;
    std::string E = "E"; //eyes position
    std::string V = "V"; //first is eye direction. second is up direction to fix the sight
    std::string F = "F"; //view angle  
    std::string R = "R"; //resolution (w, h)
    std::string M = "M"; //material  (r g b) (Ka Kd Ks) (exp::specularity) Reflect
    std::string S = "S"; //sphere data
    std::string T = "T"; //triangle data
    std::string L = "L"; // light position

    // initialize variable to catch input
    point3 eye_position;
    vec3 eye_dir;
    vec3 eye_dir_up;
    double view_angle = 0;
    int resolution[2] = {0};
    ///material
    std::shared_ptr<material> material_data;
    ///sphere data
    std::shared_ptr<sphere> sphere_data;
    ///triangle data
    std::shared_ptr<triangle> triangel_data;

    point3 light_postiton;

    //other default
    const int max_hit_depth = 10;
    hittable_list world;

    while(input>>read){
        if(read == E){         //eyes position
            double x = 0, y = 0, z = 0;
            input>>x>>y>>z;
            eye_position = point3(x, y, z);
            //printf("eye position:%lf, %lf, %lf\n", eye_position.x(), eye_position.y(), eye_position.z());
        }else if(read == V){   //first is eye direction. second is up direction to fix the sight
            double x1 = 0, y1 = 0, z1 = 0, x2 = 0, y2 = 0,  z2 = 0;
            input>>x1>>y1>>z1;
            eye_dir = point3(x1, y1, z1);
            //printf("eye direction: %lf, %lf, %lf\n", eye_dir.e[0], eye_dir.e[1], eye_dir.e[2]);
            input>>x2>>y2>>z2;
            eye_dir_up = point3(x2, y2,z2);
            //printf("eye up: %lf, %lf, %lf\n", eye_dir_up.e[0], eye_dir_up.e[1], eye_dir_up.e[2]);
        }else if(read == F){   //view angle 
            input>>view_angle;
            //printf("view angle(V): %lf\n", view_angle);
            
        }else if(read == R){   //resolution (w, h)
            input>>resolution[0]>>resolution[1];
            //printf("resolution(w ,h): (%d, %d)\n", resolution[0], resolution[1]);
            
        }else if(read == M){   //material  (r g b) (Ka Kd Ks) (exp::specularity) Reflect
            double r = 0, g = 0, b = 0, ka = 0, kd = 0, ks = 0, exp = 0, reflact = 0;
            input>>r>>g>>b>>ka>>kd>>ks>>exp>>reflact;
            material_data = std::make_shared<material>(point3(r, g, b), ka, kd, ks, exp, reflact);
            //material_data->check_material_data();

        }else if(read == S){   //sphere data
            double x = 0, y = 0, z = 0, radius = 0;
            input>>x>>y>>z>>radius;
            sphere_data = std::make_shared<sphere> (point3(x, y, z), radius, material_data);
            //printf("sphere data: center(%lf, %lf, %lf), radius: %lf\n", sphere_data->center.x(), sphere_data->center.y(), sphere_data->center.z(), sphere_data->radius);
            //printf("sphere material:");
            //material_data->check_material_data();
            world.add(sphere_data);

        }else if(read == T){   //triangle data
            double x1 = 0, y1 = 0, z1 = 0, x2 = 0, y2 = 0, z2 = 0, x3 = 0, y3 = 0, z3 = 0;
            input>>x1>>y1>>z1>>x2>>y2>>z2>>x3>>y3>>z3;
            triangel_data = std::make_shared<triangle> (point3(x1, y1, z1), point3(x2, y2, z2), point3(x3, y3, z3), material_data);
            //printf("triangle data: \npoint1(%lf, %lf, %lf),\npoint2(%lf, %lf, %lf),\npoint2(%lf, %lf, %lf)\n", triangel_data->point1.x(), triangel_data->point1.y(), triangel_data->point1.z(),  triangel_data->point2.x(), triangel_data->point2.y(), triangel_data->point2.z(),  triangel_data->point3.x(), triangel_data->point3.y(), triangel_data->point3.z());
            //printf("triangle material:");
            //material_data->check_material_data();
            world.add(triangel_data);

        }else if(read == L){   //light position
            double x = 0, y = 0, z = 0;
            input>>x>>y>>z;
            light_postiton = point3(x, y, z);
        }else{ std::cerr <<"input file error"<< std::flush;}
    }
    input.close();

    //-----------------------------------------------------------read data end-------------------------------------------------------------------------------------

    //-----------------------------------------------------------read obj file-------------------------------------------------------------------------
    std::cout<<"Load OBJ file!\n";
   
    std::vector<V3> vertices; //v 資料索引
    std::vector<V2> uvs;      //vt 資料索引
    std::vector<V3> normals;  //vn 資料索引
    std::vector<int> vertIdxs, uvIdxs, nrmIdxs; //f 資料索引

    FILE* file = fopen(FILENAME, "r");
    if (!file){
        printf("File open failded");
    }else{
        printf("Now loading\n");
    
        while (true){
            char lineHeader[128];
            
            int res = fscanf(file, "%s", lineHeader);// read the first word of the line

            if (res == EOF){
                printf("Read End!\n");
                break;
            }

            if (strcmp(lineHeader, "v") == 0){
                V3 vert; //頂點
                fscanf(file, "%f %f %f\n", &vert.x, &vert.y, &vert.z );
                vertices.push_back(vert); //儲存頂點
            }else if(strcmp(lineHeader, "vt") == 0){
                V2 uv; //uv
                fscanf(file, "%f %f %f\n", &uv.x, &uv.y );
                uvs.push_back(uv); //儲存uv
            }else if(strcmp(lineHeader, "vn") == 0){
                V3 normal; //法線
                fscanf(file, "%f %f %f\n", &normal.x, &normal.y, &normal.z);
                normals.push_back(normal); //儲存法線資料
            }else if(strcmp(lineHeader, "f") == 0){
                int vertIndex[3], uvIndex[3], normalIndex[3]; //uv
                int surface = fscanf(file, "%d//%d %d//%d %d//%d\n", 
                                    &vertIndex[0], &normalIndex[0],
                                    &vertIndex[1], &normalIndex[1],
                                    &vertIndex[2], &normalIndex[2]
                                );
                if (surface != 6){
                    printf("interface read failded\n");
                }

                vertIdxs.push_back(vertIndex[0]);
                vertIdxs.push_back(vertIndex[1]);
                vertIdxs.push_back(vertIndex[2]);

                // uvIdxs.push_back(uvIndex[0]);
                // uvIdxs.push_back(uvIndex[1]);
                // uvIdxs.push_back(uvIndex[2]);

                // nrmIdxs.push_back(normalIndex[0]);
                // nrmIdxs.push_back(normalIndex[1]);
                // nrmIdxs.push_back(normalIndex[2]);

            }else{
                //不需要的資料就跳過
                char stupidBuffer[1000];
                fgets(stupidBuffer, 1000, file);
            }
        }
        fclose(file);
        printf("file close\n");

        /*
        printf("vertices[0]: x:%lf, y:%lf, z:%lf\n", vertices[0].x, vertices[0].y, vertices[0].z);
        printf("vertices[1]: x:%lf, y:%lf, z:%lf\n", vertices[1].x, vertices[1].y, vertices[1].z);
        printf("vertices[2]: x:%lf, y:%lf, z:%lf\n\n", vertices[2].x, vertices[2].y, vertices[2].z);

        printf("vertIdxs[0]: %d\n", vertIdxs[0]);
        printf("vertIdxs[1]: %d\n", vertIdxs[1]);
        printf("vertIdxs[2]: %d\n\n", vertIdxs[2]);
        */
        int total_vert_num = vertices.size();
        //printf("total_vert_num:%d\n", total_vert_num);
        // ===================以下做面的資料重組==============================
        for(int i = 0; i < vertIdxs.size(); i +=3){//vertIdxs.size(); i++){
            // printf("vertIdxs:%d, now i:%d \n", vertIdxs.size(), i);
            
            int vertIndex1 = vertIdxs[i];
            if(vertIndex1<0){
                vertIndex1 = total_vert_num + vertIndex1 + 1; 
                //printf("vertIndex1:%d\n", vertIndex1);
            }
            int vertIndex2 = vertIdxs[i+1];
            if(vertIndex2<0){
                vertIndex2 = total_vert_num + vertIndex2 + 1; 
                //printf("vertIndex2:%d\n", vertIndex2);
            }
            int vertIndex3 = vertIdxs[i+2];
            if(vertIndex3<0){
                vertIndex3 = total_vert_num + vertIndex3 + 1; 
                //printf("vertIndex3:%d\n", vertIndex3);
            }

            //printf("debug:%d\n", vertIndex);
            V3 vert1 = vertices[vertIndex1 - 1];
            V3 vert2 = vertices[vertIndex2 - 1];
            V3 vert3 = vertices[vertIndex3 - 1];
            // printf("vert: x:%lf, y:%lf, z:%lf\n", vert1.x, vert1.y, vert1.z);
            // printf("vert: x:%lf, y:%lf, z:%lf\n", vert2.x, vert2.y, vert2.z);
            // printf("vert: x:%lf, y:%lf, z:%lf\n", vert3.x, vert3.y, vert3.z);

            //int uvIndex1   = uvIdxs[i];
            //int uvIndex2   = uvIdxs[i];
            //int uvIndex3   = uvIdxs[i];

            //V2 uv1 = uvs[uvIndex1 - 1];
            //V2 uv2 = uvs[uvIndex2 - 1];
            //V2 uv3 = uvs[uvIndex3 - 1];

            //int normalIndex1 = nrmIdxs[i];
            //int normalIndex2 = nrmIdxs[i];
            //int normalIndex3 = nrmIdxs[i];

            //V3 nrm1 = normals[normalIndex1 - 1];
            //V3 nrm2 = normals[normalIndex2 - 1];
            //V3 nrm3 = normals[normalIndex3 - 1];

            triangel_data = std::make_shared<triangle> (point3(vert1.x, vert1.y, vert1.z), point3(vert2.x, vert2.y, vert2.z), point3(vert3.x, vert3.y, vert3.z), material_data);
            //triangel_data = std::make_shared<triangle> (point3(vert1.x*20, -vert1.y*20, vert1.z*0), point3(vert2.x*-20, vert2.y*-20, vert2.z*30), point3(vert3.x*20, vert3.y*-20, vert3.z*-30), material_data);
            world.add(triangel_data);

        }
        printf("load obj file done!\n");
    }

    double st = omp_get_wtime();
    world.make_kd_tree();
    st = omp_get_wtime()-st;
    printf("make kd-tree done, %f(s)!\n", st);
    //------------------------------------------------------------camera----------------------------------------------------------------------
    // eye_position;  (lookfrom)
    // eye_dir;       (lookat)
    // eye_dir_up;
    // view_angle
    // resolution[2] => (w, h);
    //vec3 vup(0,1,0);
    double focal_length = 20;  //(dist_to_focus)
    auto aperture = 1.5/2.2;
    double viewport_screen_ratio = resolution[0]/resolution[1];
    // double viewport_height = focal_length*tan(view_angle*pi/180)*2;
    // double viewport_width = viewport_height*viewport_screen_ratio;
    
    camera cam(eye_position, eye_dir, eye_dir_up, view_angle, viewport_screen_ratio, aperture, focal_length);
    //--------------------------------------------------------------VIEWPORT---------------------------------------------------------------------------------------
    //resolution(w, h)
    
    //printf("Windows Property: viewport(W:%lf, H:%lf), focal length: %lf, resolution:(W:%d, H:%d)\n", viewport_width, viewport_height, focal_length, resolution[0], resolution[1]);

    // vec3 horizotal_vector = viewport_width * unit_vector(cross(eye_dir, eye_dir_up));
    // vec3 vertical_vector = viewport_height * unit_vector(-eye_dir_up);
    // printf("horizotal_vector: (x:%lf, y:%lf, z:%lf)\n", horizotal_vector.x(), horizotal_vector.y(), horizotal_vector.z());
    // printf("vertical_vector: (x:%lf, y:%lf, z:%lf)\n", vertical_vector.x(), vertical_vector.y(), vertical_vector.z());
    // point3 top_left_corner = eye_position + focal_length*unit_vector(eye_dir) + viewport_height/2 * unit_vector(eye_dir_up) - viewport_width/2 * unit_vector(horizotal_vector);
    // printf("top_left_corner:( %lf, %lf, %lf)\n", top_left_corner.x(), top_left_corner.y(), top_left_corner.z());

    //---------------------------------------------------------------OUT PUT PPM--------------------------------------------------------------------------------------
    ColorImage image;
    Pixel pixel = {0, 0, 0};
    image.init(resolution[0], resolution[1]);


    //--------------------------------------------------------------Antialiasing--------------------------------------------------------------------------------------
    int samples_per_pixel = 1;

    //

    st = omp_get_wtime();
    #pragma omp parallel for private(pixel)
    for(int ij = 0; ij < resolution[0]*resolution[1]; ij++){
        double temp_color[3] = {};

        //for(int k = 0; k < samples_per_pixel; k++){
            // auto u = double(i+random_double()) / (resolution[0]-1);
            // auto v = double(j+random_double()) / (resolution[1]-1);
            auto u = double(ij/resolution[1]) / (resolution[0]-1);
            auto v = double(ij%resolution[1]) / (resolution[1]-1);

            //printf("%lf, %lf\n", u, v);
            //ray ray(eye_position, (top_left_corner + u*horizotal_vector + v*vertical_vector) - eye_position);
            ray ray = cam.get_ray(u, v);
            //std::cout<<ray.at(1).x()<<"\n";
            //printf("(%lf, %lf, %lf)\n", ray.Eye_direction.x(), ray.Eye_direction.y(), ray.Eye_direction.z());
            color pixel_color = ray_color(ray, world, max_hit_depth, light_postiton);

            temp_color[0] = pixel_color.x();
            temp_color[1] = pixel_color.y();
            temp_color[2] = pixel_color.z();

        //}
        
        // pixel.R = static_cast<int>(255.999 *temp_color[0]/samples_per_pixel);
        // pixel.G = static_cast<int>(255.999 *temp_color[1]/samples_per_pixel);
        // pixel.B = static_cast<int>(255.999 *temp_color[2]/samples_per_pixel);

        pixel.R = static_cast<int>(255.0*temp_color[0]);
        pixel.G = static_cast<int>(255.0*temp_color[1]);
        pixel.B = static_cast<int>(255.0*temp_color[2]);

        image.writePixel(ij/resolution[1], ij%resolution[1], pixel);
    }
    st = omp_get_wtime()-st;
    printf("time = %f(s)\n", st);

    char output_name[15] = "output_pic.ppm" ;
    image.outputPPM(output_name);
    std::cerr << "\nDone.\n";


    return 0;
}