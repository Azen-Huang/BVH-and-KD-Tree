#include "color.h"
#include "vec3.h"
#include "2ppm.h"
#include "ray.h"
#include "matrix3.h"
#include "material.h"
#include "camera.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <memory>
#include <ctime>
using namespace std;
#define FILENAME "box.obj"
#define INF 114514.0

struct V3{  //暫存頂點或法線
    float x;
    float y;
    float z;
};
struct V2{  // 暫存vt順序
    float x;
    float y;
};

#define PI 3.1415926535897932385
vec3 E; //視線眼睛的座標
vec3 UL,UR,LL,LR; //2D視窗的3D座標

int image_width,image_height; //image的長寬

point3 Center; //圓的圓心座標
float Radius; //圓的半徑

vector<point3> centers, t1, t2, t3;
vector<float> radius;


point3 light_pos; //燈的位置
vector<float> R,G,B; //各個物體的RGB顏色
vector<float> ka,kd,ks; //各個物體的材質參數 圓：0 三角：1（根據讀檔輸入順序會不同）
vector<float> Exp,Ref; //Exp 算phong method的參數，Ref各個物體反射顏色參數
int angle; //視線能看到的角度
vec3 D;  //視線的方向向量
vec3 U;  //視線向上的方向向量

string read_txt = "hw3_input.txt";

void readfile(){
    ifstream source;                    // build a read-Stream
    source.open(read_txt, ios_base::in);
    for(string line; getline(source, line); )   //read stream line by line
    {
        stringstream in(line);      //make a stream for the line itself
        string type;
        in >> type;                  //and read the first whitespace-separated token
        if(type == "E")       
        {
            float x, y, z;
            in >> x >> y >> z;       //now read the whitespace-separated floats
            E = vec3(x, y, z);
        }
        else if(type == "O"){
            float ulx, uly, ulz, urx, ury, urz, llx, lly, llz, lrx, lry, lrz;
            in >> ulx >> uly >> ulz >> urx >> ury >> urz >> llx >> lly >> llz >> lrx >> lry >> lrz;
            UL = vec3(ulx, uly, ulz);
            UR = vec3(urx, ury, urz);
            LL = vec3(llx, lly, llz);
            LR = vec3(lrx, lry, lrz);
        }
        else if(type == "R"){
            in >> image_width >> image_height;
        }
        else if(type == "S"){
            float x, y, z, r;
            in >> x >> y >> z >> r;
            Center = point3(x, y, z);
            Radius = r;
            centers.push_back(Center);
            radius.push_back(Radius);
        }
        else if(type == "T"){
            float ulx, uly, ulz, urx, ury, urz, llx, lly, llz;
            vec3 T1, T2, T3;
            in >> ulx >> uly >> ulz >> urx >> ury >> urz >> llx >> lly >> llz;
            T1 = vec3(ulx, uly, ulz);
            T2 = vec3(urx, ury, urz);
            T3 = vec3(llx, lly, llz);
            t1.push_back(T1);
            t2.push_back(T2);
            t3.push_back(T3);
        }
        else if(type == "M"){
            float r, g, b, kka, kkd, kks, ex, re;
            in >> r >> g >> b >> kka >> kkd >> kks >> ex >> re;
            R.push_back(r);
            G.push_back(g);
            B.push_back(b);
            ka.push_back(kka);
            kd.push_back(kkd);
            ks.push_back(kks);
            Exp.push_back(ex);
            Ref.push_back(re);
        }
        else if(type == "F"){
            in >> angle;
        }
        else if(type == "V"){
            float dx, dy, dz, ux, uy, uz;
            in >> dx >> dy >> dz >> ux >> uy >> uz;
            D = vec3(dx, dy, dz);
            U = vec3(ux, uy, uz);
        }
        else if(type == "L"){
            float x,y,z;
            in >> x >> y >> z;
            light_pos = point3(x, y, z);
        }
    }
}

double random_double() {
    // Returns a random real in [0,1).
    return rand() / (RAND_MAX + 1.0);
}

float hit_sphere(const point3& center, float radius, const ray& r) {
    vec3 oc = r.origin() - center;
    auto a = dot(r.direction(), r.direction());
    auto b = 2.0 * dot(oc, r.direction());
    auto c = dot(oc, oc) - radius*radius;
    auto discriminant = b*b - 4*a*c;
    if (discriminant < 0) {
        return -1.0;
    } else {
        return (-b - sqrt(discriminant) ) / (2.0*a);
    }
}

float hit_triangle(const point3& t1,const point3& t2,const point3& t3,const ray& r){
    auto v1 = vec3(t1.x(), t1.y(), t1.z());
    auto v2 = vec3(t2.x() - t1.x(),t2.y() - t1.y(),t2.z() - t1.z());
    auto v3 = vec3(t3.x() - t1.x(),t3.y() - t1.y(),t3.z() - t1.z());
    auto o = r.origin();
    auto d = r.direction();
    auto m = matrix3(
        v2.x(), v3.x(), -d.x(),
        v2.y(), v3.y(), -d.y(),
        v2.z(), v3.z(), -d.z());

    auto inver_m = matrix3::inverse(m);

    float s1 = (((o.x() - v1.x()) * inver_m.M[0][0]) + ((o.y() - v1.y()) * inver_m.M[0][1]) + ((o.z() - v1.z()) * inver_m.M[0][2]));
    float s2 = (((o.x() - v1.x()) * inver_m.M[1][0]) + ((o.y() - v1.y()) * inver_m.M[1][1]) + ((o.z() - v1.z()) * inver_m.M[1][2]));
    float  t = (((o.x() - v1.x()) * inver_m.M[2][0]) + ((o.y() - v1.y()) * inver_m.M[2][1]) + ((o.z() - v1.z()) * inver_m.M[2][2]));
    
    if(t > 0.0 && (s1 + s2 <= 1.0) && s1 <= 1.0 && s1 >= 0.0 && s2 <= 1.0 && s2 >= 0.0)
        return t;
    else
        return -1;
}

float Diffuse_Compont(const vec3& N,const vec3& L){
    return dot(L, N);
}

float Phong_model(const vec3& N,const vec3& H,float e){
    //cout << dot(N, L) << endl;
    return pow(dot(N, H),e);
}

//------------------------------------------------BVH------------------------------------------------------
typedef struct Triangle {
    point3 p1, p2, p3;   // 三点
    point3 center;       // 中心
    Triangle(point3 a, point3 b, point3 c)
    {
        p1 = a, p2 = b, p3 = c;
        center = (p1 + p2 + p3) / 3;
    }
} Triangle;

std::vector<Triangle> triangles;

//按照三角形中心排序 -- 比較函數
bool cmpx(const Triangle& t1, const Triangle& t2) {
    return t1.center.x() < t2.center.x();
}
bool cmpy(const Triangle& t1, const Triangle& t2) {
    return t1.center.y() < t2.center.y();
}
bool cmpz(const Triangle& t1, const Triangle& t2) {
    return t1.center.z() < t2.center.z();
}

// BVH 樹節點
struct BVHNode {
    BVHNode* left = NULL;       // 左右子樹
    BVHNode* right = NULL;
    int n, index;               // 葉子節點訊息              
    vec3 AA, BB;                // AABB碰撞盒
};
//其中 n 不為 0 時表示是葉子節點，n 是存儲的三角形的個數，而 triangles 數組中 [index, index + n -1] 範圍的三角形都屬於該節點
// SAH 優化建構 BVH
int treecount = 0;

// n代表上面的n
BVHNode *buildBVHwithSAH(std::vector<Triangle> &triangles, int l, int r, int n)
{
    if (l > r) return 0;
    if(treecount % 200000 == 0)
        cout << "build " << treecount << endl;
    treecount++;
    BVHNode *node = new BVHNode();
    node->AA = vec3(1145141919, 1145141919, 1145141919);
    node->BB = vec3(-1145141919, -1145141919, -1145141919);

    // 計算 AABB
    for (int i = l; i <= r; i++) {
        // 最小點 AA
        double minx = min(triangles[i].p1.x(), min(triangles[i].p2.x(), triangles[i].p3.x()));
        double miny = min(triangles[i].p1.y(), min(triangles[i].p2.y(), triangles[i].p3.y()));
        double minz = min(triangles[i].p1.z(), min(triangles[i].p2.z(), triangles[i].p3.z()));
        node->AA.e[0] = min(node->AA.x(), minx);
        node->AA.e[1] = min(node->AA.y(), miny);
        node->AA.e[2] = min(node->AA.z(), minz);
        // 最大點 BB
        double maxx = max(triangles[i].p1.x(), max(triangles[i].p2.x(), triangles[i].p3.x()));
        double maxy = max(triangles[i].p1.y(), max(triangles[i].p2.y(), triangles[i].p3.y()));
        double maxz = max(triangles[i].p1.z(), max(triangles[i].p2.z(), triangles[i].p3.z()));
        node->BB.e[0] = max(node->BB.x(), maxx);
        node->BB.e[1] = max(node->BB.x(), maxy);
        node->BB.e[2] = max(node->BB.x(), maxz);
    }

    // 不多於 n 個三角形 返回葉子節點
    if ((r - l + 1) <= n) {
        node->n = r - l + 1;
        node->index = l;
        return node;
    }

    // 否则遞迴建樹
    double Cost = INF;
    int Axis = 0;
    int Split = (l + r) / 2;
    for (int axis = 0; axis < 3; axis++) {
        // 分别按 x，y，z 軸排序
        if (axis == 0) std::sort(&triangles[0] + l, &triangles[0] + r + 1, cmpx);
        if (axis == 1) std::sort(&triangles[0] + l, &triangles[0] + r + 1, cmpy);
        if (axis == 2) std::sort(&triangles[0] + l, &triangles[0] + r + 1, cmpz);

        // leftMax[i]: [l, i] 中最大的 xyz 值
        // leftMin[i]: [l, i] 中最小的 xyz 值
        std::vector<vec3> leftMax(r - l + 1, vec3(-INF, -INF, -INF));
        std::vector<vec3> leftMin(r - l + 1, vec3(INF, INF, INF));
        // 计算前缀 注意 i-l 以对齐到下标 0
        for (int i = l; i <= r; i++) {
            Triangle& t = triangles[i];
            int bias = (i == l) ? 0 : 1;  // 第一个元素特殊处理

            leftMax[i - l].e[0] = max(leftMax[i - l - bias].x(), max(t.p1.x(), max(t.p2.x(), t.p3.x())));
            leftMax[i - l].e[1] = max(leftMax[i - l - bias].y(), max(t.p1.y(), max(t.p2.y(), t.p3.y())));
            leftMax[i - l].e[2] = max(leftMax[i - l - bias].z(), max(t.p1.z(), max(t.p2.z(), t.p3.z())));

            leftMin[i - l].e[0] = min(leftMin[i - l - bias].x(), min(t.p1.x(), min(t.p2.x(), t.p3.x())));
            leftMin[i - l].e[1] = min(leftMin[i - l - bias].y(), min(t.p1.y(), min(t.p2.y(), t.p3.y())));
            leftMin[i - l].e[2] = min(leftMin[i - l - bias].z(), min(t.p1.z(), min(t.p2.z(), t.p3.z())));
        }

        // rightMax[i]: [i, r] 中最大的 xyz 值
        // rightMin[i]: [i, r] 中最小的 xyz 值
        std::vector<vec3> rightMax(r - l + 1, vec3(-INF, -INF, -INF));
        std::vector<vec3> rightMin(r - l + 1, vec3(INF, INF, INF));
        // 计算后缀 注意 i-l 以对齐到下标 0
        for (int i = r; i >= l; i--) {
            Triangle& t = triangles[i];
            int bias = (i == r) ? 0 : 1;  // 第一个元素特殊处理

            rightMax[i - l].e[0] = max(rightMax[i - l + bias].x(), max(t.p1.x(), max(t.p2.x(), t.p3.x())));
            rightMax[i - l].e[1] = max(rightMax[i - l + bias].y(), max(t.p1.y(), max(t.p2.y(), t.p3.y())));
            rightMax[i - l].e[2] = max(rightMax[i - l + bias].z(), max(t.p1.z(), max(t.p2.z(), t.p3.z())));

            rightMin[i - l].e[0] = min(rightMin[i - l + bias].x(), min(t.p1.x(), min(t.p2.x(), t.p3.x())));
            rightMin[i - l].e[1] = min(rightMin[i - l + bias].y(), min(t.p1.y(), min(t.p2.y(), t.p3.y())));
            rightMin[i - l].e[2] = min(rightMin[i - l + bias].z(), min(t.p1.z(), min(t.p2.z(), t.p3.z())));
        }

        // 遍历寻找分割
        float cost = INF;
        int split = l;
        for (int i = l; i <= r-1; i++) {
            float lenx, leny, lenz;
            // 左侧 [l, i]
            vec3 leftAA = leftMin[i - l];
            vec3 leftBB = leftMax[i - l];
            lenx = leftBB.x() - leftAA.x();
            leny = leftBB.y() - leftAA.y();
            lenz = leftBB.z() - leftAA.z();
            float leftS = 2.0 * ((lenx * leny) + (lenx * lenz) + (leny * lenz));
            float leftCost = leftS * (i - l + 1);

            // 右侧 [i+1, r]
            vec3 rightAA = rightMin[i + 1 - l];
            vec3 rightBB = rightMax[i + 1 - l];
            lenx = rightBB.x() - rightAA.x();
            leny = rightBB.y() - rightAA.y();
            lenz = rightBB.z() - rightAA.z();
            float rightS = 2.0 * ((lenx * leny) + (lenx * lenz) + (leny * lenz));
            float rightCost = rightS * (r - i);

            // 记录每个分割的最小答案
            float totalCost = leftCost + rightCost;
            if (totalCost < cost) {
                cost = totalCost;
                split = i;
            }
        }
        // 记录每个轴的最佳答案
        if (cost < Cost) {
            Cost = cost;
            Axis = axis;
            Split = split;
        }
    }
    
    // 按最佳轴分割
    if (Axis == 0) std::sort(&triangles[0] + l, &triangles[0] + r + 1, cmpx);
    if (Axis == 1) std::sort(&triangles[0] + l, &triangles[0] + r + 1, cmpy);
    if (Axis == 2) std::sort(&triangles[0] + l, &triangles[0] + r + 1, cmpz);

    // 递归
    node->left = buildBVHwithSAH(triangles, l, Split, n);
    node->right = buildBVHwithSAH(triangles, Split+1, r, n);

    return node;
}

float hitAABB(const ray &r, vec3 AA, vec3 BB) {
    // 1.0 / direction
    vec3 invdir = vec3(1.0 / r.dir.x(), 1.0 / r.dir.y(), 1.0 / r.dir.z());

    vec3 in = (BB - r.orig) * invdir;
    vec3 out = (AA - r.orig) * invdir;

    vec3 tmax = max(in, out);
    vec3 tmin = min(in, out);

    float t1 = min(tmax.x(), min(tmax.y(), tmax.z()));
    float t0 = max(tmin.x(), max(tmin.y(), tmin.z()));

    return (t1 >= t0) ? ((t0 > 0.0) ? (t0) : (t1)) : (-1);
}

// 求交結果
struct HitResult {
    Triangle* triangle = NULL;
    float distance = INF;
};

float hitTriangle(const point3& p1,const point3& p2,const point3& p3,const ray& r){
    vec3 S = r.origin();        // 射线起点
    vec3 d = r.direction();    // 射线方向
    vec3 outward_normal = cross(p2 - p1, p3 - p1);
    bool front_face = dot(r.direction(), outward_normal) < 0; //判斷法線指向外為1, 內為0
    vec3 N = unit_vector(front_face ? outward_normal :-outward_normal);  // 對法線加+-號

    //vec3 N = normalize();    // 法向量
    if (dot(N, d) > 0.0f) N = -N;   // 获取正确的法向量

    // 如果视线和三角形平行
    if (fabs(dot(N, d)) < 0.00001f) return INF;

    // 距离
    float t = (dot(N, p1) - dot(S, N)) / dot(d, N);
    if (t < 0.0005f) return INF;    // 如果三角形在光线背面

    // 交点计算
    vec3 P = S + d * t;

    // 判断交点是否在三角形中
    vec3 c1 = cross(p2 - p1, P - p1);
    vec3 c2 = cross(p3 - p2, P - p2);
    vec3 c3 = cross(p1 - p3, P - p3);
    if (dot(c1, N) > 0 && dot(c2, N) > 0 && dot(c3, N) > 0) return t;
    if (dot(c1, N) < 0 && dot(c2, N) < 0 && dot(c3, N) < 0) return t;

    return INF;
}

HitResult hitTriangleArray(const ray& ray, vector<Triangle>& triangles, int l, int r) {
    HitResult res;
    for (int i = l; i <= r; i++) {
        float d = hitTriangle(triangles[i].p1,triangles[i].p2,triangles[i].p3, ray);
        if (d < INF && d < res.distance) {
            res.distance = d;
            res.triangle = &triangles[i];
        }
    }
    return res;
}

// 在 BVH 上遍历求交
HitResult hitBVH(const ray& ray, std::vector<Triangle>& triangles, BVHNode* root) {
    if (root == NULL) return HitResult();

    // 如果是葉子 暴力搜索
    if (root->n > 0) {
        //cout << root->n << " " << root->n + root->index - 1 << endl;
        return hitTriangleArray(ray, triangles, root->n, root->n + root->index - 1);
    }

    // 和左右子樹 AABB 求交
    float d1 = INF, d2 = INF;
    if (root->left) d1 = hitAABB(ray, root->left->AA, root->left->BB);
    if (root->right) d2 = hitAABB(ray, root->right->AA, root->right->BB);

    // 遞迴結果
    HitResult r1, r2;
    if (d1>0) r1 = hitBVH(ray, triangles, root->left);
    if (d2>0) r2 = hitBVH(ray, triangles, root->right);
    
    return r1.distance < r2.distance ? r1 : r2;
}

int hitcount = 0;
color bvh_ray_color(const ray &r, BVHNode *root)
{
    float s_t;
    float color_r = 0;
    float color_g = 0;
    float color_b = 0;
    for (int i = 0; i < centers.size(); i++)
    {
        cout << "wtf" << endl;
        s_t = hit_sphere(centers[i], radius[i], r);
        if (s_t > 0.0){
            vec3 N = unit_vector(r.at(s_t) - centers[i]);
            vec3 L = unit_vector(light_pos - r.at(s_t));
            vec3 V = unit_vector(-r.direction());
            vec3 H = unit_vector((L+V) / 2);
            float diffuse = Diffuse_Compont(N, L);
            float specular = Phong_model(N, H,Exp[0]);
            cout << diffuse << " " << specular << " " << endl;
            color_r = R[0] * (ka[0] + kd[0] * diffuse + ks[0] * specular);
            color_g = G[0] * (ka[0] + kd[0] * diffuse + ks[0] * specular);
            color_b = B[0] * (ka[0] + kd[0] * diffuse + ks[0] * specular);
            if(color_r < 0)
                color_r = 0;
            if (color_g < 0)
                color_g = 0;
            if(color_b < 0)
               color_b = 0;

            if(color_r > 1)
                color_r = 1;
            if (color_g > 1)
                color_g = 1;
            if(color_b > 1)
                color_b = 1;
            //return color(r,g,b);
            auto ref_vec = reflect(r.at(s_t), N);
            ray rl(r.at(s_t), ref_vec);
            for (int i = 0; i < t1.size(); i++)
            {
                float t = hit_triangle(t1[i], t2[i], t3[i], rl);
                if (t > 0.0){
                    color_r *= 1-Ref[0];
                    color_g *= 1-Ref[0];
                    color_b *= 1-Ref[0];
                }
            }
        }
        
    }
    
    //float t = hit_triangle(t1[i], t2[i], t3[i], r);
    HitResult res;
    res = hitBVH(r, triangles, root);
    if (res.triangle != NULL)
    {
        //cout << "hit count" << hitcount << endl;
        //hitcount++;
        // cout << "runtriangle" << endl;
        point3 tri1 = res.triangle->p1;
        point3 tri2 = res.triangle->p2;
        point3 tri3 = res.triangle->p3;
        auto v1 = vec3(tri2.x() - tri1.x(), tri2.y() - tri1.y(), tri2.z() - tri1.z());
        auto v2 = vec3(tri3.x() - tri1.x(),tri3.y() - tri1.y(),tri3.z() - tri1.z());

        vec3 outward_normal = cross(v1, v2);
        bool front_face = dot(r.dir, outward_normal) < 0; //判斷法線指向外為1, 內為0
        vec3 N = unit_vector(front_face ? outward_normal :-outward_normal);  // 對法線加+-號

        vec3 L = unit_vector(light_pos - r.at(res.distance));
        vec3 V = unit_vector(-r.direction());
        vec3 H = unit_vector((L+V) / 2);
        float diffuse = Diffuse_Compont(N, L);
        float specular = Phong_model(N, H,Exp[0]);
        //cout << diffuse << " " << specular << " " << endl;
        color_r = R[0] * (ka[0] + kd[0] * diffuse + ks[0] * specular);
        color_g = G[0] * (ka[0] + kd[0] * diffuse + ks[0] * specular);
        color_b = B[0] * (ka[0] + kd[0] * diffuse + ks[0] * specular);
        if(color_r < 0)
            color_r = 0;
        if (color_g < 0)
            color_g = 0;
        if(color_b < 0)
            color_b = 0;

        if(color_r > 1)
            color_r = 1;
        if (color_g > 1)
            color_g = 1;
        if(color_b > 1)
            color_b = 1;
        //return color(color_r,color_g,color_b);
        //draw shadow
        ray l(light_pos, r.at(res.distance) - light_pos);
        for (int i = 0; i < centers.size(); i++)
        {
            s_t = hit_sphere(centers[i], radius[i], l);
            if (s_t > 0.0)

                return color(0, 0, 0);
        }

        //draw reflect
        auto ref_vec = reflect(r.at(res.distance)-E , N);
        ray rl(r.at(res.distance), ref_vec);
        //cout << dot(r.at(t), N) << " " << dot(ref_vec, N) << " " << N.length() << endl;
        for (int i = 0; i < centers.size(); i++)
        {
            s_t = hit_sphere(centers[i], radius[i], rl);
            if (s_t > 0.0){
                color_r *= 1-Ref[0];
                color_g *= 1-Ref[0];
                color_b *= 1-Ref[0];
            }
        }
        
    }
    
    return color(color_r, color_g, color_b);
}
//-----------------------------------------------BVH END---------------------------------------------------

int main()
{
    readfile();
    // Camera

    //get_observe();//** O(UL UR LL LR)平面的長寬
    //-----------------------------------------------------------read obj file-------------------------------------------------------------------------
    std::cout<<"Load OBJ file!\n";
   
    std::vector<V3> vertices; //v 資料索引
    std::vector<V2> uvs;      //vt 資料索引
    std::vector<V3> normals;  //vn 資料索引
    std::vector<int> vertIdxs, uvIdxs, nrmIdxs; //f 資料索引

    FILE* file = fopen(FILENAME, "r");
    if (!file){
        printf("File open failded\n");
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
                fscanf(file, "%f %f\n", &uv.x, &uv.y);
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

            }else{
                //不需要的資料就跳過
                char stupidBuffer[1000];
                fgets(stupidBuffer, 1000, file);
            }
        }
        fclose(file);
        printf("file close\n");

        int total_vert_num = vertices.size();
        //printf("total_vert_num:%d,%f\n", total_vert_num,double(vertIdxs.size()/3));
        // 871306
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

            //triangel_data = std::make_shared<Triangle> (point3(vert1.x, vert1.y, vert1.z), point3(vert2.x, vert2.y, vert2.z), point3(vert3.x, vert3.y, vert3.z));
            //triangel_data = std::make_shared<triangle> (point3(vert1.x*20, -vert1.y*20, vert1.z*0), point3(vert2.x*-20, vert2.y*-20, vert2.z*30), point3(vert3.x*20, vert3.y*-20, vert3.z*-30), material_data);
            //world.add(triangel_data);
            triangles.push_back(Triangle(point3(vert1.x, vert1.y, vert1.z), point3(vert2.x, vert2.y, vert2.z), point3(vert3.x, vert3.y, vert3.z)));
        }
        printf("load obj file done!\n\n");
    }
    //--------------------------------------------build bvh----------------------------------------
    time_t st, et;
    st = time(NULL);
    
    cout << "Start Build BVH" << endl;
    BVHNode *root = buildBVHwithSAH(triangles, 0, triangles.size() - 1, 1); // box.obj
    //BVHNode *root = buildBVHwithSAH(triangles, 0, triangles.size() - 1, 20); //dragon.obj
    et = time(NULL) - st;
    cout << "triangles have :" << triangles.size() << endl;
    cout << "Build BVH complete, cost time: " << et << " \nnode count: " << treecount << endl;
    //---------------------------------------------------------------------------------------------
    int per = 1;

    ColorImage image;
	Pixel p={0,0,0};

	image.init(image_width, image_height);	

    double focal_length = 20;  //(dist_to_focus)
    auto aperture = 1.5/2.2;
    double viewport_screen_ratio = image_width / image_height;
    camera cam(E, D, U, angle, viewport_screen_ratio, aperture, focal_length);
    st = time(NULL);
    for (int j = 0; j < image_height; ++j)
    {
        //if(j > 55)
            //cout << j << " / " << image_height << endl;
        for (int i = 0; i < image_width; ++i)
        {
            color pixel_color;
            for (int s = 0; s < per;s++){
                auto u = double(i) / (image_width-1);
                auto v = double(j) / (image_height-1);
                ray r = cam.get_ray(u, v);

                //pixel_color += ray_color(r);
                pixel_color += bvh_ray_color(r,root);
            }
            p.R = static_cast<int>(255.999 * pixel_color.x()/per);
            p.G = static_cast<int>(255.999 * pixel_color.y()/per);
            p.B = static_cast<int>(255.999 * pixel_color.z()/per);
            image.writePixel(i, j, p);
            //if(j >= 62)
            //    cout << "   " << i << " / " << image_width << endl;
        }
    }
        et = time(NULL) - st;
    cout << "Ray Tracing complete, cost time: " << et/60 << endl;
    char output_name[15] = "output.ppm";
    image.outputPPM(output_name);
    std::cerr << "\nDone.\n";
    
    return 0;
}