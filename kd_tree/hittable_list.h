#ifndef HITTABLE_LIST_H
#define HITTABLE_LIST_H

#include "hittable.h"

#include <memory>
#include <vector>

using std::shared_ptr;
using std::make_shared;

struct List{
    int index;
    struct List* pnext;
};

struct Node{
    struct List list;
    int axis;
    double v;
    struct Node* left;
    struct Node* right;
    double x_min;
    double x_max;
    double y_min;
    double y_max;
    double z_min;
    double z_max;
};

class hittable_list : public hittable {
    public:
        hittable_list() {}
        hittable_list(shared_ptr<hittable> object) { add(object); }

        void clear() { 
            objects.clear(); 
            clear_kd_tree();
        }
        void add(shared_ptr<hittable> object) { objects.push_back(object); }

        virtual bool hit(
            const ray& r, double t_min, double t_max, hit_record& rec) const override;
        virtual char get_kind(){ return '0';};
        virtual point3 get_center(){return point3(0, 0, 0);};
        virtual void get_outer_bound(point3* p){};
        
        bool hit_node(const struct Node *node, const ray& r, double t_min, double t_max, hit_record& rec) const;
        void make_kd_tree();
        void kd_tree(struct Node *node);
        void clear_kd_tree();
        void clear_tree(struct Node* node);
        void list_add(struct List *pnode, int index);
        bool check(const ray& r, struct Node *node, double t_min, double t_max) const;
        void free_list(struct List* pnode);

    public:
        std::vector<shared_ptr<hittable>> objects;
        struct Node root;
};

bool hittable_list::check(const ray& r, struct Node *node, double t_min, double t_max) const {
    double t;
    if(r.Eye_position.x()>=node->x_min && r.Eye_position.x()<=node->x_max && r.Eye_position.y()>=node->y_min && r.Eye_position.y()<=node->y_max && r.Eye_position.z()>=node->z_min && r.Eye_position.z()<=node->z_max){
        return true;
    }
    if(dot(r.Eye_direction, vec3(1, 0, 0))>0){
        t = (node->x_min-r.Eye_position.x())/r.Eye_direction.x();
        if(t>t_min && t<t_max && r.at(t).y()>=node->y_min && r.at(t).y()<=node->y_max && r.at(t).z()>=node->z_min && r.at(t).z()<=node->z_max){
            return true;
        }
    }else{
        t = (node->x_max-r.Eye_position.x())/r.Eye_direction.x();
        if(t>t_min && t<t_max && r.at(t).y()>=node->y_min && r.at(t).y()<=node->y_max && r.at(t).z()>=node->z_min && r.at(t).z()<=node->z_max){
            return true;
        }
    }
    if(dot(r.Eye_direction, vec3(0, 1, 0))>0){
        t = (node->y_min-r.Eye_position.y())/r.Eye_direction.y();
        if(t>t_min && t<t_max && r.at(t).x()>=node->x_min && r.at(t).x()<=node->x_max && r.at(t).z()>=node->z_min && r.at(t).z()<=node->z_max){
            return true;
        }
    }else{
        t = (node->y_max-r.Eye_position.y())/r.Eye_direction.y();
        if(t>t_min && t<t_max && r.at(t).x()>=node->x_min && r.at(t).x()<=node->x_max && r.at(t).z()>=node->z_min && r.at(t).z()<=node->z_max){
            return true;
        }
    }
    if(dot(r.Eye_direction, vec3(0, 0, 1))>0){
        t = (node->z_min-r.Eye_position.z())/r.Eye_direction.z();
        if(t>t_min && t<t_max && r.at(t).x()>=node->x_min && r.at(t).x()<=node->x_max && r.at(t).y()>=node->y_min && r.at(t).y()<=node->y_max){
            return true;
        }
    }else{
        t = (node->z_max-r.Eye_position.z())/r.Eye_direction.z();
        if(t>t_min && t<t_max && r.at(t).x()>=node->x_min && r.at(t).x()<=node->x_max && r.at(t).y()>=node->y_min && r.at(t).y()<=node->y_max){
            return true;
        }
    }
    return false;
}

bool hittable_list::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    bool hit_anything = false;
    auto closest_so_far = t_max;

    hit_anything = hit_node(&root, r, t_min, closest_so_far, rec);
/*
    if(hit_anything){
        closest_so_far = rec.t*0.99;
        while(hit_node(&root, r, t_min, closest_so_far, rec)){
            closest_so_far = rec.t*0.99;
        }
    }
*/
    return hit_anything;
}

bool hittable_list::hit_node(const struct Node *node, const ray& r, double t_min, double t_max, hit_record& rec)const{
    hit_record temp_rec;
    bool hit_anything = false;
    auto closest_so_far = t_max;
    struct List *pnode;

    if(node->left==NULL && node->right==NULL){
        pnode = node->list.pnext;
        while(pnode){
            if (objects[pnode->index]->hit(r, t_min, closest_so_far, temp_rec)) {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
            }
            pnode = pnode->pnext;
        }
    }else{
        if(hit_anything) closest_so_far = rec.t;
        if(check(r, node->left, t_min, closest_so_far)){
            hit_anything |= hit_node(node->left, r, t_min, closest_so_far, rec);
        }
        if(hit_anything) closest_so_far = rec.t;
        if(check(r, node->right, t_min, closest_so_far)){
            hit_anything |= hit_node(node->right, r, t_min, closest_so_far, rec);
        }
    }

    return hit_anything;
}

void hittable_list::free_list(struct List *pnode){
    struct List *pfree=NULL;
    while(pnode->pnext){
        pfree = pnode->pnext;
        pnode->pnext = pfree->pnext;
        free(pfree);
    }
    pnode->index = 0;
    return ;
}

void hittable_list::kd_tree(struct Node *node){
    double max_x=-1000000, min_x=1000000, mean_x=0, max_y=-1000000, min_y=1000000, mean_y=0, max_z=-1000000, min_z=1000000, mean_z=0;
    int n = 0;
    point3 center, p[6];
    char kind;
    struct List *pnode=NULL, *pnew=NULL;

    pnode = node->list.pnext;
    while(pnode) {
        center = objects[pnode->index]->get_center();
        if(center.x()>max_x) max_x = center.x();
        if(center.x()<min_x) min_x = center.x();
        if(center.y()>max_y) max_y = center.y();
        if(center.y()<min_y) min_y = center.y();
        if(center.z()>max_z) max_z = center.z();
        if(center.z()<min_z) min_z = center.z();
        mean_x += center.x();
        mean_y += center.y();
        mean_z += center.z();
        n++;
        pnode = pnode->pnext;
    }

    mean_x = mean_x/n;
    mean_y = mean_y/n;
    mean_z = mean_z/n;

    node->left = (struct Node*)malloc(sizeof(struct Node));
    node->right = (struct Node*)malloc(sizeof(struct Node));
    node->left->left = NULL;
    node->left->right = NULL;
    node->left->list.index = 0;
    node->left->list.pnext = NULL;
    node->left->x_min = node->x_min;
    node->left->x_max = node->x_max;
    node->left->y_min = node->y_min;
    node->left->y_max = node->y_max;
    node->left->z_min = node->z_min;
    node->left->z_max = node->z_max;

    node->right->left = NULL;
    node->right->right = NULL;
    node->right->list.index = 0;
    node->right->list.pnext = NULL;
    node->right->x_min = node->x_min;
    node->right->x_max = node->x_max;
    node->right->y_min = node->y_min;
    node->right->y_max = node->y_max;
    node->right->z_min = node->z_min;
    node->right->z_max = node->z_max;

    while(true){
        if(max_x-min_x>=max_y-min_y && max_x-min_x>=max_z-min_z){
            node->axis = 0;
            node->left->x_max = mean_x;
            node->right->x_min = mean_x;
            node->v = mean_x;
            pnode = node->list.pnext;
            while(pnode) {
                kind = objects[pnode->index]->get_kind();
                objects[pnode->index]->get_outer_bound(&p[0]);
                if(kind=='s'){
                    if(p[0].x()<mean_x || p[1].x()<mean_x || p[2].x()<mean_x || p[3].x()<mean_x || p[4].x()<mean_x || p[5].x()<mean_x){
                        list_add(&node->left->list, pnode->index);
                    }
                    if(p[0].x()>=mean_x || p[1].x()>=mean_x || p[2].x()>=mean_x || p[3].x()>=mean_x || p[4].x()>=mean_x || p[5].x()>=mean_x){
                        list_add(&node->right->list, pnode->index);
                    }
                }else if(kind=='t'){
                    if(p[0].x()<mean_x || p[1].x()<mean_x || p[2].x()<mean_x){
                        list_add(&node->left->list, pnode->index);
                    }
                    if(p[0].x()>=mean_x || p[1].x()>=mean_x || p[2].x()>=mean_x){
                        list_add(&node->right->list, pnode->index);
                    }
                }
                pnode = pnode->pnext;
            }
        }else if(max_y-min_y>=max_x-min_x && max_y-min_y>=max_z-min_z){
            node->axis = 1;
            node->left->y_max = mean_y;
            node->right->y_min = mean_y;
            node->v = mean_y;
            pnode = node->list.pnext;
            while(pnode) {
                kind = objects[pnode->index]->get_kind();
                objects[pnode->index]->get_outer_bound(&p[0]);
                if(kind=='s'){
                    if(p[0].y()<mean_y || p[1].y()<mean_y || p[2].y()<mean_y || p[3].y()<mean_y || p[4].y()<mean_y || p[5].y()<mean_y){
                        list_add(&node->left->list, pnode->index);
                    }
                    if(p[0].y()>=mean_y || p[1].y()>=mean_y || p[2].y()>=mean_y || p[3].y()>=mean_y || p[4].y()>=mean_y || p[5].y()>=mean_y){
                        list_add(&node->right->list, pnode->index);
                    }
                }else if(kind=='t'){
                    if(p[0].y()<mean_y || p[1].y()<mean_y || p[2].y()<mean_y){
                        list_add(&node->left->list, pnode->index);
                    }
                    if(p[0].y()>=mean_y || p[1].y()>=mean_y || p[2].y()>=mean_y){
                        list_add(&node->right->list, pnode->index);
                    }
                }
                pnode = pnode->pnext;
            }
        }else if(max_z-min_z>=max_x-min_x && max_z-min_z>=max_y-min_y){
            node->axis = 2;
            node->left->z_max = mean_z;
            node->right->z_min = mean_z;
            node->v = mean_z;
            pnode = node->list.pnext;
            while(pnode) {
                kind = objects[pnode->index]->get_kind();
                objects[pnode->index]->get_outer_bound(&p[0]);
                if(kind=='s'){
                    if(p[0].z()<mean_z || p[1].z()<mean_z || p[2].z()<mean_z || p[3].z()<mean_z || p[4].z()<mean_z || p[5].z()<mean_z){
                        list_add(&node->left->list, pnode->index);
                    }
                    if(p[0].z()>=mean_z || p[1].z()>=mean_z || p[2].z()>=mean_z || p[3].z()>=mean_z || p[4].z()>=mean_z || p[5].z()>=mean_z){
                        list_add(&node->right->list, pnode->index);
                    }
                }else if(kind=='t'){
                    if(p[0].z()<mean_z || p[1].z()<mean_z || p[2].z()<mean_z){
                        list_add(&node->left->list, pnode->index);
                    }
                    if(p[0].z()>=mean_z || p[1].z()>=mean_z || p[2].z()>=mean_z){
                        list_add(&node->right->list, pnode->index);
                    }
                }
                pnode = pnode->pnext;
            }
        }
        if(node->left->list.index==node->list.index || node->right->list.index==node->list.index){
            free_list(&node->left->list);
            free_list(&node->right->list);
            if(node->axis==0) max_x = min_x;
            else if(node->axis==1) max_y = min_y;
            else if(node->axis==2) max_z = min_z;

            if(max_x-min_x==0 && max_y-min_y==0 && max_z-min_z==0){
                free(node->left);
                free(node->right);
                node->left = NULL;
                node->right = NULL;
                return ;
            }
        }else{
            break;
        }
    }

    //printf("%d->%d, %d\n", node->list.index, node->left->list.index, node->right->list.index);
    free_list(&node->list);
    //printf("%f, %f, %f, %f, %f, %f\n", node->x_min, node->x_max, node->y_min, node->y_max, node->z_min, node->z_max);
    if(node->left->list.index>32) kd_tree(node->left);
    if(node->right->list.index>32) kd_tree(node->right);

    return ;
}

void hittable_list::make_kd_tree(){
    root.left = NULL;
    root.right = NULL;
    root.list.index = 0;
    root.list.pnext = NULL;
    root.x_min = -infinity;
    root.x_max = infinity;
    root.y_min = -infinity;
    root.y_max = infinity;
    root.z_min = -infinity;
    root.z_max = infinity;
    if(objects.size()==0) return ;
    for(int i=objects.size()-1;i>=0;i--){
        list_add(&root.list, i);
    }
    //if(root.list.index<32) return ;
    kd_tree(&root);
    return ;
}

void hittable_list::clear_kd_tree(){
    if(root.left!=NULL){
        clear_tree(root.left);
    }
    if(root.right!=NULL){
        clear_tree(root.right);
    }
    struct List* pnode=NULL;
    while(root.list.pnext){
        pnode = root.list.pnext;
        root.list.pnext = pnode->pnext;
        free(pnode);
    }
    root.left = NULL;
    root.right = NULL;
    return ;
}

void hittable_list::clear_tree(struct Node* node){
    if(node->left!=NULL){
        clear_tree(node->left);
    }
    if(node->right!=NULL){
        clear_tree(node->right);
    }
    struct List *pnode=NULL;
    while(node->list.pnext){
        pnode = node->list.pnext;
        node->list.pnext = pnode->pnext;
        free(pnode);
    }
    free(node);
    return ;
}

void hittable_list::list_add(struct List *pnode, int index){
    struct List *pnew=NULL;
    pnode->index++;
    pnew = (struct List*)malloc(sizeof(struct List));
    pnew->pnext = pnode->pnext;
    pnew->index = index;
    pnode->pnext = pnew;
    return ;
}

#endif