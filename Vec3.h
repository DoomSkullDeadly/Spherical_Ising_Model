#ifndef SPHERICAL_ISING_MODEL_VEC3_H
#define SPHERICAL_ISING_MODEL_VEC3_H


typedef struct {
    double x;
    double y;
    double z;
}Vec3;


void vec_add(Vec3*, Vec3);

void vec_diff(Vec3*, Vec3);

void vec_mult(Vec3*, double);

void vec_div(Vec3*, double);

double vec_dot_prod(Vec3 v1, Vec3 v2);

double vec_mag(Vec3);


#endif //SPHERICAL_ISING_MODEL_VEC3_H