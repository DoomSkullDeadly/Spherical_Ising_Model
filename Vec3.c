#include "Vec3.h"
#include <math.h>


void vec_add(Vec3* v1, Vec3 v2) {
    v1->x += v2.x;
    v1->y += v2.y;
    v1->z += v2.z;
}


void vec_diff(Vec3* v1, Vec3 v2) {
    v1->x -= v2.x;
    v1->y -= v2.y;
    v1->z -= v2.z;
}


void vec_mult(Vec3* v1, double s) {
    v1->x *= s;
    v1->y *= s;
    v1->z *= s;
}


void vec_div(Vec3* v1, double s) {
    v1->x /= s;
    v1->y /= s;
    v1->z /= s;
}


double vec_dot_prod(Vec3 v1, Vec3 v2) {
    double dp = 0;
    dp += v1.x * v2.x;
    dp += v1.y * v2.y;
    dp += v1.z * v2.z;
    return dp;
}


double vec_mag(Vec3 v) {
    double sum = 0;
    sum += pow(v.x, 2);
    sum += pow(v.y, 2);
    sum += pow(v.z, 2);
    return sqrt(sum);
}