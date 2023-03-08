#ifndef SPHERICAL_ISING_MODEL_MODEL_H
#define SPHERICAL_ISING_MODEL_MODEL_H

#include "Vec3.h"


typedef struct {
    Vec3 coord;
    Vec3 moment;
} Dipole;


typedef struct {
    Vec3 coord;
    unsigned char spin;
    Vec3 B;
    int close_6[6];
} Point;


typedef struct {
    int n_points;
    int n_dipoles;
    double energy;
    double mag;
    int evolve_steps;
    int step;
    int delta_checks;
    double T;
    Vec3 B;
    int output;
    int randomise;
    int field_type;
    Dipole* dipoles;
    Point* points;
} Model;


void distribute_points(Model*);

double energy(Model*);

double norm_mag(Model*);

double point_distance_squared(Point p1, Point p2);

void nns(Model*);

void randomise(Model*);

void evolve(Model*);

void output(Model*);

void B_from_dipoles(Model*);

void set_evolve(Model*);

#endif //SPHERICAL_ISING_MODEL_MODEL_H
