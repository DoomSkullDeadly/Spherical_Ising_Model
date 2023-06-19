#ifndef SPHERICAL_ISING_MODEL_MODEL_H
#define SPHERICAL_ISING_MODEL_MODEL_H

#include "Vec3.h"
#include <pthread.h>


typedef struct {
    Vec3 coord;
    Vec3 moment;
} Dipole;


typedef struct {
    Vec3 coord;
    float spin;
    Vec3 B;
    double direct_B;
    double dipole_B;
    int n_nns; // I fucking hate needing to have this, but fuck me ig
    int* nns;
} Point;


typedef struct {
    int thread_count;
    int n_points;
    int n_dipoles;
    double energy;
    double mag;
    int evolve_steps;
    int step;
    int delta_checks;
    double T;
    Vec3 B;
    int Ti;
    int Bi;
    int output;
    int randomise;
    double radius; // for spherical
    double point_spacing; // used for all, recalculates radius at n points to abide to this
    int length, width, height; // for rectangular and cuboidal lattices
    int lattice_type; // 1 = rectangular, 2 = cuboidal, 3 = spherical
    int field_type; // 1 = linear, 2 = dipole, 3 = both
    int working;
    Dipole* dipoles;
    Point* points;
    double** mags;
    double** energies;
} Model;


typedef struct {
    int thread_count;
    int to_run;
    pthread_t* tid;
    Model* models;
} Settings;


void distribute_points(Model*);

double energy(Model*);

double norm_mag(Model*);

double point_distance_squared(Point p1, Point p2);

void nns(Model*);

void randomise(Model*);

void evolve(Model*);

void *thread_evolve(void*);

void output(Model*);

void B_from_dipoles(Model*);

void set_evolve(Model*);

void free_Points(Model*);

void var_T(Model*, double, double, double, int);

void var_B(Model*, double, double, double, int);

void var_T_B(Model*, Settings*, double, double, double, double, double, int);

int find_thread(Settings*);

void copy_model(Model*, Model*);

void precalc_B(Model*);

Model* create_model();


#endif //SPHERICAL_ISING_MODEL_MODEL_H
