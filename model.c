#include "model.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#define J 1.
#define mu_b 1.
#define mu_0 1.
#define k_b 1.


void distribute_points(Model* model) {
    double start = (-1. + 1. / (model->n_points - 1.));
    double increment = (2. - 2. / (model->n_points - 1.)) / (model->n_points - 1.);

    for (int i = 0; i < model->n_points; ++i) {
        double s = start + i * increment;
        double X = s * (0.1 + 1.2 * model->n_points);
        double Y = M_PI / 2. * copysign(1, s) * (1. - sqrt(1. - fabs(s)));
        model->points[i].coord.x = cos(X) * cos(Y);
        model->points[i].coord.y = sin(X) * cos(Y);
        model->points[i].coord.z = sin(Y);
    }
}


void randomise(Model* model) {
    for (int i = 0; i < model->n_points; ++i) {
        model->points[i].spin = rand() % 2;
    }
}


double point_distance_squared(Point p1, Point p2) {
    double sum = 0;
    sum += pow(p1.coord.x - p2.coord.x, 2);
    sum += pow(p1.coord.y - p2.coord.y, 2);
    sum += pow(p1.coord.z - p2.coord.z, 2);
    return sum;
}


void nns(Model* model) {
    int indexes[model->n_points];
    double distances[model->n_points];

    for (int point = 0; point < model->n_points; ++point) {  // looping over all points
        printf("%i\n", point);
        for (int i = 0; i < model->n_points; ++i) {  // setting up index and distances arrays
            indexes[i] = i;
            distances[i] = point_distance_squared(model->points[point], model->points[i]);
        }

        for (int i = 1; i < model->n_points; ++i) {  // insertion sort to get all indexes of points sorted by distance
            int key = indexes[i];
            int j = i-1;
            while (j >= 0 && distances[key] < distances[indexes[j]]) {
                indexes[j+1] = indexes[j];
                j -= 1;
            }
            indexes[j+1] = key;
        }
        for (int k = 0; k < 6; ++k) {  // closest 6 points are stored for each point, excluding the point itself [0]
            model->points[point].close_6[k] = indexes[k+1];
        }
    }
}


double norm_mag(Model* model) {
    double M = 0;
    for (int point = 0; point < model->n_points; ++point) {
        M += model->points[point].spin;
    }
    M *= 2. / model->n_points;
    return M;
}


double energy(Model* model) {
    double E = 0;
    for (int point = 0; point < model->n_points; ++point) {
        double neighbour_sum = 0;
        for (int i = 0; i < 6; ++i) {
            neighbour_sum -= (double)model->points[point].spin-.5 * (double)model->points[model->points[point].close_6[i]].spin-.5;
        }
        E -= J * neighbour_sum;

        double p_dot_b = 0;
        if (model->field_type == 1) {
            p_dot_b = vec_dot_prod(model->points[point].coord, model->B);
        }
        else if (model->field_type == 2) {
            p_dot_b = vec_dot_prod(model->points[point].coord, model->points[point].B);
        }

        E -= mu_b * p_dot_b / vec_mag(model->points[point].coord) * (double)model->points[point].spin-.5;
    }
    return E;
}


void evolve(Model* model) {
    int running = 1;
    if (model->output) {
        output(model);
    }
    double d_E[model->delta_checks];
    while (running) {
        model->step++;
        for (int i = 0; i < model->n_points; ++i) {  // single step is considered as n points tested
            int point = rand() % model->n_points;
            int current_spin = model->points[point].spin;
            model->points[point].spin = current_spin ? 0 : 1;  // sets spin to opposite

            double new_E = energy(model);
            double delta_E = new_E - model->energy;

            if (delta_E > 0 && (float) (rand() % 100000) / 100000 > exp(-delta_E / (k_b * model->T))) { // set back to original
                model->points[point].spin = current_spin;
            }
            else {
                model->energy = new_E;
            }
        }

        if (model->step == // if the model is set to run for a specified amount of steps it'll stop when those are reached
            model->evolve_steps && model->evolve_steps) {
            running = 0;
        }
        else if (!model->evolve_steps) {
            d_E[model->step % model->delta_checks] = model->energy; // sets current pos to current energy
            if (fabs(d_E[(model->step+1) % model->delta_checks] - model->energy) < fabs(model->energy * 0.001)) {
                running = 0; // if the change in energy over the last n steps is less than 0.1% of the current energy it will stop running the model. This seems to work reliably.
            }
        }
        if (model->output) {
            output(model);
        }
    }
}


void set_evolve(Model* model) { // performs evolution once
    model->energy = energy(model);
    model->mag = norm_mag(model);
    printf("E = %g, M = %g\n", model->energy, model->mag);

    evolve(model);

    model->energy = energy(model);
    model->mag = norm_mag(model);
    printf("E = %g, M = %g\n", model->energy, model->mag);
}


void output(Model* model) {
    FILE *file;
    char buf[12+(int)(model->step/10)];
    sprintf(buf, "output/%i.txt", model->step);
    file = fopen(buf, "w");

    for (int i = 0; i < model->n_points; ++i) {
        fprintf(file, "%g %g %g %i\n", model->points[i].coord.x, model->points[i].coord.y, model->points[i].coord.z, model->points[i].spin);
    }
    fclose(file);
}


void B_from_dipoles(Model* model) {
    for (int point = 0; point < model->n_points; ++point) {
        printf("%i\n", point);
        model->points[point].B.x = 0.;
        model->points[point].B.y = 0.;
        model->points[point].B.z = 0.;
        for (int i = 0; i < model->n_dipoles; ++i) {
            Vec3 r = {0, 0, 0};
            vec_add(&r, model->points[point].coord);
            vec_diff(&r, model->dipoles[i].coord);
            double mag_r = vec_mag(r);
            double m_dot_r = vec_dot_prod(model->dipoles[i].moment, r);

            vec_add(&model->points[point].B, r);
            vec_mult(&model->points[point].B, 3 * m_dot_r);
            vec_div(&model->points[point].B, pow(mag_r, 5));

            Vec3 half2 = {0, 0, 0};
            vec_add(&half2, model->dipoles[i].moment);
            vec_div(&half2, pow(mag_r, 3));
            vec_diff(&model->points[point].B, half2);
        }
        vec_mult(&model->points[point].B, mu_0/(4 * M_PI));
    }
}