#include "model.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#define J 1.
#define mu_b 1.
#define mu_0 1.
#define k_b 1.


void distribute_points(Model* model) {
    if (model->lattice_type == 1) { // rectangular lattice

        model->n_points = model->length * model->width;
        model->points = (Point*) realloc(model->points, model->n_points * sizeof(Point));
        if (model->points == NULL) {
            printf("Error occurred allocating memory!\n");
            exit(0);
        }

        for (int y = 0; y < model->width; ++y) {
            for (int x = 0; x < model->length; ++x) {
                int index = model->length * y + x;
                model->points[index].coord.x = x * model->point_spacing;
                model->points[index].coord.y = y * model->point_spacing;
                model->points[index].coord.z = 0;

                model->points[index].n_nns = 4;
                model->points[index].nns = (int*) malloc(4 * sizeof(int)); //TODO: free memory cause realloc fucks this
                if (model->points[index].nns == NULL) {
                    printf("Error occurred allocating memory!\n");
                    exit(0);
                }
            }
        }
    }

    else if (model->lattice_type == 2) { // cuboidal lattice
        int i = 0;
        for (int z = 0; z < model->height; ++z) {
            for (int y = 0; y < model->width; ++y) {
                for (int x = 0; x < model->length; ++x) {
                    if (((z != 0 || z != model->height-1) && (x == 0 || y == 0 || x == model->length-1 || y == model->width-1)) || (z == 0 || z == model->height-1)) {

                        model->points = (Point*) realloc(model->points, (i+1) * sizeof(Point));
                        if (model->points == NULL) {
                            printf("Error occurred allocating memory!\n");
                            exit(0);
                        }

                        model->points[i].coord.x = x * model->point_spacing;
                        model->points[i].coord.y = y * model->point_spacing;
                        model->points[i].coord.z = z * model->point_spacing;

                        int cond_x = (x == 0 || x == model->length-1) ? 1 : 0;
                        int cond_y = (y == 0 || y == model->width-1) ? 1 : 0;
                        int cond_z = (z == 0 || z == model->height-1) ? 1 : 0;
                        int pnns = (cond_x + cond_y + cond_z == 3) ? 3 : 4;  // checks if a corner point
                        model->points[i].n_nns = pnns;
                        model->points[i].nns = (int*) malloc(pnns * sizeof(int)); //TODO: free memory cause realloc fucks it pt.2
                        if (model->points[i].nns == NULL) {
                            printf("Error occurred allocating memory!\n");
                            exit(0);
                        }

                        i++;
                    }
                }
            }
        }
        model->n_points = i;
    }

    else if (model->lattice_type == 3) { // spherical lattice

        model->points = (Point*) realloc(model->points, model->n_points * sizeof(Point));
        if (model->points == NULL) {
            printf("Error occurred allocating memory!\n");
            exit(0);
        }

        double start = (-1. + 1. / (model->n_points - 1.));  //TODO: credit alg
        double increment = (2. - 2. / (model->n_points - 1.)) / (model->n_points - 1.);

        for (int i = 0; i < model->n_points; ++i) {
            double s = start + i * increment;
            double X = s * (0.1 + 1.2 * model->n_points);
            double Y = M_PI / 2. * copysign(1, s) * (1. - sqrt(1. - fabs(s)));
            model->points[i].coord.x = cos(X) * cos(Y);
            model->points[i].coord.y = sin(X) * cos(Y);
            model->points[i].coord.z = sin(Y);

            model->points[i].n_nns = 6;
            model->points[i].nns = (int*) realloc(model->points[i].nns, 6 * sizeof(int));
            if (model->points[i].nns == NULL) {
                printf("Error occurred allocating memory!\n");
                exit(0);
            }
        }
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

    if (model->lattice_type == 1) {
        for (int i = 0; i < model->n_points; ++i) {
            model->points[i].nns[0] = (model->n_points + i - 1) % model->n_points;
            model->points[i].nns[1] = (model->n_points + i + 1) % model->n_points;
            model->points[i].nns[2] = (model->n_points + i - model->length) % model->n_points;
            model->points[i].nns[3] = (model->n_points + i + model->length) % model->n_points;
        }
    }

    else if (model->lattice_type != 1) {
        int indexes[model->n_points];
        double distances[model->n_points];
        int pnns;

        for (int point = 0; point < model->n_points; ++point) {  // looping over all points
            printf("%i\n", point);
            for (int i = 0; i < model->n_points; ++i) {  // setting up index and distances arrays
                indexes[i] = i;
                distances[i] = point_distance_squared(model->points[point], model->points[i]);
            }

            for (int i = 1;
                 i < model->n_points; ++i) {  // insertion sort to get all indexes of points sorted by distance
                int key = indexes[i];
                int j = i - 1;
                while (j >= 0 && distances[key] < distances[indexes[j]]) {
                    indexes[j + 1] = indexes[j];
                    j -= 1;
                }
                indexes[j + 1] = key;
            }

            for (int k = 0; k < model->points[point].n_nns; ++k) {  // closest points are stored for each point, excluding the point itself [0]
                model->points[point].nns[k] = indexes[k + 1];  // number of points depends on the lattice type, and point position.
            }
        }
    }
}


double norm_mag(Model* model) {
    double M = 0;
    for (int point = 0; point < model->n_points; ++point) {
        M += model->points[point].spin -.5;
    }
    M *= 2. / model->n_points;
    return M;
}


double energy(Model* model) {
    double E = 0;
    for (int point = 0; point < model->n_points; ++point) {
        double neighbour_sum = 0;

        for (int i = 0; i < model->points[point].n_nns; ++i) {
            neighbour_sum += ((double)model->points[point].spin -.5) * ((double)model->points[model->points[point].nns[i]].spin -.5);
        }
        E -= J * neighbour_sum;

        double B; // the following could all be done by using nns to find plane and thus normal to point to use with dot prod but that would probably be more computationally taxing.
        if (model->lattice_type == 1) {
            B = model->B.z + model->points[point].B.z;
            E -= mu_b * B * ((double)model->points[point].spin - .5);
        }

        else if (model->lattice_type == 2) {
            int cond_x = (model->points[point].coord.x == 0 || model->points[point].coord.x == model->length-1) ? 1 : 0;
            int cond_y = (model->points[point].coord.y == 0 || model->points[point].coord.y == model->width-1) ? 1 : 0;
            int cond_z = (model->points[point].coord.z == 0 || model->points[point].coord.z == model->height-1) ? 1 : 0;
            int cond_sum = cond_x + cond_y + cond_z;
            Vec3 condition = {copysign(1, model->points[point].coord.x-cond_x) * cond_x / sqrt(cond_sum),
                              copysign(1, model->points[point].coord.y-cond_y) * cond_y / sqrt(cond_sum),
                              copysign(1, model->points[point].coord.z-cond_z) * cond_z / sqrt(cond_sum)};
            B = vec_dot_prod(condition, model->points[point].B);
            B += vec_dot_prod(condition, model->B);
            E -= mu_b * B * (double)model->points[point].spin - .5;
        }

        else if (model->lattice_type == 3) {
            B = vec_dot_prod(model->points[point].coord, model->points[point].B);
            B += vec_dot_prod(model->points[point].coord, model->B);
            E -= mu_b * B / vec_mag(model->points[point].coord) * (double)model->points[point].spin - .5;
        }
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

            if (delta_E > 0 && (double) (rand() % 100000) / 100000 > exp(-delta_E / (k_b * model->T))) { // set back to original
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