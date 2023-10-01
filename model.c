#include "model.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

#define J 1.
#define mu_b 1.
#define mu_0 1.
#define k_b 1.
#define spin_mag 1.


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
                model->points[index].nns = (int*) malloc(4 * sizeof(int));
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
                        int p_nns = (cond_x + cond_y + cond_z == 3) ? 3 : 4;  // checks if a corner point
                        model->points[i].n_nns = p_nns;
                        model->points[i].nns = (int*) malloc(p_nns * sizeof(int));
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
        model->points[i].spin = rand() % 2 ? spin_mag : -spin_mag; // NOLINT
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

        for (int point = 0; point < model->n_points; ++point) {  // looping over all points
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
        M += model->points[point].spin;
    }
    M /= model->n_points / spin_mag;
    return fabs(M);
}


double energy(Model* model) {
    double E = 0;
    for (int point = 0; point < model->n_points; ++point) {
        double neighbour_sum = 0;

        for (int i = 0; i < model->points[point].n_nns; ++i) {
            neighbour_sum += ((double)model->points[point].spin) * ((double)model->points[model->points[point].nns[i]].spin);
        }
        E -= J * neighbour_sum / 2;

        double B = 0;
        if (model->field_type != 2) {
            B += model->points[point].direct_B;
        }
        if (model->field_type != 1) {
            B += model->points[point].dipole_B;
        }

        if (model->lattice_type == 3) {
            E -= mu_b * B / vec_mag(model->points[point].coord) * (double)model->points[point].spin;
        }
        else {
            E -= mu_b * B * (double)model->points[point].spin;
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
            int point = rand() % model->n_points; // NOLINT
            float current_spin = model->points[point].spin;
            model->points[point].spin *= -1;  // sets spin to opposite

            double new_E = energy(model);
            double delta_E = new_E - model->energy;

            if (delta_E > 0 && (double) (rand() % 100000) / 100000 > exp(-delta_E / (k_b * model->T))) { // set back to original // NOLINT
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

    evolve(model);

    model->energy = energy(model);
    model->mag = norm_mag(model);
}


void* thread_evolve(void *m) {
    Model* model = (Model*)m;
    set_evolve(model);

    model->mags[model->Ti][model->Bi] += model->mag;
    model->energies[model->Ti][model->Bi] += model->energy;
    model->m_sq[model->Ti][model->Bi] += model->mag * model->mag;
    model->e_sq[model->Ti][model->Bi] += model->energy * model->energy;

    for (int i = 0; i < model->n_points; ++i) {
        free(model->points[i].nns);
    }
    free(model->points);
    free(model->dipoles);

    model->working = 0;
    return NULL;
}


void output(Model* model) {
    FILE *file;
    char buf[12+(int)(model->step/10)];
    sprintf(buf, "output/%i.txt", model->step);
    file = fopen(buf, "w");

    for (int i = 0; i < model->n_points; ++i) {
        fprintf(file, "%g %g %g %g\n", model->points[i].coord.x, model->points[i].coord.y, model->points[i].coord.z, model->points[i].spin);
    }
    fclose(file);
}


void B_from_dipoles(Model* model) {
    for (int point = 0; point < model->n_points; ++point) {
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


void precalc_B(Model* model) {
    for (int point = 0; point < model->n_points; ++point) { // the following could all be done by using nns to find plane and thus normal to point to use with dot prod but that would probably be more computationally taxing.
        if (model->lattice_type == 1) {
            model->points[point].direct_B = model->B.z;
            model->points[point].dipole_B = model->points[point].B.z;
        }

        else if (model->lattice_type == 2) {
            int cond_x = (model->points[point].coord.x == 0 || model->points[point].coord.x == model->length-1) ? 1 : 0;
            int cond_y = (model->points[point].coord.y == 0 || model->points[point].coord.y == model->width-1) ? 1 : 0;
            int cond_z = (model->points[point].coord.z == 0 || model->points[point].coord.z == model->height-1) ? 1 : 0;
            int cond_sum = cond_x + cond_y + cond_z;
            Vec3 condition = {copysign(1, model->points[point].coord.x-cond_x) * cond_x / sqrt(cond_sum),
                              copysign(1, model->points[point].coord.y-cond_y) * cond_y / sqrt(cond_sum),
                              copysign(1, model->points[point].coord.z-cond_z) * cond_z / sqrt(cond_sum)};

            model->points[point].direct_B = vec_dot_prod(condition, model->B);
            model->points[point].dipole_B = vec_dot_prod(condition, model->points[point].B);
        }

        else if (model->lattice_type == 3) {
            model->points[point].direct_B = vec_dot_prod(model->points[point].coord, model->B);
            model->points[point].dipole_B = vec_dot_prod(model->points[point].coord, model->points[point].B);
        }
    }
}


void free_Points(Model* model) {
    for (int i = 0; i < model->n_points; ++i) {
        free(model->points[i].nns);
    }
    free(model->points);
}


void var_T(Model* model, double start, double end, double increment, int repeats) {
    distribute_points(model);
    nns(model);
    B_from_dipoles(model);
    precalc_B(model);
    int arr_length = (int)((end - start) / increment) + 1;
    double* mags = (double*)calloc(arr_length, sizeof(double));
    double* energies = (double*)calloc(arr_length, sizeof(double));
    for (int i = 0; i < repeats; ++i) {
        for (int j = 0; j < arr_length ; ++j) {
            model->T = start + j * increment;
            model->step = 0;
            randomise(model);
            set_evolve(model);
            mags[j] += model->mag;
            energies[j] += model->energy;
        }
    }
    printf("T\tM\tE\n");
    for (int i = 0; i < arr_length; ++i) {
        mags[i] /= repeats;
        energies[i] /= repeats * model->n_points;
//        printf("%g\t%g\t%g\n", start + i * increment, mags[i], energies[i]);
    }

    FILE *file;
    char buf[50];
    sprintf(buf, "output/MvT-B=%g.txt", model->B.z);
    file = fopen(buf, "w");

    fprintf(file, "T\tM\tE\n");
    for (int i = 0; i < arr_length; ++i) {
        fprintf(file, "%g\t%g\t%g\n", start + i * increment, mags[i], energies[i]);
    }
    fclose(file);
}


void var_B(Model* model, double start, double end, double increment, int repeats) {
    distribute_points(model);
    nns(model);
    B_from_dipoles(model);
    precalc_B(model);
    int arr_length = (int)((end - start) / increment) + 1;
    double* mags = (double*)calloc(arr_length, sizeof(double));
    double* energies = (double*)calloc(arr_length, sizeof(double));
    for (int i = 0; i < repeats; ++i) {
        for (int j = 0; j < arr_length ; ++j) {
            model->B.z = start + j * increment;
            model->step = 0;
            randomise(model);
            set_evolve(model);
            mags[j] += model->mag;
            energies[j] += model->energy;
        }
    }
    printf("B\tM\tE\n");
    for (int i = 0; i < arr_length; ++i) {
        mags[i] /= repeats;
        energies[i] /= repeats * model->n_points;
//        printf("%g\t%g\t%g\n", start + i * increment, mags[i], energies[i]);
    }

    FILE *file;
    char buf[50];
    sprintf(buf, "output/MvB-T=%g.txt", model->T);
    file = fopen(buf, "w");

    fprintf(file, "B\tM\tE\n");
    for (int i = 0; i < arr_length; ++i) {
        fprintf(file, "%g\t%g\t%g\n", start + i * increment, mags[i], energies[i]);
    }
    fclose(file);
}


void var_T_B(Model* model, Settings* settings, double start_B, double end_B, double start_T, double end_T, double increment, int repeats) {
    distribute_points(model);
    nns(model);
    B_from_dipoles(model);
    precalc_B(model);
    model->norm = 1. / (repeats * model->n_points);
    int length_B = (int)((end_B - start_B) / increment) + 1;
    int length_T = (int)((end_T - start_T) / increment) + 1;
    model->mags = (double**)calloc(length_T, sizeof(double*));
    model->energies = (double**)calloc(length_T, sizeof(double*));
    model->m_sq = (double**)calloc(length_T, sizeof(double*));
    model->e_sq = (double**)calloc(length_T, sizeof(double*));
    for (int i = 0; i < length_T; ++i) {
        model->mags[i] = (double*)calloc(length_B, sizeof(double));
        model->energies[i] = (double*)calloc(length_B, sizeof(double));
        model->m_sq[i] = (double*)calloc(length_B, sizeof(double));
        model->e_sq[i] = (double*)calloc(length_B, sizeof(double));
    }
    settings->tid = (pthread_t*) malloc(model->thread_count * sizeof(pthread_t));
    settings->models = (Model*) malloc(model->thread_count * sizeof(Model));
    for (int i = 0; i < model->thread_count; ++i) {
        copy_model(model, &settings->models[i]);
        settings->models[i].working = 0;
    }

    for (int i = 0; i < repeats; ++i) {
        printf("repeat: %i\n", i);
        for (int T = 0; T < length_T; ++T) {
            for (int B = 0; B < length_B; ++B) {
//                printf("%g %g\n", T, B);
                int run = 1;
                int free_thread;
                while (run) {
                    int result = find_thread(settings);
                    if (result != -1) {
                        free_thread = result;
                        run = 0;
                    }
                }
//                printf("free thread: %i\n", free_thread);
                copy_model(model, &settings->models[free_thread]);

                settings->models[free_thread].T = start_T + T * increment;
                settings->models[free_thread].B.z = start_B + B * increment;
                settings->models[free_thread].step = 0;

                randomise(&settings->models[free_thread]);
                precalc_B(&settings->models[free_thread]);

                settings->models[free_thread].Ti = T;
                settings->models[free_thread].Bi = B;
                settings->to_run = free_thread;
                settings->models[free_thread].working = 1;

                pthread_create(&settings->tid[free_thread], NULL, thread_evolve, &settings->models[free_thread]);
                pthread_detach(settings->tid[free_thread]);
            }
        }
    }
    printf("Simulation finished\n");

    int wait = 1;
    while (wait) {
        int done = 0;
        for (int i = 0; i < model->thread_count; ++i) {
            if (!settings->models[i].working) {++done;}
        }
        if (done == model->thread_count) {wait=0;}
    }

    for (int i = 0; i < length_T; ++i) {
        for (int j = 0; j < length_B; ++j) {
            model->mags[i][j] /= repeats;
            model->energies[i][j] *= model->norm;
            model->m_sq[i][j] *= model->norm;
            model->e_sq[i][j] *= model->norm;
        }
    }

    FILE *file;
    file = fopen("output/MvTvB.txt", "w");
    fprintf(file, "T\tB\tM\tE\tS\tC\n");
    for (int i = 0; i < length_T; ++i) {
        double T = start_T + i * increment;
        double beta = 1 / (k_b * T);
        for (int j = 0; j < length_B; ++j) {
            double S = (1 / T) * (model->m_sq[i][j] - (pow(model->mags[i][j], 2) * model->n_points));
            double C = (1 / (T*T)) * (model->e_sq[i][j] - (pow(model->energies[i][j], 2) * model->n_points));
            fprintf(file, "%g\t%g\t%g\t%g\t%g\t%g\n",
                    start_T + i * increment, start_B + j * increment, model->mags[i][j], model->energies[i][j], S, C);
        }
    }
    fclose(file);
}


int find_thread(Settings* settings) {
    for (int i = 0; i < settings->thread_count; ++i) {
        if (!settings->models[i].working) {
            return i;
        }
    }
    return -1;
}


void copy_model(Model* old, Model* new) { // couldn't be asked to figure out how to deepcopy, so did it manually
    new->thread_count = old->thread_count;
    new->n_points = old->n_points;
    new->lattice_type = old->lattice_type;
    new->length = old->length;
    new->step = old->step;
    new->B.x = old->B.x;
    new->B.y = old->B.y;
    new->B.z = old->B.z;
    new->width = old->width;
    new->point_spacing = old->point_spacing;
    new->height = old->height;
    new->T = old->T;
    new->delta_checks = old->delta_checks;
    new->evolve_steps = old->evolve_steps;
    new->n_dipoles = old->n_dipoles;
    new->dipoles = (Dipole*) malloc(old->n_dipoles * sizeof(Dipole));
    new->output = old->output;
    new->field_type = old->field_type;
    new->randomise = old->field_type;
    new->radius = old->radius;
    new->working = old->working;
    new->points = (Point*) malloc(old->n_points * sizeof(Point));
    for (int i = 0; i < old->n_points; ++i) {
        new->points[i].coord.x = old->points[i].coord.x;
        new->points[i].coord.y = old->points[i].coord.y;
        new->points[i].coord.z = old->points[i].coord.z;
        new->points[i].n_nns = old->points[i].n_nns;
        new->points[i].nns = (int*) malloc(old->points[i].n_nns * sizeof(int));
        for (int nns = 0; nns < new->points[i].n_nns; ++nns) {
            new->points[i].nns[nns] = old->points[i].nns[nns];
        }
        new->points[i].B.x = old->points[i].B.x;
        new->points[i].B.y = old->points[i].B.y;
        new->points[i].B.z = old->points[i].B.z;
        new->points[i].direct_B = old->points[i].direct_B;
        new->points[i].dipole_B = old->points[i].dipole_B;
    }
    for (int i = 0; i < new->n_dipoles; ++i) {
        new->dipoles[i].moment.x = old->dipoles[i].moment.x;
        new->dipoles[i].moment.y = old->dipoles[i].moment.y;
        new->dipoles[i].moment.z = old->dipoles[i].moment.z;
        new->dipoles[i].coord.x = old->dipoles[i].coord.x;
        new->dipoles[i].coord.y = old->dipoles[i].coord.y;
        new->dipoles[i].coord.z = old->dipoles[i].coord.z;
    }
    new->mags = old->mags;
    new->energies = old->energies;
    new->m_sq = old->m_sq;
    new->e_sq = old->e_sq;
}