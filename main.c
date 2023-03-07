#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Vec3.h"

// other todos:
// TODO: cubic lattice
// TODO: functions of T and B
// TODO: input validation
// TODO: rewrite menu system to use functions rather than overcomplicate main()

#define J 1.
#define mu_b 1.
#define mu_0 1.
#define k_b 1.


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


typedef struct { // TODO: settings should be their own struct, such as evolve_steps, delta_checks, output, randomise
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

// TODO: rewrite functions into external file to keep main for menu only
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


int main() {
    srand(time(NULL));
    Model model = {1000, 0, 0, 0, 0, 0, 20, 3., 20., 0., 0., 1, 1, 2};
    // TODO: initial params declared earlier as constants, allow external file for user set default settings
    int running = 1;
    double lower, upper, increment;
    // TODO: rewrite menu system to accommodate other changes
    while (running) {
        printf("Model Parameters:\n");
        printf("N points\tTemperature\tMagnetic Field\tEvolution Steps\tDelta Checks\n");
        printf("%i\t\t%gK\t\t%gT\t\t%i\t\t%i\n\n", model.n_points, model.T, model.B, model.evolve_steps, model.delta_checks);
        printf("(1) Run Model\n(2) Settings\n(3) Magnetisation vs Temperature\n(4) Magnetisation vs Magnetic Field\n(5) Lattice From File\n(6) Create Video\n(0) Exit\n");
        int input;
        scanf("%i", &input);

        switch (input) {
            default:
                running = 0;
                break;

            case 1:
                model.points = (Point*) malloc(model.n_points * sizeof(Point));
                if (model.points == NULL) {
                    printf("Error occurred allocating memory!\n");
                    exit(0);
                }
                distribute_points(&model);
                if (model.randomise) {
                    randomise(&model);
                }
                nns(&model);
                set_evolve(&model);
                break;

            case 2: // editing model settings
                printf("Edit Parameters:\n");
                int settings = 1;
                while (settings) {
                    printf("(1) Model Parameters\n");
                    printf("(2) Evolution Settings\n");
                    printf("\n(0) Exit\n");
                    int setting;
                    scanf("%i", &setting);

                    switch (setting) {
                        default:
                            settings = 0;
                            break;

                        case 1: // settings for model parameters
                            printf("");
                            int settings1 = 1;
                            while (settings1) {
                                printf("(1) Number of points:\t%i\n", model.n_points);
                                printf("(2) Temperature:\t%g\n", model.T);
                                printf("(3) Magnetic Field Type:\t%s\n", model.field_type == 1 ? "Direct" : "Dipole");
                                printf("(4) Edit %s\n", model.field_type == 1 ? "direct field" : "dipoles");
                                printf("\n(0) Exit\n");
                                int param = 0;
                                scanf("%i", &param);

                                switch (param) { // switch parameters to new user inputs
                                    default:
                                        settings1 = 0;
                                        break;

                                    case 1:
                                        printf("Enter new number of lattice points:");
                                        scanf("%i", &model.n_points);
                                        model.points = (Point*) realloc(model.points, model.n_points * sizeof(Point));
                                        if (model.points == NULL) {
                                            printf("Error occurred allocating memory!\n");
                                            exit(0);
                                        }
                                        break;

                                    case 2:
                                        printf("Enter new temperature:");
                                        scanf("%lf", &model.T);
                                        break;

                                    case 3:
                                        printf("Type of field to be used:\n");
                                        printf("(1) Direct\n(2) Dipole\n");
                                        scanf("%i", &model.field_type);
                                        break;

                                    case 4:
                                        switch (model.field_type) {
                                            default:
                                                printf("Invalid field type selected: %i\n", model.field_type);
                                                break;

                                            case 1:
                                                printf("Enter magnetic field (x y z):\n");
                                                scanf("%lf %lf %lf", &model.B.x, &model.B.y, &model.B.z);
                                                break;

                                            case 2:
                                                printf("");
                                                int edit_dipoles = 1;
                                                while (edit_dipoles) {
                                                    printf("Current dipoles:\n");
                                                    for (int i = 0; i < model.n_dipoles; i++) {
                                                        printf("Dipole %i:\t(%g, %g, %g)J/T\t(%g, %g, %g)\n", i+1, model.dipoles[i].moment.x, model.dipoles[i].moment.y, model.dipoles[i].moment.z, model.dipoles[i].coord.x, model.dipoles[i].coord.y, model.dipoles[i].coord.z);
                                                    }
                                                    printf("(1) Add Dipole\n");
                                                    printf("(2) Delete Dipole\n");
                                                    printf("(3) Edit Dipole\n");
                                                    printf("\n(0) Exit\n");

                                                    int dipole_choice = 0;
                                                    scanf("%i", &dipole_choice);
                                                    switch (dipole_choice) {
                                                        default:
                                                            edit_dipoles = 0;
                                                            break;

                                                        case 1: // add dipole
                                                            model.n_dipoles++;
                                                            model.dipoles = (Dipole*) realloc(model.dipoles, model.n_dipoles * sizeof(Dipole));
                                                            if (model.dipoles == NULL) {
                                                                printf("Error occurred allocating memory!\n");
                                                                exit(0);
                                                            }
                                                            printf("Enter dipole moment vector x y z:\n");
                                                            scanf("%lf %lf %lf", &model.dipoles[model.n_dipoles-1].moment.x, &model.dipoles[model.n_dipoles-1].moment.y, &model.dipoles[model.n_dipoles-1].moment.z);
                                                            printf("Enter dipole position vector x y z:\n");
                                                            scanf("%lf %lf %lf", &model.dipoles[model.n_dipoles-1].coord.x, &model.dipoles[model.n_dipoles-1].coord.y, &model.dipoles[model.n_dipoles-1].coord.z);
                                                            break;

                                                        case 2: // delete dipole
                                                            printf("Enter ID of dipole to delete:");
                                                            int to_yeet;
                                                            scanf("%i", &to_yeet);
                                                            if (0 < to_yeet <= model.n_dipoles) {
//                                                                Dipole tmp[model.n_dipoles-1]; // set up temp array
//                                                                for (int i = 0; i < to_yeet; ++i) { // copy array
//                                                                    tmp[i] = model.dipoles[i];
//                                                                }
//                                                                for (int i = to_yeet; i < model.n_dipoles; ++i) {
//                                                                    tmp[i] = model.dipoles[i+1];
//                                                                }
                                                                for (int i = to_yeet-1; i < model.n_dipoles-1; ++i) {
                                                                    model.dipoles[i] = model.dipoles[i+1];
                                                                }
                                                                for (int i = 0; i < model.n_dipoles; i++) {
                                                                    printf("Dipole %i:\t(%g, %g, %g)J/T\t(%g, %g, %g)\n", i+1, model.dipoles[i].moment.x, model.dipoles[i].moment.y, model.dipoles[i].moment.z, model.dipoles[i].coord.x, model.dipoles[i].coord.y, model.dipoles[i].coord.z);
                                                                }

                                                                model.n_dipoles--;
                                                                model.dipoles = (Dipole*) realloc(model.dipoles, model.n_dipoles * sizeof(Dipole));
                                                                if (model.dipoles == NULL) {
                                                                    printf("Error occurred allocating memory!\n");
                                                                    exit(0);
                                                                }
                                                            }
                                                            break;

                                                        case 3:
                                                            printf("Enter ID of dipole to edit:");
                                                            int to_edit;
                                                            scanf("%i", &to_edit);
                                                            if (0 < to_edit <= model.n_dipoles) {
                                                                printf("Enter dipole moment vector x y z:\n");
                                                                scanf("%lf %lf %lf", &model.dipoles[to_edit-1].moment.x, &model.dipoles[to_edit-1].moment.y, &model.dipoles[to_edit-1].moment.z);
                                                                printf("Enter dipole position vector x y z:\n");
                                                                scanf("%lf %lf %lf", &model.dipoles[to_edit-1].coord.x, &model.dipoles[to_edit-1].coord.y, &model.dipoles[to_edit-1].coord.z);
                                                            }
                                                            break;
                                                    }
                                                }
                                        }
                                }
                            }
                            break;

                        case 2: // settings for evolution of the model
                            printf("");
                            int settings2 = 1;
                            while (settings2) {
                                int param = 0;
                                printf("(1) Evolution steps:\t%i\n", model.evolve_steps);
                                printf("(2) Starting step:\t%i\n", model.step);
                                printf("(3) Delta checks:\t%i\n", model.delta_checks);
                                printf("(4) Random at start:\t%s\n", model.randomise ? "Yes" : "No");
                                printf("(5) Output steps:\t%s\n", model.output ? "Yes" : "No");
                                printf("\n(0) Exit\n");
                                scanf("%i", &param);

                                switch (param) {
                                    default:
                                        settings2 = 0;
                                        break;

                                    case 1:
                                        printf("Enter number of evolution steps (0=auto):");
                                        scanf("%i", &model.evolve_steps);
                                        break;

                                    case 2:
                                        printf("Starting step of model (if file input):");
                                        scanf("%i", &model.step);
                                        break;

                                    case 3:
                                        printf("Number of delta checks (auto stop):");
                                        scanf("%i", &model.delta_checks);
                                        break;

                                    case 4:
                                        printf("Randomise model at start?\n");
                                        printf("(1) Yes\n(0) No\n");
                                        scanf("%i", &model.randomise);
                                        break;

                                    case 5:
                                        printf("Output evolution steps?\n");
                                        printf("(1) Yes\n(0) No\n");
                                        scanf("%i", &model.output);
                                        break;
                                }
                            }
                            break;
                    }

                }
                break;

            case 3:  // temporary for testing
                model.points = (Point*) malloc(model.n_points * sizeof(Point));
                if (model.points == NULL) {
                    printf("Error occurred allocating memory!\n");
                    exit(0);
                }
                distribute_points(&model);
                randomise(&model);
                nns(&model);

                model.dipoles = (Dipole*) malloc(model.n_dipoles * sizeof(Dipole));
                if (model.dipoles == NULL) {
                    printf("Error occurred allocating memory!\n");
                    exit(0);
                }
                model.dipoles[0].moment.x = 0;
                model.dipoles[0].moment.y = 0;
                model.dipoles[0].moment.z = 1000;
                model.dipoles[0].coord.x = 1;
                model.dipoles[0].coord.y = 0;
                model.dipoles[0].coord.z = 0;
                B_from_dipoles(&model);

                set_evolve(&model);
                break;
        }
    }

    return 0;
}


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