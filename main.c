#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define J 1
#define mu_b 1


typedef struct {
    double x;
    double y;
    double z;
} Vec3;


typedef struct {
    Vec3 coord;
    unsigned char spin;
    int close_6[6];
} Point;


typedef struct {
    int n_points;
    double energy;
    double mag;
    int evolve_steps;
    int step;
    int delta_checks;
    double T;
    Vec3 B;
    int output;
    int randomise;
    Point *points;
} Model;


void distribute_points(Model*);

double energy(Model*);

double norm_mag(Model*);

double point_distance(Point, Point);

void nns(Model*);

void randomise(Model*);

void evolve(Model*);

void output(Model*);

void set_evolve(Model*);

double dot_prod(Vec3, Vec3);

double vec_mag(Vec3);


int main() {
    srand(time(NULL));
    Model model = {100, 0, 0, 0, 0, 20, 1., 2., 0, 0, 0, 1};

    int running = 1;
    double lower, upper, increment;

    while (running) {
        printf("Model Parameters:\n");
        printf("N points\tTemperature\tMagnetic Field\tEvolution Steps\tDelta Checks\n");
        printf("%i\t\t%gK\t\t%gT\t\t%i\t\t%i\n\n", model.n_points, model.T, model.B, model.evolve_steps, model.delta_checks);
        printf("(1) Run Model\n(2) Edit Parameters\n(3) Magnetisation vs Temperature\n(4) Magnetisation vs Magnetic Field\n(5) Lattice From File\n(6) Create Video\n(0) Exit\n");
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
                set_evolve(&model);
                break;

            case 2:  // temporary for testing
                model.points = (Point*) malloc(model.n_points * sizeof(Point));
                if (model.points == NULL) {
                    printf("Error occurred allocating memory!\n");
                    exit(0);
                }
                distribute_points(&model);
                nns(&model);
                randomise(&model);
                output(&model);
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


double point_distance(Point p1, Point p2) {
    double sum = 0;
    sum += pow(p1.coord.x - p2.coord.x, 2);
    sum += pow(p1.coord.y - p2.coord.y, 2);
    sum += pow(p1.coord.z - p2.coord.z, 2);
    return sum;
}


void nns(Model* model) {
    int indexes[model->n_points];
    for (int i = 0; i < model->n_points; ++i) {  // setting up index array
        indexes[i] = i;
    }
    for (int point = 0; point < model->n_points; ++point) {  // looping over all points
        for (int i = 1; i < model->n_points; ++i) {  // insertion sort to get all indexes of points sorted by distance
            int key = indexes[i];
            int j = i-1;
            while (j >= 0 && point_distance(model->points[point], model->points[key]) < point_distance(model->points[point], model->points[indexes[j]])) {
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
        double p_dot_b = dot_prod(model->points[point].coord, model->B);
        E -= mu_b * p_dot_b / vec_mag(model->points[point].coord) * (double)model->points[point].spin-.5;
    }
    return E;
}


void evolve(Model*) {}


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


double dot_prod(Vec3 v1, Vec3 v2) {
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