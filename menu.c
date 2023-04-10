#include <stdio.h>
#include <stdlib.h>
#include "menu.h"
#include "model.h"


void settings_main(Model *model) {
    int running = 1;
    while (running) {
        printf("(1) Model Parameters\n");
        printf("(2) Evolution Settings\n");
        printf("(3) Import/Export Settings\n");
        printf("\n(0) Exit\n");
        int choice;
        scanf("%i", &choice);

        switch (choice) {
            default:
                running = 0;
                break;

            case 1:
                settings_param(model);
                break;

            case 2:
                settings_evo(model);
                break;
        }
    }
}


void settings_param(Model *model) {
    int running = 1;
    while (running) {
        printf("(1) Number of points:\t%i\n", model->n_points);
        printf("(2) Temperature:\t%g\n", model->T);
        printf("(3) Magnetic Field Type:\t%s\n", model->field_type == 1 ? "Linear" : model->field_type == 2 ? "Dipole" : "Linear and Dipole");
        printf("(4) Edit %s\n", model->field_type == 1 ? "linear field" : model->field_type == 2 ? "dipoles" : "linear and dipoles");
        printf("(5) Lattice Type: %s\n", model->lattice_type == 1 ? "Rectangular" : model->lattice_type == 2 ? "Cuboidal" : "Spherical");
        printf("(6) Edit Lattice Parameters\n");
        printf("\n(0) Exit\n");
        int choice = 0;
        scanf("%i", &choice);

        switch (choice) {
            default:
                running = 0;
                break;

            case 1:
                printf("Enter new number of lattice points:");
                scanf("%i", &model->n_points);
                model->points = (Point*) realloc(model->points, model->n_points * sizeof(Point));
                if (model->points == NULL) {
                    printf("Error occurred allocating memory!\n");
                    exit(0);
                }
                break;

            case 2:
                printf("Enter new temperature:");
                scanf("%lf", &model->T);
                break;

            case 3:
                printf("Type of field to be used:\n");
                printf("(1) Direct\n(2) Dipole\n");
                scanf("%i", &model->field_type);
                break;

            case 4:
                switch (model->field_type) {
                    default:
                        printf("Invalid field type selected: %i\n", model->field_type);
                        break;

                    case 1:
                        printf("Enter magnetic field (x y z):\n");
                        scanf("%lf %lf %lf", &model->B.x, &model->B.y, &model->B.z);
                        break;

                    case 2:
                        settings_dipoles(model);
                }
                break;

            case 5:
                printf("Enter lattice type:");
                printf("(1) Rectangular\n(2) Cuboidal\n(3) Spherical");
                scanf("%i", &model->lattice_type);
                break;

            case 6:
                break;
        }
    }
}


void settings_dipoles(Model *model) {
    int running = 1;
    while (running) {
        printf("Current dipoles:\n");
        for (int i = 0; i < model->n_dipoles; i++) {
            printf("Dipole %i:\t(%g, %g, %g)J/T\t(%g, %g, %g)\n", i+1, model->dipoles[i].moment.x, model->dipoles[i].moment.y, model->dipoles[i].moment.z, model->dipoles[i].coord.x, model->dipoles[i].coord.y, model->dipoles[i].coord.z);
        }
        printf("(1) Add Dipole\n");
        printf("(2) Delete Dipole\n");
        printf("(3) Edit Dipole\n");
        printf("\n(0) Exit\n");

        int choice = 0;
        scanf("%i", &choice);
        switch (choice) {
            default:
                running = 0;
                break;

            case 1: // add dipole
                model->n_dipoles++;
                model->dipoles = (Dipole*) realloc(model->dipoles, model->n_dipoles * sizeof(Dipole));
                if (model->dipoles == NULL) {
                    printf("Error occurred allocating memory!\n");
                    exit(0);
                }
                printf("Enter dipole moment vector x y z:\n");
                scanf("%lf %lf %lf", &model->dipoles[model->n_dipoles-1].moment.x, &model->dipoles[model->n_dipoles-1].moment.y, &model->dipoles[model->n_dipoles-1].moment.z);
                printf("Enter dipole position vector x y z:\n");
                scanf("%lf %lf %lf", &model->dipoles[model->n_dipoles-1].coord.x, &model->dipoles[model->n_dipoles-1].coord.y, &model->dipoles[model->n_dipoles-1].coord.z);
                break;

            case 2: // delete dipole
                printf("Enter ID of dipole to delete:");
                int to_yeet;
                scanf("%i", &to_yeet);
                if (0 < to_yeet <= model->n_dipoles) {
                    for (int i = to_yeet-1; i < model->n_dipoles-1; ++i) {
                        model->dipoles[i] = model->dipoles[i+1];
                    }
                    for (int i = 0; i < model->n_dipoles; i++) {
                        printf("Dipole %i:\t(%g, %g, %g)J/T\t(%g, %g, %g)\n", i+1, model->dipoles[i].moment.x, model->dipoles[i].moment.y, model->dipoles[i].moment.z, model->dipoles[i].coord.x, model->dipoles[i].coord.y, model->dipoles[i].coord.z);
                    }

                    model->n_dipoles--;
                    model->dipoles = (Dipole*) realloc(model->dipoles, model->n_dipoles * sizeof(Dipole));
                    if (model->dipoles == NULL) {
                        printf("Error occurred allocating memory!\n");
                        exit(0);
                    }
                }
                break;

            case 3: // edit dipole
                printf("Enter ID of dipole to edit:");
                int to_edit;
                scanf("%i", &to_edit);
                if (0 < to_edit <= model->n_dipoles) {
                    printf("Enter dipole moment vector x y z:\n");
                    scanf("%lf %lf %lf", &model->dipoles[to_edit-1].moment.x, &model->dipoles[to_edit-1].moment.y, &model->dipoles[to_edit-1].moment.z);
                    printf("Enter dipole position vector x y z:\n");
                    scanf("%lf %lf %lf", &model->dipoles[to_edit-1].coord.x, &model->dipoles[to_edit-1].coord.y, &model->dipoles[to_edit-1].coord.z);
                }
                break;
        }
    }
}


void settings_evo(Model *model) {
    int running = 1;
    while (running) {
        int choice = 0;
        printf("(1) Evolution steps:\t%i\n", model->evolve_steps);
        printf("(2) Starting step:\t%i\n", model->step);
        printf("(3) Delta checks:\t%i\n", model->delta_checks);
        printf("(4) Random at start:\t%s\n", model->randomise ? "Yes" : "No");
        printf("(5) Output steps:\t%s\n", model->output ? "Yes" : "No");
        printf("\n(0) Exit\n");
        scanf("%i", &choice);

        switch (choice) {
            default:
                running = 0;
                break;

            case 1:
                printf("Enter number of evolution steps (0=auto):");
                scanf("%i", &model->evolve_steps);
                break;

            case 2:
                printf("Starting step of model (if file input):");
                scanf("%i", &model->step);
                break;

            case 3:
                printf("Number of delta checks (auto stop):");
                scanf("%i", &model->delta_checks);
                break;

            case 4:
                printf("Randomise model at start?\n");
                printf("(1) Yes\n(0) No\n");
                scanf("%i", &model->randomise);
                break;

            case 5:
                printf("Output evolution steps?\n");
                printf("(1) Yes\n(0) No\n");
                scanf("%i", &model->output);
                break;
        }
    }
}