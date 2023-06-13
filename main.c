#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "model.h"
#include "menu.h"


// other todos:
// TODO: input validation


int main() {
    srand(time(NULL)); // NOLINT
    Model model = {12,
                   1600,
                   0,
                   0,
                   0,
                   25,
                   0,
                   20,
                   5.,
                   0., 0., 0.,
                   0, 0,
                   0,
                   1,
                   1,
                   1,
                   10,
                   10,
                   10,
                   1,
                   1,
                   0};
    // TODO: initial params declared earlier as constants, allow external file for user set default settings
    int running = 1;
    while (running) {
        printf("Model Parameters:\n");
        printf("N points\tTemperature\tMagnetic Field\tEvolution Steps\tDelta Checks\n");
        printf("%i\t\t%gK\t\t%gT\t\t%i\t\t%i\n\n", model.n_points, model.T, model.B, model.evolve_steps, model.delta_checks);
        printf("(1) Run Model\n(2) Settings\n(3) Magnetisation vs Temperature\n(4) Magnetisation vs Magnetic Field\n(5) Varying B and T\n(6) Lattice From File\n(7) Create Video\n(0) Exit\n");
        int input;
        scanf("%i", &input);
//        input = 1;  // to debug so don't have to faf with inputs
        switch (input) {
            default:
                running = 0;
                break;

            case 1:
                printf("Distributing points\n");
                distribute_points(&model);
                if (model.randomise) {
                    printf("Randomising points\n");
                    randomise(&model);
                }
                printf("Finding nearest neighbours\n");
                nns(&model);
                printf("Calculating magnetic field\n");
                B_from_dipoles(&model);
                printf("Evolving model\n");
                set_evolve(&model);
                free_Points(&model);
                break;

            case 2: // editing model settings
                settings_main(&model);
                break;

            case 3: {
                double start, end, increment;
                int repeats;

                printf("Enter start, end, increment for T:\n");
                scanf("%lf %lf %lf", &start, &end, &increment); //NOLINT
                printf("Enter number of repeats:\n");
                scanf("%i", &repeats); //NOLINT

                var_T(&model, start, end, increment, repeats);
                break;
            }

            case 4: {
                double start, end, increment;
                int repeats;

                printf("Enter start, end, increment for B:\n");
                scanf("%lf %lf %lf", &start, &end, &increment); //NOLINT
                printf("Enter number of repeats:\n");
                scanf("%i", &repeats); //NOLINT

                var_B(&model, start, end, increment, repeats);
                break;
            }

            case 5: {
                double start_B, end_B, start_T, end_T, increment;
                int repeats;

                printf("Enter start and end for T:\n");
                scanf("%lf %lf", &start_T, &end_T); //NOLINT
                printf("Enter start and end for B:\n");
                scanf("%lf %lf", &start_B, &end_B); //NOLINT
                printf("Enter increment for T and B:\n");
                scanf("%lf", &increment); // NOLINT
                printf("Enter number of repeats:\n");
                scanf("%i", &repeats); //NOLINT

                Settings settings = {12, 0};

                var_T_B(&model, &settings, start_B, end_B, start_T, end_T, increment, repeats);
                break;
            }

            case 6:
                break;
        }
    }

    return 0;
}
