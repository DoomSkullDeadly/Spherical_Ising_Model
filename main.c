#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "model.h"
#include "menu.h"


// other todos:
// TODO: functions of T and B
// TODO: input validation


int main() {
    srand(time(NULL));
    Model model = {400,
                   0,
                   0,
                   0,
                   0,
                   0,
                   20,
                   3.,
                   0., 0., 20.,
                   0,
                   1,
                   3,
                   1,
                   20,
                   20,
                   20,
                   1,
                   1};
    // TODO: initial params declared earlier as constants, allow external file for user set default settings
    int running = 1;
    while (running) {
        printf("Model Parameters:\n");
        printf("N points\tTemperature\tMagnetic Field\tEvolution Steps\tDelta Checks\n");
        printf("%i\t\t%gK\t\t%gT\t\t%i\t\t%i\n\n", model.n_points, model.T, model.B, model.evolve_steps, model.delta_checks);
        printf("(1) Run Model\n(2) Settings\n(3) Magnetisation vs Temperature\n(4) Magnetisation vs Magnetic Field\n(5) Lattice From File\n(6) Create Video\n(0) Exit\n");
        int input;
        scanf("%i", &input);
//        input = 1;  // to debug so don't have to faf with inputs
        switch (input) {
            default:
                running = 0;
                break;

            case 1:
                printf("Allocating memory for points\n");
                model.points = (Point*) malloc(model.n_points * sizeof(Point));
                if (model.points == NULL) {
                    printf("Error occurred allocating memory!\n");
                    exit(0);
                }
                printf("Distributing points\n");
                distribute_points(&model);
                if (model.randomise) {
                    printf("Randomising points\n");
                    randomise(&model);
                }
                printf("Finding nearest neighbours\n");
                nns(&model);
                printf("Evolving model\n");
                set_evolve(&model);
                break;

            case 2: // editing model settings
                settings_main(&model);
                break;

            case 3:
                break;

            case 4:
                break;

            case 5:
                break;

            case 6:
                break;
        }
    }

    return 0;
}
