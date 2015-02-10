#ifndef ISING_H
#define ISING_H

#include <sys/stat.h>
#include <stdbool.h>

#define MAGPERSPINOUTPUT "output/Temp_Mag_BField" // Unterordner + Name der Datei in der die Magnetisierung pro Spin fuer jeden Sweepschritt am Ende der Simulation gespeichert wird (für m(T) und m(B))
#define ENERGYPERMAG "output/Step_Energy_Mag" // Unterordner + Name der Datei in die pro Step Energie und Magnetisierung geschrieben wird (für E(step), m(step), E(m))
#define MT_MAX 4294967295  // 2^32 - 1, hoechster von mt_random() generierter Wert

typedef struct
{
    bool **cluster;
    double addProbability;
    
} Cluster_Vars;

typedef struct
{
    void *spins;
    int N, dim, steps;
    double J, kB;     // kann scheinbar nicht mehr const sein :/
    double sweep_init, sweep_end, sweep_step, C;
    char sweepPar, matMode, sweepMode, hystMode, clusterMode, filmMode;
    
    Cluster_Vars cv;
  
} Parameters;


void runSweep(Parameters *parameters);


static void usage(char* progname) // typical usage-function
{    
    printf("\nUsage: %s [N] [dim] (steps) (sweepVar) (S_i) (S_e) (S_s) (C) (sweepMode) (hystMode) (clusterMode) (filmMode)\n", progname);
    printf("Where values in [] are required, while values in () are optional \n\n");
    printf("  - N: length of Matrix in one dimesion ( int value > 0 ) \n");
    printf("  - dim: Dimension of Matrix ( int value = 2 or 3 ) \n");
    printf("  - steps: Number of MonteCarlo steps ( int value > 0 ) \n");
    printf("  - sweepPar: sweep parameter for montecarlo-simulations ( magnetic field  or temperatures )( B or T )\n");
    printf("  - S_i: initial value of sweep parameter for the first montecarlo-simulation ( double value ) \n");
    printf("  - S_e: value of sweep parameter where the Simualtion will end ( double value ) \n");
    printf("  - S_s: step size of sweep parameter ( double value ) \n");
    printf("  - C: is the value of extern constant parameter over the sweep ( double value ) \n");
    printf("  - matMode: fill method for matrix start config ( r, p or n ) \n");
    printf("  - sweepMode: sweep over one lattice (otherwise there is a new lattice for every step) ( y or n ) \n");
    printf("  - hystMode: additional simulation back from S_e to S_i ( y or n ) \n");
    printf("  - clusterMode: cluster update as Montecarlo scheme ( y or n ) \n");
    printf("  - filmMode: additional lattice saving ( y or n ) \n");
    printf("\n");
    printf("Example: %s 25 2 1000 T 1.50 4.00 0.10 0.00 r n n n n \n", progname);
    exit(EXIT_FAILURE);
}


void getParameters(int argc, char **argv, Parameters *para)
{
    //zunaechst Standardwerte
    para->steps = 1000;
    para->J = 1.0; para->kB = 1.0; 
    para->sweep_init = 1.5; para->sweep_end = 4.0; para->sweep_step = 0.1; para->C = 0.0000;
    para->sweepPar = 'T'; para->matMode = 'r'; para->sweepMode = 'n'; para->hystMode = 'n'; para->clusterMode = 'n'; para->filmMode = 'n';   

    if(argc < 3 
      || sscanf(argv[1], "%d", &para->N) != 1 || para->N < 1  
      || sscanf(argv[2], "%d", &para->dim) != 1 || (para->dim != 2 && para->dim != 3)
      || (argc > 3 && sscanf(argv[3], "%d", &para->steps) != 1) || para->steps < 1
      || (argc > 4 && sscanf(argv[4], "%c", &para->sweepPar) != 1) || (para->sweepPar != 'B' && para->sweepPar != 'T')
      || (argc > 5 && sscanf(argv[5], "%lf", &para->sweep_init) != 1)
      || (argc > 6 && sscanf(argv[6], "%lf", &para->sweep_end) != 1)
      || (argc > 7 && sscanf(argv[7], "%lf", &para->sweep_step) != 1)
      || (argc > 8 && sscanf(argv[8], "%lf", &para->C) != 1)
      || (argc > 9 && sscanf(argv[9], "%c", &para->matMode) != 1) || (para->matMode != 'r' && para->matMode != 'p' && para->matMode != 'n')
      || (argc > 10 && sscanf(argv[10], "%c", &para->sweepMode) != 1) || (para->sweepMode != 'y' && para->sweepMode != 'n')
      || (argc > 11 && sscanf(argv[11], "%c", &para->hystMode) != 1) || (para->hystMode != 'y' && para->hystMode != 'n')
      || (argc > 12 && sscanf(argv[12], "%c", &para->clusterMode) != 1) || (para->clusterMode != 'y' && para->clusterMode != 'n')
      || (argc > 13 && sscanf(argv[13], "%c", &para->filmMode) != 1) || (para->filmMode != 'y' && para->filmMode != 'n'))
        usage(argv[0]);    
        
    printf("\nattempting to execute %s %d %d %d %c %f %f %f %f %c %c %c %c %c \n\n", argv[0], para->N, para->dim, para->steps, para->sweepPar, para->sweep_init, para->sweep_end, para->sweep_step, para->C, para->matMode, para->sweepMode, para->hystMode, para->clusterMode, para->filmMode); 
}


void initialize(Parameters *para)
{
    FILE *fp;
    int i;
    
    mkdir("output", 0777); // creates a directory like mkdir does
    if (para->filmMode == 'y') mkdir("film", 0777);
    
    if ((fp = fopen(MAGPERSPINOUTPUT, "w")) == NULL) // Dateiinhalt löschen
        fprintf(stderr, "Konnte nicht in Datei %s schreiben \n", MAGPERSPINOUTPUT);
    else
        fclose(fp);
    
    srand(time(NULL)); // muss weiterhin gemacht werden, da der Twister mit rand() initialisiert wird
    mt_init();
    
    for (i = 0; i < 600000; ++i) mt_random();// erstmal den Twister ordentlich aufwaermen!
}   


#endif
