/*
 *
 *  Created by Samuel Melvin Flaxman on 2/24/15.
 *  Copyright 2015 Samuel Melvin Flaxman. All rights reserved.
 *
 */
 
 //To do: get all variables into main function and use malloc; pass variables properly between functions
 //CUTOFF_DEFAULT issues?
 

const char *version = "1.0.3";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>

#include "MT/dSFMT.h"

// code for using Mersenne Twister RNG
dsfmt_t dsfmt;
#define seedRand(s) dsfmt_init_gen_rand(&dsfmt, s)
#define	randU() dsfmt_genrand_close_open(&dsfmt)
#define	randI() (unsigned)dsfmt_genrand_uint32(&dsfmt)


// parameters
#define nDEMES 2
#define TOTAL_N_DEFAULT 10000 // population size
#define TOTAL_GENERATIONS_DEFAULT 10000 // total length of time to run
#define nGENOMIC_REGIONS_DEFAULT 2 // number of large scale regions of potential divergence to simulate
#define nWINDOWS_PER_REGION_DEFAULT 60 // number of windows in each genome
#define WINDOW_SIZE_BASES_DEFAULT 20000 // window size in number of bases
#define nSITES_PER_WINDOW_DEFAULT 20 // maximum number of polymorphic sites within a window
#define MIGRATION_RATE_DEFAULT 0.05 // fraction of each deme migrating to other deme per generation
#define MEAN_S_BG_DEFAULT 0.1 // selection coefficient mean against deleterious mutations ("background")
#define MEAN_S_DIVERGENT_DEFAULT 0.01 // mean selection coefficient for divergently selected mutations
#define MEAN_S_BENEFICIAL_DEFAUT 0.01 // mean selection coefficient for universally beneficial mutations
#define PROB_MUT_BG_DEFAULT 0.0 // probability a new mutation is deleterious
#define PROB_MUT_DIVERGENT_DEFAULT 0.0 // probability a new mutation is divergently selected
#define PROB_MUT_BENEFICIAL_DEFAULT 0.0 // probability a new mutation is universally adaptive
#define MU_HIGH_DEFAULT 1E-06 // "high" mutation rate per site per individual per meiosis
#define MU_LOW_DEFAULT 1E-06 // "low" mutation rate per site per individual per meiosis
#define DURATION_ALLOPATRY_DEFAULT 0 // number of generations of allopatry at beginning
#define DETERMINISTIC_DEFAULT 1 // use RNG seed in file or not (Boolean)'
#define HIGH_RECOMB_PROB_PER_KB_DEFAULT 0.001 // recombination probability per 1000 bases per meiosis
#define LOW_RECOMB_PROB_PER_KB_DEFAULT 0.0001
#define INITIAL_CONDITION_DEFAULT 0 // for setting up alleles and frequencies in initialization
#define TS_SAMP_FREQ 1000
#define CUTOFF_DEFAULT 0.3


// useful flags
#define ALLELE_TYPE_NEUTRAL 0
#define ALLELE_TYPE_DIVERGENT 1
#define ALLELE_TYPE_BG 2

#define INIT_NO_INITIAL_DIVERGENCE 0
#define INIT_NEUTRAL_STANDING_HWE 1

#define MODEL_TYPE_NEUTRAL_ONLY 0
#define MODEL_TYPE_SELECTION 1

#define REGION_FLAG_START 1
#define REGION_FLAG_END 2
#define REGION_FLAG_MIDDLE 0

#define WINDOW_FLAG_START 1
#define WINDOW_FLAG_END 2
#define WINDOW_FLAG_MIDDLE 0

#define HIGH_MUT_FLAG 1
#define LOW_MUT_FLAG 0


// function declarations
short int addAnewMutation(int locus, short int parentalAllele);
double boxMuller(double mu, double sd);
void buildFitnessArray( double *selectionCoefficients, int *demeIndexesStart, int deme, int numInDeme );
void calculateMetricsAndStats(void);
void finalPrints(void);
void initializePopulation(void);
void makeOffspring(int momi, int dadi, short int *offspringGTpt);
void migration(int *demeIndexes);
void openDataRecordingFiles(void);
int pickParent( double *fitnessArray, int prevChosen, int nInDeme );
long int Poisson(double mm);
void reproduction(int *demeIndexes);
void RNGsetup(void);
void usage(char *s);
void useShortTestValuesOfParameters(void);
int compare(const void *p, const void *q);
double Cutoff(double viabilityArray);
double walker(double* dum, int index);

// global variables for command line options
int TOTAL_N = TOTAL_N_DEFAULT;
long int TOTAL_GENERATIONS = TOTAL_GENERATIONS_DEFAULT;
int nGENOMIC_REGIONS = nGENOMIC_REGIONS_DEFAULT;
int nWINDOWS_PER_REGION = nWINDOWS_PER_REGION_DEFAULT;
int WINDOW_SIZE_BASES = WINDOW_SIZE_BASES_DEFAULT;
int nSITES_PER_WINDOW = nSITES_PER_WINDOW_DEFAULT;
double MIGRATION_RATE = MIGRATION_RATE_DEFAULT;
double MEAN_S_BENEFICIAL = MEAN_S_BENEFICIAL_DEFAUT;
double MEAN_S_BG = MEAN_S_BG_DEFAULT;
double MEAN_S_DIVERGENT = MEAN_S_DIVERGENT_DEFAULT;
double PROB_MUT_BENEFICIAL = PROB_MUT_BENEFICIAL_DEFAULT;
double PROB_MUT_BG = PROB_MUT_BG_DEFAULT;
double PROB_MUT_DIVERGENT = PROB_MUT_DIVERGENT_DEFAULT;
double MU_HIGH = MU_HIGH_DEFAULT;
double MU_LOW = MU_LOW_DEFAULT;
int DURATION_ALLOPATRY = DURATION_ALLOPATRY_DEFAULT;
int DETERMINISTIC = DETERMINISTIC_DEFAULT;
double HIGH_RECOMB_PROB_PER_KB = HIGH_RECOMB_PROB_PER_KB_DEFAULT;
double LOW_RECOMB_PROB_PER_KB = LOW_RECOMB_PROB_PER_KB_DEFAULT;
int INITIAL_CONDITION = INITIAL_CONDITION_DEFAULT;
int MODEL_TYPE = MODEL_TYPE_NEUTRAL_ONLY;
int MONITOR_PROGRESS = 0;
double CUTOFF = CUTOFF_DEFAULT;

// other global variables
int demeCounts[nDEMES];
int totalSitesInGenome;
int nInDeme;
FILE *alleleFrequenciesOverTime, *mutationLog, *fixationLog, *piOverTime;
FILE *DaOverTime, *DxyOverTime, *variableLociOverTime, *FSToverTime;
long int t = 0;
int initializer=0, cinitializer=0;
double* selectionCoefficients;
const int cutoff;

// global pointers for malloc() calls
double *genomeMap;
short int *genotypes;
int *isVariableSite;
short int *fixedAllele;
double *alleleFrequenciesByDeme;
double *alleleFrequencies;
double *alleleSelectionCoeffs;
int *demeLocations;
int *regionFlags;
int *windowFlags;
int *regionMembership;
int *windowMembership;
int *mutationRateClass;

int main(int argc, char *argv[])
{
    int ch;
    char *progname = argv[0];

    // read in optional command line arguments
    while ((ch = getopt(argc, argv, "A:C:D:I:L:m:M:N:p:P:q:r:R:s:S:t:T:u:U:w:W:z:Z?")) != -1) {
        switch (ch) {
            case 'A':
                DURATION_ALLOPATRY = atoi(optarg);
                break;
            case 'C':
                CUTOFF = strtod(optarg, (char**)NULL);
                break;
            case 'D':
                DETERMINISTIC = atoi(optarg);
                break;
            case 'I':
                INITIAL_CONDITION = atoi(optarg);
                break;
            case 'L':
                nSITES_PER_WINDOW = atoi(optarg);
                break;
            case 'm':
                MIGRATION_RATE = strtod(optarg, (char **)NULL);
                break;
            case 'M':
                MODEL_TYPE = atoi(optarg);
                break;
            case 'N':
                TOTAL_N = atoi(optarg);
                break;
            case 'p':
                PROB_MUT_BENEFICIAL = strtod(optarg, (char **)NULL);
                break;
            case 'P':
                PROB_MUT_BG = strtod(optarg, (char **)NULL);
                break;
            case 'q':
                PROB_MUT_DIVERGENT = strtod(optarg, (char **)NULL);
                break;
            case 'Q':
                HIGH_RECOMB_PROB_PER_KB = strtod(optarg, (char **)NULL);
                break;
            case 'r':
                LOW_RECOMB_PROB_PER_KB = strtod(optarg, (char **)NULL);
                break;
            case 'R':
                nGENOMIC_REGIONS = atoi(optarg);
                break;
            case 's':
                MEAN_S_BENEFICIAL = strtod(optarg, (char **)NULL);
                break;
            case 'S':
                MEAN_S_BG = strtod(optarg, (char **)NULL);
                break;
            case 't':
                MEAN_S_DIVERGENT = strtod(optarg, (char **)NULL);
                break;
            case 'T':
                TOTAL_GENERATIONS = atoi(optarg);
                break;
            case 'u':
                MU_LOW = strtod(optarg, (char **)NULL);
                break;
            case 'U':
                MU_HIGH = strtod(optarg, (char **)NULL);
                break;
            case 'w':
                WINDOW_SIZE_BASES = atoi(optarg);
                break;
            case 'W':
                nWINDOWS_PER_REGION = atoi(optarg);
                break;
            case 'z':
                MONITOR_PROGRESS = atoi(optarg);
                break;
            case 'Z':
                useShortTestValuesOfParameters();
                break;
            case '?':
            default:
                usage(progname);
                exit(-1);
        }
    }
    
    int *demeIndexes;
    demeIndexes = (int *) malloc( (sizeof(int) * TOTAL_N * nDEMES));
    
    // set up
    RNGsetup();
    initializePopulation();
    openDataRecordingFiles();
    
    
    
    // main operations
    for ( t = 1; t <= TOTAL_GENERATIONS; t++ ) {
        
        if ( t > DURATION_ALLOPATRY )
            migration( demeIndexes );
        
        reproduction( demeIndexes );
        
        calculateMetricsAndStats();
        
        if ( MONITOR_PROGRESS ) {
            if ( t % MONITOR_PROGRESS == 0 )
                fprintf(stdout, "%li\n", t);
        }
    }
    
    
    finalPrints();
    
    // free blocks from malloc
    free(genomeMap);
    free(genotypes);
    free(isVariableSite);
    free(fixedAllele);
    free(alleleFrequenciesByDeme);
    free(alleleSelectionCoeffs);
    free(demeLocations);
    free(demeIndexes);
    free(regionFlags);
    free(windowFlags);
    free(regionMembership);
    free(windowMembership);
    free(mutationRateClass);
    free(alleleFrequencies);
    
    // close files
    fclose(mutationLog);
    fclose(alleleFrequenciesOverTime);
    fclose(fixationLog);
    fclose(piOverTime);
    fclose(DaOverTime);
    fclose(DxyOverTime);
    fclose(variableLociOverTime);
    
    return 0;
}


short int addAnewMutation(int locus, short int parentalAllele)
{
    short int newAllele;
    
    isVariableSite[locus] = 1;
    fixedAllele[locus] = -1;
    
    if ( parentalAllele )
        newAllele = 0;
    else
        newAllele = 1;
    
    if ( MODEL_TYPE == MODEL_TYPE_SELECTION ) {
        fprintf(stderr, "\nError in addAnewMutation(): routine for adding selected mutations not yet built!\n");
        exit(-1);
    }
    
    fprintf(mutationLog, "%li,%i,%i,%i,%i\n", t, locus, newAllele, regionMembership[locus], mutationRateClass[locus]);
    
    return newAllele;
}



double boxMuller(double mu, double sd)
{
  //  /* boxmuller.c           Implements the Polar form of the Box-Muller
	 /*Transformation
	 
	 (c) Copyright 1994, Everett F. Carter Jr.
	 Permission is granted by the author to use
	 this software for any application provided this
	 copyright notice is preserved.
	 [ smf's note: accessed online 12.03.06 at
	 http:www.taygeta.com/random/boxmuller.html ]*/
	 
    
    double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;
	
	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * ( randU() ) - 1.0;
			x2 = 2.0 * ( randU() ) - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );
		
		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}
	
	return ( mu + y1 * sd );
}


void buildFitnessArray(int *demeIndexesStart, int deme, int numInDeme )
{
double siteTypeProbabilities[3], dum;
int i;
    siteTypeProbabilities[0] = PROB_MUT_BENEFICIAL;
    siteTypeProbabilities[1] = siteTypeProbabilities[0] + PROB_MUT_BG;
    siteTypeProbabilities[2] = siteTypeProbabilities[1] + PROB_MUT_DIVERGENT;
    
    // make a check to see if the cumulative probabilites exceed 1
    
	if (siteTypeProbabilities[0]+siteTypeProbabilities[1]+siteTypeProbabilities[2] > 1){
		printf("error: cumulative site type probabilities exceed 1! (buildFitnessArray, line 372)");
		exit(1);
	}

    // choose which sites are under selection
if (initializer = 0){
    for ( i = 0; i < totalSitesInGenome; i++ ) {
        dum = randU();
        if ( dum < siteTypeProbabilities[0] ) {
            // for a selected site i
            // do once, at beginning of simulation
            selectionCoefficients[i] = log( 1.0 - randU() ) * -MEAN_S_BENEFICIAL;
        }
        else if ( dum < siteTypeProbabilities[1] ){
            selectionCoefficients[i] = log( 1.0 - randU() ) * MEAN_S_BG;
        }
        else if ( dum < siteTypeProbabilities[2] ) {
            selectionCoefficients[i] = log( 1.0 - randU() ) * -MEAN_S_DIVERGENT;
        }
        else {
            selectionCoefficients[i] = 0.0;
        }
     initializer++;
}
}

void calculateMetricsAndStats(void)
{
    int i, j, k, locus, n, deme, alleleSum, numVariableInWindow;
    short int *gtpt;
    double *dpt1, *dpt2, p1, p2, pi, kbinv, *dpt3, pg;
    double dum, ninv, Dxy, FST;
    double pitotal, piwithin, n1frac, n2frac;
    int totalAlleleCount, demeAlleleCount;
    _Bool atLeastOneVariable;
    
    dpt1 = alleleFrequencies;
    ninv = 1.0 / ((double) (2 * TOTAL_N));
    kbinv = 1.0 / ((double) nSITES_PER_WINDOW);
    
    // allele frequencies and by deme:
    for ( i = 0; i < totalSitesInGenome; i++ ) {
        // go sites 1 by 1
        dpt2 = alleleFrequenciesByDeme + i;  // pointer to allele frequency for site i in deme 0
        
        if ( isVariableSite[i] ) { // if is variable site in population

            gtpt = genotypes + (2 * i); // pointer to first allele of this site in first individual
            totalAlleleCount = 0;
            
            for ( deme = 0; deme < nDEMES; deme++ ) {
                demeAlleleCount = 0;
                n = demeCounts[deme];
                for ( j = 0; j < n; j++ ) {
                    alleleSum = (*gtpt) + (*(gtpt + 1));
                    totalAlleleCount += alleleSum;
                    demeAlleleCount += alleleSum;
                    gtpt += 2 * totalSitesInGenome;
                }
                
                *dpt2 = ((double) demeAlleleCount) / ((double) (2 * n)); // 2 for diploid
                dpt2 += totalSitesInGenome;
            }
            
            if ( totalAlleleCount == 2 * TOTAL_N ) { // 2 for diploid
                isVariableSite[i] = 0;
                fixedAllele[i] = 1;
                *dpt1 = 1.0;
                fprintf(fixationLog, "%li,%i,1,%i\n", t, i, regionMembership[i]);
            }
            else if ( totalAlleleCount == 0 ) {
                isVariableSite[i] = 0;
                fixedAllele[i] = 0;
                fprintf(fixationLog, "%li,%i,0,%i\n", t, i, regionMembership[i]);
                *dpt1 = 0.0;
            }
            else
                *dpt1 = ((double) totalAlleleCount) * ninv;
            
        }
        else {
          if ( fixedAllele[i] )
                dum = 1.0;
            else
                dum = 0.0;
            
            *dpt1 = dum;
            for ( j = 0; j < nDEMES; j++ ) {
                *dpt2 = dum;
                dpt2 += totalSitesInGenome;
            }
        }
        
        dpt1++;
    }
    
    if ( t % TS_SAMP_FREQ == 0 ) {
        
        if ( nDEMES > 2 ) {
            fprintf(stderr, "\nError in calculateMetricsAndStats(): methods assume only two demes!\n");
            exit(-1);
        }
        
        fprintf(alleleFrequenciesOverTime, "%li", t);
        for ( i = 0; i < totalSitesInGenome; i++ ) {
            fprintf(alleleFrequenciesOverTime, ",%E", alleleFrequencies[i]);
        }
        dpt2 = alleleFrequenciesByDeme;
        for ( j = 0; j < nDEMES; j++ ) {
            for ( i = 0; i < totalSitesInGenome; i++ ) {
                fprintf(alleleFrequenciesOverTime, ",%E", *dpt2);
                dpt2++;
            }
        }
        fprintf(alleleFrequenciesOverTime,"\n");
        
    
        // Dxy, pi, Da
        fprintf(piOverTime, "%li", t);
        fprintf(DxyOverTime, "%li", t);
        fprintf(DaOverTime, "%li", t);
        fprintf(FSToverTime, "%li", t);
        
        dpt1 = alleleFrequenciesByDeme; // deme 1
        dpt2 = alleleFrequenciesByDeme + totalSitesInGenome; // deme 2
        dpt3 = alleleFrequencies; // global
        locus = 0;
        
        n1frac = ((double) demeCounts[0]) / ((double) TOTAL_N);
        n2frac = ((double) demeCounts[1]) / ((double) TOTAL_N);
        
        for ( i = 0; i < nGENOMIC_REGIONS; i++ ) {
            for ( j = 0; j < nWINDOWS_PER_REGION; j++ ) {
                pi = 0.0;
                Dxy = 0.0;
                piwithin = 0.0;
                pitotal = 0.0;
                atLeastOneVariable = 0;
                for ( k = 0; k < nSITES_PER_WINDOW; k++ ) {
                    if ( isVariableSite[locus] ) {
                        atLeastOneVariable = 1;
                        p1 = *dpt1;
                        p2 = *dpt2;
                        pg = *dpt3;
                        
                        // pi metric
                        pi += (p1 * (1.0 - p1)) + (p2 * (1.0 - p2));
                        
                        // Dxy
                        Dxy += (p1 * (1.0 - p2)) + (p2 * (1.0 - p1));
                        
                        // FST
                        piwithin += (n1frac * (p1 * (1.0 - p1))) + (n2frac * (p2 * (1.0 - p2)));
                        pitotal += pg * (1.0 - pg);
                        
                        if ( piwithin > pitotal ) {
                            fprintf(stderr, "\nWarning: piwithin = %E > pitotal = %E\n\t p1 = %E, p2 = %E, pg = %E\n", piwithin, pitotal, p1, p2, pg);
                        }
                        
                    }
                    dpt1++;
                    dpt2++;
                    dpt3++;
                    locus++;
                }
                pi *= kbinv;
                fprintf(piOverTime, ",%E", pi);
                
                Dxy *= kbinv;
                fprintf(DxyOverTime, ",%E", Dxy);
                
                fprintf(DaOverTime, ",%E", (Dxy - pi));
                
                if ( atLeastOneVariable ) {
                    FST = (pitotal - piwithin)/pitotal;
                    if ( FST < 0.0 ) {
                        fprintf(stderr, "\nWarning: FST = %E < 0.0\n", FST);
                    }
                    fprintf(FSToverTime, ",%E", FST);
                }
                else {
                    fprintf(FSToverTime, ",NA");
                }
                
            }
        }
        fprintf(piOverTime, "\n");
        fprintf(DxyOverTime, "\n");
        fprintf(DaOverTime, "\n");
        fprintf(FSToverTime, "\n");
        
        i = 0;
        for ( locus = 0; locus < totalSitesInGenome; locus++ ) {
            if ( isVariableSite[locus] )
                i++;
        }
        fprintf(variableLociOverTime, "%li,%i", t, i);
        
        locus = 0;
        for ( i = 0; i < nGENOMIC_REGIONS; i++ ) {
            for ( j = 0; j < nWINDOWS_PER_REGION; j++ ) {
                numVariableInWindow = 0;
                for ( k = 0; k < nSITES_PER_WINDOW; k++ ) {
                    if ( isVariableSite[locus] )
                        numVariableInWindow++;
                    locus++;
                }
                fprintf(variableLociOverTime, ",%i", numVariableInWindow);
            }
        }
        fprintf(variableLociOverTime, "\n");

    }
    
    
}


void finalPrints(void)
{
    FILE *paramPrints;
    
    paramPrints = fopen("parameters.m", "w");
    
    fprintf(paramPrints, "codeVersion = '%s';\n", version);


    
    fprintf(paramPrints, "nDEMES = %i;\nTOTAL_N = %i;\nTOTAL_GENERATIONS = %li;\n", nDEMES, TOTAL_N, TOTAL_GENERATIONS);
    fprintf(paramPrints, "nGENOMIC_REGIONS = %i;\nnWINDOWS_PER_REGION = %i;\nWINDOW_SIZE_BASES = %i;\n", nGENOMIC_REGIONS, nWINDOWS_PER_REGION, WINDOW_SIZE_BASES);
    fprintf(paramPrints, "nSITES_PER_WINDOW = %i;\nMIGRATION_RATE = %E;\nMEAN_S_BG = %E;\n", nSITES_PER_WINDOW, MIGRATION_RATE, MEAN_S_BG);
    fprintf(paramPrints, "MEAN_S_DIVERGENT = %E;\nMEAN_S_BENEFICIAL = %E;\n PROB_MUT_BG = %E;\n", MEAN_S_DIVERGENT, MEAN_S_BENEFICIAL, PROB_MUT_BG);
    fprintf(paramPrints, "PROB_MUT_DIVERGENT = %E;\nPROB_MUT_BENEFICIAL = %E;\nMU_HIGH = %E;\n", PROB_MUT_DIVERGENT, PROB_MUT_BENEFICIAL, MU_HIGH);
    fprintf(paramPrints, "MU_LOW = %E;\nDURATION_ALLOPATRY = %i;\nDETERMINISTIC = %i;\n", MU_LOW, DURATION_ALLOPATRY, DETERMINISTIC);
    fprintf(paramPrints, "HIGH_RECOMB_PROB_PER_KB = %E;\nLOW_RECOMB_PROB_PER_KB = %E;\nINITIAL_CONDITION = %i;\nTS_SAMP_FREQ = %i;\n", HIGH_RECOMB_PROB_PER_KB, LOW_RECOMB_PROB_PER_KB, INITIAL_CONDITION, TS_SAMP_FREQ);
    fprintf(paramPrints, "MODEL_TYPE = %i;\ntotalSitesInGenome = %i;\n", MODEL_TYPE, totalSitesInGenome);
    
    
    
    fclose(paramPrints);
    
}


void initializePopulation(void)
{
    int i, j, nPerDeme, count, sitesPerRegion;
    int currentRegion, currentWindow;
    short int *sipt;
    double highMapDistanceBetweenLoci, lowMapDistanceBetweenLoci, *dpt1, *dpt2;
    FILE *printMap;
    
    // set up genomic map
    totalSitesInGenome = nSITES_PER_WINDOW * nGENOMIC_REGIONS * nWINDOWS_PER_REGION;
    sitesPerRegion = nSITES_PER_WINDOW * nWINDOWS_PER_REGION;
    
    highMapDistanceBetweenLoci = ((double) WINDOW_SIZE_BASES) * ((HIGH_RECOMB_PROB_PER_KB / 1000.0) * 100.0) / ((double) nSITES_PER_WINDOW);
    lowMapDistanceBetweenLoci = ((double) WINDOW_SIZE_BASES) * ((LOW_RECOMB_PROB_PER_KB / 1000.0) * 100.0) / ((double) nSITES_PER_WINDOW);
    genomeMap = (double *) malloc( (sizeof(double) * totalSitesInGenome) );
    genomeMap[0] = 0.0;
    
    regionFlags = (int *) malloc( (sizeof(int) * totalSitesInGenome) );
    windowFlags = (int *) malloc( (sizeof(int) * totalSitesInGenome) );
    regionMembership = (int *) malloc( (sizeof(int) * totalSitesInGenome) );
    windowMembership = (int *) malloc( (sizeof(int) * totalSitesInGenome) );
    mutationRateClass = (int *) malloc( (sizeof(int) * totalSitesInGenome) );
    
    currentRegion = 0;
    currentWindow = 0;
    for ( i = 0; i < totalSitesInGenome; i++ ) {
        if ( ((i + 1) % sitesPerRegion) == 1 ) {
            regionFlags[i] = REGION_FLAG_START;
        }
        else if ( ((i + 1) % sitesPerRegion) == 0 ) {
            regionFlags[i] = REGION_FLAG_END;
        }
        else
            regionFlags[i] = REGION_FLAG_MIDDLE;
        
        if ( ((i + 1) % nSITES_PER_WINDOW) == 1 ) {
            windowFlags[i] = WINDOW_FLAG_START;
        }
        else if ( ((i + 1) % nSITES_PER_WINDOW) == 0 ) {
            windowFlags[i] = WINDOW_FLAG_END;
        }
        else
            windowFlags[i] = WINDOW_FLAG_MIDDLE;
        
        if ( (i > 0) && (i % nSITES_PER_WINDOW == 0) )
            currentWindow++;
        
        windowMembership[i] = currentWindow;
        
        if ( (i > 0) && (i % sitesPerRegion == 0) )
            currentRegion++;
        
        regionMembership[i] = currentRegion;
        
        if ( (currentRegion % 2) == 1 ) {
            mutationRateClass[i] = HIGH_MUT_FLAG;
            if ( i > 0 )
                genomeMap[i] = genomeMap[(i-1)] + highMapDistanceBetweenLoci;
        }
        else {
            mutationRateClass[i] = LOW_MUT_FLAG;
            if ( i > 0 )
                genomeMap[i] = genomeMap[(i-1)] + lowMapDistanceBetweenLoci;
        }
        

    }
    
    printMap = fopen("genomeDetails.m","w");
    fprintf(printMap, "genomeMap = [");
    for ( i = 0; i < totalSitesInGenome; i++ )
        fprintf(printMap, "%f ", genomeMap[i]);
    fprintf(printMap, "];\n\n");
    
    fprintf(printMap, "regionMembership = [");
    for ( i = 0; i < totalSitesInGenome; i++ ) {
        fprintf(printMap, "%i ", regionMembership[i]);
    }
    fprintf(printMap, "];\n\nwindowMembership = [");
    for ( i = 0; i < totalSitesInGenome; i++ ) {
        fprintf(printMap, "%i ", windowMembership[i]);
    }
    fprintf(printMap, "];\n\n");
    
    fprintf(printMap, "regionFlags = [");
    for ( i = 0; i < totalSitesInGenome; i++ ) {
        fprintf(printMap, "%i ", regionFlags[i]);
    }
    fprintf(printMap, "];\n\nwindowFlags = [");
    for ( i = 0; i < totalSitesInGenome; i++ ) {
        fprintf(printMap, "%i ", windowFlags[i]);
    }
    fprintf(printMap, "];\n");
    
    
    fclose(printMap);
    //exit(0);
    
    // set up genotypes
    genotypes = (short int *) malloc( (sizeof(short int) * 2 * TOTAL_N * totalSitesInGenome ) ); // 2 for diploidy
    isVariableSite = (int *) malloc ( (sizeof(int) * totalSitesInGenome) );
    fixedAllele = (short int *) malloc ( (sizeof(short int) * totalSitesInGenome) );
    alleleFrequenciesByDeme = (double *) malloc ( (sizeof(double) * totalSitesInGenome * nDEMES) );
    alleleFrequencies = (double *) malloc ( (sizeof(double) * totalSitesInGenome) );
    alleleSelectionCoeffs = (double *) malloc ( (sizeof(double) * totalSitesInGenome * 4) ); // two alleles, and separate effects in each deme
    demeLocations = (int *) malloc ( (sizeof(int) * TOTAL_N) );
    
    // genotype setup depends upon choices about initial conditions
    sipt = genotypes;
    dpt1 = alleleFrequenciesByDeme;
    dpt2 = alleleSelectionCoeffs;
    if ( INITIAL_CONDITION == INIT_NO_INITIAL_DIVERGENCE ) {
        // everybody is homogeneous for one fixed allele at every damn locus
        for ( i = 0; i < totalSitesInGenome; i++ ) {
            isVariableSite[i] = 0;
            fixedAllele[i] = 0;
            for ( j = 0; j < 2; j++ ) {
                *dpt1 = 1.0; // frequency of "0" = ancestral allele
                dpt1++;
            }
            for ( j = 0; j < 4; j++ ) {
                *dpt2 = 0.0; // neutral
                dpt2++;
            }
        }
        // set up genotypes
        for ( i = 0; i < TOTAL_N; i++ ) {
            for ( j = 0; j < (2 * totalSitesInGenome); j++ ) { // 2 for diploid
                *sipt = 0;
                sipt++;
            }
        }

    
        // set up deme locations
        nPerDeme = TOTAL_N / nDEMES;
        j = 0;
        for ( i = 0; i < (nDEMES-1); i++ ) {
            demeCounts[i] = nPerDeme;
            j += nPerDeme;
        }
        demeCounts[(nDEMES-1)] = TOTAL_N - j;
        
        count = 0;
        for ( i = 0; i < nDEMES; i++ ) {
            for ( j = 0; j < demeCounts[i]; j++ ) {
                demeLocations[count] = i;
                count++;
            }
        }
        if ( count != TOTAL_N ) {
            fprintf(stderr, "\nError in initializePopulation(): deme count (= %i) != TOTAL_N (= %i)\n", count, TOTAL_N);
            exit(-1);
        }
        
    }
    else if ( INITIAL_CONDITION == INIT_NEUTRAL_STANDING_HWE ) {
        fprintf(stderr, "\nError in initializePopulation() !\n Procedures for INIT_NEUTRAL_STANDING_HWE undefined!\n");
        exit(-1);
    }
    else {
        fprintf(stderr, "\nError in initializePopulation() !\n INITIAL_CONDITION = %i undefined!\n", INITIAL_CONDITION);
        exit(-1);
    }
    
}


void makeOffspring(int momi, int dadi, short int *offspringGTpt)
{
    int i, j, csome, parent;
    double nextCrossover;
    short int *sipt, *parentpt, newAllele;
    static long int highMutCountdown = -1;
    static long int lowMutCountdown = -1;
    _Bool addMutation;
    
    if ( highMutCountdown == -1 ) // initialization flag
        highMutCountdown = (long int) ( (log(1.0 - randU())) * (-1.0/MU_HIGH) );
    if ( lowMutCountdown == -1 )  // initialization flag
        lowMutCountdown = (long int) ( (log(1.0 - randU())) * (-1.0/MU_LOW) );

    
    for ( i = 0; i < 2; i++ ) { // 2 is for two parents
        
        // set offspring genotype pointer to right spot
        sipt = offspringGTpt + i;
        
        // get the right parent index
        if ( i == 0 )
            parent = momi;
        else
            parent = dadi;
        
        parentpt = genotypes + (2 * parent * totalSitesInGenome);
        
        csome = 0; // which haplotype is used from the parent
        for ( j = 0; j < totalSitesInGenome; j++ ) {
            if ( regionFlags[j] == REGION_FLAG_START ) {
                // pick a starting chromosome at random
                if ( randU() < 0.5 )
                    csome = 1;
                else
                    csome = 0;
                
                // next crossover
                nextCrossover = genomeMap[j] + (log(1.0 - randU())) * (-50.0); // crossovers are i.i.d. with a mean of 50 cM between them
            }
            else {
                if ( genomeMap[j] > nextCrossover ) {
                    do {
                        // actuate crossover; do-while allows for multiple crossovers in short spans (unlikely but possible)
                        if ( csome )
                            csome = 0;
                        else
                            csome = 1;
                        
                        nextCrossover += (log(1.0 - randU())) * (-50.0);
                        
                    } while ( genomeMap[j] > nextCrossover );
                }
            }
            
            addMutation = 0;
            if ( mutationRateClass[j] == HIGH_MUT_FLAG ) {
                if ( highMutCountdown < 1 ) {
                    addMutation = 1;
                    highMutCountdown = (long int) ( (log(1.0 - randU())) * (-1.0/MU_HIGH) );
                }
                else
                    highMutCountdown--;
            }
            else {
                if ( lowMutCountdown < 1 ) {
                    addMutation = 1;
                    lowMutCountdown = (long int) ( (log(1.0 - randU())) * (-1.0/MU_LOW) );
                }
                else
                    lowMutCountdown--;
            }
            
            // now, we know which chromosome to use, and if we need to add a mutation or not!
            newAllele = *(parentpt + csome); // get allele from parent
            
            if ( addMutation )
                newAllele = addAnewMutation(j, newAllele); // actuate mutation and update globals accordingly
            
            *sipt = newAllele; // put allele in offspring's genome
            
            sipt += 2;
            parentpt += 2;
            
        }
    }
    
    
}


void migration(int *demeIndexes)
{
    int i, j, nmig = 0;
    FILE *finalLocations;
    
    for ( i = 0; i < nDEMES; i++ )
        demeCounts[i] = 0;
    
    if ( nDEMES != 2 ) {
        fprintf(stderr, "\nError in migration()!  Not set up for more than 2 demes!\n");
        exit(-1);
    }
        
    for ( i = 0; i < TOTAL_N; i++ ) {
        if ( randU() < MIGRATION_RATE ) {
            nmig++;
            if ( demeLocations[i] == 1 )
                demeLocations[i] = 0;
            else
                demeLocations[i] = 1;
        }
        if ( demeLocations[i] == 1 ) {
            *(demeIndexes + TOTAL_N + demeCounts[1]) = i;
            demeCounts[1] = demeCounts[1] + 1;
        }
        else {
            *(demeIndexes + demeCounts[0]) = i;
            demeCounts[0] = demeCounts[0] + 1;
        }
    }
    
    if ( (demeCounts[0] + demeCounts[1]) != TOTAL_N ) {
        fprintf(stderr, "\nError in migration()!  Counts don't add up!\n");
        exit(-1);
    }
    
    if ( t == TOTAL_GENERATIONS ) {
        
        finalLocations = fopen("finalLocations.m", "w");
        
        //fprintf(stderr, "\nmig count = %i, demeCounts[0] = %i, demeCounts[1] = %i\n", nmig, demeCounts[0], demeCounts[1]);
        // array of counts of individuals in each deme
        fprintf(finalLocations, "demeCounts = [");
        for (i = 0; i < nDEMES; i++ )
            fprintf(finalLocations, "%i ", demeCounts[i]);
        fprintf(finalLocations, "];\n");
        
        // array of deme location of each individual:
        fprintf(finalLocations, "demeLocations = [");
        for (i = 0; i < TOTAL_N; i++ )
            fprintf(finalLocations, "%i ", demeLocations[i]);
        fprintf(finalLocations, "];\n");
        
        
        // arrays of indexes of individuals in each deme
        for ( i = 0; i < nDEMES; i++ ) {
            fprintf(finalLocations, "deme%iindexes = [", i);
            for ( j = 0; j < demeCounts[i]; j++ )
                fprintf(finalLocations, "%i ", *(demeIndexes + (i * TOTAL_N) + j));
            fprintf( finalLocations, "] + 1;\n");
        }
        fclose(finalLocations);
    }

}



void openDataRecordingFiles(void)
{
    int i, j, locus;
    
    // mutation log:
    mutationLog = fopen("MutationLog.csv", "w");
    fprintf(mutationLog, "Time,Locus,NewAllele,RegionMembership,MutationRateClass\n");
    
    // fixation log:
    fixationLog = fopen("FixationLog.csv", "w");
    fprintf(mutationLog, "Time,Locus,FixedAllele,RegionMembership\n");
    
    // allele freq's over time:
    alleleFrequenciesOverTime = fopen("AlleleFreqTimeSeries.csv", "w");
    fprintf(alleleFrequenciesOverTime, "Time");
    for ( i = 0; i < totalSitesInGenome; i++ ) {
        fprintf(alleleFrequenciesOverTime, ",Locus%i", i);
    }
    for ( i = 0; i < nDEMES; i++ ) {
        for ( j = 0; j < totalSitesInGenome; j++ ) {
            fprintf(alleleFrequenciesOverTime, ",Deme%iLocus%i", i, j);
        }
    }
    fprintf(alleleFrequenciesOverTime, "\n");
    
    // pi values over time: windows
    piOverTime = fopen("PiOverTime.csv", "w");
    DxyOverTime = fopen("DxyOverTime.csv", "w");
    DaOverTime = fopen("DaOverTime.csv", "w");
    FSToverTime = fopen("FSToverTime.csv", "w");
    
    fprintf(piOverTime, "Time");
    fprintf(DxyOverTime, "Time");
    fprintf(DaOverTime, "Time");
    fprintf(FSToverTime, "Time");
    
    for ( i = 0; i < nGENOMIC_REGIONS; i++ ) {
        for ( j = 0; j < nWINDOWS_PER_REGION; j++ ) {
            fprintf(piOverTime, ",Region%iWindow%i", i, j);
            fprintf(DxyOverTime, ",Region%iWindow%i", i, j);
            fprintf(DaOverTime, ",Region%iWindow%i", i, j);
            fprintf(FSToverTime, ",Region%iWindow%i", i, j);
        }
    }
    fprintf(piOverTime, "\n");
    fprintf(DxyOverTime, "\n");
    fprintf(DaOverTime, "\n");
    fprintf(FSToverTime, "\n");
    
    variableLociOverTime = fopen("VariableSitesOverTime.csv", "w");
    fprintf(variableLociOverTime, "time,totalVariable");
    for ( i = 0; i < nGENOMIC_REGIONS; i++ ) {
        for ( j = 0; j < nWINDOWS_PER_REGION; j++ ) {
            fprintf(variableLociOverTime, ",Region%iWindow%i", i, j);
        }
    }
    fprintf(variableLociOverTime, "\n");
    
}



int pickParent( double *fitnessArray, int prevChosen, int nInDeme )
{
    int i, j, n, count, dadi, momi, *ipt;
    int aliveOnes[nInDeme], nAlive=0;
    double fitnessArray[TOTAL_N];
    short int *offspringGenotyes, *sipt;
    int *offspringDemeLocations;
    offspringGenotyes = (short int *) malloc( (sizeof(short int) * 2 * totalSitesInGenome * TOTAL_N) );
    offspringDemeLocations = (int *) malloc( (sizeof(int) * TOTAL_N) );

//	momi = *(ipt + momi); // extract the right index from the overall index array; "momi" and "dadi" are "local" indexes in the consecutive fitness array
//        dadi = *(ipt + dadi);

	sipt = offspringGenotyes;

    int luckyOne;
    if ( MODEL_TYPE == MODEL_TYPE_NEUTRAL_ONLY ) {
        do {
            luckyOne = randI() % nInDeme;
        } while ( luckyOne == prevChosen );
    }
    else if ( MODEL_TYPE == MODEL_TYPE_SELECTION ) {
        //fprintf(stderr, "\nError in pickParent()! Procedures not in place for MODEL_TYPE_SELECTION\n");
        //exit(-1);

	int s_total=0;
	int s_sum;
	int individual=0;
	int i=0, j=0;
	double viabilityArray[nInDeme];
for (i; i<sitesInGenome; i++){
	s_sum += selectionCoefficients[i];
}

for(individual = 0; individual < nInDeme; individual++){
	s_total=0, i=0;
	for (j; j<sitesInGenome; j++)
		{
			for (i=0;i<2;i++){
				if (*sipt = 1){
					s_total += selectionCoefficients[j];
				}
				sipt++;
			}
		} //builds individual selection total for comparison to overall total
viabilityArray[individual]=s_total; 

}
cutoff = Cutoff(selectionCoefficients, viabilityArray);
i=0;
int deadTracker;

for (i; i<nInDeme; i++){
	if (viabilityArray[i]>cutoff){
		//aliveTracker[deadTracker] = i;
		deadTracker++;
	} //Scans through to find values less than cutoff; make 0 or 1 in living ones array instead?
}



luckyOne = randI() % deadTracker; //not nInDeme! Choose from array of live/dead
//aliveTracker currently has nInDeme elements; is the current modulus operation OK??
}

    return aliveOnes[luckyOne];
}

double Cutoff(double viabilityArray){
    double cutoffArray[nInDeme];
    memcpy(cutoffArray, viabilityArray, sizeof(double));
    qsort(cutoffArray, nInDeme, sizeof(double), compare);
    int deadMarker = CUTOFF*nInDeme, j=0;
    exit(1);

    /*for (i=0; i<nInDeme; i++){
        if (cutoffArray[i] != 0){
            alive_svalues[j] = cutoffArray[i];
            j++;
        }
    } IS THIS LOOP NECESSARY??? */ 
    return cutoffArray[deadMarker];
}

int compare(const void *p, const void *q)
{
    //obtained from wikipedia at https://en.wikipedia.org/wiki/Qsort
    int ret;
    int x = *(const int *)p;
    int y = *(const int *)q;
    
    if (x == y)
        ret = 0;
    else if (x < y)
        ret = -1;
    else
        ret = 1;

    return ret;
}

long int Poisson(double mm)
{
    /* ==================================================
     * Returns a Poisson distributed non-negative integer.
     * NOTE: use mm > 0
     * code modified from http://www.cs.wm.edu/~va/software/park/rvgs.c
     * accessed 4/12/12
     * ==================================================
     */

    double tt = 0.0;
    long int events = 0;
    
    while (tt < mm) {  // cumulative "time" between events; mm is the expected events per unit of time
        tt += -(log(1.0 - randU())); /* adding an amount of "time" from a random exponential draw with mean 1.0 unit
                                      hence, this winds up treating mm as the expected number of units of time to achieve
                                      mm events.  the number of times we have to increment tt to reach mm is the actual number of
                                      events that occurred in mm units of time + 1. x counts the number of increments. s*/
        events++;  // one more "event"
    }
    
    return (events-1); // corrects for counting one past the last "event" in the while loop
}





//############################ INSTANCE OF FITNESS ARRAY!!!!! MODIFY ONCE FITNESS ARRAY BUILDER COMPLETE ####################################
void reproduction(int *demeIndexes)
{
    int i, j, n, count, dadi, momi, *ipt;
    short int *offspringGenotyes, *sipt;
    int *offspringDemeLocations;
    
    count = 0;
    
    offspringGenotyes = (short int *) malloc( (sizeof(short int) * 2 * totalSitesInGenome * TOTAL_N) );
    offspringDemeLocations = (int *) malloc( (sizeof(int) * TOTAL_N) );
    fitnessArray = (double *) malloc( (sizeof(double)*totalSitesInGenome));
    
    sipt = offspringGenotyes;
    ipt = demeIndexes;
    
    for ( i = 0; i < nDEMES; i++ ) {
        
        n = demeCounts[i];
        
        // set up fitness array for choosing parents
        if ( MODEL_TYPE == MODEL_TYPE_SELECTION )
            buildFitnessArray( fitnessArray, ipt, i, n );// buildFitnessArray uses addres
        else if ( MODEL_TYPE != MODEL_TYPE_NEUTRAL_ONLY ) {
            fprintf(stderr, "\nError in reproduction()!  MODEL_TYPE (= %i) not recognized!\n", MODEL_TYPE);
            exit(-1);
        }
        
        // choose parents n times and make n offspring
        for ( j = 0; j < n; j++ ) {
            momi = pickParent( &selectionCoefficients[0], -1, n ); // -1 = flag for parent not chosen yet
            dadi = pickParent( &selectionCoefficients[0], momi, n );
            momi = *(ipt + momi); // extract the right index from the overall index array; "momi" and "dadi" are "local" indexes in the consecutive fitness array
            dadi = *(ipt + dadi);
            
            // /*  // test check:
//            if ( demeLocations[momi] != i || demeLocations[dadi] != i ) {
//                fprintf(stderr, "\nError in reproduction():\n\tLocations of parents not in right deme!\n\tdadi = %i, momi = %i\n", dadi, momi);
//                fprintf(stderr, "i = %i, demeLocations[momi] = %i, demeLocations[dadi] = %i\n\n", i, demeLocations[momi], demeLocations[dadi]);
//                exit(-1);
//            }
//            // */ // end test check
            
            makeOffspring(momi, dadi, sipt);
            sipt += 2 * totalSitesInGenome; // 2 for diploid; advance to next spot to start storing offspring genotype
            
            offspringDemeLocations[count] = i;
            count++;
        }
        
        ipt += TOTAL_N; // advance deme index pointer to next deme's set
    }
    
    if ( count != TOTAL_N ) {
        fprintf(stderr, "\nError in reproduction()!  Made wrong number of offspring!\n\tcount (=%i) != TOTAL_N (= %i)\n", count, TOTAL_N);
        exit(-1);
    }
    
    free(genotypes);
    genotypes = offspringGenotyes;
    
    free(demeLocations);
    demeLocations = offspringDemeLocations;
    
    //   // test check:
//    count = 0;
//    for ( i = 0; i < nDEMES; i++ ) {
//        for ( j = 0; j < demeCounts[i]; j++ ) {
//            if ( demeLocations[count] != i ) {
//                fprintf(stderr, "\nError in reproduction():\n\tOffspring count by location not consistent.\n");
//                fprintf(stderr, "deme = %i, demeCounts[%i] = %i, j = %i, count = %i, demeLocations[count] = %i\n", i, i, demeCounts[i], j, count, demeLocations[count]);
//                exit(-1);
//            }
//            count++;
//        }
//    }
    //  // end test check
    
//}



void RNGsetup(void)
{
    int seed, rcount;
    FILE *fpt;
    int stime;
    long ltime;
    FILE *rseed;
    long int i;
    
    
    
    if (DETERMINISTIC) {
        
        fpt = fopen("RnumSeed.txt","r");
        if (fpt == NULL) {
            perror("Can't open RnumSeed.txt");
            exit(-1);
        }
        
        rcount = fscanf(fpt,"%i",&seed);
        if ( rcount ) {
            seedRand(seed);	// fixed random number seed
            fclose(fpt);
        }
        else {
            fprintf(stderr, "\n\n\tError! nothing read from file! Exiting!\n\n");
            exit(-1);
        }
        //fprintf(stderr, "\n\nSeed = %i\n\n",seed);
    }
    else {
        /* use calendar time to seed random number generator.  Code adopted from Schildt's textbook */
        
        /* get the calendar time */
        ltime=time(NULL);
        stime=(unsigned) ltime/2;
        
        // generate and store random number seed
        rseed = fopen("RnumSeed.txt","w");
        fprintf(rseed,"%i\n",stime);
        fclose(rseed);
        seedRand(stime);	// get random number seed (system time)
    }
    
    // warm up
    for ( i = 0; i < 1000000; i++ )
        randU();
}
}



void usage(char *s)
{
    fprintf(stderr,  "\nUsage: %s [options]\n\n", s); // info print
	
    
    
	fprintf(stderr,  "\t[-?]\n\t\tPrint this list of options to screen and exit the program.\n\n");
    
}

void useShortTestValuesOfParameters(void)
{
    TOTAL_N = 100;
    TOTAL_GENERATIONS = 10000;
    nGENOMIC_REGIONS = 2;
    nWINDOWS_PER_REGION = 2;
    WINDOW_SIZE_BASES = 2000;
    nSITES_PER_WINDOW = 10;
}


double walker(double* dum, int index){
int i=0;
	for (i, i++, i<index){
		walker++;
	}
	return walker;
}
