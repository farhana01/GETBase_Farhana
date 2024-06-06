#include "main.h"

/**
 * Run the program: initialize a population, evolve the population for the specified number of mating events using
 * the specified fitness function, and output results.  This will be repeated runs times.
 * The command line arguments are explained below:
 * 1. Random Number Seed
 * 2. Fitness Function (0 -> Epidemic Length, 1 -> Profile Matching, 2 -> Epidemic Spread, 3 -> Epidemic Severity)
 * 3. Variants (0.0 -> No Variants, (0.0, 1.0] -> Variants with this Variant Probability)
 * 4. Initial Run Number
 * 5. Number of Runs
 * 6. Population Size
 * 7. Number of Generations/Mating Events
 * 8. Tournament Size
 * 9. Crossover Rate
 * 10. Mutation Rate
 * 11. Maximum Number of
 * Mutations (MNM) (The number of mutations will be in the range [1, MNM])
 * 12. Number of Sample Epidemics
 * 13. Number of States in the SDAs
 * If Profile Matching (aka the first argument is set to 1)
 * 14. The Path to the Profile from Working Directory (likely ./)
 * 15. The Profile Number
 * If using Variants (aka the second argument is more than 0.0 and at most 1.0)
 * 14. The Number of 1s in the Initial Variant String
 * 15. The Minimum Number of Edits for Derived Variants, inclusive
 * 16. The Maximum Number of Edits for Derived Variants, inclusive
 * 17. Coupled Immunity and Infectivity (0 -> Uncoupled, 1 -> Coupled)
 * If using Variants and Uncoupled (aka the 17th argument is 0)
 * 18. The Maximum Absolute Difference Between a Parent and Child Variant's Alpha (for Uncoupled)
 * 19. Immunity Strength (-1 -> infinite immunity, 0 -> no immunity, anything else -> immunity strength)
 * Note: The use of variants using profile matching fitness is not implemented/possible.
 *
 * @param numCmdLineArgs Number of Command Line Arguments
 * @param cmdLineArgs Those Arguments, Explained Above
 * @return Hopefully a Zero :)
 */

int main(int numCmdLineArgs, char *cmdLineArgs[]) {
    char filename[250];
    fstream runStats, expStats, readMe; // For File Output

    // Get the command line arguments and place them in their corresponding variables.
    getArgs(cmdLineArgs);

    // Create Directory for Output
    if (ctrlFitnessFctn == 0) {
        if (newVarProb > 0.0) {
            if (varCoupled) {
                sprintf(pathToOut,
                        "%sOutput - ELVarC w %.4f%%VarP, %03dPS, %06dMevs, %dTS, %03d%%CrR, %03d%%MuR, %dMNM,"
                        " %02dSEpis, %02dSt, %02dInitB, %02d-%02dEdits/", outRoot, newVarProb * 100, popsize,
                        generations, tournSize, (int) (crossoverRate * 100), (int) (mutationRate * 100),
                        maxMuts, numSampEpis, SDANumStates, initOneBits, minEdits, maxEdits);
            } else {
                sprintf(pathToOut,
                        "%sOutput - ELVarU w %.4f%%VarP, %03dPS, %06dMevs, %dTS, %03d%%CrR, %03d%%MuR, %dMNM,"
                        " %02dSEpis, %02dSt, %02dInitB, %02d-%02dEdits, %.2f%%DA/", outRoot, newVarProb * 100, popsize,
                        generations, tournSize, (int) (crossoverRate * 100), (int) (mutationRate * 100),
                        maxMuts, numSampEpis, SDANumStates, initOneBits, minEdits, maxEdits, varAlphaDelta);
            }
        } else {
            sprintf(pathToOut, "%sOutput - EL w %03dPS, %06dMevs, %dTS, %03d%%CrR, %03d%%MuR, %dMNM,"
                               " %02dSEpis, %02dSt/", outRoot, popsize, generations, tournSize,
                    (int) (crossoverRate * 100), (int) (mutationRate * 100), maxMuts, numSampEpis, SDANumStates);
        }
    } else if (ctrlFitnessFctn == 1) { // TODO: Not Yet Implemented.
        sprintf(pathToOut, "%sOutput - PM%d w %03dPS, %06dMevs, %dTS, %03d%%CrR, %03d%%MuR, %dMNM,"
                           " %02dSEpis, %02dSt/", outRoot, profileNum, popsize, generations, tournSize,
                (int) (crossoverRate * 100), (int) (mutationRate * 100), maxMuts, numSampEpis, SDANumStates);
    } else if (ctrlFitnessFctn == 2) {
        if (newVarProb > 0.0) {
            if (varCoupled) {
                sprintf(pathToOut,
                        "%sOutput - ESprVarC w %.4f%%VarP, %03dPS, %06dMevs, %dTS, %03d%%CrR, %03d%%MuR, %dMNM,"
                        " %02dSEpis, %02dSt, %02dInitB, %02d-%02dEdits/", outRoot, newVarProb * 100, popsize,
                        generations, tournSize, (int) (crossoverRate * 100), (int) (mutationRate * 100), maxMuts,
                        numSampEpis, SDANumStates, initOneBits, minEdits, maxEdits);
            } else {
                sprintf(pathToOut,
                        "%sOutput - ESprVarU w %.4f%%VarP, %03dPS, %06dMevs, %dTS, %03d%%CrR, %03d%%MuR, %dMNM,"
                        " %02dSEpis, %02dSt, %02dInitB, %02d-%02dEdits, %.2f%%DA, %02dImStr, %03d%%YgPop/", outRoot, newVarProb * 100, popsize,
                        generations, tournSize, (int) (crossoverRate * 100), (int) (mutationRate * 100), maxMuts,
                        numSampEpis, SDANumStates, initOneBits, minEdits, maxEdits, varAlphaDelta, immuStr, (int)(youngPopulation * 100));
            }
        } else {
            sprintf(pathToOut, "%sOutput - ESpr w %03dPS, %06dMevs, %dTS, %03d%%CrR, %03d%%MuR, %dMNM,"
                               " %02dSEpis, %02dSt/", outRoot, popsize, generations, tournSize,
                    (int) (crossoverRate * 100), (int) (mutationRate * 100), maxMuts, numSampEpis, SDANumStates);
        }
    } else {
        if (varCoupled) {
            sprintf(pathToOut,
                    "%sOutput - ESevVarC w %.4f%%VarP, %03dPS, %06dMevs, %dTS, %03d%%CrR, %03d%%MuR, %dMNM,"
                    " %02dSEpis, %02dSt, %02dInitB, %02d-%02dEdits/", outRoot, newVarProb * 100, popsize,
                    generations, tournSize, (int) (crossoverRate * 100), (int) (mutationRate * 100), maxMuts,
                    numSampEpis, SDANumStates, initOneBits, minEdits, maxEdits);
        } else {
            sprintf(pathToOut,
                    "%sOutput - ESevVarU w %.4f%%VarP, %03dPS, %06dMevs, %dTS, %03d%%CrR, %03d%%MuR, %dMNM,"
                    " %02dSEpis, %02dSt, %02dInitB, %02d-%02dEdits, %.2f%%DA, %02dImStr, %03d%%YgPop/", outRoot, newVarProb * 100, popsize,
                    generations, tournSize, (int) (crossoverRate * 100), (int) (mutationRate * 100), maxMuts,
                    numSampEpis, SDANumStates, initOneBits, minEdits, maxEdits, varAlphaDelta, immuStr, (int)(youngPopulation * 100));
        }
    }
    mkdir(pathToOut, 0777);

    // Determine the Location of the Different Output Files
    sprintf(filename, "%sbest%02d.dat", pathToOut, initRunNum);
    expStats.open(filename, ios::out);
    sprintf(filename, "%sreadme.dat", pathToOut);
    readMe.open(filename, ios::out);

    // Generate the Readme File
    createReadMe(readMe);
    readMe.close();

    // Okay, Let's Get Started!
    initAlg();
    cmdLineIntro(cout);

    for (int run = initRunNum; run < initRunNum + runs; run++) {
        sprintf(filename, "%srun%02d.dat", pathToOut, run);
        runStats.open(filename, ios::out);
        if (verbose) cmdLineRun(run, cout);
        initPop(); // Initialization
        report(runStats); // Initial Report
        for (int mev = 1; mev <= generations; mev++) { // Evolution
            matingEvent();
            if (mev % reportEvery == 0) { // Time to report
                if (verbose) {
                    cout << left << setw(5) << run;
                    cout << left << setw(4) << mev / reportEvery;
                }
                report(runStats);
            }
        }
        runStats.close();
        reportBest(expStats);
        cout << "Done run " << run << ".  " << runs - (run - initRunNum + 1) << " more to go. " << endl;
    }
    expStats.close();
    return (0);
}