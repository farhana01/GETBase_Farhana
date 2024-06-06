#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <string>
#include <sys/stat.h>
#include <bitset>
#include <numeric>
#include <algorithm>
#include <cassert>

#include "Graph.h"
#include "SDA.h"

// Generic Controls
int seed = 9107819;         // Random Number Seed
const bool verbose = true;  // Print Information to User?
const int numNodes = 256;

// Genetic Algorithm Controls 
int popsize = 50;
int initRunNum = 1;
int runs = 30;
int generations = 10000;                                    // Number of Mating Events
const int numReports = 100;                                 // Reporting Interval
const int reportEvery = (int) (generations / numReports);   // Report Every
int tournSize = 7;                                          // Tournament Size
double crossoverRate = 1.0;
double mutationRate = 1.0;
int maxMuts = 2;                                            // Maximum Number of Mutations

// Epidemic Controls
const int profLen = 16 + 1;     // Profile Length (Length of File and Time Step 0)
char *pathToProfile;
int profileNum = 1;
int numSampEpis = 30;           // Number of Sample Epidemics
const double alpha = 0.5;       // Probability of Infection
const int minEpiLen = 3;        // Minimum Epidemic Length
const int shortEpiRetrys = 5;   // Retries for Short Epidemics

// Self-Driving Automata (SDA) Controls
int SDANumStates = 12;
const int SDAOutLen = numNodes * (numNodes - 1) / 2;
const int SDANumChars = 2;
const int SDAMaxRespLen = 2;

// Variant Controls
int initOneBits = 32;
const int maxNumVars = numNodes;
const int maxEpiLen = 128;
int fadingImmunity;
int immuStr;
double youngPopulation;

// Genetic Algorithm Variables
SDA *SDAPop;                // Stores the Population of SDAs
vector<double> fits;        // Fitness of each SDA
bool *dead;                 // Which SDAs Failed the Necrotic Filter (true)
int ctrlFitnessFctn;        // Program Control for Which Fitness Function to Use
double globalWorstFit;

// Epidemic Variables
double *bestEpiVal;                    // Best Epidemic for each SDA
double targetProfile[profLen + 1];  // Profile dictionary

// Epidemic Variants Variables
int bestVarCount;                       // The Number of Variants in the Best Epidemic
int bestVarParents[maxNumVars];         // The Parents of Each Variant in the Best Epidemic
int bestVarStarts[maxNumVars];           // The Variant Lengths (start, end) for Variants in the Best Epidemic
vector<int> bestVarProfs[maxNumVars];   // The Variant Profiles for the Variants in the Best Epidemic
vector<vector<int>> bestVarDNAs(maxNumVars,
                                vector<int>(DNALen)); // The Variant DNA Strings for the Variants in the Best Epidemic
double bestVarAlphas[maxNumVars];
int bestVarSeverity[DNALen];
double newVarProb;                     // Probability of Generating a new Variant
bool varCoupled;
double varAlphaDelta;
int minEdits;                           // Minimum Number of Edits to New Variant String
int maxEdits;                           // Maximum Number of Edits to New Variant String

// Other Variables
char outRoot[50] = "./Output/";     // The Root Output Directory
char pathToOut[200];                // Variable to Store Path to Output for Opening File Streams
const static bool diagFill = true;  // Should the SDA Fill the UTAM Diagonally
Graph network(numNodes);
vector<int> SDAOutput;
vector<int> epiProfile;

// Method Declarations
int getArgs(char *args[]);                              // Get Command Line Arguments
void createReadMe(ostream &outStrm);                    // Make readme.dat
void initAlg();                                         // Initialize the Algorithm
void cmdLineIntro(ostream &outStrm);                    // Print Intro to User
void cmdLineRun(int run, ostream &outStrm);             // Print Column Headers to User
void initPop();                                         // Initialize the Population
bool necroticFilter();

double fitness(int idx, bool final);                    // Calculate the Fitness
double epiLenFitness(int idx, bool final);
double epiSpreadFitness(int idx, bool final);
double profileMatchingFitness(int idx, bool final);
double epiSeverityFitness(int idx, bool final);
void report(ostream &outStrm);
vector<double> calcStats(const vector<int> &goodIdxs,
                         bool biggerBetter);
void matingEvent();
vector<int> tournSelect(int size, bool decreasing);
bool compareFitness(int popIdx1, int popIdx2);
template<class T1, class T2>
int printIdxsOfVector(T1 &outp, vector<T2> vec,
                      const vector<int> &idxs,
                      const string &msg,
                      const string &sep,
                      bool newline);

void reportBest(ostream &outStrm);

int getArgs(char *args[]) {
    seed = stoi(args[1]);
    ctrlFitnessFctn = stoi(args[2]);
    newVarProb = stod(args[3]);
    // Ensure that variants are not used in combination with profile matching fitness.
    assert(ctrlFitnessFctn != 1 || newVarProb == 0.0);
    initRunNum = stoi(args[4]);
    runs = stoi(args[5]);
    popsize = stoi(args[6]);
    generations = stoi(args[7]);
    tournSize = stoi(args[8]);
    crossoverRate = stod(args[9]);
    mutationRate = stod(args[10]);
    maxMuts = stoi(args[11]);
    numSampEpis = stoi(args[12]);
    SDANumStates = stoi(args[13]);
    if (ctrlFitnessFctn == 1) { // Profile Matching
        pathToProfile = args[14];
        profileNum = stoi(args[15]);
    } else if (newVarProb > 0.0) { // Variants
        initOneBits = stoi(args[14]);
        minEdits = stoi(args[15]);
        maxEdits = stoi(args[16]);
        varCoupled = stoi(args[17]) == 1;
        if (!varCoupled) {
            varAlphaDelta = stod(args[18]);
        }
    }
    immuStr = stoi(args[19]);
    if (immuStr > 0) { // fading immunity
        fadingImmunity = 1;
    } else if (immuStr == 0) { // static immunity
        fadingImmunity = 0;
    } else { // no immunity
        fadingImmunity = -1;
    }
    youngPopulation = stod(args[20]);
    cout << "Arguments Captured!" << endl;
    return 0;
}

void createReadMe(ostream &outStrm) {
    outStrm << "This file contains the info about the files in this folder and the parameter settings for "
               "this experiment." << endl;
    outStrm << "Graph Evolution Tool." << endl;
    outStrm << "Fitness Function: ";
    if (ctrlFitnessFctn == 0) {
        outStrm << "Epidemic Length" << endl;
    } else if (ctrlFitnessFctn == 1) {
        outStrm << "Profile Matching" << endl;
        outStrm << "Profile: ";
        for (double i: targetProfile) {
            outStrm << i << " ";
        }
        outStrm << endl;
    } else if (ctrlFitnessFctn == 2) {
        outStrm << "Epidemic Spread" << endl;
    } else if (ctrlFitnessFctn == 3) {
        outStrm << "Epidemic Severity" << endl;
    }

    outStrm << "Epidemic Model Used: ";
    if (newVarProb > 0.0) {
        outStrm << "SIR with Variants" << endl;
    } else {
        outStrm << "SIR" << endl;
    }

    outStrm << "Representation: SDAs Filling an Adjacency Matrix" << endl;
    outStrm << endl;
    outStrm << "The parameter settings are as follows: " << endl;
    outStrm << "Number of sample epidemics: " << numSampEpis << endl;
    outStrm << "Alpha: " << alpha << endl;
    outStrm << "Minimum epidemic length: " << minEpiLen << endl;
    outStrm << "Re-tries for short epidemics: " << shortEpiRetrys << endl;
    outStrm << "Runs: " << runs << endl;
    outStrm << "Mating events: " << generations << endl;
    outStrm << "Population size: " << popsize << endl;
    outStrm << "Number of vertices: " << numNodes << endl;
    outStrm << "Maximum number of mutations: " << maxMuts << endl;
    outStrm << "Tournament size: " << tournSize << endl;
    outStrm << "Number of States: " << SDANumStates << endl;
    outStrm << "Variant Probability: " << newVarProb << endl;
    outStrm << "Edit Range: " << minEdits << " to " << maxEdits << endl;
    outStrm << endl;
    outStrm << "The file descriptions are as follows: " << endl;
    outStrm << "best##.dat -> the best fitness and it's associated data for one or more runs" << endl;
    outStrm << "run##.dat -> population statistics for each run" << endl;
}

void initAlg() {
    srand48(seed);

    SDAPop = new SDA[popsize];
    fits.reserve(popsize);
    bestEpiVal = new double[popsize];
    dead = new bool[popsize];
    for (int idx = 0; idx < popsize; ++idx) {
        SDAPop[idx] = SDA(SDANumStates, SDANumChars, SDAMaxRespLen, SDAOutLen);
        fits.push_back(-1);
        dead[idx] = false;
        bestEpiVal[idx] = -1;
    }
    epiProfile.reserve(maxEpiLen);
    for (int idx = 0; idx < maxEpiLen; ++idx) {
        epiProfile.push_back(-1);
    }
    SDAOutput.reserve(SDAOutLen);
    for (int i = 0; i < SDAOutLen; ++i) {
        SDAOutput.push_back(-1);
    }

    if (ctrlFitnessFctn == 1) {
        fstream inp;
        char buf[20];
        inp.open(pathToProfile, ios::in);      //open input file
        for (int i = 0; i < profLen; i++) {
            targetProfile[i] = 0;  //pre-fill missing values
        }
        targetProfile[0] = 1;  //put in patient zero
        for (int i = 1; i < profLen; i++) {      //loop over input values
            inp.getline(buf, 19);  //read in the number
            targetProfile[i] = strtod(buf, nullptr);    //translate the number
        }
        inp.close();
    }
}

void cmdLineIntro(ostream &outStrm) {
    outStrm << "Graph Evolution Tool." << endl;
    outStrm << "Fitness Function: ";

    if (ctrlFitnessFctn == 0) {
        outStrm << "Epidemic Length" << endl;
    } else if (ctrlFitnessFctn == 1) {
        outStrm << "Profile Matching" << endl;
        outStrm << "Profile: ";
        for (double i: targetProfile) {
            outStrm << i << " ";
        }
        outStrm << endl;
    } else if (ctrlFitnessFctn == 2) {
        outStrm << "Epidemic Spread" << endl;
    } else if (ctrlFitnessFctn == 3) {
        outStrm << "Epidemic Severity" << endl;
    }

    outStrm << "Epidemic Model Used: ";
    if (newVarProb > 0.0) {
        outStrm << "SIR with Variants" << endl;
    } else {
        outStrm << "SIR" << endl;
    }
    outStrm << "Representation: SDAs Filling an Adjacency Matrix" << endl;
    outStrm << "Check readme.dat for more information about parameters and output.";
    outStrm << endl;
}

void cmdLineRun(int run, ostream &outStrm) {
    outStrm << endl << "Beginning Run " << run << "." << endl;
    outStrm << left << setw(5) << "Run";
    outStrm << left << setw(4) << "RI";
    outStrm << left << setw(10) << "Mean";
    outStrm << left << setw(12) << "95% CI";
    outStrm << left << setw(10) << "SD";
    outStrm << left << setw(12) << "Best Fit";
    outStrm << left << setw(10) << "Best Epi";
    outStrm << left << setw(8) << "# Dead";
    outStrm << endl;
    outStrm << left << setw(5) << run;
    outStrm << left << setw(4) << "0";
}

void initPop() {
    // Generate the Initial Population
    for (int idx = 0; idx < popsize; idx++) {
        dead[idx] = false;
        do {
            SDAPop[idx].randomize();
            SDAPop[idx].fillOutput(SDAOutput);
        } while (necroticFilter());
        fits[idx] = fitness(idx, false);
    }

    // Determine the global worst fitness
    globalWorstFit = fits[0];
    if (ctrlFitnessFctn == 1) { // Minimizing
        for (double i: fits) {
            if (i > globalWorstFit) {
                globalWorstFit = i;
            }
        }
    } else { // Maximizing
        for (double i: fits) {
            if (i < globalWorstFit) {
                globalWorstFit = i;
            }
        }
    }
}

bool necroticFilter() {
    int count[2] = {0, 0};
    int bounds[2] = {1 * numNodes, 6 * numNodes};

    for (int val: SDAOutput) {
        count[val]++;
    }
    if (count[1] < bounds[0] || count[1] > bounds[1]) {
        return true; // DEAD
    } else {
        return false;
    }
}

double epiLenFitness(int idx, bool final) {
    int epiLen;
    int sum = 0;
    int totInf;
    bestEpiVal[idx] = 0;

    if (newVarProb == 0.0) {
        for (int epi = 0; epi < (final ? 10 * numSampEpis : numSampEpis); ++epi) {
            int epiCnt = 0;
            do {
                epiLen = network.SIR(0, alpha, epiProfile, totInf);
                epiCnt += 1;
            } while (epiLen < minEpiLen && epiCnt < shortEpiRetrys);
            sum += epiLen;
            if (epiLen > bestEpiVal[idx]) {
                bestEpiVal[idx] = epiLen;
                bestVarProfs[0] = epiProfile;
            }
        }
    } else {
        int numVars = 0;
        vector<int> varProfs[maxNumVars];
        //initialize and filling all zeros
        vector<vector<int>> varDNAs(maxNumVars, vector<int>(DNALen, 0));
        int varParents[maxNumVars];
        int varStarts[maxNumVars];
        double varAlphas[maxNumVars];
        int varInfSeverity[DNALen];
        for (int var = 0; var < maxNumVars; ++var) {
            varAlphas[var] = -1;
        }
        varAlphas[0] = alpha;
        // creating vector of granted immunity for each node
        vector<double> immunityGranted(numNodes,1.0);
        int oldNodeCount = 0;
        int youngNodeCount = numNodes; // default everyone is young;
        if(youngPopulation<1){
            youngNodeCount = numNodes*youngPopulation;
            oldNodeCount = numNodes - youngNodeCount;
            for(int i=0; i<oldNodeCount;i++){
                immunityGranted[i] = 0.5;
            }
        }

        for (int epi = 0; epi < (final ? 10 * numSampEpis : numSampEpis); ++epi) {
            int epiCnt = 0;
            do {
                epiLen = network.SIRwithVariants(0, varAlphas, varCoupled, newVarProb, numVars, maxNumVars, maxEpiLen,
                                                 varProfs, varDNAs, varParents, varStarts, varInfSeverity, initOneBits,
                                                 minEdits, maxEdits, varAlphaDelta, totInf, fadingImmunity, immuStr, immunityGranted);
                epiCnt += 1;
            } while (epiLen < minEpiLen && epiCnt < shortEpiRetrys);
            sum += epiLen;
            if (epiLen > bestEpiVal[idx]) {
                bestEpiVal[idx] = epiLen;
                if (final) {
                    bestVarCount = numVars;
                    for (int i = 0; i < DNALen; ++i) {
                        bestVarSeverity[i] = varInfSeverity[i];
                    }
                    for (int var = 0; var <= bestVarCount; ++var) {
                        bestVarStarts[var] = varStarts[var];
                        bestVarProfs[var] = varProfs[var];
                        bestVarParents[var] = varParents[var];
                        bestVarDNAs[var] = varDNAs[var];
                        bestVarAlphas[var] = varAlphas[var];
                    }
                }
            }
        }
    }
    return (double) sum / (final ? 10 * numSampEpis : numSampEpis);
}

double epiSpreadFitness(int idx, bool final) {
    int epiLen;
    int sum = 0;
    int totInf;

    bestEpiVal[idx] = 0;

    if (newVarProb == 0.0) {
        for (int epi = 0; epi < (final ? 10 * numSampEpis : numSampEpis); ++epi) {
            int epiCnt = 0;
            do {
                epiLen = network.SIR(0, alpha, epiProfile, totInf);
                epiCnt += 1;
            } while (epiLen < minEpiLen && epiCnt < shortEpiRetrys);
            sum += totInf;
            if (totInf > bestEpiVal[idx]) {
                bestEpiVal[idx] = totInf;
                bestVarProfs[0] = epiProfile;
            }
        }
    } else {
        int numVars = 0;
        vector<int> varProfs[maxNumVars];
        vector<vector<int>> varDNAs(maxNumVars, vector<int>(DNALen, 0));
        int varParents[maxNumVars];
        int varStarts[maxNumVars];
        double varAlphas[maxNumVars];
        int varInfSeverity[DNALen];
        for (int var = 0; var < maxNumVars; ++var) {
            varAlphas[var] = -1;
        }
        varAlphas[0] = alpha;

        // creating vector of granted immunity for each node
        vector<double> immunityGranted(numNodes,1.0);
        int oldNodeCount = 0;
        int youngNodeCount = numNodes; // default everyone is young;
        if(youngPopulation<1){
            youngNodeCount = numNodes*youngPopulation;
            oldNodeCount = numNodes - youngNodeCount;
            for(int i=0; i<oldNodeCount;i++){
                immunityGranted[i] = 0.5;
            }
        }

        for (int epi = 0; epi < (final ? 10 * numSampEpis : numSampEpis); ++epi) {
            int epiCnt = 0;
            do {
                epiLen = network.SIRwithVariants(0, varAlphas, varCoupled, newVarProb, numVars, maxNumVars, maxEpiLen,
                                                 varProfs, varDNAs, varParents, varStarts, varInfSeverity, initOneBits,
                                                 minEdits, maxEdits, varAlphaDelta, totInf, fadingImmunity, immuStr, immunityGranted);
                epiCnt += 1;
            } while (epiLen < minEpiLen && epiCnt < shortEpiRetrys);
            sum += totInf;
            if (totInf > bestEpiVal[idx]) {
                bestEpiVal[idx] = totInf;
                if (final) {
                    bestVarCount = numVars;
                    for (int i = 0; i < DNALen; ++i) {
                        bestVarSeverity[i] = varInfSeverity[i];
                    }
                    for (int var = 0; var <= bestVarCount; ++var) {
                        bestVarStarts[var] = varStarts[var];
                        bestVarProfs[var] = varProfs[var];
                        bestVarParents[var] = varParents[var];
                        bestVarDNAs[var] = varDNAs[var];
                        bestVarAlphas[var] = varAlphas[var];
                    }
                }
            }
        }
    }
    return (double) sum / (final ? 10 * numSampEpis : numSampEpis);
}

double epiSeverityFitness(int idx, bool final) {
    int epiLen;
    int oneSum;
    int multiSum = 0.0;
    int totInf;
    bestEpiVal[idx] = 0;

    if (newVarProb == 0) {
        cout << "Error in main file: epiSeverityFitness(...): Must have variants with epidemic severity fitness!"
             << endl;
    } else {
        int numVars = 0;
        vector<int> varProfs[maxNumVars];
        vector<vector<int>> varDNAs(maxNumVars, vector<int>(DNALen, 0));
        int varParents[maxNumVars];
        int varStarts[maxNumVars];
        double varAlphas[maxNumVars];
        int varInfSeverity[DNALen];
        for (int var = 0; var < maxNumVars; ++var) {
            varAlphas[var] = -1;
        }
        varAlphas[0] = alpha;

        // creating vector of granted immunity for each node
        vector<double> immunityGranted(numNodes,1.0);
        int oldNodeCount = 0;
        int youngNodeCount = numNodes; // default everyone is young;
        if(youngPopulation<1){
            youngNodeCount = numNodes*youngPopulation;
            oldNodeCount = numNodes - youngNodeCount;
            for(int i=0; i<oldNodeCount;i++){
                immunityGranted[i] = 0.5;
            }
        }

        for (int epi = 0; epi < (final ? 10 * numSampEpis : numSampEpis); ++epi) {
            int epiCnt = 0;
            do {
                epiLen = network.SIRwithVariants(0, varAlphas, varCoupled, newVarProb, numVars, maxNumVars, maxEpiLen,
                                                 varProfs, varDNAs, varParents, varStarts, varInfSeverity, initOneBits,
                                                 minEdits, maxEdits, varAlphaDelta, totInf, fadingImmunity, immuStr, immunityGranted);
                epiCnt += 1;
            } while (epiLen < minEpiLen && epiCnt < shortEpiRetrys);
            oneSum = 0;
            for (int i = 0; i < DNALen; i++) {
                oneSum += varInfSeverity[i] * i;
            }
            multiSum += oneSum;
            if (oneSum > bestEpiVal[idx]) {
                bestEpiVal[idx] = oneSum;
                if (final) {
                    bestVarCount = numVars;
                    for (int i = 0; i < DNALen; ++i) {
                        bestVarSeverity[i] = varInfSeverity[i];
                    }
                    for (int var = 0; var <= bestVarCount; ++var) {
                        bestVarStarts[var] = varStarts[var];
                        bestVarProfs[var] = varProfs[var];
                        bestVarParents[var] = varParents[var];
                        bestVarDNAs[var] = varDNAs[var];
                        bestVarAlphas[var] = varAlphas[var];
                    }
                }
            }
        }
    }
    return (double) multiSum / (final ? 10 * numSampEpis : numSampEpis);
}

double profileMatchingFitness(int idx, bool final) {
    int epiLen;
    double multiSum = 0.0;
    double oneSum;
    int totInf;
    bestEpiVal[idx] = MAXFLOAT;

    if (newVarProb == 0.0) {
        for (int epi = 0; epi < (final ? 10 * numSampEpis : numSampEpis); ++epi) {
            int epiCnt = 0;
            do {
                epiLen = network.SIR(0, alpha, epiProfile, totInf);
                epiCnt += 1;
            } while (epiLen < minEpiLen && epiCnt < shortEpiRetrys);
            oneSum = 0.0;
            for (int day = 0; day < max(epiLen, profLen); day++) {
                if (day >= epiLen) {
                    oneSum += pow(targetProfile[day] - 0.0, 2);
                } else if (day >= profLen) {
                    oneSum += pow(0.0 - epiProfile[day], 2);
                } else {
                    oneSum += pow(targetProfile[day] - epiProfile[day], 2);
                }
            }
            oneSum = sqrt(oneSum / max(epiLen, profLen));
            multiSum += oneSum;
            if (oneSum < bestEpiVal[idx]) {
                bestEpiVal[idx] = oneSum;
                bestVarProfs[0] = epiProfile;
            }
            bestVarProfs[0] = epiProfile;
        }
    } else {
        cout << "Error in main file: profileMatchingFitness(...): Cannot have variants with profile matching fitness!"
             << endl;
        return -1.0;
    }
    return multiSum / (final ? 10 * numSampEpis : numSampEpis);
}

double fitness(int idx, bool final) {//compute the fitness
    SDAPop[idx].fillOutput(SDAOutput);
    if (necroticFilter()) {
        dead[idx] = true;
        return globalWorstFit;
    }
    network.fill(SDAOutput, diagFill);

    if (ctrlFitnessFctn == 0) { // Epidemic Length Fitness
        return epiLenFitness(idx, final);
    } else if (ctrlFitnessFctn == 1) { // Profile Matching Fitness
        return profileMatchingFitness(idx, final);
    } else if (ctrlFitnessFctn == 2) { // Epidemic Spread Fitness
        return epiSpreadFitness(idx, final);
    } else if (ctrlFitnessFctn == 3) { // Epidemic Severity Fitness
        return epiSeverityFitness(idx, final);
    } else {
        return -1.0;
    }
}

void report(ostream &outStrm) {//make a statistical report
    vector<int> goodIdxs;
    int deaths = 0;

    for (int i = 0; i < popsize; i++) {
        if (!dead[i]) {
            goodIdxs.push_back(i);
        } else {
            deaths++;
        }
    }

//    printIdxsOfVector(cout, fits, goodIdxs, "Fits: ", "\t", true);

    vector<double> stats = calcStats(goodIdxs, ctrlFitnessFctn != 1);
    double mean = stats[0];
    double stdDev = stats[1];
    double CI95 = stats[2];
    double bestVal = stats[3];
    double worstVal = stats[4];

    // Find the Best and Update the Fitness of Dead Members of the Population
    int bestIdx = 0;
    if (ctrlFitnessFctn == 1) { // Minimizing
        for (int i = 1; i < popsize; i++) {
            if (fits[i] < fits[bestIdx]) {
                bestIdx = i; //find best fitness
            }
        }
    } else { // Maximizing
        for (int i = 1; i < popsize; i++) {
            if (fits[i] > fits[bestIdx]) {
                bestIdx = i;
            }
        }
    }
    if (worstVal != globalWorstFit) {
        for (int i = 0; i < popsize; i++) {
            if (dead[i]) { // Update to Worst Fitness
                fits[i] = worstVal;
            }
        }
        globalWorstFit = worstVal;
    }

    outStrm << left << setw(10) << mean;
    outStrm << left << setw(12) << CI95;
    outStrm << left << setw(10) << stdDev;
    outStrm << left << setw(12) << bestVal;
    outStrm << left << setw(10) << bestEpiVal[bestIdx];
    outStrm << left << setw(8) << deaths << endl;

    if (verbose) {
        cout << left << setw(10) << mean;
        cout << left << setw(12) << CI95;
        cout << left << setw(10) << stdDev;
        cout << left << setw(12) << bestVal;
        cout << left << setw(10) << bestEpiVal[bestIdx];
        cout << left << setw(8) << deaths << endl;
    }
}

vector<double> calcStats(const vector<int> &goodIdxs, bool biggerBetter) {
    double sum = 0.0;
    double bestVal = (biggerBetter ? 0.0 : MAXFLOAT);
    double worstVal = (biggerBetter ? MAXFLOAT : 0.0);

    for (int idx: goodIdxs) {
        sum += fits[idx];
        if ((biggerBetter && fits[idx] > bestVal) || (!biggerBetter && fits[idx] < bestVal)) {
            bestVal = fits[idx];
        }
        if ((biggerBetter && fits[idx] < worstVal) || (!biggerBetter && fits[idx] > worstVal)) {
            worstVal = fits[idx];
        }
    }

    double mean = sum / (double) goodIdxs.size();
    double stdDevSum = 0.0;
    for (int idx: goodIdxs) {
        stdDevSum += pow(fits[idx] - mean, 2);
    }
    double stdDev = sqrt(stdDevSum / ((double) goodIdxs.size() - 1.0));
    double CI95 = 1.96 * (stdDev / sqrt((double) goodIdxs.size()));

    return {mean, stdDev, CI95, bestVal, worstVal}; // {mean, stdDev, 95CI, best, worst}
}

void matingEvent() {
    int numMuts;
    vector<int> tournIdxs;
    // Selection
    tournIdxs = tournSelect(tournSize, ctrlFitnessFctn == 1);

    // Copy the Parents -> Children
    SDAPop[tournIdxs[0]].copy(SDAPop[tournIdxs[tournSize - 2]]);
    SDAPop[tournIdxs[1]].copy(SDAPop[tournIdxs[tournSize - 1]]);

    // Crossover
    if (drand48() < crossoverRate) SDAPop[tournIdxs[0]].twoPointCrossover(SDAPop[tournIdxs[1]]);

    // Mutation
    if (drand48() < mutationRate) {
        numMuts = (int) lrand48() % maxMuts + 1;
        SDAPop[tournIdxs[0]].mutate(numMuts);
        numMuts = (int) lrand48() % maxMuts + 1;
        SDAPop[tournIdxs[1]].mutate(numMuts);
    }

    // Reset dead SDAs
    SDAPop[tournIdxs[0]].fillOutput(SDAOutput);
    dead[tournIdxs[0]] = necroticFilter();
    SDAPop[tournIdxs[1]].fillOutput(SDAOutput);
    dead[tournIdxs[1]] = necroticFilter();

    if (!dead[tournIdxs[0]]) {
        fits[tournIdxs[0]] = fitness(tournIdxs[0], false);
    } else {
        fits[tournIdxs[0]] = globalWorstFit;
    }

    if (!dead[tournIdxs[1]]) {
        fits[tournIdxs[1]] = fitness(tournIdxs[1], false);
    } else {
        fits[tournIdxs[1]] = globalWorstFit;
    }
}

vector<int> tournSelect(int size, bool decreasing) {
    vector<int> tournIdxs;
    int idxToAdd;

    tournIdxs.reserve(size);
    if (size == popsize) {
        for (int idx = 0; idx < size; idx++) {
            tournIdxs.push_back(idx);
        }
    } else {
        do {
            idxToAdd = (int) lrand48() % popsize;
            if (count(tournIdxs.begin(), tournIdxs.end(), idxToAdd) == 0) {
                tournIdxs.push_back(idxToAdd);
            }
        } while (tournIdxs.size() < tournSize);
    }

    sort(tournIdxs.begin(), tournIdxs.end(), compareFitness);
    if (decreasing) {
        reverse(tournIdxs.begin(), tournIdxs.end());
    }
    return tournIdxs;
}

bool compareFitness(int popIdx1, int popIdx2) {
    return fits[popIdx1] < fits[popIdx2];
}

template<class T1, class T2>
int printIdxsOfVector(T1 &outp, vector<T2> vec, const vector<int> &idxs, const string &msg, const string &sep,
                      bool newline) {
    outp << msg;
    for (auto idx: idxs) {
        outp << vec[idx] << sep;
    }
    if (newline) outp << "\n";
    return 0;
}

void reportBest(ostream &outStrm) {//report the best graph
    int bestIdx;
    Graph G(numNodes);

    bestIdx = 0;
    if (ctrlFitnessFctn == 1) { // Minimizing
        for (int i = 1; i < popsize; i++) {
            if (fits[i] < fits[bestIdx]) {
                bestIdx = i;
            }
        }
    } else { // Maximizing
        for (int i = 1; i < popsize; i++) {
            if (fits[i] > fits[bestIdx]) {
                bestIdx = i;
            }
        }
    }

    outStrm << "Best Fitness: " << fits[bestIdx] << endl;

    // Print the Profiles
    fitness(bestIdx, true);
    outStrm << "Epidemic Profile" << endl;
    if (newVarProb > 0.0) {
        for (int i = 0; i <= bestVarCount; i++) {
            if (i == 0) {
                outStrm << "NA-->" << "V" << left << setw(2) << i << "\t";
            } else {
                outStrm << "V" << left << setw(2) << bestVarParents[i] << "->V" << left << setw(2) << i << "\t";
            }
            outStrm << "[" << left << setw(3) << bestVarStarts[i] << "-";
            outStrm << left << setw(3) << bestVarStarts[i] + bestVarProfs[i].size() << "]:\t";
            for (int j = 0; j < bestVarStarts[i]; j++) {
                outStrm << "\t";
            }
            for (int j: bestVarProfs[i]) {
                outStrm << j << "\t";
            }
            outStrm << endl;
        }
        for (int i = 0; i <= bestVarCount; i++) {
            outStrm << "V" << i << "\t" << left << setw(10) << bestVarAlphas[i];
            for (int j = 0; j < bestVarDNAs[i].size(); j++) {
                outStrm << bestVarDNAs[i][j];
            }
            outStrm << endl;
        }
        outStrm << "Severity Histogram: ";
        for (int i = 0; i < DNALen; i++) {
            outStrm << bestVarSeverity[i] << "\t";
        }
        outStrm << endl;
    } else {
        outStrm << left << setw(4) << "V0" << " ";
        outStrm << "[" << left << setw(2) << bestVarStarts[0] << " ";
        outStrm << bestVarProfs[0].size() << "]: ";
        for (int j: bestVarProfs[0]) {
            if (j > 0) {
                outStrm << j << " ";
            }
        }
        outStrm << endl;
    }

    // Write the SDA
    outStrm << "Self-Driving Automata" << endl;
    SDAPop[bestIdx].printSDA(outStrm);

    SDAPop[bestIdx].fillOutput(SDAOutput);
    G.fill(SDAOutput, diagFill);
    outStrm << "Graph" << endl;
    G.print(outStrm);
    outStrm << endl;
}

