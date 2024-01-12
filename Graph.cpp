#include "Graph.h"

#include <random>
#include <algorithm>

Graph::Graph() {
    numNodes = 0;
    numEdges = 0;
    totWeight = 0;
    maxWeight = 0;
}

Graph::Graph(int numNodes) {
    this->numNodes = numNodes;
    numEdges = 0;
    totWeight = 0;
    maxWeight = 0;

    numInfNeighs.reserve(numNodes);
    potStrains.reserve(numNodes);
    state.reserve(numNodes);
    adjM.reserve(numNodes);
    immunity.reserve(numNodes);

    vector<int> emptyRow(numNodes);
    vector<int>emptyBS(DNALen,0);
    for (int node = 0; node < numNodes; ++node) {
        numInfNeighs.push_back(0);
        state.push_back(0);
        adjM.push_back(emptyRow);
        potStrains.push_back((vector<int>) {});
        potStrains[node].reserve(numNodes);
        immunity.emplace_back(emptyBS);
    }
}

int Graph::fill(const vector<int> &weights, bool diag) {
    int nn = quadForm(1, -1, (-1) * (int) weights.size() * 2);
    if (nn != numNodes) {
        cout << "ERROR!  numNodes does not match length of vector provided to fill(...) method!!" << endl;
        numNodes = nn;
    }

    int idx = 0;
    if (diag) {
        int col;
        for (int iter = 1; iter < numNodes; ++iter) {
            for (int row = 0; row < numNodes - iter; ++row) {
                col = row + iter;
                adjM[row][col] = weights[idx];
                adjM[col][row] = weights[idx];
                if (adjM[row][col] > 0) numEdges += 1;
                totWeight += weights[idx];
                if (weights[idx] > maxWeight) maxWeight = weights[idx];
                idx += 1;
            }
        }
        if (idx != weights.size()) {
            cout << "ERROR! Diag fill completed incorrectly!" << endl;
        }
    } else {
        for (int row = 0; row < numNodes; row++) {
            for (int col = row + 1; col < numNodes; col++) {
                adjM[row][col] = weights[idx];
                adjM[col][row] = weights[idx];
                if (adjM[row][col] > 0) numEdges += 1;
                totWeight += weights[idx];
                idx += 1;
                if (weights[idx] > maxWeight) maxWeight = weights[idx];
            }
        }
    }
    return 0;
}

int Graph::quadForm(int A, int B, int C) {
    return ((-1) * B + (int) sqrt(B * B - 4 * A * C)) / (2 * A);
}

/**
 * This method simulates one epidemic on this Graph (network) by initially infecting p0 and allowing the epidemic
 * to spread along the edges of the network.  Each susceptible individual has likelihood alpha of being infected by
 * each adjacent neighbour that is infected every time step.  The epidemic continues until there are no more
 * infected nodes.  Infected individuals remain infected for one time step before becoming removed.
 *
 * @param p0            Patient zero.
 * @param alpha         Likelihood of infection spreading from am infected to and adjacent susceptible node.
 * @param epiProfile    The epidemic profile generated by this epidemic.
 * @param totInf        Total number of infected individuals over the course of the epidemic.
 * @return              The length of the epidemic.
 */
int Graph::SIR(int p0, double alpha, vector<int> &epiProfile, int &totInf) {
    int curInf, epiLen;

    // Set all to susceptible
    for (int i = 0; i < numNodes; i++) {
        state[i] = 0;
    }

    epiLen = 0;
    totInf = 0;
    state[p0] = 1;
    curInf = 1;
    epiProfile[0] = 1;
    totInf += curInf;

    while (curInf > 0) { // While there are still infected people
        for (int i = 0; i < numNodes; ++i) {
            numInfNeighs[i] = 0;
        }

        // Determine number of infected neighbours for each node
        for (int node = 0; node < numNodes; ++node) {
            if (state[node] == 1) {
                for (int neigh = 0; neigh < numNodes; neigh++) {
                    if (neigh != node && adjM[node][neigh] > 0) {
                        numInfNeighs[neigh] += adjM[node][neigh];
                    }
                }
            }
        }

        // Determine which susceptible nodes get infected
        for (int node = 0; node < numNodes; node++) {
            if (state[node] == 0 && numInfNeighs[node] > 0) {
                if (infect(numInfNeighs[node], alpha)) {
                    state[node] = 3;
                }
            }
        }

        // Update the state of each node
        curInf = 0;
        for (int node = 0; node < numNodes; ++node) {
            switch (state[node]) {
                case 0: // Susceptible
                    break;
                case 1: // Infected
                    state[node] = 2; // -> Removed
                    break;
                case 2: // Removed
                    break;
                case 3: // Newly Infected
                    state[node] = 1; // -> Infected
                    curInf += 1;
                    break;
            }
        }

        epiLen += 1; // Next day
        totInf += curInf;
        epiProfile[epiLen] = curInf;
    }
    return epiLen;
}

bool Graph::infect(int numInfNeighs, double alpha) {
    double beta = 1.0 - exp(numInfNeighs * log(1 - alpha));
    if (drand48() < beta) return true;
    return false;
}

int Graph::SIRwithVariants(int p0, double varAlphas[], bool coupled, double newVarProb, int &varCnt, int maxVars,
                           int maxLen, vector<int> varProfs[], vector<vector<int>> &varDNAs, int varParents[],
                           int varStarts[], int varInfSeverity[], int initBits, int minEdits, int maxEdits,
                           double alphaDelta, int &totInf, int fadingImmunity, int immuStr) {
    int curInf;
    int epiLen = 0;
    int curVarInf[maxVars];
    totInf = 0;
    //int immuStrength = fadingImmunity ? immuStr : 1;
    int immuStrength;
    if(fadingImmunity>0){ // fading immunity
        immuStrength = immuStr;
    }
    else if(fadingImmunity==0){ // static immunity
        immuStrength = 1;
    }
    else{
        immuStrength = 0; // no immunity
    }

    // Stores all the indices of the variant strings for generating initial/new varDNAs
    vector<int> randIdxVector(DNALen);
    for (int idx = 0; idx < DNALen; ++idx) {
        randIdxVector[idx] = idx;
        varInfSeverity[idx] = 0;
    }

    for (int node = 0; node < numNodes; ++node) {
        state[node] = -1; // Susceptible
        std::fill(immunity[node].begin(), immunity[node].end(),0); // Zeros in Immunity Strings
        std::fill(varDNAs[node].begin(), varDNAs[node].end(),0); // Zeros in Variant Strings
    }

    shuffle(randIdxVector.begin(), randIdxVector.end(), std::mt19937(std::random_device()()));

    //creating first variant with random 1's
    for (int idx = 0; idx < initBits; ++idx) {
        varDNAs[0][randIdxVector[idx]]=1;
    }


    for (int var = 0; var < maxVars; ++var) {
        curVarInf[var] = 0;
        varProfs[var] = (vector<int>) {};
        varProfs[var].reserve(maxLen);
        varParents[var] = -2;
        varStarts[var] = -1;
    }

    varCnt = 0;
    state[p0] = 0; // Infected with variant 0 (initial variant)
    curInf = 1;
    curVarInf[0] = 0;
    varParents[0] = -1; // Initial variant
    varProfs[0] = {1};
    varStarts[0] = 0;
    if(fadingImmunity>0) {
        immunityUpdate(immunity[0], varDNAs[0], immuStrength + 1);
    }
    else{
        immunityUpdate(immunity[0], varDNAs[0], immuStrength);
    }
    //counting 1's
    varInfSeverity[count(varDNAs[0].begin(), varDNAs[0].end(), 1)]++;

    while (curInf > 0 && epiLen < maxLen) {
        if (epiLen > maxLen) {
            cout << "ERROR!!" << endl;
        }
        epiLen += 1;

        // Determine infected neighbours
        for (int node = 0; node < numNodes; ++node) {
            potStrains[node] = (vector<int>) {};
        }
        for (int from = 0; from < numNodes; ++from) {
            if (state[from] >= 0) { // Infected
                for (int to = 0; to < numNodes; ++to) {
                    if(from != to){
                        if (adjM[from][to] > 0 && state[to] == -1) {
                            potStrains[to].push_back(state[from]);
                        }
                    }
                }
            }
        }

        // Determine which susceptible nodes get infected
        for (int node = 0; node < numNodes; ++node) {
            if (!potStrains[node].empty()) {
                state[node] = variantInfect(immunity[node], varAlphas, potStrains[node], varDNAs, varInfSeverity,
                                            maxVars, coupled);

                if (state[node] < -1) { // Infected
                    int curStrainID = state[node] + maxVars + 1;
                    immunityUpdate(immunity[node],varDNAs[curStrainID], immuStrength+2);
                    if (newVarProb > 0 && drand48() < newVarProb) { // New Variant
                        varCnt += 1;
                        if (varCnt == maxVars) {
                            cout << "ERROR!! Variant count exceeds the maximum!" << endl;
                        }
                        newVariant(varDNAs[curStrainID], varAlphas[curStrainID], varDNAs[varCnt], varAlphas[varCnt],
                                   randIdxVector, minEdits, maxEdits, alphaDelta, coupled);
                        state[node] = varCnt - maxVars - 1;
                        if(fadingImmunity>0){
                            immunityUpdate(immunity[node], varDNAs[varCnt], immuStrength + 2);
                        }
                        else{
                            immunityUpdate(immunity[node], varDNAs[varCnt], immuStrength);
                        }

                        varStarts[varCnt] = epiLen;
                        varParents[varCnt] = curStrainID;
                        //varProfs[varCnt] = {1};
                    }
                }
            }
        }

        // Update
        totInf += curInf;
        curInf = 0;
        for (int node = 0; node < numNodes; ++node) {
            if (state[node] < -1) { // Newly Infected
                state[node] = state[node] + 1 + maxVars; // -> Current Strain
                curVarInf[state[node]] += 1;
                curInf += 1;
            } else if (state[node] >= 0) { // Infected
                state[node] = -1; // -> Susceptible
            }
            // Otherwise Susceptible
        }
        for (int var = 0; var <= varCnt; ++var) {
            if (curVarInf[var] > 0) {
                varProfs[var].push_back(curVarInf[var]);
                int variantStart = 0;
                int variantEnd = varProfs[var].size();
                if(variantEnd-variantStart>immuStrength){
                    variantStart = variantEnd - immuStrength-1;
                }
                int sum = 0;
                for(int i= variantStart;i<variantEnd;i++){
                    sum+=varProfs[var][i];
                }
                if(sum>numNodes){
                    cout << "Error sum>numNodes" << endl;
                }

                curVarInf[var] = 0;
            }
        }
        //Decrease immunity
        if(fadingImmunity){
            for (int node = 0; node < numNodes; ++node) {
                if (state[node] == -1) {
                    decreaseImmunity(immunity[node]);
                }
            }

        }

    }
    return epiLen;
}

int Graph::variantInfect(vector<int> &immStr, double varAlphas[], vector<int> &potVars,
                         vector<vector<int>> &varStrs, int varInfSeverity[], int maxVars, bool coupled) const {
    vector<pair<double, int>> alphaIdxPairs;
    vector<double> uncoupledFullyImmuneProbs;
    uncoupledFullyImmuneProbs.reserve(numNodes);
    alphaIdxPairs.reserve(numNodes);

    for (int varIdx: potVars) {
        int varOnes = (int) count(varStrs[varIdx].begin(), varStrs[varIdx].end(), 1);
        int badOnes=0;
        for(int i=0;i<varStrs[varIdx].size();i++){
            if(varStrs[varIdx][i]==1 && immStr[i]==0){
                badOnes++;
            }
        }
        if (badOnes == 0) {
            // Do nothing as they cannot become infected as they are immune.
            if (!coupled){
                uncoupledFullyImmuneProbs.push_back(varAlphas[varIdx]);
            }
        } else {
            if (coupled) {
                alphaIdxPairs.emplace_back(((double) badOnes / (double) varOnes) * varAlphas[varIdx], varIdx);
            } else alphaIdxPairs.emplace_back(varAlphas[varIdx], varIdx);
        }
    }
    sort(alphaIdxPairs.begin(), alphaIdxPairs.end(), compareSeverity);
    for (auto &alphaIdxPair: alphaIdxPairs) {
        if (infect(1, alphaIdxPair.first)) {
            int badOnes = 0;
            for(int i=0;i<varStrs[alphaIdxPair.second].size();i++){
                if(varStrs[alphaIdxPair.second][i]==1 && immStr[i]==0){
                    badOnes++;
                }
            }
            varInfSeverity[badOnes]++;
            return alphaIdxPair.second - maxVars - 1;
        }
    }
    if (!coupled){
        for(double prob: uncoupledFullyImmuneProbs){
            if (infect(1, prob)){
                varInfSeverity[0]++;
                return -1;
            }
        }
    }
    return -1;
}

bool Graph::compareSeverity(pair<double, int> severity1, pair<double, int> severity2) {
    return severity1.first > severity2.first;
}

int
Graph::newVariant(vector<int> &origVar, const double &origVarAlpha, vector<int> &newVar, double &newVarAlpha,
                  vector<int> &rndIdxVec, int minEdits, int maxEdits, double alphaDelta, bool coupled) {
    shuffle(rndIdxVec.begin(), rndIdxVec.end(), std::mt19937(std::random_device()()));
    int numEdits = (int) lrand48() % (maxEdits - minEdits + 1) + minEdits;
    newVar = origVar;
    for (int idx = 0; idx < numEdits; ++idx) {
        vectorFlip(newVar,rndIdxVec[idx]);
    }
    if (coupled) {
        newVarAlpha = origVarAlpha;
    } else {
        double fullChangeRange;
        // Find value in range [0.0, 2 * alphaDelta] where alphaDelta is the limit alpha can change by.
        do {
            fullChangeRange = drand48();
        } while (fullChangeRange > alphaDelta * 2);
        // newVarAlpha is in range [origVarAlpha - alphaDelta, origVarAlpha + alphaDelta].
        newVarAlpha = origVarAlpha + (fullChangeRange - alphaDelta);
        if (newVarAlpha < 0.1) newVarAlpha = 0.1;
        if (newVarAlpha > 0.9) newVarAlpha = 0.9;
    }
    return 0;
}

int Graph::print(ostream &out) {
    out << "Nodes: " << numNodes << endl;
    out << "Edges: " << numEdges << endl;
    out << "Tot Weight: " << totWeight << endl;
    out << "Max Weight: " << maxWeight << endl;
    out << "W Hist: ";

    for (int v: weightHist()) {
        out << v << " ";
    }
    out << endl;

    for (int from = 0; from < numNodes; ++from) {
        for (int to = 0; to < numNodes; ++to) {
            if (adjM[from][to] > 0) { // Weight > 0
                for (int cnt = 0; cnt < adjM[from][to]; ++cnt) { // One print for each weight of edge
                    out << to << " ";
                }
            }
        }
        out << endl;
    }
    return 0;
}

vector<int> Graph::weightHist() {
    vector<int> rtn(maxWeight + 1);

    for (int weight = 0; weight <= maxWeight; ++weight) {
        rtn[weight] = 0;
    }

    for (int from = 0; from < numNodes; ++from) {
        for (int to = from + 1; to < numNodes; to++) {
            rtn[adjM[from][to]]++;
        }
    }
    return rtn;
}

vector<int> Graph::fill(const string filename) {
    vector<int> vals;

    numNodes = 0;
    numEdges = 0;
    totWeight = 0;

    ifstream infile(filename);
    string line;
    infile >> numNodes;
    getline(infile, line);
    getline(infile, line);

    adjM.reserve(numNodes);
    vector<int> row(numNodes);
    for (int i = 0; i < numNodes; ++i) {
        adjM.push_back(row);
    }

    int from = 0;
    while (getline(infile, line)) {
        stringstream ss(line);

        int to;
        while (ss >> to) {
            if (adjM[from][to] == 0) {
                numEdges++;
            }
            adjM[from][to]++;
//            adjM[to][from]++;
            totWeight++;
        }
        from++;
    }

    int val1, val2;
    for (int row = 0; row < numNodes; ++row) {
        for (int col = row + 1; col < numNodes; ++col) {
            val1 = adjM[row][col];
            val2 = adjM[col][row];
            if (val1 != val2) {
                totWeight += abs(val1 - val2);
                adjM[row][col] = max(val1, val2);
                adjM[col][row] = max(val1, val2);
            }
        }
    }

    vals.reserve(numNodes * (numNodes - 1) / 2);
    int col;
    int idx = 0;
    for (int iter = 1; iter < numNodes; ++iter) {
        for (int row = 0; row < numNodes - iter; ++row) {
            col = row + iter;
            vals.push_back(adjM[row][col]);
            idx++;
        }
    }
    if (idx != numNodes * (numNodes - 1) / 2) {
        cout << "ERROR!!! Made weight list from graph wrong" << endl;
    }
    return vals;
}

void Graph::vectorFlip(vector<int>&v, int pos) {
    if (v[pos] != 0 && v[pos] != 1) {
        printf("Unexpected value!");
    }
    else{
        v[pos] = v[pos] == 0 ? 1 : 0;
    }
}

//Immunity update function
void Graph::immunityUpdate(vector<int>&immunityStr, vector<int>&variantStr, int immuStrength){
    for(int i=0; i<immunityStr.size();i++){
        if(variantStr[i]==1){
            immunityStr[i]=immuStrength;
        }
    }
}

void Graph::decreaseImmunity(vector<int>&immuString){
    for(int i=0;i<immuString.size();i++){
        if(immuString[i]>0) immuString[i]--;
    }

}