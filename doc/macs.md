# Manual for `macs`
## Acquire
```bash
git clone https://github.com/gchen98/macs
cd macs
make
# or non debug mode
g++ -O3 -Wall -o marcs simulator.cpp simulator.h datastructures.cpp algorithm.cpp
mv macs msformatter ~/.local/bin # or elsewhere
```
## Generaral form of the command
```bash
macs <samplesize> <region in base pairs> [options]
```
## Options
| Switch | Arugument | Descritpion|
| -: | --- | -------------- |
| -s | *random seed* | |
| -d | | enable debugging messages|
| -i | *iterations* | |
| -h | *history*| number of previous base pairs to retain |
| -t | *mu* | mutation rate per site per 4N generations
| -F | *inputfname* $[0|1]$ | Tab delimited frequency distribution file where first column indicate range of SNP allele frequencies from previous row to current row and last column is desired bin frequency. Second parameter is 1 if SNPs with derived allele freq $> 1.0$ should have alleles flipped, 0 otherwise. |
| -r | *r* | recombination rate per site per 4N generations |
| -c | *f* $\lambda$ | $f$ = ratio of gene conversion rate to crossover rate.  tracklen $\lambda$ is mean length of tract in base pairs. |
| -R | *inputfname* | Tab delimited file where first two columns indicate range of base pair positions scaled to the unit interval and last column is ratio with respect to base line recombination rate. |
| -T | | Print each local tree in Newick format to standard out |
| -G | $\alpha$ | Assign growth rate alpha across populations where $\alpha=-\log\frac{N_p}{N_r}$ |
| -I | $n\,n_1\,n_2$ .. mig-rate | Assign all elements of the migration matrix for $n$ populations.  Values in matrix set to $\frac{\textrm{mig-rate}}{n-1}$. |
| -m | $i\,j\,m$ | Assign $i,j$-th element of migration matrix to $m$. |
| -ma | $m_{11}$ .. $m_{12}$ .. $m_{nn}$ | Assign values to all elements of migration matrix for $n$ populations. |
| -n | $i$ size | Pop $i$ has size set to size$\times N_0$
| -g | $i\,\alpha$ | If used must appear after `-m` option. |

<!-- originaly -M, I changed it to -m, as there is no -M switch -->

The following options modify parameters at time $t$.

| Switch | Arugument | Descritpion|
| -: | --- | -------------- |
| -eG | $t\,\alpha$ | Assign growth rate for all pops at time $t$. |
| -eg | $t\,i\,\alpha$ | Assign growth rate $\alpha$ of pop $i$ at time $t$. |
| -eM | $t\,m$ | Assign migrate rate $m$ for all elements of migration matrix at time $t$. |
| -em | $t\,i\,j\,m_{ij}$ | Assign migration rate for $i,j$-th element of migration matrix at time $t$. |
| -ema | $t\,n\,m_{11}$ .. $m_{12}$ .. $m_{nn}$ | Assign migration rates within the migration matrix for $n$ populations at time $t$ |
| -eN | $t$ size | New pop sizes at time $t$ for all pops where new size $= \mathrm{size}\times N_0$ |
| -en | $t\,i\,\mathrm{size}_i$ | New pop size of pop $i$ will be set to $\mathrm{size}_i\times N_0$ at time $t$. |
| -es | $t\,i\,p$ | Split two populations.  At time $t$, a proportion $p$ of chromosomes from pop $i$ will migrate to a population $i+1$. |
| -ej | $t\,i\,j$ | Join two populations.  At time $t$ all chromosomes migrate from pop $i$ to pop $j$. |

## Internal struct for scenario setup
```c++
Configuration::Configuration(){
    bDebug = false;
    bVariableRecomb = false;
    bSNPAscertainment = false;
    bFlipAlleles = false;
    bNewickFormat = false;
    bMigrationChangeEventDefined = false;
    dTheta = 0; // mutation rate parameter theta
    dRecombRateRAcrossSites = 0.0; // recombination rate parameter r
    dSeqLength = 0.;
    dGeneConvRatio = 0.;
    iGeneConvTract = 1;
    iTotalPops = 1; // total populations declared
    dBasesToTrack = 1;
    iRandomSeed = time(NULL);
    iIterations = 1;
    pAlleleFreqBinPtrSet = NULL;
}
```

## Examples
### Example 1
```bash
macs 100 1e6 -T -t .001 -r .001 -h 1e2 \
	-R example_input/hotspot.txt \
	-F example_input/ascertainment.txt 0 2>trees.txt | 
		./msformatter > haplotypes.txt
```

Simulate 100 sequences on a region $10^6$ base pairs long.
The per base pair mutation and recombination rate scaled at $4N$ is .001.
The $h$ parameter approximates the Markovian order by instructing the program to include all geneologies from the current point to $10^2$ base pairs to the left if one considers simulation proceeding from the left end to the right end of the entire sequence.
Any branches unique to the geneology that is beyond $10^2$ base pairs is pruned from the ARG.
`-T` tells `macs` to output the local trees in Newick format similar to `ms` output.

```config
# hotspot.txt
0	0.3	0.570688491664533
0.3	0.5	1.45639998785967
```

The option `-R` instructs the program to read in a variable recombination file.
The first line of the file `hotspot.txt` says that from position 0 to .3 (unit scaled on entire sequence to be simulated), the cM to Mb ratio is 0.57.
If there are any coordinates not covered by ranges in the file, the ratio defaults to 1.

```config
# ascertainment.txt
0.01 0
0.05 0.01
```

The option `-F` instructs the program to read in a SNP ascertainment file.
Entering `1` after the filename instructs the program to assume any derived allele frequency (DAF) > .5 to have a DAF of 1 - DAF.
This might be useful for scenarios where one is interested only in the minor allele frequency and the identity of the derived allele is unknown.
The flag `0` disables this behavior.
The first line of the file `ascertainment.txt` instructs the program to completely filter out all SNPs with DAF of range 0 to 0.01.
The second line says to filter SNPs with DAF from 0.01 to 0.05 to the point where SNPs in this DAF range comprise 1% of the SNPs output.

### Example 2

For very large sample sizes and/or sequence lengths, we recommend storing output in a text file to postprocess later (e.g. import into a database)

```bash
macs 10000 1e9 -t .001 -r .001 -h 1e2 \
	-F example_input/ascertainment.txt 0 2>trees.txt 1> sites.txt
```

### Example 3
This example was used in Hicky's paper to create a $N_e=100$ population.
```bash
macs 4000 100000000 -t 0.00001 -r 0.000004 \
    -eN    0.03   1.75 \
    -eN    0.06   2.00 \
    -eN    0.13   3.50 \
    -eN    0.25   5.00 \
    -eN    0.50   7.00 \
    -eN    0.75   8.20 \
    -eN    1.00   8.50 \
    -eN    1.25   9.00 \
    -eN    1.50  10.00 \
    -eN    1.75  11.00 \
    -eN    2.00  12.75 \
    -eN    2.25  13.00 \
    -eN    2.50  12.00 \
    -eN    5.00  20.00 \
    -eN    7.50  25.00 \
    -eN   10.00  30.00 \
    -eN   12.50  32.00 \
    -eN   15.00  35.00 \
    -eN   17.50  38.00 \
    -eN   20.00  40.00 \
    -eN   22.50  42.00 \
    -eN   25.00  45.00 \
    -eN   50.00  54.56 \
    -eN  100.00  73.67 \
    -eN  150.00  92.78 \
    -eN  200.00 111.90 \
    -eN  250.00 131.01 \
    -eN  500.00 226.58 \
    -eN 1000.00 417.72 \
    -eN 1500.00 608.86 \
    -eN 2000.00 800.00 2>debug.txt | msformatter > haplotypes.txt
```

### Example 4
For a pop of $N_e=1000$
```bash
macs 4000 100000000 -t 0.0001 -r 0.00004 \
       -eN   0.50  2.00 \
       -eN   0.75  2.50 \
       -eN   1.00  3.00 \
       -eN   1.25  3.20 \
       -eN   1.50  3.50 \
       -eN   1.75  3.80 \
       -eN   2.00  4.00 \
       -eN   2.25  4.20 \
       -eN   2.50  4.50 \
       -eN   5.00  5.46 \
       -eN  10.00  7.37 \
       -eN  15.00  9.28 \
       -eN  20.00 11.19 \
       -eN  25.00 13.10 \
       -eN  50.00 22.66 \
       -eN 100.00 41.77 \
       -eN 150.00 60.89 \
       -eN 200.00 80.00 \
       2>debug.txt |
    ./msformatter > haplotypes.txt
```
<!-- 

void Simulator::readInputParameters(CommandArguments arguments){
    unsigned int popId;

	int iSampleSize;
	double dDefaultPopSize,dDefaultGrowthAlpha,dDefaultMigrationRate,popSize;
	bool bAcceptFullMigrMatrix;

    unsigned int iTotalArgs = arguments.size();
    dDefaultPopSize = 1.0;
    dDefaultGrowthAlpha =0.0;

    pConfig=new Configuration();
    pConfig->iTotalPops = 1;

    EventPtrList * pEventList = new EventPtrList;

    dDefaultMigrationRate = 0.0;
    for (unsigned int i=0;i<pConfig->iTotalPops;++i){
        vector<double> newRow;
        for (unsigned int j=0;j<pConfig->iTotalPops;++j)
            newRow.push_back(dDefaultMigrationRate);
        pConfig->dMigrationMatrix.push_back(newRow);
    }

    if( arguments[0].size()!=2 ){
        cerr<< "You must enter a value for the sample size and seq length."<<endl;
        exit(1);
    }
    iSampleSize = atoi(arguments[0][0].data());
    Population newPop;
    newPop.setChrSampled(iSampleSize);
    newPop.setPopSize(dDefaultPopSize);
    newPop.setGrowthAlpha(dDefaultGrowthAlpha);

    pConfig->pPopList.push_back(newPop);
    pConfig->iSampleSize = iSampleSize;
    cerr<<"INPUT: Sample size is now "<<pConfig->iSampleSize<<endl;
	if( iSampleSize<= 0) {
	    cerr<<"First argument error. Sample size needs to be greater than 0.\n";
	    printUsage();
    }
    pConfig->dSeqLength = atof(arguments[0][1].data());
    cerr<<"INPUT: Seq length is now "<<pConfig->dSeqLength<<endl;
    set<float> eventTimes;
    for (unsigned int iCurrentArg = 1;iCurrentArg<iTotalArgs;++iCurrentArg){
        try{
            if(arguments[iCurrentArg][0][0] != '-' ) {
                cerr<<"At argument "<<iCurrentArg<<
                ", argument needs to be prefixed with a -"<<endl;
                cerr<<"You entered "<<arguments[iCurrentArg][0][0]<<endl;
                printUsage();
            }
            double dTime;
            char chType;
            EventPtr wrapper;
            short unsigned int iNoMigrPops;
            short unsigned int iMigrPops;
            short unsigned int iTotalCells;
            short unsigned int iSubOption;
            int iRunningSample=0;
            string command;
            const char * filename;
            ifstream inFile;
            bool  flipAllele;
            switch (arguments[iCurrentArg][0][1] ){
                case 'T' :
                    pConfig->bNewickFormat = true;
                    // example:
                    // (2:1.766,(4:0.505,(3:0.222,(1:0.163,5:0.163):0.059):0.283):1.261);
                    break;
                case 'd' :
                    pConfig->bDebug = true;
                    break;
                case 'h' :
                    if (arguments[iCurrentArg].size()!=2) {
                        cerr<<"For flag "<<arguments[iCurrentArg][0][1]<<
                        ", you must enter a single integer for retaining the number of previous trees\n";
                        exit(1);
                    }
                    pConfig->dBasesToTrack = atof(arguments[iCurrentArg][1].data());
                    cerr<<"INPUT: Base pairs to track is "<<pConfig->dBasesToTrack<<endl;
                    break;
                case 's' :
                    if (arguments[iCurrentArg].size()!=2) {
                        cerr<<"For flag "<<arguments[iCurrentArg][0][1]<<
                        ", you must enter a single integer for the random seed\n";
                        exit(1);
                    }
                    pConfig->iRandomSeed = atoi(arguments[iCurrentArg][1].data());
                    cerr<<"INPUT: Random seed used is "<<pConfig->iRandomSeed<<endl;
                    break;
                case 't' :  // set mutation parameter
                    if (arguments[iCurrentArg].size()!=2) {
                        cerr<<"For flag "<<arguments[iCurrentArg][0][1]<<
                        ", you must enter a single float value for the mutation parameter\n";
                        exit(1);
                    }
                    pConfig->dTheta = pConfig->dSeqLength * atof(arguments[iCurrentArg][1].data());
                    cerr<<"INPUT: Scaled mutation rate is now "<<pConfig->dTheta<<endl;
                    break;
                case 'F':
                    if (arguments[iCurrentArg].size()!=3){
                        cerr<<"For the SNP ascertainment feature you must enter the filename of the SNP "<<
                        "frequency list and whether to flip the alleles."<<endl;
                        exit(1);
                    }
                    filename = arguments[iCurrentArg][1].data();
                    flipAllele = pConfig->bFlipAlleles = atoi(arguments[iCurrentArg][2].data());
                    inFile.open(filename);
                    if (!inFile)
                        throw "SNP freq input file not found\n";
                    if (inFile.is_open()) {
                        pConfig->pAlleleFreqBinPtrSet = new AlleleFreqBinPtrSet;
                        string line;
                        int total=0;
                        double lastStart=0.;
                        double cumFreq=0.;
                        double maxFreq = flipAllele?0.5:1.;
                        while(getline(inFile,line)){
                            istringstream inStr(line);
                            double start = lastStart,end,freq;
                            inStr>>end>>freq;
                            cumFreq+=freq;
                            if (end>maxFreq) end = maxFreq;
                            if (start>=end) throw "The freq range entered is incorrect.";
                            if (pConfig->bDebug){
                                cerr<<"Frequency bin from "<<start<<" to "<<end<<" with freq "<<freq<<endl;
                            }
                            AlleleFreqBinPtr bin = AlleleFreqBinPtr(new AlleleFreqBin(start,end,freq));
                            pConfig->pAlleleFreqBinPtrSet->insert(bin);
                            lastStart = end;
                            ++total;
                        }
                        cerr<<"INPUT: Accepted "<<total<<" freqency bins"<<endl;
                        inFile.close();
                        pConfig->bSNPAscertainment = true;
                        if (cumFreq>1.0) throw "The total frequency entered exceeds one";
                        if (cumFreq<1.0){
                            if (lastStart==maxFreq){
                            }else{
                                AlleleFreqBinPtr bin = AlleleFreqBinPtr(new AlleleFreqBin(lastStart,maxFreq,1-cumFreq));
                                pConfig->pAlleleFreqBinPtrSet->insert(bin);
                                if (pConfig->bDebug){
                                    cerr<<"Added frequency bin from "<<lastStart<<" to "<<maxFreq<<" with freq "<<1-cumFreq<<endl;
                                }
                            }
                        }
                    }
                    break;
                case 'r' :
                    if (arguments[iCurrentArg].size()!=2) {
                        cerr<<"For flag "<<arguments[iCurrentArg][0][1]<<
                        ", you must enter a single float value for the recombination parameter\n";
                        exit(1);
                    }
                    pConfig->dRecombRateRAcrossSites = pConfig->dSeqLength * atof(arguments[iCurrentArg][1].data());
                    cerr<<"INPUT: Scaled recombination rate is now "<<pConfig->dRecombRateRAcrossSites<<endl;
                    break;
                case 'c' :
                    if (arguments[iCurrentArg].size()!=3) {
                        cerr<<"For flag "<<arguments[iCurrentArg][0][1]<<
                        ", you must enter the conversion to xover ratio followed by the mean tract length in bp.\n";
                        exit(1);
                    }
                    pConfig->dGeneConvRatio = atof(arguments[iCurrentArg][1].data());
                    pConfig->iGeneConvTract = atoi(arguments[iCurrentArg][2].data());
                    if (pConfig->dGeneConvRatio<0||pConfig->iGeneConvTract<0){
                        cerr<<"The gene conversion parameters must be positive\n";
                        exit(1);
                    }

                    cerr<<"INPUT: Gene conversion ratio is now "<<pConfig->dGeneConvRatio<<endl;
                    cerr<<"INPUT: Gene conversion tract length is now "<<pConfig->iGeneConvTract<<endl;
                    break;
                case 'R':
                    if (arguments[iCurrentArg].size()!=2){
                        cerr<<"For the hotspot feature you must enter the filename of the hotspot list."<<endl;
                        exit(1);
                    }
                    filename = arguments[iCurrentArg][1].data();
                    inFile.open(filename);
                    if (!inFile)
                        throw "HotSpot input file not found\n";
                    if (inFile.is_open()) {
                        pConfig->pHotSpotBinPtrList = new HotSpotBinPtrList;
                        string line;
                        int total=0;
                        while(getline(inFile,line)){
                            istringstream inStr(line);
                            double start,end,ratio;
                            inStr>>start>>end>>ratio;
                            if (pConfig->bDebug){
                                cerr<<"Hot spot from "<<start<<" to "<<end<<" with rate "<<ratio<<endl;
                            }
                            HotSpotBinPtr bin(new HotSpotBin(start,end,ratio));
                            pConfig->pHotSpotBinPtrList->push_back(bin);
                            ++total;
                        }
                        cerr<<"INPUT: Accepted "<<total<<" hotspots"<<endl;
                        inFile.close();
                        pConfig->bVariableRecomb = true;
                    }
                    break;
                case 'i' :
                    if (arguments[iCurrentArg].size()!=2) {
                        cerr<<"For flag "<<arguments[iCurrentArg][0][1]<<
                        ", you must enter a single int value for the number of iterations\n";
                        exit(1);
                    }
                    pConfig->iIterations = atoi(arguments[iCurrentArg][1].data());
                    cerr<<"INPUT: Iterations is now "<<pConfig->iIterations<<endl;
                    break;
                case 'I' :
                    if (arguments[iCurrentArg].size()<2) {
                        cerr<<"For flag "<<arguments[iCurrentArg][0][1]<<
                        ", the first parameter needs to the number of population islands\n";
                        exit(1);
                    }
                    pConfig->iTotalPops = atoi( arguments[iCurrentArg][1].data());
                    iNoMigrPops=2+pConfig->iTotalPops;
                    iMigrPops=3+pConfig->iTotalPops;

                    if (arguments[iCurrentArg].size()==iNoMigrPops||
                    arguments[iCurrentArg].size()==iMigrPops) {
                        pConfig->pPopList.clear();
                        for(unsigned int i=0; i<pConfig->iTotalPops; ++i) {
                            Population newPop;
                            const char * arg = arguments[iCurrentArg][2+i].data();
                            newPop.setChrSampled(atoi(arg));
                            iRunningSample+=newPop.getChrSampled();
                            cerr<<"INPUT: Setting chr sampled for pop "<<(i+1)<<" to "<<newPop.getChrSampled()<<endl;
                            newPop.setPopSize(dDefaultPopSize) ;
                            newPop.setGrowthAlpha(dDefaultGrowthAlpha);
                            newPop.setLastTime(0);
                            pConfig->pPopList.push_back(newPop);
                        }
                        if (arguments[iCurrentArg].size()==iMigrPops){
                            pConfig->dGlobalMigration = atof(arguments[iCurrentArg][iMigrPops-1].data());
                        }else{
                            pConfig->dGlobalMigration = dDefaultMigrationRate;
                        }
                        cerr<<"INPUT: Global migration rate to "<<pConfig->dGlobalMigration<<endl;

                    }else{
                        cerr<<"For flag "<<arguments[iCurrentArg][0][1]<<
                        ", the number of island sample sizes entered does not match the first parameter\n";
                        exit(1);
                    }
                    if (iRunningSample!=iSampleSize){
                        throw "The number of chromosomes entered in the -I option doesn't match the total sample size";
                    }

                    // Allocate migration rate matrix
                   if (pConfig->bDebug){
                     cerr<<"Constructing migration matrix of dimension "<<pConfig->iTotalPops<<endl;
                   }
                    pConfig->dMigrationMatrix.clear();
                    for (int i=0;i<pConfig->iTotalPops;++i){
                        vector<double> newRow;
                        for (int j=0;j<pConfig->iTotalPops;++j)
                            newRow.push_back(pConfig->dGlobalMigration/
                            (pConfig->iTotalPops-1));
                        pConfig->dMigrationMatrix.push_back(newRow);
                    }
                    break;
                case 'm' :
                    if( pConfig->iTotalPops < 2 ) {
                        cerr<<"You must use -I option first (i.e. specify more than one population)."<<endl;
                        exit(1);
                    }
                    if (arguments[iCurrentArg][0][2]=='a') {
                        iTotalCells = pConfig->iTotalPops * pConfig->iTotalPops + 1;
                        if (arguments[iCurrentArg].size()!=iTotalCells){
                            cerr<<"For flag "<<arguments[iCurrentArg][0]<<
                            ", the number of matrix cells does not match the total populations squared\n";
                            exit(1);
                        }
                        iSubOption = 0;
                        for(int pop1 = 0; pop1 <pConfig->iTotalPops; ++pop1){
                            for(int pop2 = 0; pop2 <pConfig->iTotalPops; ++pop2){
                                pConfig->dMigrationMatrix[pop1][pop2]=
                                atof( arguments[iCurrentArg][++iSubOption].data() ) ;
                            }
                        }
                        // set the diagonals as the sum of the off-diagonals
                        for(int pop1 = 0; pop1 < pConfig->iTotalPops; ++pop1) {
                            pConfig->dMigrationMatrix[pop1][pop1] = 0.0 ;
                            for(int pop2 = 0; pop2 < pConfig->iTotalPops; ++pop2){
                                if( pop1 != pop2 )
                                    pConfig->dMigrationMatrix[pop1][pop1] +=
                                    pConfig->dMigrationMatrix[pop1][pop2];
                            }
                        }
                    } else {
//                    // lets the user enter the entire migration by specified element
                        if (arguments[iCurrentArg].size()!=4){
                            cerr<<"For flag "<<arguments[iCurrentArg][0]<<
                            ", you need the source pop, dest pop, and the migration rate.\n";
                            exit(1);
                        }else{
                            int i = atoi( arguments[iCurrentArg][1].data() ) -1;
                            int j = atoi( arguments[iCurrentArg][2].data() ) -1;
                            double mij = atof( arguments[iCurrentArg][3].data() );
                            pConfig->dMigrationMatrix[i][i] += mij -
                            pConfig->dMigrationMatrix[i][j];
                            pConfig->dMigrationMatrix[i][j] = mij;
                        }
                    }
                    break;
                case 'n' :
//                    // specify population size for each population
                    if( pConfig->iTotalPops < 2 ) {
                        cerr<<"You must use -I option first (i.e. specify more than one population)."<<endl;
                        printUsage();
                    }
                    if (arguments[iCurrentArg].size()!=3){
                        cerr<<"For flag "<<arguments[iCurrentArg][0]<<
                        ", you need to specify the pop ID and the population size.\n";
                        exit(1);
                    }else{
                        popId = atoi( arguments[iCurrentArg][1].data() ) -1;
                        popSize = atof( arguments[iCurrentArg][2].data() );
                        if (popId>pConfig->iTotalPops||popId<0){
                            cerr<<"Invalid pop ID"<<endl;
                            exit(1);
                        }
                        pConfig->pPopList[popId].setPopSize(popSize) ;
                        cerr<<"INPUT: Pop "<<arguments[iCurrentArg][1]<<" has size: "<<popSize<<endl;
                    }
                    break;
                case 'g' :
//                    // specify growth rates
                    if( pConfig->iTotalPops < 2 ) {
                        cerr<<"You must use -I option first (i.e. specify more than one population)."<<endl;
                        printUsage();
                    }
                    if (arguments[iCurrentArg].size()!=3){
                        cerr<<"For flag "<<arguments[iCurrentArg][0]<<
                        ", you need to specify the pop ID and the population growth rate.\n";
                        exit(1);
                    }else{
                        popId = atoi( arguments[iCurrentArg][1].data() ) -1;
                        dDefaultGrowthAlpha = atof( arguments[iCurrentArg][2].data() );
                        if (popId>pConfig->iTotalPops||popId<0){
                            cerr<<"Invalid pop ID"<<endl;
                            exit(1);
                        }
                        pConfig->pPopList[popId].setGrowthAlpha(dDefaultGrowthAlpha);
                        cerr<<"INPUT: Pop "<<arguments[iCurrentArg][1].data()<<" has growth rate: "<<dDefaultGrowthAlpha<<endl;
                    }
                    break;
                case 'G' :
//                    // specify growth rates across all populations
                    if (arguments[iCurrentArg].size()!=2){
                        cerr<<"For flag "<<arguments[iCurrentArg][0]<<
                        ", you need to specify a single growth rate for all populations.\n";
                        exit(1);
                    }else{
                        float g = atof(arguments[iCurrentArg][1].data());
                        if (g<0) throw "Global growth rate must be positive";
                        dDefaultGrowthAlpha = atof( arguments[iCurrentArg][1].data() );
//                    cerr<<"INPUT: Growth rate for all pop "<<dDefaultGrowthAlpha<<endl;
                        for(int i=0; i<pConfig->iTotalPops; ++i){
                           pConfig->pPopList[i].setGrowthAlpha(dDefaultGrowthAlpha);
                           cerr<<"INPUT: Growth rate for pop "<<i<<" is "<<
                           pConfig->pPopList[i].getGrowthAlpha()<<endl;
                        }
                    }
                    break;
                case 'e' :
                    // these are events.  Be sure the times are unique
                    chType = arguments[iCurrentArg][0][2];
                    if ((arguments[iCurrentArg][0][3])=='a') bAcceptFullMigrMatrix = true;
                    else bAcceptFullMigrMatrix = false;
                    if (arguments[iCurrentArg].size()<2){
                        cerr<<"For event flags, you need to specify at least a time after "<<
                        arguments[iCurrentArg][0]<<endl;
                        exit(1);
                    }
                    dTime = atof(arguments[iCurrentArg][1].data());
                    cerr<<"INPUT: At time "<<dTime<<": ";
                    if (eventTimes.find(dTime)==eventTimes.end()){
                      eventTimes.insert(dTime);
                    }else{
                      cerr<<"Error, this event is redundant with a previous time.  Please increment it slightly from "<<dTime<<" to prevent unpredictable results\n";
                      throw "Invalid input";
                    }
                    int iPop1,iPop2;
                    double dProportion;
                    switch(chType){
                        case 'N': // global population size
                            if (arguments[iCurrentArg].size()!=3){
                                cerr<<"For flag "<<arguments[iCurrentArg][0]<<
                                ", you need to specify a single pop size for all populations.\n";
                                exit(1);
                            }else{
                                //int iType = Event::GLOBAL_POPSIZE;
                                wrapper = EventPtr(new GenericEvent(
                                Event::GLOBAL_POPSIZE,dTime,
                                atof(arguments[iCurrentArg][2].data())));
                                cerr<<"Global pop size is "<<
                                atof( arguments[iCurrentArg][2].data() )<<endl;
                            }
                            break;
                        case 'G': // global growth rate
                            if (arguments[iCurrentArg].size()!=3){
                                cerr<<"For flag "<<arguments[iCurrentArg][0]<<
                                ", you need to specify a single growth rate for all populations.\n";
                                exit(1);
                            }else{
                                float g = atof(arguments[iCurrentArg][2].data());
                                //if (g<0) throw "Global event growth rate must be positive";
                                //int iType = Event::GLOBAL_POPGROWTH;
                                wrapper = EventPtr(new GenericEvent(
                                Event::GLOBAL_POPGROWTH,dTime,
                                g));
                                cerr<<"Global growth rate is "<<
                                g<<endl;
                            }
                            break;
                        case 'M': // global migration rate
                            pConfig->bMigrationChangeEventDefined = true;
                            if (arguments[iCurrentArg].size()!=3){
                                cerr<<"For flag "<<arguments[iCurrentArg][0]<<
                                ", you need to specify a single migration rate for all populations.\n";
                                exit(1);
                            }else{
                                //iType = Event::GLOBAL_MIGRATIONRATE;
                                wrapper = EventPtr(new GenericEvent(
                                Event::GLOBAL_MIGRATIONRATE,dTime,
                                atof(arguments[iCurrentArg][2].data())));
                                cerr<<"Global migration rate is "<<
                                atof( arguments[iCurrentArg][2].data() )<<endl;
                            }
                            break;
                        case 'n' :  // subpopulation size
                            if (arguments[iCurrentArg].size()!=4){
                                cerr<<"For flag "<<arguments[iCurrentArg][0]<<
                                ", you need to specify pop id followed by the new size.\n";
                                exit(1);
                            }else{
                                //iType = Event::POPSIZE;
                                wrapper = EventPtr(new PopSizeChangeEvent(
                                Event::POPSIZE,dTime,atoi( arguments[iCurrentArg][2].data() ) -1,
                                atof( arguments[iCurrentArg][3].data() )));
                                cerr<<"For population "<<arguments[iCurrentArg][2]<<
                                ", pop size is now "<<atof( arguments[iCurrentArg][3].data() )<<endl;
                            }
                            break;
                        case 'g' :  // subpopulation growth
                            if (arguments[iCurrentArg].size()!=4){
                                cerr<<"For flag "<<arguments[iCurrentArg][0]<<
                                ", you need to specify pop id followed by the new growth rate.\n";
                                exit(1);
                            }else{
                                //iType = Event::GROWTH;
                                wrapper = EventPtr(new PopSizeChangeEvent(
                                Event::GROWTH,dTime,atoi( arguments[iCurrentArg][2].data() ) -1,
                                atof( arguments[iCurrentArg][3].data() )));
                                cerr<<"For population "<<arguments[iCurrentArg][2]<<
                                ", pop growth rate is now "<<atof( arguments[iCurrentArg][3].data() )<<endl;
                            }
                            break;
                        case 's' :  // split
                            pConfig->bMigrationChangeEventDefined = true;
                            if (arguments[iCurrentArg].size()!=4){
                                cerr<<"For flag "<<arguments[iCurrentArg][0]<<
                                ", you need to specify pop id followed by the proportion of the split.\n";
                                exit(1);
                            }else{
                                //iType = Event::POPSPLIT;
                                iPop1 = atoi( arguments[iCurrentArg][2].data() )-1;
                                dProportion = atof( arguments[iCurrentArg][3].data() );
                                if (iPop1<0||iPop1>=pConfig->iTotalPops||dProportion<0||
                                dProportion>=1){
                                    cerr<<"Bad values in parameters for pop IDs and/or proportion in pop split\n";
                                    printUsage();
                                }
                                wrapper = EventPtr(new PopSizeChangeEvent(
                                Event::POPSPLIT,dTime,iPop1,dProportion));
                                cerr<<"Population "<<arguments[iCurrentArg][2]<<
                                " splits at proportion "<<dProportion<<endl;
                            }
                            break;
                        case 'j':   // move lineages from pop1 to pop2

                            if (arguments[iCurrentArg].size()!=4){
                                cerr<<"For flag "<<arguments[iCurrentArg][0]<<
                                ", you need to specify source pop id followed by the destination pop id.\n";
                                exit(1);
                            }else{
                                //iType = Event::POPJOIN;
                                iPop1 = atoi( arguments[iCurrentArg][2].data() ) -1;
                                iPop2 = atoi( arguments[iCurrentArg][3].data() ) -1;
                                if (iPop1<0||iPop2<0){
                                  cerr<<"Bad values in parameters for pop IDs pop join\n";
                                  printUsage();
                                }
                                if (iPop1>=pConfig->iTotalPops||
                                iPop2>=pConfig->iTotalPops){
                                  cerr<<"WARNING: The pop IDs used in pop join is greater than the number specified in -I.  You must have a split event before this join event.\n";
                                }
                                wrapper = EventPtr(new PopJoinEvent(
                                Event::POPJOIN,dTime,iPop1,iPop2));
                                cerr<<"Population "<<
                                arguments[iCurrentArg][2].data()<<
                                " will merge with "<<
                                arguments[iCurrentArg][3].data()<<endl;
                            }
                            break;
                        case 'm':
                            pConfig->bMigrationChangeEventDefined = true;
                            if (bAcceptFullMigrMatrix){ // the -ema iTotalPops
//                            //<matrix element list>
                                if (arguments[iCurrentArg].size()<3){
                                    cerr<<"For flag "<<arguments[iCurrentArg][0]<<
                                    ", you need to at least specify the total number of populations.\n";
                                    exit(1);
                                }
                                int iTotalPops = atoi(arguments[iCurrentArg][2].data());
                                iTotalCells = iTotalPops * iTotalPops + 3;
                                if (arguments[iCurrentArg].size()!=iTotalCells){
                                    cerr<<"For flag "<<arguments[iCurrentArg][0]<<
                                    ", the number of cells do not match the number of pops specified squared.\n";
                                    exit(1);
                                }
                                //iType = Event::MIGRATION_MATRIX_RATE;
                                iSubOption = 2;
                                MatrixDouble dMigrationMatrix;
                                for (int i=0;i<iTotalPops;++i){
                                    vector<double> newRow;
                                    for(int j=0; j<iTotalPops; ++j) {
                                        if (i==j){
                                            newRow.push_back(0.0);
                                            ++iSubOption;
                                        }
                                        else{
                                            newRow.push_back(atof(arguments[iCurrentArg][++iSubOption].data()));
                                        }
                                    }
                                    dMigrationMatrix.push_back(newRow);
                                }
                                for(int i=0; i< iTotalPops; ++i) {
                                    for(int j=0; j<iTotalPops; ++j) {
                                        if (i!=j) dMigrationMatrix[i][i] +=
                                        dMigrationMatrix[i][j];
                                    }
                                }
                                wrapper = EventPtr(new
                                MigrationRateMatrixEvent(
                                Event::MIGRATION_MATRIX_RATE,dTime,
                                dMigrationMatrix));
                                cerr<<"Full migration matrix provided by the user\n";

                            }else{
                                 // the -em t i j x option specify just
                                 //part of the migration matrix
                                if (arguments[iCurrentArg].size()!=5){
                                    cerr<<"For flag "<<arguments[iCurrentArg][0]<<
                                    ", you must specify the source pop, dest pop, and the migration rate.\n";
                                    exit(1);
                                }else{
                                //iType = Event::MIGRATION_RATE;
                                cerr<<"Mig rate of source pop "<<arguments[iCurrentArg][2] <<" to dest pop "<<arguments[iCurrentArg][3]<<" set to "<<arguments[iCurrentArg][4]<<".\n";
                                wrapper = EventPtr(new MigrationRateEvent(
                                Event::MIGRATION_RATE,dTime,
                                atoi( arguments[iCurrentArg][2].data() ) -1,
                                atoi( arguments[iCurrentArg][3].data() ) -1,atof( arguments[iCurrentArg][4].data() ) ));
                                }
                            }
                            break;
                        default:
                            cerr<<"Invalid suboption, you entered"<<chType<<endl;
                            break;
                    }
                    pEventList->push_back(wrapper);
                    break;
                default:
                    cerr<<"Invalid option, you entered "<<arguments[iCurrentArg][0][1]<<endl;
                    printUsage();
            }
        }catch(const out_of_range & e){
            cerr<<"There were too many arguments.\n";
            printUsage();
        }
    }

    // Final sanity checks for the program before we begin:

    if (pConfig->iGeneConvTract>pConfig->dBasesToTrack){
        cerr<<"Warning: the gene conversion tract (-c 2nd parameter) cannot be "<<
        "longer than the length of sequence (-h parameter) to retain. ";
        pConfig->dBasesToTrack=2.0*pConfig->iGeneConvTract;
        cerr<<"The -h parameter is now revised to the recommend value of 2*tractlen = "
        <<pConfig->dBasesToTrack<<endl;
    }


    pEventList->sort(byEventTime());
    pConfig->pEventList = pEventList;
}
 -->
