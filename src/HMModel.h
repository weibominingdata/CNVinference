#pragma once

#include <map>
#include <vector>
#include <string>

struct ReadDepthDataItem
{
	std::string chr;
	int startPos;
	int endPos;
	int windowSize;
	double count;        // in original scale
	double hFuntionGC;   // in original scale
	double logMap;       // in log scale
	int state;
};

struct ReadDepthData
{
	char * fileName;
	std::vector<ReadDepthDataItem> data;
	double largestReadCount;
	double medianReadCount;
	double medianlogmap;
	double medianhgc;
	void loadData(char * f);
	std::string chr;
};

class HMModel
{
public:
	HMModel(void);
	~HMModel(void);
	HMModel(const HMModel & m);


	// read depth data
	void loadReadDepthData(char * filename);
	void setReadDepthVariable();

	// read depth variable
	ReadDepthData inferData;
//	void calculateMuAndPhiWithAutoRegression();
//	void calculateMuAndPhi(bool init=false); // fill in the mu and phi matrix

	void calculateMuAndPhiAllStatesCombined(bool init=false);
	void calculateMuAndPhiWithAutoRegressionAllStatesCombined();

	void startFromCoefficient();

	double **mu;
	double *phi;

	double mixtureProportion;     // proportion of mixture component
	double mixtureProportionNormal; // proportion of mixture compoent of state 2


public:
//	void loadData(char * sn, char * in, std::string chr);
	void writeResult(void);
	void inferAndEstimation(int rounds);

	bool USINGMAPPABILITY;
	bool USINGAUTOREGRESSION;
	bool USINGMIXTURECOMPONENT;
	bool REESTIMATETRANSITION;

private:
	void computAlpha(void);
	void computBeta(void);
	void computLikelihood(void);
	void computGamma(void);
	void reEstimation(bool transitionReestimate=true);//Baum-Welsh
	void findBestPath(bool viterbi=true);//Viterbi Algorithm
	void printVariable(void);
	void doOneRoundInference();


private:


private:
	int nSTATES;
	int nLength;
	int normalStates; // index of normal states
	double **pAlpha; // now store the log value
	double **pBeta;  // now store the log value
	double **pGamma; // now store the log posterial value
	double *pPi;   // the initial probability
	double **pTranTbl; 
	double **pEmissTbl;
	double ***pTran;
	std::vector<double> cLikelihood;   // log likelihood
	int * inferenceResults;
	double largestReadCount;      // largest read count
	double medianReadCount;    // median of read count
	double medianLogmap;
	double medianHgc;
	std::vector<double> median;


	int nITRATION;


//	void fillEmissTblItem(int site,int state, int pos);
	void fillEmissTblItem(int site,int state);
	void fillEmissionTbl(void);
	void fillTranContinous(void); 
	void fillTranDiscrete(void);
	double lamda;   // state duration parameter. assume now it is not state specific


public:

	


//  SNP arrary data, not use for this time   
	std::map<int, double> cLRR;
	std::map<int, double> cBAF;
	std::map<int, double> cPFB;
	std::vector<std::string> cName;
	std::vector<int> cPos;

	char * snpdataname;
	char * infodataname;
	//void setVariables(void);
	double ***pKexi;
	double **pGa;

private:
	std::string chrSymbol;
	void setFileName(char * sn, char * in);
//	double calLRR(double r, int z, double *, double *, double *);
//	double calBAF(double b, int z, double pb, double *, double **, double **);
private:
	double *pircn;
	double *murcn;
	double *sdrcn;
	double *pir;
	double *mur;
	double *sdr;
	double *pib;
	double **mub;
	double **sdb;

private:
	// for debug
	void writeKeyValue(int index);

};
