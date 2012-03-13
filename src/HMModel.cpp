#include "HMModel.h"

// starting of my header files
#include <iostream>
#include <fstream>
#include <cmath>
#include <assert.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <algorithm>

// math tools, NB regression
#include "MathTools.h"

using namespace std;
// ending of my header files


void ReadDepthData::loadData(char *f)
{
	int lastStartPos = 0;
	largestReadCount = 0;
	vector<double> readcounts;
	vector<double> logmap;
	vector<double> hgc;
	string str;
	ifstream fin(f);
	fin >> str >> str >> str >> str >> str >> str >> str 
		>> str >> str >> str >> str >> str >> str >> str;
	string temp;
	
	while(!fin.eof())
	{
		ReadDepthDataItem item;
		fin >> item.chr;
		if (fin.eof())
			break;
		fin >> item.startPos >> item.endPos >> item.windowSize;
		fin >> temp >> temp >> temp;
		fin >> item.count;
		fin >> temp >> temp >> temp;
		fin	>> item.hFuntionGC >> item.logMap >> item.state;
		if (item.count > largestReadCount)
			largestReadCount = item.count;
		readcounts.push_back(item.count);
		logmap.push_back(item.logMap);
		hgc.push_back(item.hFuntionGC);
		data.push_back(item);
	}

	chr = data[0].chr;
	
	cout << "Finish Reading Files" << endl;

	sort(readcounts.begin(),readcounts.end());
	medianReadCount = readcounts[int(readcounts.size()/2)];

	sort(logmap.begin(),logmap.end());
	medianlogmap = logmap[int(logmap.size()/2)];

	sort(hgc.begin(),hgc.end());
	medianhgc = hgc[int(hgc.size()/2)];
	//int lastStartPos = 0;
	//string str;
	//ifstream fin(f);
	//fin >> str >> str >> str >> str >> str >> str >> str;
	//string temp;
	//
	//while(!fin.eof())
	//{
	//	ReadDepthDataItem item;
	//	fin >> item.chr;
	//	if (fin.eof())
	//		break;
	//	fin >> item.startPos >> item.endPos >> item.state >> item.logMap >> item.hFuntionGC
	//		>> item.count;
	//	item.windowSize = 500;
	//	data.push_back(item);
	//}

	//chr = data[0].chr;
}



////////////////////////////////////////////////////
// an auxiliary data structure for data processing 
// store LRR and BAF value
struct LRRBAF
{
	LRRBAF() {LRR=BAF=-10.0;}
	double LRR;
	double BAF;
};


HMModel::HMModel(void)
: nSTATES(0)
, nLength(0)
, pTranTbl(NULL)
, pEmissTbl(NULL)
, pAlpha(NULL)
, pBeta(NULL)
//, nOBSERVE(0)
, pGamma(NULL)
, pKexi(NULL)
, pGa(NULL)
, pPi(NULL)
, pircn(NULL)
, murcn(NULL)
, sdrcn(NULL)
, pir(0)
, mur(NULL)
, sdr(NULL)
, pib(NULL)
, mub(NULL)
, sdb(NULL)
, chrSymbol("1")
, lamda(0.00001)
, nITRATION(0)
, largestReadCount(0)
, medianReadCount(0)
, medianLogmap(0)
, medianHgc(0)
, mixtureProportion(0.1)
, mixtureProportionNormal(0.1)
, normalStates(2)
, normalSelfTran(0.995)
, otherSelfTran(0.95)
, USINGMAPPABILITY(true)
, USINGAUTOREGRESSION(true)
, USINGMIXTURECOMPONENT(true)
, REESTIMATETRANSITION(true)
, REESTIMATEINIT(true)
{
}

void HMModel::loadReadDepthData(char * filename)
{
	setFileName(filename, filename);
	inferData.loadData(filename);
	chrSymbol = inferData.chr;
	largestReadCount = inferData.largestReadCount;
	medianReadCount = inferData.medianReadCount;
	setReadDepthVariable();
	// load median count
	//ifstream fin("TrioMedian_chr11");
	//double temp;
	//while(!fin.eof())
	//{
	//	fin >> temp;
	//	median.push_back(temp);
	//}
	//fin.close();
	
	//calculateMuAndPhi(true);
	//calculateMuAndPhiWithAutoRegression(); // auto regression part

	//calculateMuAndPhiAllStatesCombined(true);
	//calculateMuAndPhiWithAutoRegressionAllStatesCombined();

	startFromCoefficient();

	fillTranDiscrete();  // for read depth data, it is discrete time
	fillEmissionTbl();
}

void HMModel::startFromCoefficient()
{
	cout << "start with coefficient" << endl;
	cout << nSTATES << " " << normalStates << endl;
	// calculate the median of the win counts

	// this version accords with the all states combined model
	//double coverageDifference = medianReadCount/220;
	//double intercept = 4.7;
	//double newintercept = intercept+log(coverageDifference);
	double delta = 0.5;
	double coefficientforgc = 0.5;

	cout << delta << " " << coefficientforgc << " " << endl;

	//cout << coverageDifference << " " << newintercept << endl;
    double newintercept = log(medianReadCount)-log(normalStates)-medianLogmap-coefficientforgc*medianHgc;
    cout << newintercept << endl;

/*	double coverageDifference = medianReadCount/220;
	double intercept = 4.7;
	double newintercept = intercept+log(coverageDifference);
	double delta = 0.5;
	double coefficientforgc = 1.0;
	cout << coverageDifference << " " << newintercept << endl;
*/
	for(int i = 0; i < nSTATES; ++i)
	{
		phi[i] = 1.0;
	}

	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			double offset = 0;
			if (j==0)
				offset = log(delta)+inferData.data[i].logMap;
			else
				offset = log(j*1.0)+inferData.data[i].logMap;
			mu[i][j] = exp(newintercept+offset+
					coefficientforgc*inferData.data[i].hFuntionGC);
		}
	}

}

void HMModel::calculateMuAndPhiAllStatesCombined(bool init)
{
	// use glmNB to fill mu matrix and phi matrix 
	// for state 0, we only need to fit phi,
	// for other state, we need to get fitted mu and re-estimated phi

	// for state 0
	// load weights, load fitted value, in original scale, estimate phi
	int maxIt = 25;
	double convR = 1e-8;
	int nCovariate = 1;          // gcContent, //mappability
	double *y = new double[nSTATES*nLength];
	double *fitted = new double[nSTATES*nLength];
	double *x = new double[nSTATES*nLength*nCovariate];
	double *prior = new double[nSTATES*nLength];
	double *weights = new double[nSTATES*nLength]; // will be changed in computation
	double *resid = new double[nSTATES*nLength]; 
	double *Xb = new double[nSTATES*nLength*nCovariate];  // used in the program
	double *offset = new double[nSTATES*nLength];         // offset, now log(state) is a offset
	double delta = 0.5; // delta for cn0

	// loaded x, y, z, fitted, weights
	// y x are not dependant with initial value
	int index = 0;
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
			y[index++] = inferData.data[i].count;
	}
	index = 0;
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			x[index++] = inferData.data[i].hFuntionGC;
		}
	}
	//for(int i = 0; i < nLength; ++i)
	//{
	//	for(int j = 0; j < nSTATES; ++j)
	//	{
	//		x[index++] = inferData.data[i].logMap;
	//	}
	//}
	//for(int i = 0; i < nLength; ++i)
	//{
	//	x[index++] = log(delta);
	//	for(int j = 1; j < nSTATES; ++j)
	//	{
	//		x[index++] = log(j*1.0);
	//	}
	//}


	// load offset
	index = 0;
	for(int i = 0; i < nLength; ++i)
	{
		offset[index++] = log(delta);
		for(int j = 1; j < nSTATES; ++j)
		{
			offset[index++] = log(j*1.0);
		}
	}
	if (USINGMAPPABILITY)
	{
		index = 0;
		for(int i = 0; i < nLength; ++i)
		{
			offset[index++] += inferData.data[i].logMap/*+log(median[i]+0.01)*/;
			for(int j = 1; j < nSTATES; ++j)
			{
				offset[index++] += inferData.data[i].logMap/*+log(median[i]+0.01)*/;
			}
		}
	}


	// load weights
	index = 0;
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			prior[index++] = exp(pGamma[i][j]);
		}
	}

	// update weights
	// given mu, readcount, overdispersion, have likelihood
	// calculate proportion
	if (!init && USINGMIXTURECOMPONENT)
	{
		index = 0;
		double w = 1;
		double l = 0;
		double proportion = 0;
		for(int i = 0; i < nLength; ++i)
		{
			for(int j = 0; j < nSTATES; ++j)
			{
				
				l = exp(MathTools::loglik_NB(1, phi[j], &mu[i][j], &inferData.data[i].count, &w));
				if (j==normalStates)
				{
					proportion = 
					(l*(1-mixtureProportionNormal))
					/(mixtureProportionNormal/largestReadCount+(1-mixtureProportionNormal)*l);
				}
				else
				{
					proportion = 
					(l*(1-mixtureProportion))
					/(mixtureProportion/largestReadCount+(1-mixtureProportion)*l);
				}
				
				prior[index++] *= proportion;
			}
		}
	}

	int dim[5];
	dim[0] = nLength*nSTATES;
	dim[1] = nCovariate;  
	dim[2] = maxIt; //
	dim[3] = 0; // false
	dim[4] = 1; // false, no offset

	int nIter = 0;
	int linkR = LOG;
	convR = 1e-8; 
	int rank = 0;
	double PHI = 0;  // not use
	double scale = 1.0;
	int de_resid = 0;
	int family = 0;
	double twologlik = 0;
	double scoreTestP = 0;
	int trace = 0;
	double beta = 0;



	int conv = MathTools::glmNB(dim, &nIter,y,prior, &linkR, offset, x, &convR, &rank,
		Xb, fitted, resid, weights, &PHI, &scale, &de_resid, &family,
		&twologlik, &scoreTestP, &trace, &beta);

	cout << "beta for gc content is " << beta << " ";
	//cout << "beta for log map is " << beta << " ";
	//double sumIntercept = 0;
	//for (int i = 0; i < index; ++i)
	//{
	//	//cout << log(fitted[i])-offset[i]-beta*x[i] << endl;
	//	sumIntercept += log(fitted[i])-offset[i]-beta*x[i];
	//}
	//double meanIntercept = sumIntercept/index;
	//cout << "mean value of intercept " << meanIntercept << " ";
	//sumIntercept = 0;
	//for (int i = 0; i < index; ++i)
	//{
	//	double temp = log(fitted[i])-offset[i]-beta*x[i];
	//	sumIntercept += (temp-meanIntercept)*(temp-meanIntercept);
	//}
	//cout << "variance value of intercept " << sumIntercept/index << " ";



	// save fitted
	index = 0;
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			mu[i][j] = fitted[index++];
		}
	}

	// calculate the mean square root of errors for fitting
	//double error = 0;
	//for(int i = 0; i < nLength; ++i)
	//{
	//	double e = 0;
	//	for(int j = 0; j < nSTATES; ++j)
	//	{
	//		e += (mu[i][j]-inferData.data[i].count)*(mu[i][j]-inferData.data[i].count);
	//	}
	//	e /= nSTATES;
	//	error += e;
	//}
	//error /= nLength;
	//error = sqrt(error);
	//cout << "MSE of fitted data is " << error << " ";

	// calculate the weighted mean square root of errors for fitting
	//error = 0;
	//for(int i = 0; i < nLength; ++i)
	//{
	//	double e = 0;
	//	for(int j = 0; j < nSTATES; ++j)
	//	{
	//		e += exp(pGamma[i][j])*(mu[i][j]-inferData.data[i].count)*(mu[i][j]-inferData.data[i].count);
	//	}
	//	error += e;
	//}
	//error /= nLength;
	//error = sqrt(error);
	//cout << "weigthed MSE of fitted data is " << error << endl;


	// fix phi

	//delete []y; y=NULL;
	delete []x; x=NULL;
	//delete []prior; prior=NULL;
	//delete []fitted; fitted=NULL;
	delete []weights; weights=NULL;
	delete []resid; resid=NULL;
	delete []Xb; Xb=NULL;
	delete []offset;

	// run phi_ml several times to get phi
	//y = new double[nLength];
	//fitted = new double[nLength];
	//prior = new double[nLength];
	//for(int j = 0; j < nSTATES; ++j)
	//{
	//	for(int i = 0; i < nLength; ++i)
	//	{
	//		y[i] = inferData.data[i].count;
	//		fitted[i] = mu[i][j];
	//		prior[i] = exp(pGamma[i][j]);
	//	}
	//	int cvPhi = MathTools::phi_ml(y, fitted, nLength, prior, maxIt,
	//		convR, &phi[j], 0, 0);
	//}
		int cvPhi = MathTools::phi_ml(y, fitted, nLength*nSTATES, prior, maxIt,
			convR, &phi[0], 0, 0);
		if (phi[0]>20)
			phi[0]=20;
		for(int j = 1; j < nSTATES; ++j)
		{
			phi[j] = phi[0];
		}


	delete []y; y=NULL;
	delete []fitted; fitted=NULL;
	delete []prior; prior=NULL;

}

void HMModel::calculateMuAndPhiWithAutoRegressionAllStatesCombined()
{
	// use glmNB to fill mu matrix and phi matrix 
	// for state 0, we only need to fit phi,
	// for other state, we need to get fitted mu and re-estimated phi

	// for state 0
	// load weights, load fitted value, in original scale, estimate phi
	int maxIt = 25;
	double convR = 1e-8;
	int nCovariate = 2;          // gcContent, autoRegression, //logmap
	double *y = new double[nSTATES*nLength];
	double *fitted = new double[nSTATES*nLength];
	double *x = new double[nSTATES*nLength*nCovariate];
	double *prior = new double[nSTATES*nLength];
	double *weights = new double[nSTATES*nLength]; // will be changed in computation
	double *resid = new double[nSTATES*nLength]; 
	double *Xb = new double[nSTATES*nLength*nCovariate];  // used in the program
	double *offset = new double[nSTATES*nLength];
	double delta = 0.5; // delta for cn0

	// loaded x, y, z, fitted, weights
	// y x are not dependant with initial value
	int index = 0;
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
			y[index++] = inferData.data[i].count;
	}
	index = 0;
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			if (i==0)
				x[index++] = 0;
			else
			{
				if (inferData.data[i-1].count == 0)
				{
					x[index++] = log(inferData.data[i-1].count+0.1)-log(mu[i-1][j]);
				}
				else
				{
					x[index++] = log(inferData.data[i-1].count)-log(mu[i-1][j]);
				}
			}   
		}
	}
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			x[index++] = inferData.data[i].hFuntionGC;
		}
	}
	//for(int i = 0; i < nLength; ++i)
	//{
	//	for(int j = 0; j < nSTATES; ++j)
	//	{
	//		x[index++] = inferData.data[i].logMap;
	//	}
	//}
	//for(int i = 0; i < nLength; ++i)
	//{
	//	x[index++] = log(delta);
	//	for(int j = 1; j < nSTATES; ++j)
	//	{
	//		x[index++] = log(j*1.0);
	//	}
	//}


	// load offset
	index = 0;
	for(int i = 0; i < nLength; ++i)
	{
		offset[index++] = log(delta);
		for(int j = 1; j < nSTATES; ++j)
		{
			offset[index++] = log(j*1.0);
		}
	}
	if (USINGMAPPABILITY)
	{
		index = 0;
		for(int i = 0; i < nLength; ++i)
		{
			offset[index++] += inferData.data[i].logMap/*+log(median[i]+0.01)*/;
			for(int j = 1; j < nSTATES; ++j)
			{
				offset[index++] += inferData.data[i].logMap/*+log(median[i]+0.01)*/;
			}
		}
	}



	// load weights
	index = 0;
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			prior[index++] = exp(pGamma[i][j]);
		}
	}

	// update weights
	// given mu, readcount, overdispersion, have likelihood
	// calculate proportion
	if (USINGMIXTURECOMPONENT)
	{
		index = 0;
		double w = 1;
		double l = 0;
		double proportion = 0;
		for(int i = 0; i < nLength; ++i)
		{
			for(int j = 0; j < nSTATES; ++j)
			{
				
				l = exp(MathTools::loglik_NB(1, phi[j], &mu[i][j], &inferData.data[i].count, &w));
				if (j==normalStates)
				{
					proportion = 
					(l*(1-mixtureProportionNormal))
					/(mixtureProportionNormal/largestReadCount+(1-mixtureProportionNormal)*l);
				}
				else
				{
					proportion = 
					(l*(1-mixtureProportion))
					/(mixtureProportion/largestReadCount+(1-mixtureProportion)*l);
				}
				prior[index++] *= proportion;
			}
		}
	}

	int dim[5];
	dim[0] = nLength*nSTATES;
	dim[1] = nCovariate;  
	dim[2] = maxIt; //
	dim[3] = 0; // false
	dim[4] = 1; // false, no offset

	int nIter = 0;
	int linkR = LOG;
	convR = 1e-8; 
	int rank = 0;
	double PHI = 0;  // not use
	double scale = 1.0;
	int de_resid = 0;
	int family = 0;
	double twologlik = 0;
	double scoreTestP = 0;
	int trace = 0;
	double beta = 0;

	int conv = MathTools::glmNB(dim, &nIter,y,prior, &linkR, offset, x, &convR, &rank,
		Xb, fitted, resid, weights, &PHI, &scale, &de_resid, &family,
		&twologlik, &scoreTestP, &trace, &beta);

	cout << "beta for gc content is " << beta << " " << endl;
	//cout << "beta for mappability is " << beta << " " << endl;


	// save fitted
	index = 0;
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			mu[i][j] = fitted[index++];
		}
	}

	// calculate the mean square root of errors for fitting
	//double error = 0;
	//for(int i = 0; i < nLength; ++i)
	//{
	//	double e = 0;
	//	for(int j = 0; j < nSTATES; ++j)
	//	{
	//		e += (mu[i][j]-inferData.data[i].count)*(mu[i][j]-inferData.data[i].count);
	//	}
	//	e /= nSTATES;
	//	error += e;
	//}
	//error /= nLength;
	//error = sqrt(error);
	//cout << "MSE of fitted data is " << error << " ";

	// calculate the weighted mean square root of errors for fitting
	//error = 0;
	//for(int i = 0; i < nLength; ++i)
	//{
	//	double e = 0;
	//	for(int j = 0; j < nSTATES; ++j)
	//	{
	//		e += exp(pGamma[i][j])*(mu[i][j]-inferData.data[i].count)*(mu[i][j]-inferData.data[i].count);
	//	}
	//	error += e;
	//}
	//error /= nLength;
	//error = sqrt(error);
	//cout << "weigthed MSE of fitted data is " << error << endl;


	//delete []y; y=NULL;
	delete []x; x=NULL;
	//delete []prior; prior=NULL;
	//delete []fitted; fitted=NULL;
	delete []weights; weights=NULL;
	delete []resid; resid=NULL;
	delete []Xb; Xb=NULL;
	delete []offset; offset = NULL;

	// run phi_ml several times to get phi
	//y = new double[nLength];
	//fitted = new double[nLength];
	//prior = new double[nLength];
	//for(int j = 0; j < nSTATES; ++j)
	//{
	//	for(int i = 0; i < nLength; ++i)
	//	{
	//		y[i] = inferData.data[i].count;
	//		fitted[i] = mu[i][j];
	//		prior[i] = exp(pGamma[i][j]);
	//	}
	//	int cvPhi = MathTools::phi_ml(y, fitted, nLength, prior, maxIt,
	//		convR, &phi[j], 0, 0);
	//}
		int cvPhi = MathTools::phi_ml(y, fitted, nLength*nSTATES, prior, maxIt,
			convR, &phi[0], 0, 0);
		if (phi[0]>20)
			phi[0]=20;
		for(int j = 1; j < nSTATES; ++j)
		{
			phi[j] = phi[0];
		}

	delete []y; y=NULL;
	delete []fitted; fitted=NULL;
	delete []prior; prior=NULL;
}

void HMModel::setTranInitValue(double **pTran)
{
	// rules
	// self-transition > transition to other states
	// self-transition normal > self-transition other states
	// transition to normal > transition to other
	// transition to similar state > transition to other states
	// added up to 1
	// detail design:
	// state normal transit to del = state normal transit to dup
	// state del transit to normal = 2 state transit to del = 4 state transit to dup
	// state dup transit to normal = 2 state transit to dup = 4 state transit to del
	int nDel = normalStates-0+1-1;
	cout << "#del states = " << nDel << endl;
	int nDup = nSTATES-normalStates-1;
	cout << "#dup states = " << nDup << endl;
	double toNormal=4, toSame=2, toDiff=1;
	double total= toNormal+toSame+toDiff;
	for(int i = 0; i < nSTATES; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			if (i==j)
			{
				if (i==normalStates)
					pTran[i][j] = normalSelfTran;
				else
					pTran[i][j] = otherSelfTran;
			}
			else
			{
				if (i==normalStates)
				{
					double p = 1-normalSelfTran;
					if (j<normalStates)
					{
						p = p/2;
						pTran[i][j] = p/nDel;
					}
					else
					{
						p = p/2;
						pTran[i][j] = p/nDup;
					}
				}
				else if (i<normalStates)
				{
					double p = 1-otherSelfTran;
					if (j < normalStates)
					{
						p = p*(toSame)/total;
						pTran[i][j]=p/(nDel-1);
					}
					else if (j == normalStates)
					{
						p = p*(toNormal)/total;
						pTran[i][j] = p;
					}
					else
					{
						p = p*(toDiff)/total;
						pTran[i][j] = p/(nDup);
					}
				}
				else
				{
					double p = 1-otherSelfTran;
					if (j < normalStates)
					{
						p = p*(toDiff)/total;
						pTran[i][j]=p/(nDel);
					}
					else if (j == normalStates)
					{
						p = p*(toNormal)/total;
						pTran[i][j] = p;
					}
					else
					{
						p = p*(toSame)/total;
						pTran[i][j] = p/(nDup-1);
					}
				}

			}
			/*			if (i==j)
			{
				if (i==normalStates)
					pTran[i][j] = normalSelfTran;
				else
					pTran[i][j] = otherSelfTran;
			}
			else
			{
				if (j==normalStates)
					pTran[i][j] = (1-otherSelfTran)*3/(nSTATES+1);
				else
				{
					if (i==normalStates)
						pTran[i][j] = (1-normalSelfTran)/(nSTATES-1);
					else
						pTran[i][j] = (1-otherSelfTran)/(nSTATES+1);
				}
			}
*/		}
	}
}



void HMModel::setReadDepthVariable()
{
	// hard coding for the real human data
	nSTATES = 7;
	normalStates = 2;
	// mouse data modification
	//nSTATES = 4;
	//normalStates = 1;

	nLength = inferData.data.size();

	mu = new double*[nLength];
	for(int i = 0; i < nLength; ++i)
	{
		mu[i] = new double[nSTATES];
	}
	// at the beginning, mu[:][i] = readcout 
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			mu[i][j] = inferData.data[i].count;	
		}
	}

	phi = new double[nSTATES];

	inferenceResults = new int[nLength];

	pTranTbl = new double*[nSTATES];
	for(int i = 0; i < nSTATES; ++i)
	{
		pTranTbl[i] = new double[nSTATES];
	}

	pTran = new double **[nLength];
	for(int i = 0; i < nLength; ++i)
	{
		pTran[i] = new double *[nSTATES];
		for(int j = 0; j < nSTATES; ++j)
		{
			pTran[i][j] = new double[nSTATES];
		}
	}
	//for(int i = 0; i < nSTATES; ++i)
	//{
	//	for(int j = 0; j < nSTATES; ++j)
	//	{
	//		if (i==j)
	//			pTranTbl[i][j] = 0.9;
	//		else
	//			pTranTbl[i][j] = 0.1/(nSTATES-1);
	//	}
	//}
	//double normalSelfTran=0.99995;
	//double normalSelfTran = 0.995;
	//double otherSelfTran=0.95;
	cout << normalSelfTran << " " << otherSelfTran << endl;
	setTranInitValue(pTranTbl);
/*	for(int i = 0; i < nSTATES; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			if (i==j)
			{
				if (i==normalStates)
					pTranTbl[i][j] = normalSelfTran;
				else
					pTranTbl[i][j] = otherSelfTran;
			}
			else
			{
				if (j==normalStates)
					pTranTbl[i][j] = (1-otherSelfTran)*3/(nSTATES+1);
				else
				{
					if (i==normalStates)
						pTranTbl[i][j] = (1-normalSelfTran)/(nSTATES-1);
					else
						pTranTbl[i][j] = (1-otherSelfTran)/(nSTATES+1);
				}
			}
		}
	}
*/



	// revise pEmissTbl nSTATES * nLength
	pEmissTbl = new double*[nLength];
	for(int i = 0; i < nLength; ++i)
	{
		pEmissTbl[i] = new double[nSTATES];
	}

	// hard coding for real data
	pPi = new double[nSTATES];
	for(int i = 0; i < nSTATES; ++i)
	{
		pPi[i] = (1-normalSelfTran)/(nSTATES-1);
	}
	pPi[normalStates] = normalSelfTran;

	pAlpha = new double*[nLength];
	for(int i = 0; i < nLength; ++i)
	{
		pAlpha[i] = new double[nSTATES];
		for(int j = 0; j < nSTATES; ++j)
		{
			pAlpha[i][j] = 0;
		}
	}

	pBeta = new double*[nLength];
	for(int i = 0; i < nLength; ++i)
	{
		pBeta[i] = new double[nSTATES];
		for(int j = 0; j < nSTATES; ++j)
		{
			pBeta[i][j] = 0;
		}
	}

	pGamma = new double*[nLength];   // store the posterior probability
	for(int i = 0; i < nLength; ++i)
		pGamma[i] = new double[nSTATES];
	for(int i=0; i < nLength; ++i)
	{
		for(int j=0; j < nSTATES; ++j)
		{
			if (j==inferData.data[i].state)
			{
				pGamma[i][j] = log(0.9);
			}
			else
			{
				pGamma[i][j] = log(0.1/(nSTATES-1));
			}
		}
	}


	//pKexi = new double **[nLength];
	//for(int i = 0; i < nLength; ++i)
	//{
	//	pKexi[i] = new double*[nSTATES];
	//	for(int j = 0; j < nSTATES; ++j)
	//	{
	//		pKexi[i][j] = new double[nSTATES];
	//	}
	//}

	//pGa = new double *[nLength];
	//for(int i = 0; i < nLength; ++i)
	//	pGa[i] = new double[nSTATES];


	// lots of new variables
	//pircn = new double [nSTATES];
	//murcn = new double [nSTATES];
	//sdrcn = new double [nSTATES];
	//pir = new double [nSTATES];
	//mur = new double [nSTATES];
	//sdr = new double [nSTATES];
	//pib = new double [nSTATES];
	//mub = new double *[nSTATES];
	//for(int i = 0; i < nSTATES; ++i)
	//	mub[i] = new double [nSTATES];
	//sdb = new double*[nSTATES];
	//for(int i = 0; i < nSTATES; ++i)
	//	sdb[i] = new double [nSTATES];

	//pir[0] = 0.01; pir[1] = 0.01; pir[2] = 0.1; pir[3] = 0.01; pir[4] = 0.01; pir[5] = 0.01;
	//mur[0] = 0.0; mur[1] = 0.0; mur[2] = -3.5; mur[3] = -0.67; mur[4] = 0.4; mur[5] = 0.68;
	//sdr[0] = 0.15; sdr[1] = 0.15; sdr[2] = 1.3; sdr[3] = 0.28; sdr[4] = 0.2; sdr[5] = 0.19;
	//pircn[0] = 0.01; pircn[1] = 0.01; pircn[2] = 0.1; pircn[3] = 0.01; pircn[4] = 0.01; pircn[5] = 0.01;
	//murcn[0] = 0.0; murcn[1] = 0.0; murcn[2] = -2.1; murcn[3] = -0.58; murcn[4] = 0.37; murcn[5] = 0.65;
	//sdrcn[0] = 0.18; sdrcn[1] = 0.18; sdrcn[2] = 2.1; sdrcn[3] = 0.35; sdrcn[4] = 0.24; sdrcn[5] = 0.3;
	//pib[0] = 0.01; pib[1] = 0.01; pib[2] = 0.5; pib[3] = 0.01; pib[4] = 0.01; pib[5] = 0.01;

	//mub[0][0] = 0.0;mub[0][1] = 0.5;mub[0][2] = 1.0;mub[0][3] = -1.0;mub[0][4] = -1.0;mub[0][5] = -1.0;
	//mub[1][0] = 0.0;mub[1][1] = -1.0;mub[1][2] = -1.0;mub[1][3] = 1.0;mub[1][4] = -1.0;mub[1][5] = -1.0;
	//mub[2][0] = 0.5;mub[2][1] = -1.0;mub[2][2] = -1.0;mub[2][3] = -1.0;mub[2][4] = -1.0;mub[2][5] = -1.0;
	//mub[3][0] = 0.0;mub[3][1] = -1.0;mub[3][2] = -1.0;mub[3][3] = 1.0;mub[3][4] = -1.0;mub[3][5] = -1.0;
	//mub[4][0] = 0.0;mub[4][1] = 0.333;mub[4][2] = 0.667;mub[4][3] = 1.0;mub[4][4] = -1.0;mub[4][5] = -1.0;
	//mub[5][0] = 0.0;mub[5][1] = 0.25;mub[5][2] = 0.5;mub[5][3] = 0.75;mub[5][4] = 1.0;mub[5][5] = -1.0;

	//sdb[0][0] = 0.016; sdb[0][1] = 0.035; sdb[0][2] = 0.016; sdb[0][3] = -1.0; sdb[0][4] = -1.0; sdb[0][5] = -1.0;
	//sdb[1][0] = 0.016; sdb[1][1] = -1.0; sdb[1][2] = -1.0; sdb[1][3] = 0.16; sdb[1][4] = -1.0; sdb[1][5] = -1.0;
	//sdb[2][0] = 0.15; sdb[2][1] = -1.0; sdb[2][2] = -1.0; sdb[2][3] = -1.0; sdb[2][4] = -1.0; sdb[2][5] = -1.0;
	//sdb[3][0] = 0.016; sdb[3][1] = -1.0; sdb[3][2] = -1.0; sdb[3][3] = 0.16; sdb[3][4] = -1.0; sdb[3][5] = -1.0;
	//sdb[4][0] = 0.016; sdb[4][1] = 0.042; sdb[4][2] = 0.042; sdb[4][3] = 0.016; sdb[4][4] = -1.0; sdb[4][5] = -1.0;
	//sdb[5][0] = 0.016; sdb[5][1] = 0.042; sdb[5][2] = 0.035; sdb[5][3] = 0.42; sdb[5][4] = 0.016; sdb[5][5] = -1.0;

}



HMModel::~HMModel(void)
{
	//for(int i = 0; i < nSTATES; ++i)
	//{
	//	delete []pTranTbl[i];
	//	pTranTbl[i] = NULL;
	//}
	//delete []pTranTbl;
	//pTranTbl = NULL;

	//delete []inferenceResults;
	//inferenceResults = NULL;

	//for(int i = 0; i < nLength; ++i)
	//{
	//	for(int j = 0; j < nSTATES; ++j)
	//	{
	//		delete []pTran[i][j];
	//		pTran[i][j] = NULL;
	//	}
	//	delete []pTran[i];
	//	pTran[i] = NULL;
	//}
	//delete []pTran;
	//pTran = NULL;
	//

	//for(int i = 0; i < nLength; ++i)
	//{
	//	delete []pEmissTbl[i];
	//	pEmissTbl[i] = NULL;
	//}
	//delete []pEmissTbl;
	//pEmissTbl = NULL;

	//for(int i = 0; i < nLength; ++i)
	//{
	//	delete []pAlpha[i];
	//	pAlpha[i] = NULL;
	//}
	//delete []pAlpha;
	//pAlpha = NULL;

	//for(int i = 0; i < nLength; ++i)
	//{
	//	delete []pBeta[i];
	//	pBeta[i] = NULL;
	//}
	//delete []pBeta;
	//pBeta = NULL;

	//for(int i = 0; i < nLength; ++i)
	//{
	//	delete []pGa[i];
	//	pGa[i] = NULL;
	//}
	//delete []pGa;
	//pGa = NULL;

	//for(int i = 0; i < nLength; ++i)
	//{
	//	for(int j = 0; j < nSTATES; ++j)
	//	{
	//		delete []pKexi[i][j];
	//		pKexi[i][j] = NULL;
	//	}
	//	delete []pKexi[i];
	//	pKexi[i] = NULL;
	//}
	//delete []pKexi;
	//pKexi = NULL;


	//for(int i = 0; i < nLength; ++i)
	//{
	//	delete []pGamma[i];
	//	pGamma[i] = NULL;
	//}
	//delete []pGamma;
	//pGamma = NULL;

	//delete []pPi;
	//pPi = NULL;

	//// lots of deconstruction 
	//delete	[]pircn;
	//pircn = NULL;
	//delete	[]murcn;
	//murcn = NULL;
	//delete	[]sdrcn;
	//sdrcn = NULL;
	//delete	[]pir;
	//pir = NULL;
	//delete	[]mur;
	//mur = NULL;
	//delete	[]sdr;
	//sdr = NULL;
	//delete	[]pib;
	//pib = NULL;
	//for(int i = 0; i < nSTATES; ++i)
	//{
	//	delete	[]mub[i];
	//	mub[i] = NULL;
	//}
	//delete	[]mub;
	//mub = NULL;
	//for(int i = 0; i < nSTATES; ++i)
	//{
	//	delete	[]sdb[i];
	//	sdb[i] = NULL;
	//}
	//delete	[]sdb;
	//sdb = NULL;
}

HMModel::HMModel(const HMModel & m)
{
}




void HMModel::inferAndEstimation(int rounds)
{
	writeKeyValue(0);
	for(int i = 0; i < rounds; ++i)
	{
		doOneRoundInference();
		reEstimation(REESTIMATETRANSITION, REESTIMATEINIT);
		writeKeyValue(i+1);
	}
	findBestPath(false);
	printVariable();
}

void HMModel::doOneRoundInference()
{
	nITRATION++;
	computAlpha();
	computBeta();
	computLikelihood();
	computGamma();
}

void HMModel::computAlpha(void)
{

	for(int i = 0; i < nSTATES; ++i)
	{
		pAlpha[0][i] = log(pPi[i]) + log(pEmissTbl[0][i]);
	}

	double *v = new double[nSTATES];
	for(int i = 1; i < nLength; ++i)
	{
		//if (i%100 == 0) cout << i << endl;
		for(int j = 0; j < nSTATES; ++j)
		{
			for(int k = 0; k < nSTATES; ++k)
			{
				v[k] = pAlpha[i-1][k]+log(pTran[i][k][j]);

			}
			pAlpha[i][j] = MathTools::logsumexp(v, nSTATES)+log(pEmissTbl[i][j]);
		}
	}
	delete []v;
}

void HMModel::computBeta(void)
{
	for(int i = 0;  i < nSTATES; ++i)
	{
		pBeta[nLength-1][i] = 0;
	}
	double *v = new double[nSTATES];
	for(int i = nLength-2; i >=0; --i)
	{
		//if (i%100 == 0) cout << i << endl;
		for(int j = 0; j < nSTATES; ++j)
		{
			for(int k = 0; k < nSTATES; ++k)
			{
				v[k] = pBeta[i+1][k] + log(pTran[i+1][j][k]) + log(pEmissTbl[i+1][k]);
			}
			pBeta[i][j] = MathTools::logsumexp(v, nSTATES);
		}
	}
	delete []v;
}

void HMModel::computGamma(void)
{
	// make sure comput Gamma is called after the corresponding computAlpha and computBeta
	for(int step = 0; step < nLength; ++step)
	{
		for(int i = 0; i < nSTATES; ++i)
		{
			pGamma[step][i] = pAlpha[step][i]+pBeta[step][i]-cLikelihood[nITRATION-1];
		}
	}
}

void HMModel::computLikelihood()
{
	double likelihood = 0;
	double *v = new double[nSTATES];
	for(int k = 0; k < nSTATES; ++k)
	{
		v[k] = pAlpha[nLength-1][k];
	}
	likelihood = MathTools::logsumexp(v, nSTATES);
	cLikelihood.push_back(likelihood);
	delete []v;
}



void HMModel::writeResult(void)
{
	
	//computGamma();  // now we have posterior probability
	//int * maxPostPIndx = new int[nLength];
	int * cn = new int[nLength];
	//double max = -10;
	for(int i = 0; i < nLength; ++i)
	{
		int index = inferenceResults[i];
		//max = -10;
		//for(int j = 0; j < nSTATES; ++j)
		//{
			//if (pGamma[i][j] > pGamma[i][index])
			//if (pAlpha[i][j]*pBeta[i][j] > max)
			//{
			//	max = pAlpha[i][j]*pBeta[i][j];
			//	index = j;
			//}
			//maxPostPIndx[i] = index;
			//if (index == 0 || index ==1)
			//{
			//	cn[i] = 2;
			//}
			//else if (index == 2)
			//{
			//	cn[i] = 0;
			//}
			//else if (index == 3)
			//{
			//	cn[i] = 1;
			//}
			//else if (index == 4)
			//{
			//	cn[i] = 3;
			//}
			//else if (index == 5)
			//{
			//	cn[i] = 4;
			//}
		//}
		// assign cn by the closet read count, not the HMM inferred states
		// purpose: how well is the fitting
/*		index=-1;
		double diff = 100000.0;

		for(int j = 0; j < nSTATES; ++j)
		{
			if (fabs(mu[i][j]-inferData.data[i].count) < diff)
			{
				diff = fabs(mu[i][j]-inferData.data[i].count);
				index = j;
		    }
		}
*/
		cn[i] = index;
	}

	string fileName("chr");
	fileName += chrSymbol;
	int sPos =(string(snpdataname).find_last_of("/")==string::npos)?0:string(snpdataname).find_last_of("/")+1;
	fileName += string(snpdataname).substr(sPos,string(snpdataname).length());
	sPos =(string(infodataname).find_last_of("/")==string::npos)?0:string(infodataname).find_last_of("/")+1;
	fileName += string(infodataname).substr(sPos,string(infodataname).length());
	string snpName = "Jingerbread_"+fileName+"_SNP.dat";
	string segName = "Jingerbread_"+fileName+"_segment.dat";
	ofstream out(snpName.c_str());
	//ofstream out("JSNP.dat");
	ofstream out1(segName.c_str());
	//out << "name\t" << "state\t" << "stateP\t" << "CN\t" << endl;
	//out1 << "chr\t" << "start\t" << "end\t" << "state\t" << "cn\t" << "sample\t" << "snp1\t" << "snp2\t" << "score\t" << "n\t" << endl;
	out << "state\t" << "stateP\t" << "CN\t" << endl;
	out1 << "chr\t" << "start\t" << "end\t" << "state\t" << "cn\t" << "sample\t" << "score\t" << "n\t" << endl;
	
	int start = 0;
	int end = 0;
	bool cnv = false;
	int cnvtype = -1;
	for(int i = 0; i < nLength; ++i)
	{
		if (cn[i] != normalStates)
		{
			if (!cnv)
			{
				cnv = true;
				start = i;
				cnvtype = cn[i];
			}
			else
			{
				if (cnvtype != cn[i])
				{
					end = i-1;
					double score = 0;
					for(int j = start; j <= end; ++j)
						score += exp(pGamma[j][inferenceResults[j]]);
					//out1 << chrSymbol << "\t" << cPos[start] << "\t" << cPos[end] << "\t"
					//	<< inferenceResults[start]+1 << "\t" << cn[start] << "\t" << "test\t"
					//	<< cName[start] << "\t" << cName[end] << "\t" << score << "\t"
					//	<< end-start+1 << endl;
					out1 << chrSymbol << "\t" << inferData.data[start].startPos << "\t" << inferData.data[end].endPos << "\t"
						<< inferenceResults[start] << "\t" << cn[start] << "\t" << "test\t"
						<< score << "\t" << end-start+1 << endl;
					// a start of a new CNV
					start = i;
					cnvtype = cn[i];
				}
			}
		}
		else
		{
			if (cnv)
			{
				end = i-1;
				// out a cnv segment
				double score = 0;
				for(int j = start; j <= end; ++j)
					score += exp(pGamma[j][inferenceResults[j]]);
				//out1 << chrSymbol << "\t" << cPos[start] << "\t" << cPos[end] << "\t"
				//	<< inferenceResults[start]+1 << "\t" << cn[start] << "\t" << "test\t"
				//	<< cName[start] << "\t" << cName[end] << "\t" << score << "\t"
				//	<< end-start+1 << endl;
				out1 << chrSymbol << "\t" << inferData.data[start].startPos << "\t" << inferData.data[end].endPos << "\t"
					<< inferenceResults[start] << "\t" << cn[start] << "\t" << "test\t"
					<< score << "\t" << end-start+1 << endl;				
			}
			cnv = false;
		}
		//out << cName[i] << "\t" << inferenceResults[i] << "\t" << exp(pGamma[i][inferenceResults[i]]) <<"\t"
		//	<< cn[i] << endl;
		out << inferenceResults[i] << "\t" << exp(pGamma[i][inferenceResults[i]]) <<"\t"
			<< cn[i] << endl;
	}

	if (cnv)
	{
		end = nLength-1;
		// out a cnv segment
		double score = 0;
		for(int j = start; j <= end; ++j)
			score += exp(pGamma[j][inferenceResults[j]]);
		//out1 << chrSymbol << "\t" << cPos[start] << "\t" << cPos[end] << "\t"
		//	<< inferenceResults[start]+1 << "\t" << cn[start] << "\t" << "test\t"
		//	<< cName[start] << "\t" << cName[end] << "\t" << score << "\t"
		//	<< end-start+1 << endl;
		out1 << chrSymbol << "\t" << inferData.data[start].startPos << "\t" << inferData.data[end].endPos << "\t"
			<< inferenceResults[start] << "\t" << cn[start] << "\t" << "test\t"
			<< score << "\t" << end-start+1 << endl;				
	}


	out.close();
	out1.close();
	//delete []maxPostPIndx;
	//maxPostPIndx = NULL;
	delete []cn;
	cn = NULL;

}



void HMModel::reEstimation(bool transitionReestimate, bool initReestimation)
{
    if (initReestimation)
    {
		cout << "initial probability re-estimated" << endl;
		// update initial probability
		for(int i = 0; i < nSTATES; ++i)
		{
			pPi[i] = exp(pAlpha[0][i] + pBeta[0][i] - cLikelihood[nITRATION-1]);
		}
    }

	if (transitionReestimate)
	{
		cout << "transition re-estimated" << endl;
		// update transition probability
		// first we need to create a temp transition matrix to store new values
		double **newTran = new double*[nSTATES];
		for(int i = 0; i < nSTATES; ++i)
		{
			newTran[i] = new double[nSTATES];
			for(int j = 0; j < nSTATES; ++j)
				newTran[i][j] = 0;
		}
		// then we create a temp cjk function
		double **c = new double *[nSTATES];
		for(int i = 0; i < nSTATES; ++i)
		{
			c[i] = new double[nSTATES];
			for(int j = 0; j < nSTATES; ++j)
				c[i][j] = 0;
		}
		double *v = new double[nLength-1];
		for(int j = 0; j < nSTATES; ++j)
		{
			for(int k = 0; k < nSTATES; ++k)
			{
				for(int i = 1; i < nLength; ++i)
				{
					v[i-1] = pAlpha[i-1][j]+log(pEmissTbl[i][k])+pBeta[i][k];
				}
				c[j][k] = MathTools::logsumexp(v, nLength-1);
			}
		}
		delete []v;
		v = new double[nSTATES];
		for(int j = 0; j < nSTATES; ++j)
		{
			for(int k = 0; k < nSTATES; ++k)
			{
				newTran[j][k] = log(pTranTbl[j][k]);
				newTran[j][k] += c[j][k];
				for(int l = 0; l < nSTATES; ++l)
				{
					v[l] = log(pTranTbl[j][l])+c[j][l];
				}
				double vsum = MathTools::logsumexp(v, nSTATES);
				newTran[j][k] -= vsum;
				newTran[j][k] = exp(newTran[j][k]);
			}
		}
		delete []v;
		// put the new value to Tran Table
		for(int i = 0; i < nSTATES; ++i)
		{
			for(int j = 0; j < nSTATES; ++j)
				pTranTbl[i][j] = newTran[i][j];
		}
		fillTranDiscrete();
		for(int i = 0; i < nSTATES; ++i)
		{
			delete []newTran[i];
			newTran[i] = NULL;
		}
		delete []newTran;
		newTran = NULL;
		for(int i = 0; i < nSTATES; ++i)
		{
			delete []c[i];
			c[i] = NULL;
		}
		delete []c;
		c = NULL;
	}
	//calculateMuAndPhi();
	//calculateMuAndPhiWithAutoRegression(); // calculate with auto regression

	calculateMuAndPhiAllStatesCombined();
	if (USINGAUTOREGRESSION)
		calculateMuAndPhiWithAutoRegressionAllStatesCombined();


	fillEmissionTbl();
	//for(int step = 0; step < nLength-1; ++step)
	//{
	//	double sum = 0;
	//	for(int i = 0; i < nSTATES; ++i)
	//	{
	//		for(int j = 0; j < nSTATES; ++j)
	//		{
	//			sum += pAlpha[step][i]*pTran[step+1][i][j]*pEmissTbl[step+1][j]*pBeta[step+1][j];
	//		}
	//	}
	//	for(int i = 0; i < nSTATES; ++i)
	//	{
	//		for(int j = 0; j < nSTATES; ++j)
	//		{
	//			pKexi[step][i][j] = pAlpha[step][i]*pTran[step+1][i][j]*pEmissTbl[step+1][j]*pBeta[step+1][j]/sum;
	//		}
	//	}
	//}
	//for(int step = 0; step < nLength-1; ++step)
	//{
	//	for(int i = 0; i < nSTATES; ++i)
	//	{
	//		pGa[step][i] = 0;
	//		for(int j = 0; j < nSTATES; ++j)
	//		{
	//			pGa[step][i] += pKexi[step][i][j];
	//		}
	//	}
	//}
	//double sum = 0;
	//for(int i = 0; i < nSTATES; ++i)
	//{
	//	pGa[nLength-1][i] = pAlpha[nLength-1][i]*pBeta[nLength-1][i];
	//	sum += pGa[nLength-1][i];
	//}
	//for(int i = 0; i < nSTATES; ++i)
	//{
	//	pGa[nLength-1][i] /= sum;
	//}
	//// re-estimate init probability
	//for(int i = 0; i < nSTATES; ++i)
	//{
	//	pPi[i] = pGa[0][i];
	//}
	//// re-estimate transition probability
	//for(int i = 0; i < nSTATES; ++i)
	//{
	//	double sumi = 0;
	//	for(int step = 0; step < nLength-1; ++step)
	//	{
	//		sumi += pGa[step][i];
	//	}
	//	double sumii = 0;
	//	double *sumij = new double[nSTATES];
	//	for(int j = 0; j < nSTATES; ++j)
	//	{
	//		sumij[j] = 0;
	//		for(int step = 0; step < nLength-1; ++step)
	//		{
	//			sumij[j] += pKexi[step][i][j];
	//		}
	//		if (j == i)
	//		{
	//			sumii = sumij[j];
	//			sumij[j] = 0;
	//		}		
	//	}
	//	for(int j = 0; j < nSTATES; ++j)
	//		pTranTbl[i][j] = sumij[j] / (sumi-sumii);
	//	delete []sumij;
	//}
	//fillTran();

}



void HMModel::fillEmissionTbl(void)
{
	//for(int i = 0; i < nSTATES; ++i)
	//{
	//	int j = 0;
	//	map<int, double>::iterator it = cLRR.begin();
	//	for(; it != cLRR.end(); ++it, ++j)
	//	{
	//		fillEmissTblItem(j, i, (*it).first);
	//	}
	//}
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			fillEmissTblItem(i,j);
		}
	}
}

void HMModel::fillEmissTblItem(int site, int state)
{
	//double th = theta[state];
	//double u = 0;
	//if (state ==0)
	//{
	//	u = beta0;
	//}
	//else
	//{
	//	u = exp(beta0+beta1*log(state*1.0)
	//		+beta2*inferData.data[site].logMap
	//	+beta3*inferData.data[site].hFuntionGC);
	//}
	//pEmissTbl[site][state] = exp(
	//	MathTools::loglik_NB(1, th, &u, &inferData.data[site].count));
	double m = mu[site][state];
	double y = inferData.data[site].count;
	double w = 1; // weight is not useful when calculating the emission probability
	double e = exp(MathTools::loglik_NB(1, phi[state], &m, &y, &w));
	//pEmissTbl[site][state] = e;
	if (state==2)
	{
		pEmissTbl[site][state] = mixtureProportionNormal/largestReadCount + (1-mixtureProportionNormal)*e;
	}
	else
	{
		pEmissTbl[site][state] = mixtureProportion/largestReadCount + (1-mixtureProportion)*e;
	}
	if (pEmissTbl[site][state] <1e-12)
		pEmissTbl[site][state] = 1e-12;
	//if (pEmissTbl[site][state] <= 0)
	//{
	//	cout << "wired" << endl;
	//	cout << site << " " << state << endl;
	//	cout << m << " " << y << " " << phi[state] << endl;
	//	cout << endl;
	//}
}

void HMModel::setFileName(char * sn, char * in)
{
	snpdataname = sn;
	infodataname = in;
}

void HMModel::printVariable(void)
{
	for(int i = 0; i < nITRATION; ++i)
		cout << "The log likelihood value for the " << i+1 << " round inference is " << cLikelihood[i] << endl;
	//ofstream ofile("emissionTbl.csv");
	//for(int i = 0; i < nSTATES; ++i)
	//{
	//	for(int j = 0; j < nLength; ++j)
	//	{
	//		ofile << pEmissTbl[i][j] << ",";
	//	}
	//	ofile << endl;
	//}
	//ofile.close();
	//
	//ofile.open("initProb.csv");
	//for(int i = 0; i < nSTATES; ++i)
	//	ofile << pPi[i] << ",";
	//ofile.close();

	//ofile.open("TranProb.csv");
	//for(int i = 0; i < nSTATES; ++i)
	//{
	//	for(int j = 0; j < nSTATES; ++j)
	//	{
	//		ofile << pTranTbl[i][j] << ",";
	//	}
	//	ofile << endl;
	//}
	//ofile.close();

	//ofile.open("Tran.csv");
	//for(int i = 0; i < nSTATES; ++i)
	//{
	//	for(int j = 0; j < nSTATES; ++j)
	//	{
	//		for(int step = 0; step < nLength; ++step)
	//		{
	//			ofile << pTran[step][i][j] << ",";
	//		}
	//		ofile << endl;
	//	}
	//}
	//ofile.close();

	//ofile.open("pGa.csv");
	//for(int i = 0; i < nSTATES; ++i)
	//{
	//	for(int j = 0; j < nLength; ++j)
	//	{
	//		ofile << pGa[j][i] << ",";
	//	}
	//	ofile << endl;
	//}
	//ofile.close();

	//ofile.open("pGamma.csv");
	//for(int i = 0; i < nSTATES; ++i)
	//{
	//	for(int j = 0; j < nLength; ++j)
	//	{
	//		ofile << pGamma[j][i] << ",";
	//	}
	//	ofile << endl;
	//}
	//ofile.close();

	//ofile.open("pAlpha.csv");
	//for(int i = 0; i < nSTATES; ++i)
	//{
	//	for(int j = 0; j < nLength; ++j)
	//	{
	//		ofile << pAlpha[j][i] << ",";
	//	}
	//	ofile << endl;
	//}
	//ofile.close();

	//ofile.open("pBeta.csv");
	//for(int i = 0; i < nSTATES; ++i)
	//{
	//	for(int j = 0; j < nLength; ++j)
	//	{
	//		ofile << pBeta[j][i] << ",";
	//	}
	//	ofile << endl;
	//}
	//ofile.close();

}


void HMModel::fillTranContinous()
{
	for(int i = 0; i < nSTATES; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
		{
			pTran[0][i][j] = pTranTbl[i][j];
		}
	}
	int step = 1;
	map<int, double>::iterator it = cLRR.begin();
	int pre = (*it).first;
	++it;
	for(; it != cLRR.end(); ++it, ++step)
	{
		int cur = (*it).first;
		for(int i = 0; i < nSTATES; ++i)
		{
			double s = 0;
			for(int j = 0; j < nSTATES; ++j)
			{
				if (i != j)
				{
					double temp = pTranTbl[i][j]*(1-exp(-lamda*(cur-pre)));
					if (temp > 1e-10)
						pTran[step][i][j] = temp;
					else
						pTran[step][i][j] = 1e-10;
					s+=pTran[step][i][j];
				}
				else
				{
					pTran[step][i][j] = 0;
				}
			}
			pTran[step][i][i] = 1-s;
		}
		pre = cur;
	}
}

void HMModel::fillTranDiscrete()
{
	for(int step = 0; step < nLength; ++step)
		for(int i = 0; i < nSTATES; ++i)
			for(int j = 0; j < nSTATES; ++j)
			{
				if (pTranTbl[i][j] > 1e-10)
					pTran[step][i][j] = pTranTbl[i][j];
				else
					pTran[step][i][j] = 1e-10;
			}

}

void HMModel::findBestPath(bool viterbi)
{
	if (viterbi)
	{
		// viterbi
		double **v = new double*[nLength];  // in log scale
		int **pathm = new int*[nLength];
		for(int i = 0; i < nLength; ++i)
		{
			v[i] = new double[nSTATES];
			pathm[i] = new int[nSTATES];
		}

		for(int i = 0; i < nSTATES; ++i)
		{
			v[0][i] = log(pPi[i])+log(pEmissTbl[0][i]);
		}

		for(int i = 1; i < nLength; ++i)
		{
			for(int z = 0; z < nSTATES; ++z)
			{
				double maxValue = -1e10;
				int index = -1;
				for(int j = 0; j < nSTATES; ++j)
				{
					double value = v[i-1][j]+log(pTran[i][j][z]);
					if (value > maxValue)
					{
						maxValue = value;
						index = j;
					}
				}
				v[i][z] = maxValue + log(pEmissTbl[i][z]);
				pathm[i-1][z] = index;
			}
		}

		ofstream fout("v.log");
		for(int i = 0; i < nLength; ++i)
		{
			for(int j = 0; j < nSTATES; ++j)
				fout << v[i][j] << " ";
			fout << endl;
		}
		fout.close();	

		fout.open("path.log");
		for(int i = 0; i < nLength; ++i)
		{
			for(int j = 0; j < nSTATES; ++j)
				fout << pathm[i][j] << " ";
			fout << endl;
		}
		fout.close();



		double maxValue = -1e10;
		int index = -1;
		for(int z = 0; z < nSTATES; ++z)
		{
			if (v[nLength-1][z] > maxValue)
			{
				maxValue = v[nLength-1][z];
				index = z;
			}
		}

		inferenceResults[nLength-1] = index;
		cout << "The log probability of the most likely path is " << maxValue << endl;

		for(int i = nLength -2; i >= 0; --i)
			inferenceResults[i] = pathm[i][inferenceResults[i+1]];



		for(int i = 0; i < nLength; ++i)
		{
			delete []v[i];
			v[i] = NULL;
			delete []pathm[i];
			pathm[i] = NULL;
		}
		delete []v;
		v = NULL;
		delete []pathm;
		pathm = NULL;
	}
	else
	{
		for(int i = 0; i < nLength; ++i)
		{

			int index = 2;
			double max = pGamma[i][index];
			for(int j = 0; j < nSTATES; ++j)
			{
				if (pGamma[i][j] > max)
				{
					index = j;
					max = pGamma[i][j];
				}
			}
			inferenceResults[i] = index;
		}
	}



	//int index = -1;
	//double max = 0;
	//for(int i = 0; i < nLength; ++i)
	//{
	//	index = -1;
	//	max = -10;
	//	for(int j = 0; j < nSTATES; ++j)
	//	{
	//		//if (pGamma[i][j] > max)
	//		if (pAlpha[i][j]*pBeta[i][j] > max)
	//		{
	//			//max = pGamma[i][j];
	//			max = pAlpha[i][j]*pBeta[i][j];
	//			index = j;
	//		}
	//	}
	//	if (index == 0 || index ==1)
	//	{
	//		inferenceResults[i] = 2;
	//	}
	//	else if (index == 2)
	//	{
	//		inferenceResults[i] = 0;
	//	}
	//	else if (index == 3)
	//	{
	//		inferenceResults[i] = 1;
	//	}
	//	else if (index == 4)
	//	{
	//		inferenceResults[i] = 3;
	//	}
	//	else if (index == 5)
	//	{
	//		inferenceResults[i] = 4;
	//	}
	//}

}

void HMModel::writeKeyValue(int index)
{
	// 
	stringstream s;
	s << index;
	string postfix = s.str()+".log";
	// alpha,beta,gamma,transition,emission
	string filename("alpha");
	filename += postfix;
	ofstream fout(filename.c_str());
	fout.precision(14);
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
			fout << pAlpha[i][j] << " ";
		fout << endl;
	}
	fout.close();

	filename = "beta";
	filename += postfix;
	fout.open(filename.c_str());
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
			fout << pBeta[i][j] << " ";
		fout << endl;
	}
	fout.close();

	filename = "gamma";
	filename += postfix;
	fout.open(filename.c_str());
	for(int i = 0; i < nLength; ++i)
	{
		double sum = 0;
		for(int j = 0; j < nSTATES; ++j)
		{
			fout << exp(pGamma[i][j]) << " ";
			sum += exp(pGamma[i][j]);
		}
		fout << sum;
		fout << endl;
	}
	fout.close();

	filename = "transition";
	filename += postfix;
	fout.open(filename.c_str());
	for(int i = 0; i < nSTATES; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
			fout << pTranTbl[i][j] << " ";
		fout << endl;
	}
	fout.close();

	filename = "emission";
	filename += postfix;
	fout.open(filename.c_str());
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
			fout << log(pEmissTbl[i][j]) << " ";
		fout << endl;
	}
	fout.close();

	filename = "mu";
	filename += postfix;
	fout.open(filename.c_str());
	for(int i = 0; i < nLength; ++i)
	{
		for(int j = 0; j < nSTATES; ++j)
			fout << mu[i][j] << " ";
		fout << endl;
	}
	fout.close();

	filename = "phi";
	filename += postfix;
	fout.open(filename.c_str());
	for(int j = 0; j < nSTATES; ++j)
		fout << phi[j] << " ";
	fout.close();

	filename = "init";
	filename += postfix;
	fout.open(filename.c_str());
	for(int j = 0; j < nSTATES; ++j)
		fout << pPi[j] << " ";
	fout.close();
}

////////////////////////////////////// codes that are not active ///////////////////


//void HMModel::calculateMuAndPhi(bool init)
//{
//	// use glmNB to fill mu matrix and phi matrix 
//	// for state 0, we only need to fit phi,
//	// for other state, we need to get fitted mu and re-estimated phi
//
//	// for state 0
//	// load weights, load fitted value, in original scale, estimate phi
//	int maxIt = 25;
//	double convR = 1e-8;
//	int nCovariate = 1;          // logMap
//	double *y = new double[nLength];
//	double *fitted = new double[nLength];
//	double *x = new double[nLength*nCovariate];
//	double *prior = new double[nLength];
//	double *weights = new double[nLength]; // will be changed in computation
//	double *resid = new double[nLength]; // not use
//	double *Xb = new double[nLength*nCovariate];  // used in the program
//
//
//	for(int i = 0; i < nLength; ++i)
//		y[i] = inferData.data[i].count;
//
//	for(int i = 0; i < nLength; ++i)
//	{
//		prior[i] = exp(pGamma[i][0]);
//	}
//
//	for(int i = 0; i < nLength; ++i)
//		x[i] = inferData.data[i].logMap;
//
//
//	int dim[5];
//	dim[0] = nLength;
//	dim[1] = nCovariate;     // 1 covariate  
//	dim[2] = maxIt; // 
//	dim[3] = 0; // false
//	dim[4] = 0; // false, no offset
//
//	int nIter = 0;
//	int linkR = LOG;
//	convR = 1e-8; 
//	int rank = 0;
//	double PHI = 0;  // not use
//	double scale = 1.0;
//	int de_resid = 0;
//	int family = 0;
//	double twologlik = 0;
//	double scoreTestP = 0;
//	int trace = 0;
//	double beta = 0;
//
//
//	int conv = MathTools::glmNB(dim, &nIter,y,prior, &linkR, NULL, x, &convR, &rank,
//		Xb, fitted, resid, weights, &PHI, &scale, &de_resid, &family,
//		&twologlik, &scoreTestP, &trace, &beta);
//	//cout << "beta for state 0 cn is " << beta << " ";
//
//
//
//	for(int i = 0; i < nLength; ++i)
//	{
//		prior[i] = exp(pGamma[i][0]);
//	}
//
//	for(int i = 0; i < nLength; ++i)
//	{
//		mu[i][0] = fitted[i];
//	}
//	
//	int cvPhi = MathTools::phi_ml(y, fitted, nLength, prior, maxIt, convR,
//		&phi[0], 0, 0);
//	
//	
//	delete []y; y=NULL;
//	delete []fitted; fitted=NULL;
//	delete []x; x=NULL;
//	delete []prior; prior=NULL;
//	delete []weights; weights=NULL;
//	delete []resid; resid=NULL;
//	delete []Xb; Xb=NULL;
//
//	//// for states other than 0,
//	//// load y, load x, glmNB to get fitted value (in original scale)
//	//// estimate phi separately
//	//
//	//create x, y, z, fitted,weights, residual
//
//	nCovariate = 2;     // log(state), logMappability
//	y = new double[(nSTATES-1)*nLength];
//	x = new double[(nSTATES-1)*nLength*nCovariate]; // 3 covariates
//	prior = new double[(nSTATES-1)*nLength];
//	fitted = new double[(nSTATES-1)*nLength];
//	weights = new double[(nSTATES-1)*nLength]; // temp value
//	resid = new double[(nSTATES-1)*nLength]; // not use
//	Xb = new double[(nSTATES-1)*nLength*nCovariate];  // not use
//
//	// loaded x, y, z, fitted, weights
//	// y x are not dependant with initial value
//	int index = 0;
//	for(int i = 0; i < nLength; ++i)
//	{
//		for(int j = 1; j < nSTATES; ++j)
//			y[index++] = inferData.data[i].count;
//	}
//	index = 0;
//	//for(int i = 0; i < nLength; ++i)
//	//{
//	//	for(int j = 1; j < nSTATES; ++j)
//	//	{
//	//		x[index++] = inferData.data[i].hFuntionGC;
//	//	}
//	//}
//	for(int i = 0; i < nLength; ++i)
//	{
//		for(int j = 1; j < nSTATES; ++j)
//		{
//			x[index++] = inferData.data[i].logMap;
//		}
//	}
//	for(int i = 0; i < nLength; ++i)
//	{
//		for(int j = 1; j < nSTATES; ++j)
//		{
//			x[index++] = log(j*1.0);
//		}
//	}
//
//	// load weights
//	index = 0;
//	for(int i = 0; i < nLength; ++i)
//	{
//		for(int j = 1; j < nSTATES; ++j)
//		{
//			prior[index++] = exp(pGamma[i][j]);
//		}
//	}
//
//	dim[0] = nLength*(nSTATES-1);
//	dim[1] = nCovariate;  
//	dim[2] = maxIt; //
//	dim[3] = 0; // false
//	dim[4] = 0; // false, no offset
//
//	nIter = 0;
//	linkR = LOG;
//	convR = 1e-8; 
//	rank = 0;
//	PHI = 0;  // not use
//	scale = 1.0;
//	de_resid = 0;
//	family = 0;
//	twologlik = 0;
//	scoreTestP = 0;
//	trace = 0;
//	beta = 0;
//
//	conv = MathTools::glmNB(dim, &nIter,y,prior, &linkR, NULL, x, &convR, &rank,
//		Xb, fitted, resid, weights, &PHI, &scale, &de_resid, &family,
//		&twologlik, &scoreTestP, &trace, &beta);
//
//	cout << "beta for states1-6 cn is " << beta << " ";
//
//
//	// save fitted
//	index = 0;
//	for(int i = 0; i < nLength; ++i)
//	{
//		for(int j = 1; j < nSTATES; ++j)
//		{
//			mu[i][j] = fitted[index++];
//		}
//	}
//
//	// calculate the mean square root of errors for fitting
//	double error = 0;
//	for(int i = 0; i < nLength; ++i)
//	{
//		double e = 0;
//		for(int j = 0; j < nSTATES; ++j)
//		{
//			e += (mu[i][j]-inferData.data[i].count)*(mu[i][j]-inferData.data[i].count);
//		}
//		e /= nSTATES;
//		error += e;
//	}
//	error /= nLength;
//	error = sqrt(error);
//	//cout << "MSE of fitted data is " << error << " ";
//
//	// calculate the weighted mean square root of errors for fitting
//	error = 0;
//	for(int i = 0; i < nLength; ++i)
//	{
//		double e = 0;
//		for(int j = 0; j < nSTATES; ++j)
//		{
//			e += exp(pGamma[i][j])*(mu[i][j]-inferData.data[i].count)*(mu[i][j]-inferData.data[i].count);
//		}
//		error += e;
//	}
//	error /= nLength;
//	error = sqrt(error);
//	//cout << "weigthed MSE of fitted data is " << error << endl;
//
//
//	delete []y; y=NULL;
//	delete []x; x=NULL;
//	delete []prior; prior=NULL;
//	delete []fitted; fitted=NULL;
//	delete []weights; weights=NULL;
//	delete []resid; resid=NULL;
//	delete []Xb; Xb=NULL;
//
//	// run phi_ml several times to get phi
//	y = new double[nLength];
//	fitted = new double[nLength];
//	prior = new double[nLength];
//	for(int j = 1; j < nSTATES; ++j)
//	{
//		for(int i = 0; i < nLength; ++i)
//		{
//			y[i] = inferData.data[i].count;
//			fitted[i] = mu[i][j];
//			prior[i] = exp(pGamma[i][j]);
//		}
//		cvPhi = MathTools::phi_ml(y, fitted, nLength, prior, maxIt,
//			convR, &phi[j], 0, 0);
//	}
//
//	delete []y; y=NULL;
//	delete []fitted; fitted=NULL;
//	delete []prior; prior=NULL;
//
//
//	//for(int i = 0; i < nSTATES; ++i)
//	//	phi[i] = 1.0/32.65;
//
////#ifdef _DEBUG
////	ofstream fout("ob.txt");
////	for(int i = 0; i < nLength; ++i)
////		fout << inferData.data[i].count << " " << inferData.data[i].hFuntionGC << " " << inferData.data[i].logMap << endl;
////	fout.close();
////	fout.open("fitted.txt");
////	fout.precision(20);
////	for(int i = 0; i < nLength; ++i)
////	{
////		for(int j = 1; j < nSTATES; ++j)
////			fout << mu[i][j] << " ";
////		fout << endl;
////	}
////	fout.close();
////	fout.open("theta.txt");
////	for(int i = 1; i < nSTATES; ++i)
////	{
////		fout << 1/phi[i] << endl;
////	}
////	fout.close();
////	fout.open("prior.txt");
////	for(int i = 0; i < nLength; ++i)
////	{
////		for(int j = 1; j < nSTATES; ++j)
////			fout << pGamma[i][j] << " ";
////		fout << endl;
////	}
////	fout.close();
////
////#endif
//
//
//}

//void HMModel::calculateMuAndPhiWithAutoRegression()
//{
//	// use glmNB to fill mu matrix and phi matrix 
//	// for state 0, we only need to fit phi,
//	// for other state, we need to get fitted mu and re-estimated phi
//
//	// for state 0
//	// load weights, load fitted value, in original scale, estimate phi
//	int maxIt = 25;
//	int nCovariate = 2;          // logMap and residual
//	double convR = 1e-8;
//	double *y = new double[nLength];
//	double *fitted = new double[nLength];
//	double *x = new double[nLength*nCovariate];  // two covariates, mappability and residual 
//	double *prior = new double[nLength];
//	double *weights = new double[nLength]; // will be changed in computation
//	double *resid = new double[nLength]; // not use
//	double *Xb = new double[nLength*nCovariate];  // 
//
//
//	for(int i = 0; i < nLength; ++i)
//		y[i] = inferData.data[i].count;
//
//	for(int i = 0; i < nLength; ++i)
//	{
//		prior[i] = exp(pGamma[i][0]);
//	}
//
//	int index = 0;
//	for(int i = 0; i < nLength; ++i)
//		x[index++] = inferData.data[i].logMap;
//	for(int i = 0; i < nLength; ++i)
//	{
//		if (i==0)
//			x[index++] = 0;
//		else
//		{
//			if (inferData.data[i-1].count == 0)
//			{
//				x[index++] = log(inferData.data[i-1].count+0.1)-log(mu[i-1][0]);
//			}
//			else
//			{
//				x[index++] = log(inferData.data[i-1].count)-log(mu[i-1][0]);
//			}
//		}
//	}
//
//	//debug
//	//ofstream fout("temp.txt");
//	//fout << "count logMap autoR weights" << endl;
//	//for(int i = 0; i < nLength; ++i)
//	//{
//	//	fout << y[i] << " " << x[i] << " " << x[i+1000] << " " << prior[i] << endl;
//	//}
//	//fout.close();
//
//	/////////////////////////////////////////////////////
//
//
//	int dim[5];
//	dim[0] = nLength;
//	dim[1] = nCovariate;     // 2 covariate  
//	dim[2] = maxIt; // 
//	dim[3] = 0; // false
//	dim[4] = 0; // false, no offset
//
//	int nIter = 0;
//	int linkR = LOG;
//	convR = 1e-8; 
//	int rank = 0;
//	double PHI = 0;  // not use
//	double scale = 1.0;
//	int de_resid = 0;
//	int family = 0;
//	double twologlik = 0;
//	double scoreTestP = 0;
//	int trace = 0;
//	double beta = 0;
//
//
//	//double *z = new double[nLength];
//
//	//int conv = MathTools::glmNB(true, dim, &nIter,y,z, &linkR, NULL, x, &convR, &rank,
//	//	Xb, fitted, resid, weights, &PHI, &scale, &de_resid, &family,
//	//	&twologlik, &scoreTestP, &trace, &beta);
//
//
//
//	int conv = MathTools::glmNB(dim, &nIter,y,prior, &linkR, NULL, x, &convR, &rank,
//		Xb, fitted, resid, weights, &PHI, &scale, &de_resid, &family,
//		&twologlik, &scoreTestP, &trace, &beta);
//
//	cout << "beta for state 0 residual is " << beta << " ";
//
//
//
//	for(int i = 0; i < nLength; ++i)
//	{
//		prior[i] = exp(pGamma[i][0]);
//	}
//
//	for(int i = 0; i < nLength; ++i)
//	{
//		mu[i][0] = fitted[i];
//	}
//	
//	int cvPhi = MathTools::phi_ml(y, fitted, nLength, prior, maxIt, convR,
//		&phi[0], 0, 0);
//	
//	
//	delete []y; y=NULL;
//	delete []fitted; fitted=NULL;
//	delete []x; x=NULL;
//	delete []prior; prior=NULL;
//	delete []weights; weights=NULL;
//	delete []resid; resid=NULL;
//	delete []Xb; Xb=NULL;
//
//	//cout << "beta of residual on state 0 is " << beta << " ";
//
//	//// for states other than 0,
//	//// load y, load x, glmNB to get fitted value (in original scale)
//	//// estimate phi separately
//	//
//	//create x, y, z, fitted,weights, residual
//
//	nCovariate = 3;     // log(state), logMappability, residual
//	y = new double[(nSTATES-1)*nLength];
//	x = new double[(nSTATES-1)*nLength*nCovariate]; // 3 covariates
//	prior = new double[(nSTATES-1)*nLength];
//	fitted = new double[(nSTATES-1)*nLength];
//	weights = new double[(nSTATES-1)*nLength]; // temp value
//	resid = new double[(nSTATES-1)*nLength]; // not use
//	Xb = new double[(nSTATES-1)*nLength*nCovariate];  // used in the program
//
//	// loaded x, y, z, fitted, weights
//	// y x are not dependant with initial value
//	index = 0;
//	for(int i = 0; i < nLength; ++i)
//	{
//		for(int j = 1; j < nSTATES; ++j)
//			y[index++] = inferData.data[i].count;
//	}
//	index = 0;
//	//for(int i = 0; i < nLength; ++i)
//	//{
//	//	for(int j = 1; j < nSTATES; ++j)
//	//	{
//	//		x[index++] = inferData.data[i].hFuntionGC;
//	//	}
//	//}
//	for(int i = 0; i < nLength; ++i)
//	{
//		for(int j = 1; j < nSTATES; ++j)
//		{
//			x[index++] = inferData.data[i].logMap;
//		}
//	}
//	for(int i = 0; i < nLength; ++i)
//	{
//		for(int j = 1; j < nSTATES; ++j)
//		{
//			x[index++] = log(j*1.0);
//		}
//	}
//	for(int i = 0; i < nLength; ++i)
//	{
//		for(int j = 1; j < nSTATES; ++j)
//		{
//			if (i==0)
//				x[index++] = 0;
//			else
//			{
//				if (inferData.data[i-1].count == 0)
//				{
//					x[index++] = log(inferData.data[i-1].count+0.1)-log(mu[i-1][j]);
//				}
//				else
//				{
//					x[index++] = log(inferData.data[i-1].count)-log(mu[i-1][j]);
//				}
//			}
//		}
//	}
//
//
//
//	// load weights
//	index = 0;
//	for(int i = 0; i < nLength; ++i)
//	{
//		for(int j = 1; j < nSTATES; ++j)
//		{
//			prior[index++] = exp(pGamma[i][j]);
//		}
//	}
//
//	dim[0] = nLength*(nSTATES-1);
//	dim[1] = nCovariate;  
//	dim[2] = maxIt; //
//	dim[3] = 0; // false
//	dim[4] = 0; // false, no offset
//
//	nIter = 0;
//	linkR = LOG;
//	convR = 1e-8; 
//	rank = 0;
//	PHI = 0;  // not use
//	scale = 1.0;
//	de_resid = 0;
//	family = 0;
//	twologlik = 0;
//	scoreTestP = 0;
//	trace = 0;
//	beta = 0;
//
//	conv = MathTools::glmNB(dim, &nIter,y,prior, &linkR, NULL, x, &convR, &rank,
//		Xb, fitted, resid, weights, &PHI, &scale, &de_resid, &family,
//		&twologlik, &scoreTestP, &trace, &beta);
//
//	cout << "beta for residual on states1-6 is " << beta << " ";
//	//cout << "beta for states1-6 cn is " << beta << " ";
//
//
//	// save fitted
//	index = 0;
//	for(int i = 0; i < nLength; ++i)
//	{
//		for(int j = 1; j < nSTATES; ++j)
//		{
//			mu[i][j] = fitted[index++];
//		}
//	}
//
//	// calculate the mean square root of errors for fitting
//	double error = 0;
//	for(int i = 0; i < nLength; ++i)
//	{
//		double e = 0;
//		for(int j = 0; j < nSTATES; ++j)
//		{
//			e += (mu[i][j]-inferData.data[i].count)*(mu[i][j]-inferData.data[i].count);
//		}
//		e /= nSTATES;
//		error += e;
//	}
//	error /= nLength;
//	error = sqrt(error);
//	cout << "MSE of fitted data is " << error << " ";
//
//	// calculate the weighted mean square root of errors for fitting
//	error = 0;
//	for(int i = 0; i < nLength; ++i)
//	{
//		double e = 0;
//		for(int j = 0; j < nSTATES; ++j)
//		{
//			e += exp(pGamma[i][j])*(mu[i][j]-inferData.data[i].count)*(mu[i][j]-inferData.data[i].count);
//		}
//		error += e;
//	}
//	error /= nLength;
//	error = sqrt(error);
//	cout << "weigthed MSE of fitted data is " << error << endl;
//
//
//	delete []y; y=NULL;
//	delete []x; x=NULL;
//	delete []prior; prior=NULL;
//	delete []fitted; fitted=NULL;
//	delete []weights; weights=NULL;
//	delete []resid; resid=NULL;
//	delete []Xb; Xb=NULL;
//
//	// run phi_ml several times to get phi
//	y = new double[nLength];
//	fitted = new double[nLength];
//	prior = new double[nLength];
//	for(int j = 1; j < nSTATES; ++j)
//	{
//		for(int i = 0; i < nLength; ++i)
//		{
//			y[i] = inferData.data[i].count;
//			fitted[i] = mu[i][j];
//			prior[i] = exp(pGamma[i][j]);
//		}
//		cvPhi = MathTools::phi_ml(y, fitted, nLength, prior, maxIt,
//			convR, &phi[j], 0, 0);
//	}
//
//	delete []y; y=NULL;
//	delete []fitted; fitted=NULL;
//	delete []prior; prior=NULL;
//}

//void HMModel::loadData(char * sn, char * in, string chr)
//{
//	//"E:\\Project with Jin-HSMM\\TestingData\\NA12878.null.null2.cnv.dat2.txt","E:\\Project with Jin-HSMM\\TestingData\\hhall.hg18.null2.pfb"
//	setFileName(sn, in);
//	chrSymbol = chr;
//	setVariables();
//	//fillTranContinous();
//	fillTranDiscrete();  // for read depth data, it is discrete time
//	fillEmissionTbl();
//}

//void HMModel::setVariables(void)
//{
//	// hard coding for the real data
//	nSTATES = 6;
//	nLength = 0;
//
//	ifstream ifile(infodataname);
//	string name;
//	vector<string> markerName;
//	vector<int> POS;
//	//map<int, int> dupPOS;  // to deal with duplicate pos
//	vector<double> PFB;
//	map<string, LRRBAF> markerValue;
//
//	int pos;
//	double pfb;
//	string curchr;
//	// deal with first row
//	ifile >> name >> name >> name >> name;
//	while(!ifile.eof())
//	{
//		// assuming each row is good
//		// markers from the same chromosome are coming in ascending pos order
//		ifile >> name >> curchr >> pos >> pfb;
//		// at this point not dealing with
//		if (curchr != chrSymbol)
//		{
//			continue;
//		}
//		//if (dupPOS.find(pos) == dupPOS.end())   // not a duplicate position
//		//{
//			markerName.push_back(name);
//			POS.push_back(pos);
//			PFB.push_back(pfb);
//			LRRBAF temp;
//			markerValue[name] = temp;
//		//}
//	}
//	// check the problem of reading the last row twice
//	if (markerName.size() > 1)
//	{
//		// markerName, POS, PFB should have the same size
//		if (POS[POS.size()-2] == POS[POS.size()-1]
//		&&  PFB[PFB.size()-2] == PFB[PFB.size()-1])
//		{
//			// last row is read twice
//			markerName.pop_back();
//			POS.pop_back();
//			PFB.pop_back();
//		}
//	}
//	ifile.clear();
//	ifile.close();
//
//
//	ifile.open(snpdataname);
//	double l;
//	double b;
//	// read the first row
//	ifile >> name >> name >> name;
//	// dealing with incomplete data
//	string temp1, temp2;
//	while (!ifile.eof())
//	{
//		ifile >> name >> temp1 >> temp2;
//		// we have a good data ?
//		if (temp1 == "NA" || temp2 == "NA")
//			continue;
//		// if the data we are interested in?
//		if (markerValue.find(name) != markerValue.end())
//		{
//			l = atof(temp1.c_str());
//			b = atof(temp2.c_str());
//			markerValue[name].LRR = l;
//			markerValue[name].BAF = b;
//		}
//	}
//	ifile.clear();
//	ifile.close();
//
//	// orgnizing the data
//	for(int i = 0; i < (int)markerName.size(); ++i)
//	{
//		if (markerValue[markerName[i]].LRR != -10.0
//			&&  markerValue[markerName[i]].BAF != -10.0)
//		{
//			// we have complete data for this marker
//			if (cLRR.find(POS[i]) == cLRR.end())
//			{
//				// not a duplicated record
//				cName.push_back(markerName[i]);
//				cPos.push_back(POS[i]);
//				cLRR[POS[i]] = markerValue[markerName[i]].LRR;
//				cBAF[POS[i]] = markerValue[markerName[i]].BAF;
//				cPFB[POS[i]] = PFB[i];
//			}
//		}
//	}
//	assert(cName.size() == cPos.size());
//	assert(cName.size() == cLRR.size());
//	nLength = cPFB.size();
//
//	inferenceResults = new int[nLength];
//
//	pTranTbl = new double*[nSTATES];
//	for(int i = 0; i < nSTATES; ++i)
//	{
//		pTranTbl[i] = new double[nSTATES];
//	}
//
//	pTran = new double **[nLength];
//	for(int i = 0; i < nLength; ++i)
//	{
//		pTran[i] = new double *[nSTATES];
//		for(int j = 0; j < nSTATES; ++j)
//		{
//			pTran[i][j] = new double[nSTATES];
//		}
//	}
//
//
//	// hard coding for real data
//	//pTranTbl[0][0] = 0.00;pTranTbl[0][1] = 0.01;pTranTbl[0][2] = 0.09;pTranTbl[0][3] = 0.80;pTranTbl[0][4] = 0.09;pTranTbl[0][5] = 0.01;
//	//pTranTbl[1][0] = 0.20;pTranTbl[1][1] = 0.00;pTranTbl[1][2] = 0.20;pTranTbl[1][3] = 0.20;pTranTbl[1][4] = 0.20;pTranTbl[1][5] = 0.20;
//	//pTranTbl[2][0] = 0.96;pTranTbl[2][1] = 0.01;pTranTbl[2][2] = 0.00;pTranTbl[2][3] = 0.01;pTranTbl[2][4] = 0.01;pTranTbl[2][5] = 0.01;
//	//pTranTbl[3][0] = 0.96;pTranTbl[3][1] = 0.01;pTranTbl[3][2] = 0.01;pTranTbl[3][3] = 0.00;pTranTbl[3][4] = 0.01;pTranTbl[3][5] = 0.01;
//	//pTranTbl[4][0] = 0.96;pTranTbl[4][1] = 0.01;pTranTbl[4][2] = 0.01;pTranTbl[4][3] = 0.01;pTranTbl[4][4] = 0.00;pTranTbl[4][5] = 0.01;
//	//pTranTbl[5][0] = 0.96;pTranTbl[5][1] = 0.01;pTranTbl[5][2] = 0.01;pTranTbl[5][3] = 0.01;pTranTbl[5][4] = 0.01;pTranTbl[5][5] = 0.00;
//	for(int i = 0; i < nSTATES; ++i)
//	{
//		for(int j = 0; j < nSTATES; ++j)
//		{
//			if (i==j)
//				pTranTbl[i][j] = 0.9;
//			else
//				pTranTbl[i][j] = 0.1/5;
//		}
//	}
//
//
//
//	// revise pEmissTbl nSTATES * nLength
//	pEmissTbl = new double*[nLength];
//	for(int i = 0; i < nLength; ++i)
//	{
//		pEmissTbl[i] = new double[nSTATES];
//	}
//
//	// hard coding for real data
//	pPi = new double[nSTATES];
//	pPi[0] = 0.995;
//	pPi[1] = 5e-5;
//	pPi[2] = 0.00045;
//	pPi[3] = 0.004;
//	pPi[4] = 0.00045;
//	pPi[5] = 5e-5;
//
//	pAlpha = new double*[nLength];
//	for(int i = 0; i < nLength; ++i)
//		pAlpha[i] = new double[nSTATES];
//
//	pBeta = new double*[nLength];
//	for(int i = 0; i < nLength; ++i)
//		pBeta[i] = new double[nSTATES];
//
//	pGamma = new double*[nLength];
//	for(int i = 0; i < nLength; ++i)
//		pGamma[i] = new double[nSTATES];
//
//	pKexi = new double **[nLength];
//	for(int i = 0; i < nLength; ++i)
//	{
//		pKexi[i] = new double*[nSTATES];
//		for(int j = 0; j < nSTATES; ++j)
//		{
//			pKexi[i][j] = new double[nSTATES];
//		}
//	}
//
//	pGa = new double *[nLength];
//	for(int i = 0; i < nLength; ++i)
//		pGa[i] = new double[nSTATES];
//
//
//	// lots of new variables
//	pircn = new double [nSTATES];
//	murcn = new double [nSTATES];
//	sdrcn = new double [nSTATES];
//	pir = new double [nSTATES];
//	mur = new double [nSTATES];
//	sdr = new double [nSTATES];
//	pib = new double [nSTATES];
//	mub = new double *[nSTATES];
//	for(int i = 0; i < nSTATES; ++i)
//		mub[i] = new double [nSTATES];
//	sdb = new double*[nSTATES];
//	for(int i = 0; i < nSTATES; ++i)
//		sdb[i] = new double [nSTATES];
//
//	pir[0] = 0.01; pir[1] = 0.01; pir[2] = 0.1; pir[3] = 0.01; pir[4] = 0.01; pir[5] = 0.01;
//	mur[0] = 0.0; mur[1] = 0.0; mur[2] = -3.5; mur[3] = -0.67; mur[4] = 0.4; mur[5] = 0.68;
//	sdr[0] = 0.15; sdr[1] = 0.15; sdr[2] = 1.3; sdr[3] = 0.28; sdr[4] = 0.2; sdr[5] = 0.19;
//	pircn[0] = 0.01; pircn[1] = 0.01; pircn[2] = 0.1; pircn[3] = 0.01; pircn[4] = 0.01; pircn[5] = 0.01;
//	murcn[0] = 0.0; murcn[1] = 0.0; murcn[2] = -2.1; murcn[3] = -0.58; murcn[4] = 0.37; murcn[5] = 0.65;
//	sdrcn[0] = 0.18; sdrcn[1] = 0.18; sdrcn[2] = 2.1; sdrcn[3] = 0.35; sdrcn[4] = 0.24; sdrcn[5] = 0.3;
//	pib[0] = 0.01; pib[1] = 0.01; pib[2] = 0.5; pib[3] = 0.01; pib[4] = 0.01; pib[5] = 0.01;
//
//	mub[0][0] = 0.0;mub[0][1] = 0.5;mub[0][2] = 1.0;mub[0][3] = -1.0;mub[0][4] = -1.0;mub[0][5] = -1.0;
//	mub[1][0] = 0.0;mub[1][1] = -1.0;mub[1][2] = -1.0;mub[1][3] = 1.0;mub[1][4] = -1.0;mub[1][5] = -1.0;
//	mub[2][0] = 0.5;mub[2][1] = -1.0;mub[2][2] = -1.0;mub[2][3] = -1.0;mub[2][4] = -1.0;mub[2][5] = -1.0;
//	mub[3][0] = 0.0;mub[3][1] = -1.0;mub[3][2] = -1.0;mub[3][3] = 1.0;mub[3][4] = -1.0;mub[3][5] = -1.0;
//	mub[4][0] = 0.0;mub[4][1] = 0.333;mub[4][2] = 0.667;mub[4][3] = 1.0;mub[4][4] = -1.0;mub[4][5] = -1.0;
//	mub[5][0] = 0.0;mub[5][1] = 0.25;mub[5][2] = 0.5;mub[5][3] = 0.75;mub[5][4] = 1.0;mub[5][5] = -1.0;
//
//	sdb[0][0] = 0.016; sdb[0][1] = 0.035; sdb[0][2] = 0.016; sdb[0][3] = -1.0; sdb[0][4] = -1.0; sdb[0][5] = -1.0;
//	sdb[1][0] = 0.016; sdb[1][1] = -1.0; sdb[1][2] = -1.0; sdb[1][3] = 0.16; sdb[1][4] = -1.0; sdb[1][5] = -1.0;
//	sdb[2][0] = 0.15; sdb[2][1] = -1.0; sdb[2][2] = -1.0; sdb[2][3] = -1.0; sdb[2][4] = -1.0; sdb[2][5] = -1.0;
//	sdb[3][0] = 0.016; sdb[3][1] = -1.0; sdb[3][2] = -1.0; sdb[3][3] = 0.16; sdb[3][4] = -1.0; sdb[3][5] = -1.0;
//	sdb[4][0] = 0.016; sdb[4][1] = 0.042; sdb[4][2] = 0.042; sdb[4][3] = 0.016; sdb[4][4] = -1.0; sdb[4][5] = -1.0;
//	sdb[5][0] = 0.016; sdb[5][1] = 0.042; sdb[5][2] = 0.035; sdb[5][3] = 0.42; sdb[5][4] = 0.016; sdb[5][5] = -1.0;
//
//}

//void HMModel::fillEmissTblItem(int site, int state, int pos)
//{
//	// hard coding for real data
//	if (cPFB[pos] > 1)
//	{
//		double temp = calLRR(cLRR[pos], state, pircn, murcn, sdrcn);
//		if (temp > 1e-10)
//			pEmissTbl[site][state] = temp;
//		else
//			pEmissTbl[site][state] = 1e-10;
//	}
//	else
//	{
//		double lrr = calLRR(cLRR[pos], state, pir, mur, sdr);
//		double baf = calBAF(cBAF[pos], state, cPFB[pos], pib, mub, sdb);
//		double temp = lrr*baf;
//		if (temp > 1e-10)
//			pEmissTbl[site][state] = temp;
//		else
//			pEmissTbl[site][state] = 1e-10;
//	}
//
//}


//double HMModel::calLRR(double r, int z, double *pir, double *mur, double *sdr )
//{
//	double v = MathTools::dnorm(r, mur[z], sdr[z]);
//	return (pir[z]+(1-pir[z])*v);
//}
//
//double HMModel::calBAF(double b, int z, double pb, double *pib, double **mub, double **sdb)
//{
//	int H[6] = {3,2,1,2,4,5}; // how many genotypes does each state have
//	int idx[6][6] = {{1,2,3,0,0,0},{1,4,0,0,0,0},{1,0,0,0,0,0},{1,4,0,0,0,0},{1,2,3,4,0,0},{1,2,3,4,5,0}}; // where each genotype information will be stored
//	double ws[6] = {0,0,0,0,0,0};
//	if (z==0)
//	{
//		for(int k = 0; k <= 2; ++k)
//			ws[k] = MathTools::dbinom(k, 2, pb);
//	}
//	else if (z == 1 || z == 3)
//	{
//		ws[0]  = 1-pb;
//		ws[1] = pb;
//	}
//	else if (z == 2)
//	{
//		ws[0] = 1;
//	}
//	else if (z==4)
//	{
//		for(int k = 0; k <= 3; ++k)
//			ws[k] = MathTools::dbinom(k, 3, pb);
//	}
//	else if (z == 5)
//	{
//		for(int k = 0; k <= 4; ++k)
//			ws[k] = MathTools::dbinom(k, 4, pb);
//	}
//	double v = 0;
//	if ((b > 0 && b < 1) ||  z== 2)
//	{
//		v = pib[z];
//		for(int hh = 0; hh < H[z]; ++hh)
//		{
//			v+= (1-pib[z])*ws[hh]*MathTools::dnorm(b, mub[z][idx[z][hh]-1], sdb[z][idx[z][hh]-1]);
//		}
//	}
//	else if (b == 0) // All A
//	{
//		v = (1-pib[z])*ws[0];
//	}
//	else if (b == 1) // All B what is pib?
//	{
//		v = (1-pib[z])*ws[H[z]-1];
//	}
//	return v;
//}


