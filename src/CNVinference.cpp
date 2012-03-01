// CNVinference.cpp : Defines the entry point for the console application.
//

#include "HMModel.h"
#include "stdlib.h"


int main(int argc, char* argv[])
{
	//if (argc != 3)
	//	return 1;
	HMModel model;
	//model.loadData("simu1_dat.txt", "simu1_inf.txt", "22");
	//model.loadData("NA12878.null.null2.cnv.dat1.txt", "hhall.hg18.null2.pfb", "22");
	//model.loadReadDepthData("data4hmm_chr17revise.txt");
	
	// current codes

	model.mixtureProportion = atof(argv[3]);
	model.mixtureProportionNormal = atof(argv[4]);

	if (argc > 5)
	{
		if (atoi(argv[5])==0)
			model.USINGMAPPABILITY = false;
		if (atoi(argv[6])==0)
			model.USINGAUTOREGRESSION = false;
		if (atoi(argv[7])==0)
			model.USINGMIXTURECOMPONENT = false;
	}
	if (argc > 8)
	{
		if (atoi(argv[8])==0)
			model.REESTIMATETRANSITION = false;
	}

	model.loadReadDepthData(argv[1]);
	model.inferAndEstimation(atoi(argv[2]));

	// for debug
	//model.loadReadDepthData("sim2_0_5.txt");
	//model.inferAndEstimation(1);

	model.writeResult();
	return 0;
}

