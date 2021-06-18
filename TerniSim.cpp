#include "headSim.h"



	//------------------------------------Global variables-------------------------//
	int32 seed = (int32)time(0);
	CRandomMersenne rg(seed);
	Pars pa;//parameter in struct
	initials in; //initial conditions
	totals to;
	Times ti;
	//------------------------------------------------------------------------------//

	int main(void)
		{
		//---------------------------Input parameters----------------------------//
			cin>>ti.maxT;
			cin>>ti.NumRuns;
			cin>>in.startNum;
			cin>>in.driver_time;
			cin>>pa.Fgamma;
			cin>>pa.Mgamma;
			cin>>pa.beta;
			cin>>pa.theta;
			cin>>pa.delta;
						
			cin>>pa.xi;
			cin>>pa.Frho;
			cin>>pa.Mrho;

			cin>>pa.ef;
			cin>>pa.em;
			cin>>pa.lamFStart;
			cin>>pa.lamFEnd;
			cin>>pa.lamMStart;
			cin>>pa.lamMEnd;
			cin>>pa.kFStart;
			cin>>pa.kFEnd;
			cin>>pa.kMStart;
			cin>>pa.kMEnd;
			cin>>pa.NumEggs;
			cin>>pa.numDsx;
		//-------------------------determine fitnesses from parameters----------------//
			pa.omegaww=1;
			pa.omegafww=(1-pa.Frho); pa.omegamww=(1-pa.Mrho); pa.omegabww=(1-pa.Frho)*(1-pa.Mrho);
			pa.omegafwd=(1-pa.Frho)*(1-pa.xi); pa.omegamwd=(1-pa.Mrho)*(1-pa.xi); pa.omegabwd=(1-pa.Frho)*(1-pa.Mrho)*(1-pa.xi);
			pa.omegafwr=(1-pa.Frho); pa.omegamwr=(1-pa.Mrho); pa.omegabwr=(1-pa.Frho)*(1-pa.Mrho);
			pa.omegadd=0; pa.omegadr=0; pa.omegarr=0;
		//----------------------------------------------------------------------------//
		//-------------------Compute time dependent adult mortaliity parameter vectors--//
			for(int i=0;i<ti.maxT+1;i++)
			{
			pa.lamF[i]=pa.lamFStart+1.0*i*(pa.lamFEnd-pa.lamFStart)/double(ti.maxT);
			pa.lamM[i]=pa.lamMStart+1.0*i*(pa.lamMEnd-pa.lamMStart)/double(ti.maxT);
			pa.kF[i]=pa.kFStart+1.0*i*(pa.kFEnd-pa.kFStart)/double(ti.maxT);
			pa.kM[i]=pa.kMStart+1.0*i*(pa.kMEnd-pa.kMStart)/double(ti.maxT);
			};
		//----------------------------------------------------------------------------//

			SetFertility();
			RunNReps(ti.NumRuns);
	return 0;};
	

	void RunNReps(int N)
	      {
		for(int i=0;i<N;i++)
		{
			clearAd();
			clearJ();
			RunMaxT();
		};
		return;};

	void RunMaxT()
		{
		int TT=0;
		int swit=0;
		int dsx,wt,dsxMale;
		while (TT<ti.maxT+1)
			{
		if(TT==0 ||TT==3||TT==7||TT==11)to.Jww[TL-2]+=in.startNum;//initiate with wildtype pupae
		//---------------------introduce gene-drive mosquitoes ----------------------------//
			if(TT==in.driver_time)
					{
					int num=int((to.MwwT+to.FwwT)*pa.numDsx*0.5);
					to.Mwd[0]+=num; to.MTot+=num;
					UpdateMate();
					};
			if(TT==in.driver_time+3)
					{
					int num=int((to.MwwT+to.FwwT)*pa.numDsx*0.5);
					to.Mwd[0]+=num; to.MTot+=num;
					UpdateMate();
					};
		//--------------------------------------------------------------------------------//
	
			OneStep(TT);
		ComputeTotals();
		dsx= to.Mwd[0]+to.Mdd[0]+to.Mdr[0]+ to.fVwd[0]+ to.mVwd[0]+ to.bVwd[0]+ to.Vdd[0]+to.Vdr[0];
		wt= to.Mww[0]+to.Mwr[0]+to.Mrr[0]+ to.Vww[0]+to.fVww[0]+ to.mVww[0]+ to.bVww[0]+ +to.Vwr[0]+to.Vrr[0];
		dsxMale=to.Mwd[0]+to.Mdd[0]+to.Mdr[0]+ to.Vdd[0]+to.Vdr[0]; 
		if(TT%7==0 || TT%7==4)cout<<TT-in.driver_time<<"   "<<wt<<"   "<<dsx<<"    "<<dsxMale<<"   "<<to.neggs<<"    "<<to.keepeggs<<"    "<<TT%7<<"    "<<to.r2<<"   "<<to.Mww[0]+to.Mwr[0]+to.Mrr[0]<<"   "<<to.Vww[0]+to.fVww[0]+ to.mVww[0]+ to.bVww[0]+ +to.Vwr[0]+to.Vrr[0]<<"    "<<to.fVwd[0]+ to.mVwd[0]+ to.bVwd[0]<<"     "<<to.Mwd[0]+to.Mdd[0]+to.Mdr[0]<<"      "<<to.Vdd[0]+to.Vdr[0]<<endl;
			TT++;
			};
        return;};
void ComputeTotals(void){

	to.MwwT=std::accumulate(to.Mww,to.Mww+mA,0); to.MwdT=std::accumulate(to.Mwd,to.Mwd+mA,0); to.MddT=std::accumulate(to.Mdd,to.Mdd+mA,0); to.MwrT=std::accumulate(to.Mwr,to.Mwr+mA,0); to.MrrT=std::accumulate(to.Mrr,to.Mrr+mA,0); to.MdrT=std::accumulate(to.Mdr,to.Mdr+mA,0);
	to.VwwT=std::accumulate(to.Vww,to.Vww+mA,0);
	to.fVwwT=std::accumulate(to.fVww,to.fVww+mA,0);
	to.mVwwT=std::accumulate(to.mVww,to.mVww+mA,0);
	to.bVwwT=std::accumulate(to.bVww,to.bVww+mA,0);
	to.fVwdT=std::accumulate(to.fVwd,to.fVwd+mA,0);
	to.mVwdT=std::accumulate(to.mVwd,to.mVwd+mA,0);
	to.bVwdT=std::accumulate(to.bVwd,to.bVwd+mA,0);
	to.VddT=std::accumulate(to.Vdd,to.Vdd+mA,0);
	to.VwrT=std::accumulate(to.Vwr,to.Vwr+mA,0);
	to.fVwrT=std::accumulate(to.fVwr,to.fVwr+mA,0);
	to.mVwrT=std::accumulate(to.mVwr,to.mVwr+mA,0);
	to.bVwrT=std::accumulate(to.bVwr,to.bVwr+mA,0);
	to.VrrT=std::accumulate(to.Vrr,to.Vrr+mA,0);
	to.VdrT=std::accumulate(to.Vdr,to.Vdr+mA,0);
		       
		
	to.FwwT= std::accumulate(to.Fwwww,to.Fwwww+mA,0) +std::accumulate(to.Fwwwd,to.Fwwwd+mA,0) +std::accumulate(to.Fwwdd,to.Fwwdd+mA,0) +std::accumulate(to.Fwwwr,to.Fwwwr+mA,0)+std::accumulate(to.Fwwrr,to.Fwwrr+mA,0)+std::accumulate(to.Fwwdr,to.Fwwdr+mA,0);
	to.fFwwT= std::accumulate(to.fFwwww,to.fFwwww+mA,0) +std::accumulate(to.fFwwwd,to.fFwwwd+mA,0) +std::accumulate(to.fFwwdd,to.fFwwdd+mA,0) +std::accumulate(to.fFwwwr,to.fFwwwr+mA,0)+std::accumulate(to.fFwwrr,to.fFwwrr+mA,0)+std::accumulate(to.fFwwdr,to.fFwwdr+mA,0);
	to.mFwwT= std::accumulate(to.mFwwww,to.mFwwww+mA,0) +std::accumulate(to.mFwwwd,to.mFwwwd+mA,0) +std::accumulate(to.mFwwdd,to.mFwwdd+mA,0) +std::accumulate(to.mFwwwr,to.mFwwwr+mA,0)+std::accumulate(to.mFwwrr,to.mFwwrr+mA,0)+std::accumulate(to.mFwwdr,to.mFwwdr+mA,0);
	to.bFwwT= std::accumulate(to.bFwwww,to.bFwwww+mA,0) +std::accumulate(to.bFwwwd,to.bFwwwd+mA,0) +std::accumulate(to.bFwwdd,to.bFwwdd+mA,0) +std::accumulate(to.bFwwwr,to.bFwwwr+mA,0)+std::accumulate(to.bFwwrr,to.bFwwrr+mA,0)+std::accumulate(to.bFwwdr,to.bFwwdr+mA,0);
	
	to.fFwdT=std::accumulate(to.fFwdww,to.fFwdww+mA,0)+std::accumulate(to.fFwdwd,to.fFwdwd+mA,0)+std::accumulate(to.fFwddd,to.fFwddd+mA,0)+std::accumulate(to.fFwdwr,to.fFwdwr+mA,0)+std::accumulate(to.fFwdrr,to.fFwdrr+mA,0)+std::accumulate(to.fFwddr,to.fFwddr+mA,0);
	to.mFwdT= std::accumulate(to.mFwdww,to.mFwdww+mA,0) +std::accumulate(to.mFwdwd,to.mFwdwd+mA,0) +std::accumulate(to.mFwddd,to.mFwddd+mA,0) +std::accumulate(to.mFwdwr,to.mFwdwr+mA,0)+std::accumulate(to.mFwdrr,to.mFwdrr+mA,0)+std::accumulate(to.mFwddr,to.mFwddr+mA,0);
	to.bFwdT= std::accumulate(to.bFwdww,to.bFwdww+mA,0) +std::accumulate(to.bFwdwd,to.bFwdwd+mA,0) +std::accumulate(to.bFwddd,to.bFwddd+mA,0) +std::accumulate(to.bFwdwr,to.bFwdwr+mA,0)+std::accumulate(to.bFwdrr,to.bFwdrr+mA,0)+std::accumulate(to.bFwddr,to.bFwddr+mA,0);
	to.FddT= std::accumulate(to.Fddww,to.Fddww+mA,0) +std::accumulate(to.Fddwd,to.Fddwd+mA,0) +std::accumulate(to.Fdddd,to.Fdddd+mA,0) +std::accumulate(to.Fddwr,to.Fddwr+mA,0)+std::accumulate(to.Fddrr,to.Fddrr+mA,0)+std::accumulate(to.Fdddr,to.Fdddr+mA,0);
	to.FwrT= std::accumulate(to.Fwrww,to.Fwrww+mA,0) +std::accumulate(to.Fwrwd,to.Fwrwd+mA,0) +std::accumulate(to.Fwrdd,to.Fwrdd+mA,0) +std::accumulate(to.Fwrwr,to.Fwrwr+mA,0)+std::accumulate(to.Fwrrr,to.Fwrrr+mA,0)+std::accumulate(to.Fwrdr,to.Fwrdr+mA,0);
	to.fFwrT= std::accumulate(to.fFwrww,to.fFwrww+mA,0) +std::accumulate(to.fFwrwd,to.fFwrwd+mA,0) +std::accumulate(to.fFwrdd,to.fFwrdd+mA,0) +std::accumulate(to.fFwrwr,to.fFwrwr+mA,0)+std::accumulate(to.fFwrrr,to.fFwrrr+mA,0)+std::accumulate(to.fFwrdr,to.fFwrdr+mA,0);
	to.mFwrT= std::accumulate(to.mFwrww,to.mFwrww+mA,0) +std::accumulate(to.mFwrwd,to.mFwrwd+mA,0) +std::accumulate(to.mFwrdd,to.mFwrdd+mA,0) +std::accumulate(to.mFwrwr,to.mFwrwr+mA,0)+std::accumulate(to.mFwrrr,to.mFwrrr+mA,0)+std::accumulate(to.mFwrdr,to.mFwrdr+mA,0);
	to.bFwrT= std::accumulate(to.bFwrww,to.bFwrww+mA,0) +std::accumulate(to.bFwrwd,to.bFwrwd+mA,0) +std::accumulate(to.bFwrdd,to.bFwrdd+mA,0) +std::accumulate(to.bFwrwr,to.bFwrwr+mA,0)+std::accumulate(to.bFwrrr,to.bFwrrr+mA,0)+std::accumulate(to.bFwrdr,to.bFwrdr+mA,0);
	to.FrrT= std::accumulate(to.Frrww,to.Frrww+mA,0) +std::accumulate(to.Frrwd,to.Frrwd+mA,0) +std::accumulate(to.Frrdd,to.Frrdd+mA,0) +std::accumulate(to.Frrwr,to.Frrwr+mA,0)+std::accumulate(to.Frrrr,to.Frrrr+mA,0)+std::accumulate(to.Frrdr,to.Frrdr+mA,0);
	to.FdrT= std::accumulate(to.Fdrww,to.Fdrww+mA,0) +std::accumulate(to.Fdrwd,to.Fdrwd+mA,0) +std::accumulate(to.Fdrdd,to.Fdrdd+mA,0) +std::accumulate(to.Fdrwr,to.Fdrwr+mA,0)+std::accumulate(to.Fdrrr,to.Fdrrr+mA,0)+std::accumulate(to.Fdrdr,to.Fdrdr+mA,0);
return;};			





	void clearJ(void){
		to.JTot=0;
		for(int i=0;i<TL;i++){
		to.Jww[i]=0; to.fJww[i]=0; to.mJww[i]=0; to.bJww[i]=0; to.fJwd[i]=0; to.mJwd[i]=0; to.bJwd[i]=0; to.Jdd[i]=0; to.Jwr[i]=0; to.fJwr[i]=0; to.mJwr[i]=0; to.bJwr[i]=0; to.Jrr[i]=0; to.Jdr[i]=0;
		to.neggs=0;to.keepeggs=0;to.r2=0;

					};
	return;};



	void clearAd(void){
		to.MTot=0;to.VTot=0;to.FTot=0;
		to.MwwT=0;to.MwdT=0;to.MddT=0;to.MwrT=0;to.MrrT=0;to.MdrT=0;
		to.VwwT=0; to.fVwwT=0; to.mVwwT=0; to.bVwwT=0; to.fVwdT=0; to.mVwdT=0; to.bVwdT=0; to.VddT=0; to.VwrT=0; to.fVwrT=0; to.mVwrT=0; to.bVwrT=0; to.VrrT=0; to.VdrT=0;

		to.FwwT=0; to.fFwwT=0; to.mFwwT=0; to.bFwwT=0; to.fFwdT=0; to.mFwdT=0; to.bFwdT=0; to.FddT=0; to.FwrT=0; to.fFwrT=0; to.mFwrT=0; to.bFwrT=0; to.FrrT=0; to.FdrT=0;

		for(int i=0;i<mA;i++){
			to.Mww[i]=0; to.Mwd[i]=0; to.Mdd[i]=0; to.Mwr[i]=0; to.Mrr[i]=0; to.Mdr[i]=0;
			to.Vww[i]=0; to.fVww[i]=0; to.mVww[i]=0; to.bVww[i]=0; to.fVwd[i]=0; to.mVwd[i]=0; to.bVwd[i]=0; to.Vdd[i]=0; to.Vwr[i]=0; to.fVwr[i]=0; to.mVwr[i]=0; to.bVwr[i]=0; to.Vrr[i]=0; to.Vdr[i]=0;

			to.Fwwww[i]=0; to.Fwwwd[i]=0; to.Fwwdd[i]=0; to.Fwwwr[i]=0; to.Fwwrr[i]=0; to.Fwwdr[i]=0;
			to.fFwwww[i]=0; to.fFwwwd[i]=0; to.fFwwdd[i]=0; to.fFwwwr[i]=0; to.fFwwrr[i]=0; to.fFwwdr[i]=0;
			to.mFwwww[i]=0; to.mFwwwd[i]=0; to.mFwwdd[i]=0; to.mFwwwr[i]=0; to.mFwwrr[i]=0; to.mFwwdr[i]=0;
			to.bFwwww[i]=0; to.bFwwwd[i]=0; to.bFwwdd[i]=0; to.bFwwwr[i]=0; to.bFwwrr[i]=0; to.bFwwdr[i]=0;
			to.fFwdww[i]=0; to.fFwdwd[i]=0; to.fFwddd[i]=0; to.fFwdwr[i]=0; to.fFwdrr[i]=0; to.fFwddr[i]=0;
			to.mFwdww[i]=0; to.mFwdwd[i]=0; to.mFwddd[i]=0; to.mFwdwr[i]=0; to.mFwdrr[i]=0; to.mFwddr[i]=0;
			to.bFwdww[i]=0; to.bFwdwd[i]=0; to.bFwddd[i]=0; to.bFwdwr[i]=0; to.bFwdrr[i]=0; to.bFwddr[i]=0;
			to.Fddww[i]=0; to.Fddwd[i]=0; to.Fdddd[i]=0; to.Fddwr[i]=0; to.Fddrr[i]=0; to.Fdddr[i]=0;
			to.Fwrww[i]=0; to.Fwrwd[i]=0; to.Fwrdd[i]=0; to.Fwrwr[i]=0; to.Fwrrr[i]=0; to.Fwrdr[i]=0;
			to.fFwrww[i]=0; to.fFwrwd[i]=0; to.fFwrdd[i]=0; to.fFwrwr[i]=0; to.fFwrrr[i]=0; to.fFwrdr[i]=0;
			to.mFwrww[i]=0; to.mFwrwd[i]=0; to.mFwrdd[i]=0; to.mFwrwr[i]=0; to.mFwrrr[i]=0; to.mFwrdr[i]=0;
			to.bFwrww[i]=0; to.bFwrwd[i]=0; to.bFwrdd[i]=0; to.bFwrwr[i]=0; to.bFwrrr[i]=0; to.bFwrdr[i]=0;
			to.Frrww[i]=0; to.Frrwd[i]=0; to.Frrdd[i]=0; to.Frrwr[i]=0; to.Frrrr[i]=0; to.Frrdr[i]=0;
			to.Fdrww[i]=0; to.Fdrwd[i]=0; to.Fdrdd[i]=0; to.Fdrwr[i]=0; to.Fdrrr[i]=0; to.Fdrdr[i]=0;
					};
	return;};
	
	void initiate(void){
		to.Jww[TL-1]+=in.startNum;
		return;};



	void OneStep(int day){
		JuvGetOlder();
		VirginsMate();
		if(day%7 ==5||day%7==2)FemaleCount();
		if(day%7 ==1||day%7==4)LayEggs();
		AdultsDie(day);
		JuvEmerge();
		UpdateMate();
	return;};
	

	void JuvEmerge(void){
		int surv,survM;
		surv=to.Jww[TL-1]; 
		if(surv>0)
		{		
			survM=random_binomial(surv,0.5);
			to.Mww[0]+=survM;
		       	to.MTot+=survM;
			to.Vww[0]+=surv-survM;
			to.VTot+=surv-survM;
		};
		surv=to.fJww[TL-1]; 
		if(surv>0)
		{		
			survM=random_binomial(surv,0.5);
			to.Mww[0]+=survM;
		       	to.MTot+=survM;
			to.fVww[0]+=surv-survM;
			to.VTot+=surv-survM;
		};
		surv=to.mJww[TL-1]; 
		if(surv>0)
		{		
			survM=random_binomial(surv,0.5);
			to.Mww[0]+=survM;
		       	to.MTot+=survM;
			to.mVww[0]+=surv-survM;
			to.VTot+=surv-survM;
		};
		surv=to.bJww[TL-1]; 
		if(surv>0)
		{		
			survM=random_binomial(surv,0.5);
			to.Mww[0]+=survM;
		       	to.MTot+=survM;
			to.bVww[0]+=surv-survM;
			to.VTot+=surv-survM;
		};
		surv=to.fJwd[TL-1]; 
		if(surv>0)
		{		
			survM=random_binomial(surv,0.5);
			to.Mwd[0]+=survM;
		       	to.MTot+=survM;
			to.fVwd[0]+=surv-survM;
			to.VTot+=surv-survM;
		};

		surv=to.mJwd[TL-1]; 
		if(surv>0)
		{		
			survM=random_binomial(surv,0.5);
			to.Mwd[0]+=survM;
		       	to.MTot+=survM;
			to.mVwd[0]+=surv-survM;
			to.VTot+=surv-survM;
		};
		surv=to.bJwd[TL-1]; 
		if(surv>0)
		{		
			survM=random_binomial(surv,0.5);
			to.Mwd[0]+=survM;
		       	to.MTot+=survM;
			to.bVwd[0]+=surv-survM;
			to.VTot+=surv-survM;
		};
		surv=to.Jdd[TL-1]; 
		if(surv>0)
		{		
			survM=random_binomial(surv,0.5);
			to.Mdd[0]+=survM;
		       	to.MTot+=survM;
			to.Vdd[0]+=surv-survM;
			to.VTot+=surv-survM;
		};
		surv=to.Jwr[TL-1]; 
		if(surv>0)
		{		
			survM=random_binomial(surv,0.5);
			to.Mwr[0]+=survM;
		       	to.MTot+=survM;
			to.Vwr[0]+=surv-survM;
			to.VTot+=surv-survM;
		};
		surv=to.fJwr[TL-1]; 
		if(surv>0)
		{		
			survM=random_binomial(surv,0.5);
			to.Mwr[0]+=survM;
		       	to.MTot+=survM;
			to.fVwr[0]+=surv-survM;
			to.VTot+=surv-survM;
		};
		surv=to.mJwr[TL-1]; 
		if(surv>0)
		{		
			survM=random_binomial(surv,0.5);
			to.Mwr[0]+=survM;
		       	to.MTot+=survM;
			to.mVwr[0]+=surv-survM;
			to.VTot+=surv-survM;
		};
		surv=to.bJwr[TL-1]; 
		if(surv>0)
		{		
			survM=random_binomial(surv,0.5);
			to.Mwr[0]+=survM;
		       	to.MTot+=survM;
			to.bVwr[0]+=surv-survM;
			to.VTot+=surv-survM;
		};
		surv=to.Jrr[TL-1]; 
		if(surv>0)
		{		
			survM=random_binomial(surv,0.5);
			to.Mrr[0]+=survM;
		       	to.MTot+=survM;
			to.Vrr[0]+=surv-survM;
			to.VTot+=surv-survM;
		};
		surv=to.Jdr[TL-1]; 
		if(surv>0)
		{		
			survM=random_binomial(surv,0.5);
			to.Mdr[0]+=survM;
		       	to.MTot+=survM;
			to.Vdr[0]+=surv-survM;
			to.VTot+=surv-survM;
		};
		return;};

	void JuvGetOlder(void){
			for(int age=TL-1;age>0;age--)
			{
				to.Jww[age]=to.Jww[age-1];
				to.fJww[age]=to.fJww[age-1];
				to.mJww[age]=to.mJww[age-1];
				to.bJww[age]=to.bJww[age-1];
				to.fJwd[age]=to.fJwd[age-1];
				to.mJwd[age]=to.mJwd[age-1];
				to.bJwd[age]=to.bJwd[age-1];
				to.Jdd[age]=to.Jdd[age-1];
				to.Jwr[age]=to.Jwr[age-1];
				to.fJwr[age]=to.fJwr[age-1];
				to.mJwr[age]=to.mJwr[age-1];
				to.bJwr[age]=to.bJwr[age-1];
				to.Jrr[age]=to.Jrr[age-1];
				to.Jdr[age]=to.Jdr[age-1];
			};
			to.Jww[0]=0; to.fJww[0]=0; to.mJww[0]=0; to.bJww[0]=0; to.fJwd[0]=0; to.mJwd[0]=0; to.bJwd[0]=0; to.Jdd[0]=0; to.Jwr[0]=0; to.fJwr[0]=0; to.mJwr[0]=0; to.bJwr[0]=0; to.Jrr[0]=0; to.Jdr[0]=0;
			for(int age=mA-1;age>0;age--)
			{
			to.Mww[age]=to.Mww[age-1]; to.Mwd[age]=to.Mwd[age-1]; to.Mdd[age]=to.Mdd[age-1]; to.Mwr[age]=to.Mwr[age-1]; to.Mrr[age]=to.Mrr[age-1]; to.Mdr[age]=to.Mdr[age-1]; 
			
			to.Vww[age]=to.Vww[age-1]; 
			to.fVww[age]=to.fVww[age-1]; 
			to.mVww[age]=to.mVww[age-1]; 
			to.bVww[age]=to.bVww[age-1]; 
			to.fVwd[age]=to.fVwd[age-1]; 
			to.mVwd[age]=to.mVwd[age-1]; 
			to.bVwd[age]=to.bVwd[age-1]; 
			to.Vdd[age]=to.Vdd[age-1]; 
			to.Vwr[age]=to.Vwr[age-1]; 
			to.fVwr[age]=to.fVwr[age-1]; 
			to.mVwr[age]=to.mVwr[age-1]; 
			to.bVwr[age]=to.bVwr[age-1]; 
			to.Vrr[age]=to.Vrr[age-1]; 
			to.Vdr[age]=to.Vdr[age-1]; 
			
			to.Fwwww[age]=to.Fwwww[age-1]; to.Fwwwd[age]=to.Fwwwd[age-1]; to.Fwwdd[age]=to.Fwwdd[age-1]; to.Fwwwr[age]=to.Fwwwr[age-1]; to.Fwwrr[age]=to.Fwwrr[age-1]; to.Fwwdr[age]=to.Fwwdr[age-1]; 
			to.fFwwww[age]=to.fFwwww[age-1]; to.fFwwwd[age]=to.fFwwwd[age-1]; to.fFwwdd[age]=to.fFwwdd[age-1]; to.fFwwwr[age]=to.fFwwwr[age-1]; to.fFwwrr[age]=to.fFwwrr[age-1]; to.fFwwdr[age]=to.fFwwdr[age-1]; 
			to.mFwwww[age]=to.mFwwww[age-1]; to.mFwwwd[age]=to.mFwwwd[age-1]; to.mFwwdd[age]=to.mFwwdd[age-1]; to.mFwwwr[age]=to.mFwwwr[age-1]; to.mFwwrr[age]=to.mFwwrr[age-1]; to.mFwwdr[age]=to.mFwwdr[age-1]; 
			to.bFwwww[age]=to.bFwwww[age-1]; to.bFwwwd[age]=to.bFwwwd[age-1]; to.bFwwdd[age]=to.bFwwdd[age-1]; to.bFwwwr[age]=to.bFwwwr[age-1]; to.bFwwrr[age]=to.bFwwrr[age-1]; to.bFwwdr[age]=to.bFwwdr[age-1]; 
			to.fFwdww[age]=to.fFwdww[age-1]; to.fFwdwd[age]=to.fFwdwd[age-1]; to.fFwddd[age]=to.fFwddd[age-1]; to.fFwdwr[age]=to.fFwdwr[age-1]; to.fFwdrr[age]=to.fFwdrr[age-1]; to.fFwddr[age]=to.fFwddr[age-1]; 
			to.mFwdww[age]=to.mFwdww[age-1]; to.mFwdwd[age]=to.mFwdwd[age-1]; to.mFwddd[age]=to.mFwddd[age-1]; to.mFwdwr[age]=to.mFwdwr[age-1]; to.mFwdrr[age]=to.mFwdrr[age-1]; to.mFwddr[age]=to.mFwddr[age-1]; 
			to.bFwdww[age]=to.bFwdww[age-1]; to.bFwdwd[age]=to.bFwdwd[age-1]; to.bFwddd[age]=to.bFwddd[age-1]; to.bFwdwr[age]=to.bFwdwr[age-1]; to.bFwdrr[age]=to.bFwdrr[age-1]; to.bFwddr[age]=to.bFwddr[age-1]; 
			to.Fddww[age]=to.Fddww[age-1]; to.Fddwd[age]=to.Fddwd[age-1]; to.Fdddd[age]=to.Fdddd[age-1]; to.Fddwr[age]=to.Fddwr[age-1]; to.Fddrr[age]=to.Fddrr[age-1]; to.Fdddr[age]=to.Fdddr[age-1]; 
			to.Fwrww[age]=to.Fwrww[age-1]; to.Fwrwd[age]=to.Fwrwd[age-1]; to.Fwrdd[age]=to.Fwrdd[age-1]; to.Fwrwr[age]=to.Fwrwr[age-1]; to.Fwrrr[age]=to.Fwrrr[age-1]; to.Fwrdr[age]=to.Fwrdr[age-1]; 
			to.fFwrww[age]=to.fFwrww[age-1]; to.fFwrwd[age]=to.fFwrwd[age-1]; to.fFwrdd[age]=to.fFwrdd[age-1]; to.fFwrwr[age]=to.fFwrwr[age-1]; to.fFwrrr[age]=to.fFwrrr[age-1]; to.fFwrdr[age]=to.fFwrdr[age-1]; 
			to.mFwrww[age]=to.mFwrww[age-1]; to.mFwrwd[age]=to.mFwrwd[age-1]; to.mFwrdd[age]=to.mFwrdd[age-1]; to.mFwrwr[age]=to.mFwrwr[age-1]; to.mFwrrr[age]=to.mFwrrr[age-1]; to.mFwrdr[age]=to.mFwrdr[age-1]; 
			to.bFwrww[age]=to.bFwrww[age-1]; to.bFwrwd[age]=to.bFwrwd[age-1]; to.bFwrdd[age]=to.bFwrdd[age-1]; to.bFwrwr[age]=to.bFwrwr[age-1]; to.bFwrrr[age]=to.bFwrrr[age-1]; to.bFwrdr[age]=to.bFwrdr[age-1]; 
			to.Frrww[age]=to.Frrww[age-1]; to.Frrwd[age]=to.Frrwd[age-1]; to.Frrdd[age]=to.Frrdd[age-1]; to.Frrwr[age]=to.Frrwr[age-1]; to.Frrrr[age]=to.Frrrr[age-1]; to.Frrdr[age]=to.Frrdr[age-1]; 
			to.Fdrww[age]=to.Fdrww[age-1]; to.Fdrwd[age]=to.Fdrwd[age-1]; to.Fdrdd[age]=to.Fdrdd[age-1]; to.Fdrwr[age]=to.Fdrwr[age-1]; to.Fdrrr[age]=to.Fdrrr[age-1]; to.Fdrdr[age]=to.Fdrdr[age-1]; 
			
			
			};
			to.Mww[0]=0; to.Mwd[0]=0; to.Mdd[0]=0; to.Mwr[0]=0; to.Mrr[0]=0; to.Mdr[0]=0;
			to.Vww[0]=0; to.fVww[0]=0; to.mVww[0]=0; to.bVww[0]=0; to.fVwd[0]=0; to.mVwd[0]=0; to.bVwd[0]=0; to.Vdd[0]=0; to.Vwr[0]=0; to.fVwr[0]=0; to.mVwr[0]=0; to.bVwr[0]=0; to.Vrr[0]=0; to.Vdr[0]=0;

			to.Fwwww[0]=0; to.Fwwwd[0]=0; to.Fwwdd[0]=0; to.Fwwwr[0]=0; to.Fwwrr[0]=0; to.Fwwdr[0]=0;
			to.fFwwww[0]=0; to.fFwwwd[0]=0; to.fFwwdd[0]=0; to.fFwwwr[0]=0; to.fFwwrr[0]=0; to.fFwwdr[0]=0;
			to.mFwwww[0]=0; to.mFwwwd[0]=0; to.mFwwdd[0]=0; to.mFwwwr[0]=0; to.mFwwrr[0]=0; to.mFwwdr[0]=0;
			to.bFwwww[0]=0; to.bFwwwd[0]=0; to.bFwwdd[0]=0; to.bFwwwr[0]=0; to.bFwwrr[0]=0; to.bFwwdr[0]=0;
			to.fFwdww[0]=0; to.fFwdwd[0]=0; to.fFwddd[0]=0; to.fFwdwr[0]=0; to.fFwdrr[0]=0; to.fFwddr[0]=0;
			to.mFwdww[0]=0; to.mFwdwd[0]=0; to.mFwddd[0]=0; to.mFwdwr[0]=0; to.mFwdrr[0]=0; to.mFwddr[0]=0;
			to.bFwdww[0]=0; to.bFwdwd[0]=0; to.bFwddd[0]=0; to.bFwdwr[0]=0; to.bFwdrr[0]=0; to.bFwddr[0]=0;
			to.Fddww[0]=0; to.Fddwd[0]=0; to.Fdddd[0]=0; to.Fddwr[0]=0; to.Fddrr[0]=0; to.Fdddr[0]=0;
			to.Fwrww[0]=0; to.Fwrwd[0]=0; to.Fwrdd[0]=0; to.Fwrwr[0]=0; to.Fwrrr[0]=0; to.Fwrdr[0]=0;
			to.fFwrww[0]=0; to.fFwrwd[0]=0; to.fFwrdd[0]=0; to.fFwrwr[0]=0; to.fFwrrr[0]=0; to.fFwrdr[0]=0;
			to.mFwrww[0]=0; to.mFwrwd[0]=0; to.mFwrdd[0]=0; to.mFwrwr[0]=0; to.mFwrrr[0]=0; to.mFwrdr[0]=0;
			to.bFwrww[0]=0; to.bFwrwd[0]=0; to.bFwrdd[0]=0; to.bFwrwr[0]=0; to.bFwrrr[0]=0; to.bFwrdr[0]=0;
			to.Frrww[0]=0; to.Frrwd[0]=0; to.Frrdd[0]=0; to.Frrwr[0]=0; to.Frrrr[0]=0; to.Frrdr[0]=0;
			to.Fdrww[0]=0; to.Fdrwd[0]=0; to.Fdrdd[0]=0; to.Fdrwr[0]=0; to.Fdrrr[0]=0; to.Fdrdr[0]=0;

	return;};

	void VirginsMate(){
		int vww=0;	int vwd=0;	int vdd=0;	int vwr=0;	int vrr=0;	int vdr=0;	
		int fww=0;	int fwd=0;	int fdd=0;	int fwr=0;	int frr=0;	int fdr=0;
int fvww=0; int mvww=0; int bvww=0; int fvwd=0; int mvwd=0; int bvwd=0; int fvwr=0; int mvwr=0; int bvwr=0;

		int* mates;
		int probs[6]={to.MwwT,to.MwdT,to.MddT,to.MwrT,to.MrrT,to.MdrT};
	
		for(int age=0;age<mA;age++)
		{
		vww=random_binomial(to.Vww[age],to.mate_rate);
		fvww=random_binomial(to.fVww[age],to.mate_rate);
		mvww=random_binomial(to.mVww[age],to.mate_rate);
		bvww=random_binomial(to.bVww[age],to.mate_rate);
		fvwd=random_binomial(to.fVwd[age],to.mate_rate);
		mvwd=random_binomial(to.mVwd[age],to.mate_rate);
		bvwd=random_binomial(to.bVwd[age],to.mate_rate);
		vdd=random_binomial(to.Vdd[age],to.mate_rate);
		vwr=random_binomial(to.Vwr[age],to.mate_rate);
		fvwr=random_binomial(to.fVwr[age],to.mate_rate);
		mvwr=random_binomial(to.mVwr[age],to.mate_rate);
		bvwr=random_binomial(to.bVwr[age],to.mate_rate);
		vrr=random_binomial(to.Vrr[age],to.mate_rate);
		vdr=random_binomial(to.Vdr[age],to.mate_rate);
		if(vww>0)
		{
		mates=random_multinom(vww,probs);
		to.Fwwww[age]+=*mates; to.Fwwwd[age]+=*(mates+1);to.Fwwdd[age]+=*(mates+2);to.Fwwwr[age]+=*(mates+3);to.Fwwrr[age]+=*(mates+4);to.Fwwdr[age]+=*(mates+5);
		delete[] mates;
		to.Vww[age]-=vww; to.VTot-=vww; to.FTot+=vww;
		};
		if(fvww>0)
		{
		mates=random_multinom(fvww,probs);
		to.fFwwww[age]+=*mates; to.fFwwwd[age]+=*(mates+1);to.fFwwdd[age]+=*(mates+2);to.fFwwwr[age]+=*(mates+3);to.fFwwrr[age]+=*(mates+4);to.fFwwdr[age]+=*(mates+5);
		delete[] mates;
		to.fVww[age]-=fvww; to.VTot-=fvww; to.FTot+=fvww;
		};
		if(mvww>0)
		{
		mates=random_multinom(mvww,probs);
		to.mFwwww[age]+=*mates; to.mFwwwd[age]+=*(mates+1);to.mFwwdd[age]+=*(mates+2);to.mFwwwr[age]+=*(mates+3);to.mFwwrr[age]+=*(mates+4);to.mFwwdr[age]+=*(mates+5);
		delete[] mates;
		to.mVww[age]-=mvww; to.VTot-=mvww; to.FTot+=mvww;
		};
		if(bvww>0)
		{
		mates=random_multinom(bvww,probs);
		to.bFwwww[age]+=*mates; to.bFwwwd[age]+=*(mates+1);to.bFwwdd[age]+=*(mates+2);to.bFwwwr[age]+=*(mates+3);to.bFwwrr[age]+=*(mates+4);to.bFwwdr[age]+=*(mates+5);
		delete[] mates;
		to.bVww[age]-=bvww; to.VTot-=bvww; to.FTot+=bvww;
		};
		if(fvwd>0)
		{
		mates=random_multinom(fvwd,probs);
		to.fFwdww[age]+=*mates; to.fFwdwd[age]+=*(mates+1);to.fFwddd[age]+=*(mates+2);to.fFwdwr[age]+=*(mates+3);to.fFwdrr[age]+=*(mates+4);to.fFwddr[age]+=*(mates+5);
		delete[] mates;
		to.fVwd[age]-=fvwd; to.VTot-=fvwd; to.FTot+=fvwd;
		};
		if(mvwd>0)
		{
		mates=random_multinom(mvwd,probs);
		to.mFwdww[age]+=*mates; to.mFwdwd[age]+=*(mates+1);to.mFwddd[age]+=*(mates+2);to.mFwdwr[age]+=*(mates+3);to.mFwdrr[age]+=*(mates+4);to.mFwddr[age]+=*(mates+5);
		delete[] mates;
		to.mVwd[age]-=mvwd; to.VTot-=mvwd; to.FTot+=mvwd;
		};
		if(bvwd>0)
		{
		mates=random_multinom(bvwd,probs);
		to.bFwdww[age]+=*mates; to.bFwdwd[age]+=*(mates+1);to.bFwddd[age]+=*(mates+2);to.bFwdwr[age]+=*(mates+3);to.bFwdrr[age]+=*(mates+4);to.bFwddr[age]+=*(mates+5);
		delete[] mates;
		to.bVwd[age]-=bvwd; to.VTot-=bvwd; to.FTot+=bvwd;
		};
		if(vdd>0)
		{
		mates=random_multinom(vdd,probs);
		to.Fddww[age]+=*mates; to.Fddwd[age]+=*(mates+1);to.Fdddd[age]+=*(mates+2);to.Fddwr[age]+=*(mates+3);to.Fddrr[age]+=*(mates+4);to.Fdddr[age]+=*(mates+5);
		delete[] mates;
		to.Vdd[age]-=vdd; to.VTot-=vdd; to.FTot+=vdd;
		};
		if(vwr>0)
		{
		mates=random_multinom(vwr,probs);
		to.Fwrww[age]+=*mates; to.Fwrwd[age]+=*(mates+1);to.Fwrdd[age]+=*(mates+2);to.Fwrwr[age]+=*(mates+3);to.Fwrrr[age]+=*(mates+4);to.Fwrdr[age]+=*(mates+5);
		delete[] mates;
		to.Vwr[age]-=vwr; to.VTot-=vwr; to.FTot+=vwr;
		};
		if(fvwr>0)
		{
		mates=random_multinom(fvwr,probs);
		to.fFwrww[age]+=*mates; to.fFwrwd[age]+=*(mates+1);to.fFwrdd[age]+=*(mates+2);to.fFwrwr[age]+=*(mates+3);to.fFwrrr[age]+=*(mates+4);to.fFwrdr[age]+=*(mates+5);
		delete[] mates;
		to.fVwr[age]-=fvwr; to.VTot-=fvwr; to.FTot+=fvwr;
		};
		if(mvwr>0)
		{
		mates=random_multinom(mvwr,probs);
		to.mFwrww[age]+=*mates; to.mFwrwd[age]+=*(mates+1);to.mFwrdd[age]+=*(mates+2);to.mFwrwr[age]+=*(mates+3);to.mFwrrr[age]+=*(mates+4);to.mFwrdr[age]+=*(mates+5);
		delete[] mates;
		to.mVwr[age]-=mvwr; to.VTot-=mvwr; to.FTot+=mvwr;
		};
		if(bvwr>0)
		{
		mates=random_multinom(bvwr,probs);
		to.bFwrww[age]+=*mates; to.bFwrwd[age]+=*(mates+1);to.bFwrdd[age]+=*(mates+2);to.bFwrwr[age]+=*(mates+3);to.bFwrrr[age]+=*(mates+4);to.bFwrdr[age]+=*(mates+5);
		delete[] mates;
		to.bVwr[age]-=bvwr; to.VTot-=bvwr; to.FTot+=bvwr;
		};
		if(vrr>0)
		{
		mates=random_multinom(vrr,probs);
		to.Frrww[age]+=*mates; to.Frrwd[age]+=*(mates+1);to.Frrdd[age]+=*(mates+2);to.Frrwr[age]+=*(mates+3);to.Frrrr[age]+=*(mates+4);to.Frrdr[age]+=*(mates+5);
		delete[] mates;
		to.Vrr[age]-=vrr; to.VTot-=vrr; to.FTot+=vrr;
		};
		if(vdr>0)
		{
		mates=random_multinom(vdr,probs);
		to.Fdrww[age]+=*mates; to.Fdrwd[age]+=*(mates+1);to.Fdrdd[age]+=*(mates+2);to.Fdrwr[age]+=*(mates+3);to.Fdrrr[age]+=*(mates+4);to.Fdrdr[age]+=*(mates+5);
		delete[] mates;
		to.Vdr[age]-=vdr; to.VTot-=vdr; to.FTot+=vdr;
		};
		};
	return;};

	void FemaleCount (){
	to.toFwwww=std::accumulate(to.Fwwww,to.Fwwww+mA,0); to.toFwwwd=std::accumulate(to.Fwwwd,to.Fwwwd+mA,0); to.toFwwdd=std::accumulate(to.Fwwdd,to.Fwwdd+mA,0); to.toFwwwr=std::accumulate(to.Fwwwr,to.Fwwwr+mA,0); to.toFwwrr=std::accumulate(to.Fwwrr,to.Fwwrr+mA,0); to.toFwwdr=std::accumulate(to.Fwwdr,to.Fwwdr+mA,0);
	to.tofFwwww=std::accumulate(to.fFwwww,to.fFwwww+mA,0); to.tofFwwwd=std::accumulate(to.fFwwwd,to.fFwwwd+mA,0); to.tofFwwdd=std::accumulate(to.fFwwdd,to.fFwwdd+mA,0); to.tofFwwwr=std::accumulate(to.fFwwwr,to.fFwwwr+mA,0); to.tofFwwrr=std::accumulate(to.fFwwrr,to.fFwwrr+mA,0); to.tofFwwdr=std::accumulate(to.fFwwdr,to.fFwwdr+mA,0);
	to.tomFwwww=std::accumulate(to.mFwwww,to.mFwwww+mA,0); to.tomFwwwd=std::accumulate(to.mFwwwd,to.mFwwwd+mA,0); to.tomFwwdd=std::accumulate(to.mFwwdd,to.mFwwdd+mA,0); to.tomFwwwr=std::accumulate(to.mFwwwr,to.mFwwwr+mA,0); to.tomFwwrr=std::accumulate(to.mFwwrr,to.mFwwrr+mA,0); to.tomFwwdr=std::accumulate(to.mFwwdr,to.mFwwdr+mA,0);
	to.tobFwwww=std::accumulate(to.bFwwww,to.bFwwww+mA,0); to.tobFwwwd=std::accumulate(to.bFwwwd,to.bFwwwd+mA,0); to.tobFwwdd=std::accumulate(to.bFwwdd,to.bFwwdd+mA,0); to.tobFwwwr=std::accumulate(to.bFwwwr,to.bFwwwr+mA,0); to.tobFwwrr=std::accumulate(to.bFwwrr,to.bFwwrr+mA,0); to.tobFwwdr=std::accumulate(to.bFwwdr,to.bFwwdr+mA,0);
	to.tofFwdww=std::accumulate(to.fFwdww,to.fFwdww+mA,0); to.tofFwdwd=std::accumulate(to.fFwdwd,to.fFwdwd+mA,0); to.tofFwddd=std::accumulate(to.fFwddd,to.fFwddd+mA,0); to.tofFwdwr=std::accumulate(to.fFwdwr,to.fFwdwr+mA,0); to.tofFwdrr=std::accumulate(to.fFwdrr,to.fFwdrr+mA,0); to.tofFwddr=std::accumulate(to.fFwddr,to.fFwddr+mA,0);
	to.tomFwdww=std::accumulate(to.mFwdww,to.mFwdww+mA,0); to.tomFwdwd=std::accumulate(to.mFwdwd,to.mFwdwd+mA,0); to.tomFwddd=std::accumulate(to.mFwddd,to.mFwddd+mA,0); to.tomFwdwr=std::accumulate(to.mFwdwr,to.mFwdwr+mA,0); to.tomFwdrr=std::accumulate(to.mFwdrr,to.mFwdrr+mA,0); to.tomFwddr=std::accumulate(to.mFwddr,to.mFwddr+mA,0);
	to.tobFwdww=std::accumulate(to.bFwdww,to.bFwdww+mA,0); to.tobFwdwd=std::accumulate(to.bFwdwd,to.bFwdwd+mA,0); to.tobFwddd=std::accumulate(to.bFwddd,to.bFwddd+mA,0); to.tobFwdwr=std::accumulate(to.bFwdwr,to.bFwdwr+mA,0); to.tobFwdrr=std::accumulate(to.bFwdrr,to.bFwdrr+mA,0); to.tobFwddr=std::accumulate(to.bFwddr,to.bFwddr+mA,0);
	to.toFddww=std::accumulate(to.Fddww,to.Fddww+mA,0); to.toFddwd=std::accumulate(to.Fddwd,to.Fddwd+mA,0); to.toFdddd=std::accumulate(to.Fdddd,to.Fdddd+mA,0); to.toFddwr=std::accumulate(to.Fddwr,to.Fddwr+mA,0); to.toFddrr=std::accumulate(to.Fddrr,to.Fddrr+mA,0); to.toFdddr=std::accumulate(to.Fdddr,to.Fdddr+mA,0);
	to.toFwrww=std::accumulate(to.Fwrww,to.Fwrww+mA,0); to.toFwrwd=std::accumulate(to.Fwrwd,to.Fwrwd+mA,0); to.toFwrdd=std::accumulate(to.Fwrdd,to.Fwrdd+mA,0); to.toFwrwr=std::accumulate(to.Fwrwr,to.Fwrwr+mA,0); to.toFwrrr=std::accumulate(to.Fwrrr,to.Fwrrr+mA,0); to.toFwrdr=std::accumulate(to.Fwrdr,to.Fwrdr+mA,0);
	to.tofFwrww=std::accumulate(to.fFwrww,to.fFwrww+mA,0); to.tofFwrwd=std::accumulate(to.fFwrwd,to.fFwrwd+mA,0); to.tofFwrdd=std::accumulate(to.fFwrdd,to.fFwrdd+mA,0); to.tofFwrwr=std::accumulate(to.fFwrwr,to.fFwrwr+mA,0); to.tofFwrrr=std::accumulate(to.fFwrrr,to.fFwrrr+mA,0); to.tofFwrdr=std::accumulate(to.fFwrdr,to.fFwrdr+mA,0);
	to.tomFwrww=std::accumulate(to.mFwrww,to.mFwrww+mA,0); to.tomFwrwd=std::accumulate(to.mFwrwd,to.mFwrwd+mA,0); to.tomFwrdd=std::accumulate(to.mFwrdd,to.mFwrdd+mA,0); to.tomFwrwr=std::accumulate(to.mFwrwr,to.mFwrwr+mA,0); to.tomFwrrr=std::accumulate(to.mFwrrr,to.mFwrrr+mA,0); to.tomFwrdr=std::accumulate(to.mFwrdr,to.mFwrdr+mA,0);
	to.tobFwrww=std::accumulate(to.bFwrww,to.bFwrww+mA,0); to.tobFwrwd=std::accumulate(to.bFwrwd,to.bFwrwd+mA,0); to.tobFwrdd=std::accumulate(to.bFwrdd,to.bFwrdd+mA,0); to.tobFwrwr=std::accumulate(to.bFwrwr,to.bFwrwr+mA,0); to.tobFwrrr=std::accumulate(to.bFwrrr,to.bFwrrr+mA,0); to.tobFwrdr=std::accumulate(to.bFwrdr,to.bFwrdr+mA,0);
	to.toFrrww=std::accumulate(to.Frrww,to.Frrww+mA,0); to.toFrrwd=std::accumulate(to.Frrwd,to.Frrwd+mA,0); to.toFrrdd=std::accumulate(to.Frrdd,to.Frrdd+mA,0); to.toFrrwr=std::accumulate(to.Frrwr,to.Frrwr+mA,0); to.toFrrrr=std::accumulate(to.Frrrr,to.Frrrr+mA,0); to.toFrrdr=std::accumulate(to.Frrdr,to.Frrdr+mA,0);
	to.toFdrww=std::accumulate(to.Fdrww,to.Fdrww+mA,0); to.toFdrwd=std::accumulate(to.Fdrwd,to.Fdrwd+mA,0); to.toFdrdd=std::accumulate(to.Fdrdd,to.Fdrdd+mA,0); to.toFdrwr=std::accumulate(to.Fdrwr,to.Fdrwr+mA,0); to.toFdrrr=std::accumulate(to.Fdrrr,to.Fdrrr+mA,0); to.toFdrdr=std::accumulate(to.Fdrdr,to.Fdrdr+mA,0);
	return;};

	void LayEggs (){
	FemaleCount();
	int i;
	double TotEggs[14]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	int* eggs;
	if(to.toFwwww>0 && pa.omegaww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegaww*random_binomial(to.toFwwww,pa.delta)),16,pa.fwwww,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFwwwd>0 && pa.omegaww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegaww*random_binomial(to.toFwwwd,pa.delta)),16,pa.fwwwd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFwwdd>0 && pa.omegaww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegaww*random_binomial(to.toFwwdd,pa.delta)),16,pa.fwwdd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFwwwr>0 && pa.omegaww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegaww*random_binomial(to.toFwwwr,pa.delta)),16,pa.fwwwr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFwwrr>0 && pa.omegaww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegaww*random_binomial(to.toFwwrr,pa.delta)),16,pa.fwwrr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFwwdr>0 && pa.omegaww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegaww*random_binomial(to.toFwwdr,pa.delta)),16,pa.fwwdr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };

	if(to.tofFwwww>0 && pa.omegafww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegafww*random_binomial(to.tofFwwww,pa.delta)),16,pa.ffwwww,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tofFwwwd>0 && pa.omegafww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegafww*random_binomial(to.tofFwwwd,pa.delta)),16,pa.ffwwwd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tofFwwdd>0 && pa.omegafww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegafww*random_binomial(to.tofFwwdd,pa.delta)),16,pa.ffwwdd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tofFwwwr>0 && pa.omegafww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegafww*random_binomial(to.tofFwwwr,pa.delta)),16,pa.ffwwwr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tofFwwrr>0 && pa.omegafww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegafww*random_binomial(to.tofFwwrr,pa.delta)),16,pa.ffwwrr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tofFwwdr>0 && pa.omegafww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegafww*random_binomial(to.tofFwwdr,pa.delta)),16,pa.ffwwdr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };

	if(to.tomFwwww>0 && pa.omegamww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegamww*random_binomial(to.tomFwwww,pa.delta)),16,pa.mfwwww,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tomFwwwd>0 && pa.omegamww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegamww*random_binomial(to.tomFwwwd,pa.delta)),16,pa.mfwwwd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tomFwwdd>0 && pa.omegamww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegamww*random_binomial(to.tomFwwdd,pa.delta)),16,pa.mfwwdd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tomFwwwr>0 && pa.omegamww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegamww*random_binomial(to.tomFwwwr,pa.delta)),16,pa.mfwwwr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tomFwwrr>0 && pa.omegamww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegamww*random_binomial(to.tomFwwrr,pa.delta)),16,pa.mfwwrr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tomFwwdr>0 && pa.omegamww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegamww*random_binomial(to.tomFwwdr,pa.delta)),16,pa.mfwwdr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };

	if(to.tobFwwww>0 && pa.omegabww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegabww*random_binomial(to.tobFwwww,pa.delta)),16,pa.bfwwww,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tobFwwwd>0 && pa.omegabww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegabww*random_binomial(to.tobFwwwd,pa.delta)),16,pa.bfwwwd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tobFwwdd>0 && pa.omegabww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegabww*random_binomial(to.tobFwwdd,pa.delta)),16,pa.bfwwdd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tobFwwwr>0 && pa.omegabww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegabww*random_binomial(to.tobFwwwr,pa.delta)),16,pa.bfwwwr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tobFwwrr>0 && pa.omegabww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegabww*random_binomial(to.tobFwwrr,pa.delta)),16,pa.bfwwrr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tobFwwdr>0 && pa.omegabww>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegabww*random_binomial(to.tobFwwdr,pa.delta)),16,pa.bfwwdr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };

	if(to.tofFwdww>0 && pa.omegafwd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegafwd*random_binomial(to.tofFwdww,pa.delta)),16,pa.ffwdww,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tofFwdwd>0 && pa.omegafwd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegafwd*random_binomial(to.tofFwdwd,pa.delta)),16,pa.ffwdwd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tofFwddd>0 && pa.omegafwd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegafwd*random_binomial(to.tofFwddd,pa.delta)),16,pa.ffwddd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tofFwdwr>0 && pa.omegafwd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegafwd*random_binomial(to.tofFwdwr,pa.delta)),16,pa.ffwdwr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tofFwdrr>0 && pa.omegafwd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegafwd*random_binomial(to.tofFwdrr,pa.delta)),16,pa.ffwdrr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tofFwddr>0 && pa.omegafwd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegafwd*random_binomial(to.tofFwddr,pa.delta)),16,pa.ffwddr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };

	if(to.tomFwdww>0 && pa.omegamwd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegamwd*random_binomial(to.tomFwdww,pa.delta)),16,pa.mfwdww,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tomFwdwd>0 && pa.omegamwd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegamwd*random_binomial(to.tomFwdwd,pa.delta)),16,pa.mfwdwd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tomFwddd>0 && pa.omegamwd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegamwd*random_binomial(to.tomFwddd,pa.delta)),16,pa.mfwddd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tomFwdwr>0 && pa.omegamwd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegamwd*random_binomial(to.tomFwdwr,pa.delta)),16,pa.mfwdwr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tomFwdrr>0 && pa.omegamwd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegamwd*random_binomial(to.tomFwdrr,pa.delta)),16,pa.mfwdrr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tomFwddr>0 && pa.omegamwd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegamwd*random_binomial(to.tomFwddr,pa.delta)),16,pa.mfwddr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };

	if(to.tobFwdww>0 && pa.omegabwd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegabwd*random_binomial(to.tobFwdww,pa.delta)),16,pa.bfwdww,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tobFwdwd>0 && pa.omegabwd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegabwd*random_binomial(to.tobFwdwd,pa.delta)),16,pa.bfwdwd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tobFwddd>0 && pa.omegabwd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegabwd*random_binomial(to.tobFwddd,pa.delta)),16,pa.bfwddd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tobFwdwr>0 && pa.omegabwd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegabwd*random_binomial(to.tobFwdwr,pa.delta)),16,pa.bfwdwr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tobFwdrr>0 && pa.omegabwd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegabwd*random_binomial(to.tobFwdrr,pa.delta)),16,pa.bfwdrr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tobFwddr>0 && pa.omegabwd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegabwd*random_binomial(to.tobFwddr,pa.delta)),16,pa.bfwddr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };

	if(to.toFddww>0 && pa.omegadd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegadd*random_binomial(to.toFddww,pa.delta)),16,pa.fddww,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFddwd>0 && pa.omegadd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegadd*random_binomial(to.toFddwd,pa.delta)),16,pa.fddwd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFdddd>0 && pa.omegadd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegadd*random_binomial(to.toFdddd,pa.delta)),16,pa.fdddd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFddwr>0 && pa.omegadd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegadd*random_binomial(to.toFddwr,pa.delta)),16,pa.fddwr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFddrr>0 && pa.omegadd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegadd*random_binomial(to.toFddrr,pa.delta)),16,pa.fddrr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFdddr>0 && pa.omegadd>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegadd*random_binomial(to.toFdddr,pa.delta)),16,pa.fdddr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };

	if(to.toFwrww>0 && pa.omegawr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegawr*random_binomial(to.toFwrww,pa.delta)),16,pa.fwrww,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFwrwd>0 && pa.omegawr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegawr*random_binomial(to.toFwrwd,pa.delta)),16,pa.fwrwd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFwrdd>0 && pa.omegawr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegawr*random_binomial(to.toFwrdd,pa.delta)),16,pa.fwrdd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFwrwr>0 && pa.omegawr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegawr*random_binomial(to.toFwrwr,pa.delta)),16,pa.fwrwr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFwrrr>0 && pa.omegawr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegawr*random_binomial(to.toFwrrr,pa.delta)),16,pa.fwrrr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFwrdr>0 && pa.omegawr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegawr*random_binomial(to.toFwrdr,pa.delta)),16,pa.fwrdr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };

	if(to.tofFwrww>0 && pa.omegafwr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegafwr*random_binomial(to.tofFwrww,pa.delta)),16,pa.ffwrww,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tofFwrwd>0 && pa.omegafwr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegafwr*random_binomial(to.tofFwrwd,pa.delta)),16,pa.ffwrwd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tofFwrdd>0 && pa.omegafwr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegafwr*random_binomial(to.tofFwrdd,pa.delta)),16,pa.ffwrdd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tofFwrwr>0 && pa.omegafwr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegafwr*random_binomial(to.tofFwrwr,pa.delta)),16,pa.ffwrwr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tofFwrrr>0 && pa.omegafwr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegafwr*random_binomial(to.tofFwrrr,pa.delta)),16,pa.ffwrrr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tofFwrdr>0 && pa.omegafwr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegafwr*random_binomial(to.tofFwrdr,pa.delta)),16,pa.ffwrdr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };

	if(to.tomFwrww>0 && pa.omegamwr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegamwr*random_binomial(to.tomFwrww,pa.delta)),16,pa.mfwrww,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tomFwrwd>0 && pa.omegamwr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegamwr*random_binomial(to.tomFwrwd,pa.delta)),16,pa.mfwrwd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tomFwrdd>0 && pa.omegamwr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegamwr*random_binomial(to.tomFwrdd,pa.delta)),16,pa.mfwrdd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tomFwrwr>0 && pa.omegamwr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegamwr*random_binomial(to.tomFwrwr,pa.delta)),16,pa.mfwrwr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tomFwrrr>0 && pa.omegamwr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegamwr*random_binomial(to.tomFwrrr,pa.delta)),16,pa.mfwrrr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tomFwrdr>0 && pa.omegamwr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegamwr*random_binomial(to.tomFwrdr,pa.delta)),16,pa.mfwrdr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };


	if(to.tobFwrww>0 && pa.omegabwr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegabwr*random_binomial(to.tobFwrww,pa.delta)),16,pa.bfwrww,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tobFwrwd>0 && pa.omegabwr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegabwr*random_binomial(to.tobFwrwd,pa.delta)),16,pa.bfwrwd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tobFwrdd>0 && pa.omegabwr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegabwr*random_binomial(to.tobFwrdd,pa.delta)),16,pa.bfwrdd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tobFwrwr>0 && pa.omegabwr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegabwr*random_binomial(to.tobFwrwr,pa.delta)),16,pa.bfwrwr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tobFwrrr>0 && pa.omegabwr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegabwr*random_binomial(to.tobFwrrr,pa.delta)),16,pa.bfwrrr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.tobFwrdr>0 && pa.omegabwr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegabwr*random_binomial(to.tobFwrdr,pa.delta)),16,pa.bfwrdr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	
	if(to.toFrrww>0 && pa.omegarr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegarr*random_binomial(to.toFrrww,pa.delta)),16,pa.frrww,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFrrwd>0 && pa.omegarr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegarr*random_binomial(to.toFrrwd,pa.delta)),16,pa.frrwd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFrrdd>0 && pa.omegarr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegarr*random_binomial(to.toFrrdd,pa.delta)),16,pa.frrdd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFrrwr>0 && pa.omegarr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegarr*random_binomial(to.toFrrwr,pa.delta)),16,pa.frrwr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFrrrr>0 && pa.omegarr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegarr*random_binomial(to.toFrrrr,pa.delta)),16,pa.frrrr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFrrdr>0 && pa.omegarr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegarr*random_binomial(to.toFrrdr,pa.delta)),16,pa.frrdr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFdrww>0 && pa.omegadr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegadr*random_binomial(to.toFdrww,pa.delta)),16,pa.fdrww,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFdrwd>0 && pa.omegadr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegadr*random_binomial(to.toFdrwd,pa.delta)),16,pa.fdrwd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFdrdd>0 && pa.omegadr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegadr*random_binomial(to.toFdrdd,pa.delta)),16,pa.fdrdd,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFdrwr>0 && pa.omegadr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegadr*random_binomial(to.toFdrwr,pa.delta)),16,pa.fdrwr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFdrrr>0 && pa.omegadr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegadr*random_binomial(to.toFdrrr,pa.delta)),16,pa.fdrrr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };
	if(to.toFdrdr>0 && pa.omegadr>0)
	{ eggs=random_multinom_var(random_poisson(pa.theta*pa.omegadr*random_binomial(to.toFdrdr,pa.delta)),16,pa.fdrdr,1);
	for(i=0;i<14;i++)TotEggs[i]+=*(eggs+i); delete[] eggs; };


		int* eggkeep;
		int eT=std::accumulate(TotEggs,TotEggs+14,0);
		to.neggs=eT;
		to.r2= TotEggs[6]+ TotEggs[7]+ TotEggs[8]+ 2*TotEggs[9]+ TotEggs[10]+ TotEggs[13];
		int ne;
		ne=min(eT,pa.NumEggs);
		to.keepeggs=ne;
		eggkeep= random_multinom_var(ne,14,TotEggs,eT);
				
//fraction of 1:ww,2:fww,3:mww,4:fwd,5:mwd,6:dd,7:wr,8:fwr,9:mwr,10:rr,11:dr,12:bww,13:bwd,14:bwr
	to.Jww[0]=*eggkeep; to.fJww[0]=*(eggkeep+1); to.mJww[0]=*(eggkeep+2); to.fJwd[0]=*(eggkeep+3); to.mJwd[0]=*(eggkeep+4); to.Jdd[0]=*(eggkeep+5); to.Jwr[0]=*(eggkeep+6); to.fJwr[0]=*(eggkeep+7); to.mJwr[0]=*(eggkeep+8); to.Jrr[0]=*(eggkeep+9); to.Jdr[0]=*(eggkeep+10); to.bJww[0]=*(eggkeep+11); to.bJwd[0]=*(eggkeep+12); to.bJwr[0]=*(eggkeep+13);

		to.JTot+=ne;
				return;};

	void AdultsDie(int day){
		int num;
		double mort;
		for(int age=0;age<mA;age++)
		{
	mort=1-std::exp(std::pow((1.0*age/pa.lamM[day]),pa.kM[day])-std::pow(((1.0+1.0*age)/pa.lamM[day]),pa.kM[day]));
		num=random_binomial(to.Mww[age],mort);to.Mww[age]-=num;to.MTot-=num;
		num=random_binomial(to.Mwr[age],mort);to.Mwr[age]-=num;to.MTot-=num;
		num=random_binomial(to.Mrr[age],mort);to.Mrr[age]-=num;to.MTot-=num;
		num=random_binomial(to.Mdr[age],mort);to.Mdr[age]-=num;to.MTot-=num;
		num=random_binomial(to.Mwd[age],mort);to.Mwd[age]-=num;to.MTot-=num;
		num=random_binomial(to.Mdd[age],mort);to.Mdd[age]-=num;to.MTot-=num;

	mort=1-std::exp(std::pow((1.0*age/pa.lamF[day]),pa.kF[day])-std::pow(((1.0+1.0*age)/pa.lamF[day]),pa.kF[day]));
		num=random_binomial(to.Vww[age],mort);to.Vww[age]-=num;to.VTot-=num;
		num=random_binomial(to.fVww[age],mort);to.fVww[age]-=num;to.VTot-=num;
		num=random_binomial(to.mVww[age],mort);to.mVww[age]-=num;to.VTot-=num;
		num=random_binomial(to.bVww[age],mort);to.bVww[age]-=num;to.VTot-=num;
		num=random_binomial(to.fVwd[age],mort);to.fVwd[age]-=num;to.VTot-=num;
		num=random_binomial(to.mVwd[age],mort);to.mVwd[age]-=num;to.VTot-=num;
		num=random_binomial(to.bVwd[age],mort);to.bVwd[age]-=num;to.VTot-=num;
		num=random_binomial(to.Vdd[age],mort);to.Vdd[age]-=num;to.VTot-=num;
		num=random_binomial(to.Vwr[age],mort);to.Vwr[age]-=num;to.VTot-=num;
		num=random_binomial(to.fVwr[age],mort);to.fVwr[age]-=num;to.VTot-=num;
		num=random_binomial(to.mVwr[age],mort);to.mVwr[age]-=num;to.VTot-=num;
		num=random_binomial(to.bVwr[age],mort);to.bVwr[age]-=num;to.VTot-=num;
		num=random_binomial(to.Vrr[age],mort);to.Vrr[age]-=num;to.VTot-=num;
		num=random_binomial(to.Vdr[age],mort);to.Vdr[age]-=num;to.VTot-=num;


		num=random_binomial(to.Fwwww[age],mort);to.Fwwww[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fwwwd[age],mort);to.Fwwwd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fwwdd[age],mort);to.Fwwdd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fwwwr[age],mort);to.Fwwwr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fwwrr[age],mort);to.Fwwrr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fwwdr[age],mort);to.Fwwdr[age]-=num;to.FTot-=num;	

		num=random_binomial(to.fFwwww[age],mort);to.fFwwww[age]-=num;to.FTot-=num;	
		num=random_binomial(to.fFwwwd[age],mort);to.fFwwwd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.fFwwdd[age],mort);to.fFwwdd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.fFwwwr[age],mort);to.fFwwwr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.fFwwrr[age],mort);to.fFwwrr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.fFwwdr[age],mort);to.fFwwdr[age]-=num;to.FTot-=num;	

		num=random_binomial(to.mFwwww[age],mort);to.mFwwww[age]-=num;to.FTot-=num;	
		num=random_binomial(to.mFwwwd[age],mort);to.mFwwwd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.mFwwdd[age],mort);to.mFwwdd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.mFwwwr[age],mort);to.mFwwwr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.mFwwrr[age],mort);to.mFwwrr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.mFwwdr[age],mort);to.mFwwdr[age]-=num;to.FTot-=num;	

		num=random_binomial(to.bFwwww[age],mort);to.bFwwww[age]-=num;to.FTot-=num;	
		num=random_binomial(to.bFwwwd[age],mort);to.bFwwwd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.bFwwdd[age],mort);to.bFwwdd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.bFwwwr[age],mort);to.bFwwwr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.bFwwrr[age],mort);to.bFwwrr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.bFwwdr[age],mort);to.bFwwdr[age]-=num;to.FTot-=num;	

		num=random_binomial(to.fFwdww[age],mort);to.fFwdww[age]-=num;to.FTot-=num;	
		num=random_binomial(to.fFwdwd[age],mort);to.fFwdwd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.fFwddd[age],mort);to.fFwddd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.fFwdwr[age],mort);to.fFwdwr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.fFwdrr[age],mort);to.fFwdrr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.fFwddr[age],mort);to.fFwddr[age]-=num;to.FTot-=num;	

		num=random_binomial(to.mFwdww[age],mort);to.mFwdww[age]-=num;to.FTot-=num;	
		num=random_binomial(to.mFwdwd[age],mort);to.mFwdwd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.mFwddd[age],mort);to.mFwddd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.mFwdwr[age],mort);to.mFwdwr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.mFwdrr[age],mort);to.mFwdrr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.mFwddr[age],mort);to.mFwddr[age]-=num;to.FTot-=num;	

		num=random_binomial(to.bFwdww[age],mort);to.bFwdww[age]-=num;to.FTot-=num;	
		num=random_binomial(to.bFwdwd[age],mort);to.bFwdwd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.bFwddd[age],mort);to.bFwddd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.bFwdwr[age],mort);to.bFwdwr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.bFwdrr[age],mort);to.bFwdrr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.bFwddr[age],mort);to.bFwddr[age]-=num;to.FTot-=num;	

		num=random_binomial(to.Fddww[age],mort);to.Fddww[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fddwd[age],mort);to.Fddwd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fdddd[age],mort);to.Fdddd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fddwr[age],mort);to.Fddwr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fddrr[age],mort);to.Fddrr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fdddr[age],mort);to.Fdddr[age]-=num;to.FTot-=num;	

		num=random_binomial(to.Fwrww[age],mort);to.Fwrww[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fwrwd[age],mort);to.Fwrwd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fwrdd[age],mort);to.Fwrdd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fwrwr[age],mort);to.Fwrwr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fwrrr[age],mort);to.Fwrrr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fwrdr[age],mort);to.Fwrdr[age]-=num;to.FTot-=num;	

		num=random_binomial(to.fFwrww[age],mort);to.fFwrww[age]-=num;to.FTot-=num;	
		num=random_binomial(to.fFwrwd[age],mort);to.fFwrwd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.fFwrdd[age],mort);to.fFwrdd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.fFwrwr[age],mort);to.fFwrwr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.fFwrrr[age],mort);to.fFwrrr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.fFwrdr[age],mort);to.fFwrdr[age]-=num;to.FTot-=num;	

		num=random_binomial(to.mFwrww[age],mort);to.mFwrww[age]-=num;to.FTot-=num;	
		num=random_binomial(to.mFwrwd[age],mort);to.mFwrwd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.mFwrdd[age],mort);to.mFwrdd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.mFwrwr[age],mort);to.mFwrwr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.mFwrrr[age],mort);to.mFwrrr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.mFwrdr[age],mort);to.mFwrdr[age]-=num;to.FTot-=num;	

		num=random_binomial(to.bFwrww[age],mort);to.bFwrww[age]-=num;to.FTot-=num;	
		num=random_binomial(to.bFwrwd[age],mort);to.bFwrwd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.bFwrdd[age],mort);to.bFwrdd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.bFwrwr[age],mort);to.bFwrwr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.bFwrrr[age],mort);to.bFwrrr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.bFwrdr[age],mort);to.bFwrdr[age]-=num;to.FTot-=num;	

		num=random_binomial(to.Frrww[age],mort);to.Frrww[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Frrwd[age],mort);to.Frrwd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Frrdd[age],mort);to.Frrdd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Frrwr[age],mort);to.Frrwr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Frrrr[age],mort);to.Frrrr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Frrdr[age],mort);to.Frrdr[age]-=num;to.FTot-=num;	

		num=random_binomial(to.Fdrww[age],mort);to.Fdrww[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fdrwd[age],mort);to.Fdrwd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fdrdd[age],mort);to.Fdrdd[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fdrwr[age],mort);to.Fdrwr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fdrrr[age],mort);to.Fdrrr[age]-=num;to.FTot-=num;	
		num=random_binomial(to.Fdrdr[age],mort);to.Fdrdr[age]-=num;to.FTot-=num;	
				};
		return;};


	void UpdateMate(void){
				to.mate_rate=to.MTot/(pa.beta+to.MTot);
		return;};
	void SetFertility()
{
//fraction of 1:ww,2:fww,3:mww,4:fwd,5:mwd,6:dd,7:wr,8:fwr,9:mwr,10:rr,11:dr,12:bww,13:bwd,14:bwr
	double Fwwww[14]={1,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double Fwwwd[14]={0,0,(1-pa.em-pa.Mgamma)*0.5,0,(1+pa.em)*0.5,0,0,0,pa.Mgamma*0.5,0,0,0,0,0};
	double Fwwdd[14]={0,0,0,0,1,0,0,0,0,0,0,0,0,0};
	double Fwwwr[14]={0.5,0,0,0,0,0,0.5,0,0,0,0,0,0,0};
	double Fwwrr[14]={0,0,0,0,0,0,1,0,0,0,0,0,0,0};
	double Fwwdr[14]={0,0,0,0,0.5,0,0.5,0,0,0,0,0,0,0};

	double fFwwww[14]={1,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double fFwwwd[14]={0,0,(1-pa.em-pa.Mgamma)*0.5,0,(1+pa.em)*0.5,0,0,0,pa.Mgamma*0.5,0,0,0,0,0};
	double fFwwdd[14]={0,0,0,0,1,0,0,0,0,0,0,0,0,0};
	double fFwwwr[14]={0.5,0,0,0,0,0,0.5,0,0,0,0,0,0,0};
	double fFwwrr[14]={0,0,0,0,0,0,1,0,0,0,0,0,0,0};
	double fFwwdr[14]={0,0,0,0,0.5,0,0.5,0,0,0,0,0,0,0};
	double mFwwww[14]={1,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double mFwwwd[14]={0,0,(1-pa.em-pa.Mgamma)*0.5,0,(1+pa.em)*0.5,0,0,0,pa.Mgamma*0.5,0,0,0,0,0};
	double mFwwdd[14]={0,0,0,0,1,0,0,0,0,0,0,0,0,0};
	double mFwwwr[14]={0.5,0,0,0,0,0,0.5,0,0,0,0,0,0,0};
	double mFwwrr[14]={0,0,0,0,0,0,1,0,0,0,0,0,0,0};
	double mFwwdr[14]={0,0,0,0,0.5,0,0.5,0,0,0,0,0,0,0};

	double bFwwww[14]={1,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double bFwwwd[14]={0,0,(1-pa.em-pa.Mgamma)*0.5,0,(1+pa.em)*0.5,0,0,0,pa.Mgamma*0.5,0,0,0,0,0};
	double bFwwdd[14]={0,0,0,0,1,0,0,0,0,0,0,0,0,0};
	double bFwwwr[14]={0.5,0,0,0,0,0,0.5,0,0,0,0,0,0,0};
	double bFwwrr[14]={0,0,0,0,0,0,1,0,0,0,0,0,0,0};
	double bFwwdr[14]={0,0,0,0,0.5,0,0.5,0,0,0,0,0,0,0};


	double fFwdww[14]={0,(1-pa.ef-pa.Fgamma)*0.5,0,(1+pa.ef)*0.5,0,0,0,pa.Fgamma*0.5,0,0,0,0,0,0};
	double fFwdwd[14]={0,0,0,0,0, (1+pa.ef)*(1+pa.em)*0.25,0,0,0, pa.Mgamma*pa.Fgamma*0.25,(1+pa.ef)*pa.Mgamma*0.25+(1+pa.em)*pa.Fgamma*0.25,(1-pa.ef-pa.Fgamma)*(1-pa.em-pa.Mgamma)*0.25,(1-pa.ef-pa.Fgamma)*(1+pa.em)*0.25+ (1+pa.ef)*(1-pa.em-pa.Mgamma)*0.25,(1-pa.ef-pa.Fgamma)*pa.Mgamma*0.25+ (1-pa.em-pa.Mgamma)*pa.Fgamma*0.25};
	double fFwddd[14]={0,0,0,0,0,(1+pa.ef)*0.5,0,0,0,0,0,0,(1-pa.ef-pa.Fgamma)*0.5,pa.Fgamma*0.5};
	double fFwdwr[14]={0,(1-pa.ef-pa.Fgamma)*0.25,0,(1+pa.ef)*0.25,0,0,0,(1-pa.ef-pa.Fgamma)*0.25+pa.Fgamma*0.25,0,pa.Fgamma*0.25,(1+pa.ef)*0.25,0,0,0};
	double fFwdrr[14]={0,0,0,0,0,0,0,(1-pa.ef-pa.Fgamma)*0.5,0,pa.Fgamma*0.5,(1+pa.ef)*0.5,0,0,0};
	double fFwddr[14]={0,0,0,0,0,(1+pa.ef)*0.25,0,0,0,pa.Fgamma*0.25,((1+pa.ef)*0.25+pa.Fgamma*0.25),0,(1-pa.ef-pa.Fgamma)*0.25,(1-pa.ef-pa.Fgamma)*0.25};
	double mFwdww[14]={0,(1-pa.ef-pa.Fgamma)*0.5,0,(1+pa.ef)*0.5,0,0,0,pa.Fgamma*0.5,0,0,0,0,0,0};
	double mFwdwd[14]={0,0,0,0,0,(1+pa.ef)*(1+pa.em)*0.25,0,0,0 , pa.Mgamma*pa.Fgamma*0.25,(1+pa.em)*pa.Fgamma*0.25+(1+pa.ef)*pa.Mgamma*0.25,(1-pa.ef-pa.Fgamma)*(1-pa.em-pa.Mgamma)*0.25,(1-pa.ef-pa.Fgamma)*(1+pa.em)*0.25+ (1+pa.ef)*(1-pa.em-pa.Mgamma)*0.25,(1-pa.ef-pa.Fgamma)*pa.Mgamma*0.25+ (1-pa.em-pa.Mgamma)*pa.Fgamma*0.25};
	double mFwddd[14]={0,0,0,0,0,(1+pa.ef)*0.5,0,0,0,0,pa.Fgamma*0.5,0,(1-pa.ef-pa.Fgamma)*0.5,0};
	double mFwdwr[14]={0,(1-pa.ef-pa.Fgamma)*0.25,0,(1+pa.ef)*0.25,0,0,0,((1-pa.ef-pa.Fgamma)*0.25+pa.Fgamma*0.25),0,pa.Fgamma*0.25,(1+pa.ef)*0.25,0,0,0};
	double mFwdrr[14]={0,0,0,0,0,0,0,(1-pa.ef-pa.Fgamma)*0.5,pa.Fgamma*0.5,0,(1+pa.ef)*0.5,0,0,0};
	double mFwddr[14]={0,0,0,0,0,(1+pa.ef)*0.25,0,0,0,pa.Fgamma*0.25,((1+pa.ef)*0.25+pa.Fgamma*0.25),0,(1-pa.ef-pa.Fgamma)*0.25,(1-pa.ef-pa.Fgamma)*0.25};
	double bFwdww[14]={0,(1-pa.ef-pa.Fgamma)*0.5,0,(1+pa.ef)*0.5,0,0,0,pa.Fgamma*0.5,0,0,0,0,0,0};
	double bFwdwd[14]={0,0,0,0,0,(1+pa.ef)*(1+pa.em)*0.25,0,0,0 , pa.Mgamma*pa.Fgamma*0.25,(1+pa.em)*pa.Fgamma*0.25+(1+pa.ef)*pa.Mgamma*0.25,(1-pa.ef-pa.Fgamma)*(1-pa.em-pa.Mgamma)*0.25,(1-pa.ef-pa.Fgamma)*(1+pa.em)*0.25+ (1+pa.ef)*(1-pa.em-pa.Mgamma)*0.25,(1-pa.ef-pa.Fgamma)*pa.Mgamma*0.25+ (1-pa.em-pa.Mgamma)*pa.Fgamma*0.25};
	double bFwddd[14]={0,0,0,0,0,(1+pa.ef)*0.5,0,0,0,0,pa.Fgamma*0.5,0,(1-pa.ef-pa.Fgamma)*0.5,0};
	double bFwdwr[14]={0,(1-pa.ef-pa.Fgamma)*0.25,0,(1+pa.ef)*0.25,0,0,0,((1-pa.ef-pa.Fgamma)*0.25+pa.Fgamma*0.25),0,pa.Fgamma*0.25,(1+pa.ef)*0.25,0,0,0};
	double bFwdrr[14]={0,0,0,0,0,0,0,(1-pa.ef-pa.Fgamma)*0.5,pa.Fgamma*0.5,0,(1+pa.ef)*0.5,0,0,0};
	double bFwddr[14]={0,0,0,0,0,(1+pa.ef)*0.25,0,0,0,pa.Fgamma*0.25,((1+pa.ef)*0.25+pa.Fgamma*0.25),0,(1-pa.ef-pa.Fgamma)*0.25,(1-pa.ef-pa.Fgamma)*0.25};
	double Fddww[14]={0,0,0,1,0,0,0,0,0,0,0,0,0,0};
	double Fddwd[14]={0,0,0,0,0,(1+pa.em)*0.5,0,0,0,0,pa.Mgamma*0.5,0,(1-pa.em-pa.Mgamma)*0.5,0};
	double Fdddd[14]={0,0,0,0,0,1,0,0,0,0,0,0,0,0};
	double Fddwr[14]={0,0,0,0.5,0,0,0,0,0,0,0.5,0,0,0};
	double Fddrr[14]={0,0,0,0,0,0,0,0,0,0,1,0,0,0};
	double Fdddr[14]={0,0,0,0,0,0.5,0,0,0,0,0.5,0,0,0};
	double Fwrww[14]={0.5,0,0,0,0,0,0.5,0,0,0,0,0,0,0};
	double Fwrwd[14]={0,0,(1-pa.em-pa.Mgamma)*0.25,0,(1+pa.em)*0.25,0,0,0,(pa.Mgamma*0.25+(1-pa.em-pa.Mgamma)*0.25),pa.Mgamma*0.25,(1+pa.em)*0.25,0,0,0};
	double Fwrdd[14]={0,0,0,0,0.5,0,0,0,0,0,0.5,0,0,0};
	double Fwrwr[14]={0.25,0,0,0,0,0,0.5,0,0,0.25,0,0,0,0};
	double Fwrrr[14]={0,0,0,0,0,0,0.5,0,0,0.5,0,0,0,0};
	double Fwrdr[14]={0,0,0,0,0.25,0,0,0,0.25,0.25,0.25,0,0,0};
	double fFwrww[14]={0.5,0,0,0,0,0,0.5,0,0,0,0,0,0,0};
	double fFwrwd[14]={0,0,(1-pa.em-pa.Mgamma)*0.25,0,(1+pa.em)*0.25,0,0,0,(pa.Mgamma*0.25+(1-pa.em-pa.Mgamma)*0.25),pa.Mgamma*0.25,(1+pa.em)*0.25,0,0,0};
	double fFwrdd[14]={0,0,0,0,0.5,0,0,0,0,0,0.5,0,0,0};
	double fFwrwr[14]={0.25,0,0,0,0,0,0.5,0,0,0.25,0,0,0,0};
	double fFwrrr[14]={0,0,0,0,0,0,0.5,0,0,0.5,0,0,0,0};
	double fFwrdr[14]={0,0,0,0,0.25,0,0,0,0.25,0.25,0.25,0,0,0};
	double mFwrww[14]={0.5,0,0,0,0,0,0.5,0,0,0,0,0,0,0};
	double mFwrwd[14]={0,0,(1-pa.em-pa.Mgamma)*0.25,0,(1+pa.em)*0.25,0,0,0,(pa.Mgamma*0.25+(1-pa.em-pa.Mgamma)*0.25),pa.Mgamma*0.25,(1+pa.em)*0.25,0,0,0};
	double mFwrdd[14]={0,0,0,0,0.5,0,0,0,0,0,0.5,0,0,0};
	double mFwrwr[14]={0.25,0,0,0,0,0,0.5,0,0,0.25,0,0,0,0};
	double mFwrrr[14]={0,0,0,0,0,0,0.5,0,0,0.5,0,0,0,0};
	double mFwrdr[14]={0,0,0,0,0.25,0,0,0,0.25,0.25,0.25,0,0,0};
	double bFwrww[14]={0.5,0,0,0,0,0,0.5,0,0,0,0,0,0,0};
	double bFwrwd[14]={0,0,(1-pa.em-pa.Mgamma)*0.25,0,(1+pa.em)*0.25,0,0,0,(pa.Mgamma*0.25+(1-pa.em-pa.Mgamma)*0.25),pa.Mgamma*0.25,(1+pa.em)*0.25,0,0,0};
	double bFwrdd[14]={0,0,0,0,0.5,0,0,0,0,0,0.5,0,0,0};
	double bFwrwr[14]={0.25,0,0,0,0,0,0.5,0,0,0.25,0,0,0,0};
	double bFwrrr[14]={0,0,0,0,0,0,0.5,0,0,0.5,0,0,0,0};
	double bFwrdr[14]={0,0,0,0,0.25,0,0,0,0.25,0.25,0.25,0,0,0};
	double Frrww[14]={0,0,0,0,0,0,1,0,0,0,0,0,0,0};
	double Frrwd[14]={0,0,0,0,0,0,0,0,(1-pa.em-pa.Mgamma)*0.5,0.5*pa.Mgamma,(1+pa.em)*0.5,0,0,0};
	double Frrdd[14]={0,0,0,0,0,0,0,0,0,0,1,0,0,0};
	double Frrwr[14]={0,0,0,0,0,0,0.5,0,0,0.5,0,0,0,0};
	double Frrrr[14]={0,0,0,0,0,0,0,0,0,1,0,0,0,0};
	double Frrdr[14]={0,0,0,0,0,0,0,0,0,0.5,0.5,0,0,0};
	double Fdrww[14]={0,0,0,0.5,0,0,0,0.5,0,0,0,0,0,0};
	double Fdrwd[14]={0,0,0,0,(1-pa.em-pa.Mgamma)*0.25,(1+pa.em)*0.25,0,0,(1-pa.em-pa.Mgamma)*0.25,pa.Mgamma*0.25,((1+pa.em)*0.25+pa.Mgamma*0.25),0,0,0};
	double Fdrdd[14]={0,0,0,0,0,0.5,0,0,0,0,0.5,0,0,0};
	double Fdrwr[14]={0,0,0,0.25,0,0,0,0.25,0,0.25,0.25,0,0,0};
	double Fdrrr[14]={0,0,0,0,0,0,0,0,0,0.5,0.5,0,0,0};
	double Fdrdr[14]={0,0,0,0,0,0.25,0,0,0,0.25,0.5,0,0,0};


	for(int k=0;k<16;k++)
	{
	pa.fwwww[k]=Fwwww[k]; pa.fwwwd[k]=Fwwwd[k]; pa.fwwdd[k]=Fwwdd[k]; pa.fwwwr[k]=Fwwwr[k]; pa.fwwrr[k]=Fwwrr[k]; pa.fwwdr[k]=Fwwdr[k]; 
	pa.ffwwww[k]=fFwwww[k]; pa.ffwwwd[k]=fFwwwd[k]; pa.ffwwdd[k]=fFwwdd[k]; pa.ffwwwr[k]=fFwwwr[k]; pa.ffwwrr[k]=fFwwrr[k]; pa.ffwwdr[k]=fFwwdr[k]; 
	pa.mfwwww[k]=mFwwww[k]; pa.mfwwwd[k]=mFwwwd[k]; pa.mfwwdd[k]=mFwwdd[k]; pa.mfwwwr[k]=mFwwwr[k]; pa.mfwwrr[k]=mFwwrr[k]; pa.mfwwdr[k]=mFwwdr[k]; 
	pa.bfwwww[k]=bFwwww[k]; pa.bfwwwd[k]=bFwwwd[k]; pa.bfwwdd[k]=bFwwdd[k]; pa.bfwwwr[k]=bFwwwr[k]; pa.bfwwrr[k]=bFwwrr[k]; pa.bfwwdr[k]=bFwwdr[k]; 
	
	pa.ffwdww[k]=fFwdww[k]; pa.ffwdwd[k]=fFwdwd[k]; pa.ffwddd[k]=fFwddd[k]; pa.ffwdwr[k]=fFwdwr[k]; pa.ffwdrr[k]=fFwdrr[k]; pa.ffwddr[k]=fFwddr[k]; 
	pa.mfwdww[k]=mFwdww[k]; pa.mfwdwd[k]=mFwdwd[k]; pa.mfwddd[k]=mFwddd[k]; pa.mfwdwr[k]=mFwdwr[k]; pa.mfwdrr[k]=mFwdrr[k]; pa.mfwddr[k]=mFwddr[k]; 
	pa.bfwdww[k]=bFwdww[k]; pa.bfwdwd[k]=bFwdwd[k]; pa.bfwddd[k]=bFwddd[k]; pa.bfwdwr[k]=bFwdwr[k]; pa.bfwdrr[k]=bFwdrr[k]; pa.bfwddr[k]=bFwddr[k]; 
	pa.fddww[k]=Fddww[k]; pa.fddwd[k]=Fddwd[k]; pa.fdddd[k]=Fdddd[k]; pa.fddwr[k]=Fddwr[k]; pa.fddrr[k]=Fddrr[k]; pa.fdddr[k]=Fdddr[k]; 
	pa.fwrww[k]=Fwrww[k]; pa.fwrwd[k]=Fwrwd[k]; pa.fwrdd[k]=Fwrdd[k]; pa.fwrwr[k]=Fwrwr[k]; pa.fwrrr[k]=Fwrrr[k]; pa.fwrdr[k]=Fwrdr[k]; 
	pa.ffwrww[k]=fFwrww[k]; pa.ffwrwd[k]=fFwrwd[k]; pa.ffwrdd[k]=fFwrdd[k]; pa.ffwrwr[k]=fFwrwr[k]; pa.ffwrrr[k]=fFwrrr[k]; pa.ffwrdr[k]=fFwrdr[k]; 
	pa.mfwrww[k]=mFwrww[k]; pa.mfwrwd[k]=mFwrwd[k]; pa.mfwrdd[k]=mFwrdd[k]; pa.mfwrwr[k]=mFwrwr[k]; pa.mfwrrr[k]=mFwrrr[k]; pa.mfwrdr[k]=mFwrdr[k]; 
	pa.bfwrww[k]=bFwrww[k]; pa.bfwrwd[k]=bFwrwd[k]; pa.bfwrdd[k]=bFwrdd[k]; pa.bfwrwr[k]=bFwrwr[k]; pa.bfwrrr[k]=bFwrrr[k]; pa.bfwrdr[k]=bFwrdr[k]; 
	pa.frrww[k]=Frrww[k]; pa.frrwd[k]=Frrwd[k]; pa.frrdd[k]=Frrdd[k]; pa.frrwr[k]=Frrwr[k]; pa.frrrr[k]=Frrrr[k]; pa.frrdr[k]=Frrdr[k]; 
	pa.fdrww[k]=Fdrww[k]; pa.fdrwd[k]=Fdrwd[k]; pa.fdrdd[k]=Fdrdd[k]; pa.fdrwr[k]=Fdrwr[k]; pa.fdrrr[k]=Fdrrr[k]; pa.fdrdr[k]=Fdrdr[k]; 
	};


	return;};


	int random_binomial(int N,double p){
		int ran;
		if(N==0){ran=0;}
		else if(p>0.999999){ran=N;}
		else if(p<0.000001){ran=0;}
		else if(N*p>10 && N*(1-p)>10) {ran=max(0,min(N,int(random_normal(N*p,sqrt(N*p*(1-p))))));}
		else if((N>20 && p<0.05) || (N>100 && N*p<10)){ran=random_poisson(N*p);}
		else if((N>20 && p>0.95) || (N>100 && N*(1-p)<10)){ran=N-random_poisson(N*(1-p));}
		else { variate_generator<mt19937, binomial_distribution<> > b(mt19937(int(rg.Random()*time(NULL))), binomial_distribution<>(N,p)); ran=b();};
		return ran;};


       double random_normal(double mu, double sigma)
{
	const double epsilon = std::numeric_limits<double>::min();

	static double z0, z1;
	static bool generate;
	generate = !generate;

	if (!generate)
	   return z1 * sigma + mu;

	double u1, u2;
	do
	 {
	   u1 = rg.Random();
	   u2 = rg.Random();
	 }
	while ( u1 <= epsilon );

	z0 = sqrt(-2.0 * log(u1)) * cos(TWOPI * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(TWOPI* u2);
	return z0 * sigma + mu;
};
double dist (double x1, double y1, double x2, double y2)
	{
		return double(sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)));
	};


	int random_poisson(double landa)
	{
		int k;
		if(landa<1e-5){k=0;}
		else if(landa>30) {k=max(0,(int)random_normal(landa,sqrt(landa)));}
		else
			{
			double p=exp(-landa);
			double g=p;
			double u=rg.Random()*0.999999999;
			k=0;
			while (u>g)
			    {
				p*=(landa/(double)(++k));
				g+=p;
			};
			};
	return k;
	};

int randNegBin(double r,double p)
{
int numsuc=0;
int numfail=0;
if(r>0.001)
{
while(double(numfail)<r)
{
if(rg.Random()>p)numsuc++; else numfail++;
};
};
return numsuc;};


	int* random_multinom(int N,int probs[6])
{
	int *ran=new int[6];
	double sum=double(probs[0]+probs[1]+probs[2]+probs[3]+probs[4]+probs[5]);
	ran[0]=random_binomial(N,probs[0]/sum);
	ran[1]=random_binomial(N-ran[0],probs[1]/(sum-probs[0]));
	ran[2]=random_binomial(N-ran[0]-ran[1],probs[2]/(sum-probs[0]-probs[1]));
	ran[3]=random_binomial(N-ran[0]-ran[1]-ran[2],probs[3]/(sum-probs[0]-probs[1]-probs[2]));
	ran[4]=random_binomial(N-ran[0]-ran[1]-ran[2]-ran[3],probs[4]/(sum-probs[0]-probs[1]-probs[2]-probs[3]));
	ran[5]=N-ran[0]-ran[1]-ran[2]-ran[3]-ran[4];
	return ran;};



int* random_multinom_var(int N,int howmany,double *relprobs,double tot)
{
	int *ran=new int[howmany];
	double sum=tot;
	int Nused=N;
	for(int pat=0;pat<howmany;pat++)
	{
		if(N>0)
		{
		ran[pat]=random_binomial(Nused,double(*(relprobs+pat))/sum);
		sum-=double(*(relprobs+pat));
		Nused-=ran[pat];
		}
		else ran[pat]=0;
	};
	return ran;};




//Below is code for the random number generator---------------------------------------------------------------------------------------------------------------------------------



void CRandomMersenne::Init0(uint32 seed) {
   // Detect computer architecture
   union {float f; uint32 i[2];} convert;
   convert.f = 1.0;
   if (convert.i[1] == 0x3FF00000) Architecture = LITTLE_ENDIAN1;
   else if (convert.i[0] == 0x3FF00000) Architecture = BIG_ENDIAN1;
   else Architecture = NONIEEE;

   // Seed generator
   mt[0]= seed;
   for (mti=1; mti < MERS_N; mti++) {
      mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
   }
}

void CRandomMersenne::RandomInit(uint32 seed) {
   // Initialize and seed
   Init0(seed);

   // Randomize some more
   for (int i = 0; i < 37; i++) BRandom();
}

uint32 CRandomMersenne::BRandom() {// Generate 32 random bits

   uint32 y;

   if (mti >= MERS_N) {
      // Generate MERS_N words at one time
      const uint32 LOWER_MASK = (1LU << MERS_R) - 1;       // Lower MERS_R bits
      const uint32 UPPER_MASK = 0xFFFFFFFF << MERS_R;      // Upper (32 - MERS_R) bits
      static const uint32 mag01[2] = {0, MERS_A};

      int kk;
      for (kk=0; kk < MERS_N-MERS_M; kk++) {
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

      for (; kk < MERS_N-1; kk++) {
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}

      y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
      mti = 0;
   }

   y = mt[mti++];

#if 1
   // Tempering (May be omitted):
   y ^=  y >> MERS_U;
   y ^= (y << MERS_S) & MERS_B;
   y ^= (y << MERS_T) & MERS_C;
   y ^=  y >> MERS_L;
#endif

   return y;
}


double CRandomMersenne::Random() {
   // Output random float number in the interval 0 <= x < 1
   union {float f; uint32 i[2];} convert; //Union allows one portion of the memory to be accessed as different data types.
   double ra=1.1;
	while(ra>=1)
	{
   uint32 r = BRandom();               // Get 32 random bits
   // The fastest way to convert random bits to floating point is as follows:
   // Set the binary exponent of a floating point number to 1+bias and set
   // the mantissa to random bits. This will give a random number in the
   // interval [1,2). Then subtract 1.0 to get a random number in the interval
   // [0,1). This procedure requires that we know how floating point numbers
   // are stored. The storing method is tested in function RandomInit and saved
   // in the variable Architecture.

   // This shortcut allows the compiler to optimize away the following switch
   // statement for the most common architectures:
#if defined(_M_IX86) || defined(_M_X64) || defined(__LITTLE_ENDIAN__)
   Architecture = LITTLE_ENDIAN1;
#elif defined(__BIG_ENDIAN__)
   Architecture = BIG_ENDIAN1;
#endif

   switch (Architecture) {
   case LITTLE_ENDIAN1:
      convert.i[0] =  r << 20;
      convert.i[1] = (r >> 12) | 0x3FF00000;
      return convert.f - 1.0;
   case BIG_ENDIAN1:
      convert.i[1] =  r << 20;
      convert.i[0] = (r >> 12) | 0x3FF00000;
      return convert.f - 1.0;
   case NONIEEE: default: ;
   }
   // This somewhat slower method works for all architectures, including
   // non-IEEE floating point representation:
   //return 0.0000005 +0.999999*(double)r * (1./((double)(uint32)(-1L)+1.));
   ra= (double)r * (1./((double)(uint32)(-1L)+1.));
	};
   return ra;
}


int CRandomMersenne::IRandom(int min, int max) {
   // Output random integer in the interval min <= x <= max
   // Relative error on frequencies < 2^-32
   if (max <= min) {
      if (max == min) return min; else return 0x80000000;
   }
   // Multiply interval with random and truncate
   int r = int((max - min + 1) * Random()) + min;
   if (r > max) r = max;
   return r;
}


