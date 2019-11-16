clear all
ActiveTrialNumber=20;
%%%%%%% Speed(cm/s) = Speed * 0.4 ;
%%%%%%% 0.4 is the speed factor
Mouse=4;
NBin=50;
NShuff=100;
AddressSetup;


Days=[51:57];
areaNames={'V1','LV','MV','PPC','A','S','M','RSC','Border'};


Day=Days(Mouse);

path=LoadPath{Day};
load(strcat(path,'\All_Sessions.mat'));
load(strcat(path,'\cellData.mat'));

for k=1:size(HitSE,1)
    Hit(HitSE(k,1):HitSE(k,2))=1;
end
for k=1:size(MissSE,1)
    Miss(MissSE(k,1):MissSE(k,2))=1;
end
for k=1:size(CRSE,1)
    CR(CRSE(k,1):CRSE(k,2))=1;
end
for k=1:size(FASE,1)
    FA(FASE(k,1):FASE(k,2))=1;
end


ActiveAnimal=ones(length(Lick),1);
for i=1:(length(Lick)-75*ActiveTrialNumber)
    if max(Lick(i:(i+75*ActiveTrialNumber)))==0
        ActiveAnimal(i:(i+75*ActiveTrialNumber))=0;
    end
end

X=cellData_Calcium;
cellCount=size(X,1);

TimeVec=1:length(Speed);
     Speed=double(round(Speed*10))/10;
ValidTime=Speed>0 & ActiveAnimal==1;

[SQ,BQ]=SpeedQuantEqual(Speed,ValidTime,NBin);
BQ(1)=min(Speed(ValidTime));
BQ=BQ([diff(BQ)>0,1==1]);


AverageCell=zeros(size(X,1),length(BQ),NShuff);
AverageCellSpeed0=zeros(size(X,1),NShuff);
normalizebinTime=ActiveAnimal==1 & Speed==0;
for r=1:NShuff
    r
    
        HitTime=TimeVec(Hit==1 & normalizebinTime);
        HitTime=HitTime(randperm(length(HitTime)));
        MissTime=TimeVec(Miss==1 & normalizebinTime);
        MissTime=MissTime(randperm(length(MissTime)));
        CRTime=TimeVec(CR==1 & normalizebinTime);
        CRTime=CRTime(randperm(length(CRTime)));
        FATime=TimeVec(FA==1 & normalizebinTime);
        FATime=FATime(randperm(length(FATime)));
        minSize=min([length(HitTime),length(MissTime),length(CRTime),length(FATime)]);
        TimeBinMatrix=[HitTime(1:minSize),MissTime(1:minSize),CRTime(1:minSize),FATime(1:minSize)];
        for cnum=1:cellCount
            AverageCellSpeed0(cnum,r)=mean(X(cnum,TimeBinMatrix));
        end
    
    for i=1:length(BQ)
        if i<length(BQ)
            binTime=ValidTime & Speed>=BQ(i) & Speed<=BQ(i+1);
        else
            binTime=ValidTime & Speed>=BQ(i);
        end
        
        HitTime=TimeVec(Hit==1 & binTime);
        HitTime=HitTime(randperm(length(HitTime)));
        MissTime=TimeVec(Miss==1 & binTime);
        MissTime=MissTime(randperm(length(MissTime)));
        CRTime=TimeVec(CR==1 & binTime);
        CRTime=CRTime(randperm(length(CRTime)));
        FATime=TimeVec(FA==1 & binTime);
        FATime=FATime(randperm(length(FATime)));

               
        minSize=min([length(HitTime),length(MissTime),length(CRTime),length(FATime)]);
        TimeBinMatrix=[HitTime(1:minSize),MissTime(1:minSize),CRTime(1:minSize),FATime(1:minSize)];
   
        
        for cnum=1:cellCount
            AverageCell(cnum,i,r)=mean(X(cnum,TimeBinMatrix));
        end
        
    end
    
end
BQ(length(BQ)+1)=max(Speed(Speed<400 & ValidTime));
save(strcat('C:\Users\Sadegh\Documents\VLMReborn\Reports2\2017_10_09_SpeedCorrelation\2018_01_23_FiringRate_balancedHMCF\Mouse',num2str(Mouse)),'AverageCell','AverageCellSpeed0','BQ')

SpeedValues=zeros(1,length(BQ)-1);
for bin=1:length(BQ)-1
    SpeedValues(bin)=median(Speed(Speed>=BQ(bin) & Speed<BQ(bin+1)));
end
%SpeedValues=(BQ(1:length(BQ)-1)+BQ(2:length(BQ)))/2;
AverageResp=mean(AverageCell,3);
pvalue=0.01;
[RHO,PVAL] = corr(AverageResp',SpeedValues','type','Spearman');
SpeedValues=SpeedValues*0.4;
PosCorrMean=zeros(1,length(SpeedValues));
PosCorrStd=zeros(1,length(SpeedValues));
PosCorrSem=zeros(1,length(SpeedValues));

NegCorrMean=zeros(1,length(SpeedValues));
NegCorrStd=zeros(1,length(SpeedValues));
NegCorrSem=zeros(1,length(SpeedValues));

for bin=1:length(BQ)-1
    NormResp=AverageResp(:,bin)./ mean(AverageCellSpeed0,2);
    PosCorrMean(bin)=mean(NormResp(RHO>0 & PVAL<pvalue));
    NegCorrMean(bin)=mean(NormResp(RHO<0 & PVAL<pvalue));
    
    PosCorrStd(bin)=sqrt(var(NormResp(RHO>0 & PVAL<pvalue)));
    NegCorrStd(bin)=sqrt(var(NormResp(RHO<0 & PVAL<pvalue)));
    
    PosCorrSem(bin)=sqrt(var(NormResp(RHO>0 & PVAL<pvalue)))/sum(RHO>0 & PVAL<pvalue);
    NegCorrSem(bin)=sqrt(var(NormResp(RHO<0 & PVAL<pvalue)))/sum(RHO<0 & PVAL<pvalue);
end
    




myfun = 'y~(b1)/(b2+x)^b3';
%myfun = 'y~(b1)/(b2+exp(b3*x))';

mdl = fitnlm(SpeedValues,PosCorrMean,myfun,[2,1,2]);
beta=mdl.Coefficients.Estimate;
PosCorrMeanfit=beta(1)./((beta(2)+SpeedValues).^beta(3));

mdl = fitnlm(SpeedValues,NegCorrMean,myfun,[1,1,1]);
beta=mdl.Coefficients.Estimate;
NegCorrMeanfit=beta(1)./((beta(2)+SpeedValues).^beta(3));

Result=[SpeedValues;PosCorrMean;PosCorrMeanfit;PosCorrStd;PosCorrSem;...
    NegCorrMean;NegCorrMeanfit;NegCorrStd;NegCorrSem];























