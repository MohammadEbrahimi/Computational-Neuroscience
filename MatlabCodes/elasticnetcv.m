function [B,intercept,ErrCurve,ErrVar,Lambda]=elasticnetcv(Set,group,trialNumber,crossVal,Repeat,doShuff,elasticRatio)
% Set=DataSet;
% trialNumber=shuffInd;
% crossVal=5;
% Repeat=10;
% doShuff=1;

dim=size(Set);
% Repeat=10;
% trialLength=20;
% doShuff=1;


MaxLambda=1;%1
MinLambda=1e-6;%6
l=30;



ILambda=zeros(1,l);
for i=2:l
   ILambda(i)=MinLambda * ( (MaxLambda/MinLambda)^((i-2)/(l-2))); 
end
error_FalseAlarm=zeros(crossVal,l);
error_Miss=zeros(crossVal,l);



FLambda=zeros(l,1);
FCV=zeros(1,1);
FMinIndx=zeros(1,1);
FDivIndx=zeros(1,1);


FalseAlarmCurve=zeros(l,1);
MissCurve=zeros(l,1);
Nr=Repeat;

ER=zeros(Nr,l);
ERVAR=zeros(1,l);
numberOfTrials=max(trialNumber);
timeIndex=1:dim(2);

for rep=1:Nr
    if doShuff==1
        sh=randperm(numberOfTrials);
        orderIndex=zeros(1,dim(2));
        n=0;
        for kk=1:numberOfTrials
            trialLength=sum(trialNumber == kk);
            if trialLength > 0
                orderIndex(n+1:n+trialLength)=timeIndex(trialNumber==kk);
                n=n+trialLength;
            end
        end
    else
        orderIndex=1:dim(2);
    end

    for iii=1:crossVal
        iii;
        %%%% cross validation
        [TrainSet TestSet TrainClass TestClass]=CVDataset(Set(:,orderIndex)',group(orderIndex)',iii,crossVal);
        
        %%%% lassogl input is N x P
        [B,FitInfo]=lassoglm(TrainSet,TrainClass,'binomial','Lambda',ILambda,'RelTol',0.001,'Alpha',elasticRatio);
        B1=[FitInfo.Intercept ; B];
        
        for iter=1:length(FitInfo.Intercept)
            
            dimTrain=size(TrainSet);
            dimTest=size(TestSet);
            TrainPredictClass =LRClassify([ones(dimTrain(1),1) TrainSet],B1(:,iter),0);
            TestPredictClass =LRClassify([ones(dimTest(1),1) TestSet],B1(:,iter),0);
            
            
            
            
            error_FalseAlarm(iii,iter)=sum(abs(TestPredictClass-TestClass)+TestPredictClass-TestClass)/(2*(length(TestClass)-sum(TestClass)));
            error_Miss(iii,iter)=sum(abs(TestPredictClass-TestClass)-TestPredictClass+TestClass)/(2*(sum(TestClass)));
            
            Terror_FalseAlarm(iii,iter)=sum(abs(TrainPredictClass-TrainClass)+TrainPredictClass-TrainClass)/(2*(length(TrainClass)-sum(TrainClass)));
            Terror_Miss(iii,iter)=sum(abs(TrainPredictClass-TrainClass)-TrainPredictClass+TrainClass)/(2*(sum(TrainClass)));
        end
    end
    
    FalseAlarmCurve=sum(error_FalseAlarm)/(crossVal);
    MissCurve=sum(error_Miss)/(crossVal);
    
    ER(rep,:)=0.5*(FalseAlarmCurve+MissCurve);
 rep
end
for i=1:size(ER,2)
    ERVAR(i)=var(ER(:,i));
end

if Nr>1
    ErrCurve=sum(ER)/Nr;
else
    ErrCurve=ER;
end

[~,IndexMin]=min(sum(ER)/Nr);
Lambda=ILambda(IndexMin);
ErrVar=ERVAR(IndexMin);

[B,FitInfo]=lassoglm(Set',group,'binomial','Lambda',Lambda);
intercept=FitInfo.Intercept ;


end


