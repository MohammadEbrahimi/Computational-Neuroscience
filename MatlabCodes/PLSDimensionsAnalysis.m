AddressSetup;
Mouse=54;
area=9;
%load(strcat('C:\Users\Sadegh\Documents\VLMReborn\Reports2\2018_01_26_NoiseCorrelation\Raw\Datasets-10to5_3\HitCR\Session', ...
 %   num2str(Mouse),'_mode1.mat'));
load(strcat(LoadPath{Mouse},'\cellData.mat'));
load(strcat(LoadPath{Mouse},'\All_Sessions.mat'));
load(strcat(LoadPath{Mouse},'\SessionLength.mat'));
load(strcat(LoadPath{Mouse},'\Datasets\Datasets-5to0.mat'));

[~,maxDim]=max(mean(Info_Val_Fisher(area,:,:),3));
W=FisherTuning{area,100,maxDim};
PLS=PLSRot{area,maxDim};
Score=repmat(W,[1,size(cellData_Raw,2)]).*(PLS' * cellData_Raw(CortexArea<=area,:));

for i=1:maxDim
    PLS(:,i)=PLS(:,i) / sqrt(PLS(:,i)'*PLS(:,i));
end
PLSProj=PLS' * cellData_Raw(CortexArea<=area,:);

nH=max(HitTrialNumber{1,1});
PLSH=zeros(maxDim,nH);
ScoreH=zeros(maxDim,nH);
TimeH=zeros(1,nH);
for i=1:nH
    index=HitDataset{1,1}(HitTrialNumber{1,1}==i);
    PLSH(:,i)=mean(PLSProj(:,index),2);
    ScoreH(:,i)=mean(Score(:,index),2);
    TimeH(i)=index(1)/10;
end



nC=max(CRTrialNumber{1,1});
PLSC=zeros(maxDim,nC);
ScoreC=zeros(maxDim,nC);
TimeC=zeros(1,nC);
for i=1:nC
    index=CRDataset{1,1}(CRTrialNumber{1,1}==i);
    PLSC(:,i)=mean(PLSProj(:,index),2);
    ScoreC(:,i)=mean(Score(:,index),2);
    TimeC(i)=index(1)/10;
end

nM=max(MissTrialNumber{1,1});
PLSM=zeros(maxDim,nM);
ScoreM=zeros(maxDim,nM);
TimeM=zeros(1,nM);
for i=1:nM
    index=MissDataset{1,1}(MissTrialNumber{1,1}==i);
    PLSM(:,i)=mean(PLSProj(:,index),2);
    ScoreM(:,i)=mean(Score(:,index),2);
    TimeM(i)=index(1)/10;
end


nF=max(FATrialNumber{1,1});
PLSF=zeros(maxDim,nF);
ScoreF=zeros(maxDim,nF);
TimeF=zeros(1,nF);
for i=1:nF
    index=FADataset{1,1}(FATrialNumber{1,1}==i);
    PLSF(:,i)=mean(PLSProj(:,index),2);
    ScoreF(:,i)=mean(Score(:,index),2);
    TimeF(i)=index(1)/10;
end

dim=3;
figure();
for i=1:dim
    subplot(1,dim,i)
    hold on;plot(TimeH,ScoreH(i,:)','b.');plot(TimeC,ScoreC(i,:)','r.');
    %plot(TimeM,ScoreM(i,:)','m.');plot(TimeF,ScoreF(i,:)','k.');
    for s=1:length(SessionLength)
        plot(sum(SessionLength(1:s))*[0.1,0.1],[-50,50],'k--')
    end
end

figure();hold on;
plot(TimeH,sum(ScoreH(1:dim,:))','b.');plot(TimeC,sum(ScoreC(1:dim,:))','r.');
tit='Hit/FA:   ';
 for s=1:length(SessionLength)
        NFA=sum(FASE(:,1)> sum(SessionLength(1:s-1)) & FASE(:,1)< sum(SessionLength(1:s)));
        NCR=sum(CRSE(:,1)> sum(SessionLength(1:s-1)) & CRSE(:,1)< sum(SessionLength(1:s)));
        NHit=sum(HitSE(:,1)> sum(SessionLength(1:s-1)) & HitSE(:,1)< sum(SessionLength(1:s)));
        NMiss=sum(MissSE(:,1)> sum(SessionLength(1:s-1)) & MissSE(:,1)< sum(SessionLength(1:s)));
        tit=strcat(tit,'  ',num2str(round(100*NHit/(NHit+NMiss))),'/',num2str(round(100*NFA/(NFA+NCR))),'  ::');
        plot(sum(SessionLength(1:s))*[0.1,0.1],[-50,50],'k--')
 end
 title(tit);
 

