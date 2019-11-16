%%%%%%%%%%%%%%%%% Validation Projection of HMCF on Hit/CR Decoder 
area=9;
ActiveMode=1;
SpeedMode=2;
for Mouse=[61,63:67]
    Mouse
    loadpath=LoadPath{Mouse-10};
    load(strcat('E:\Reports2\2018_03_02_FisherInstant\ROI_Bin\Datasets--5to0\HitCR_Val\Session',num2str(Mouse),'_mode1'));
    load(strcat(loadpath,'/cellData.mat'));
    load(strcat(loadpath,'/cellData_ZS.mat'));
    load(strcat(loadpath,'/SessionLength.mat'));
    load(strcat(loadpath,'/Datasets/Datasets--5to0.mat'));
    cellEvents=cellData_ROI_bin;

HitProj=zeros(max(HitTrialNumber{SpeedMode,1}),25,100);
CRProj=zeros(max(CRTrialNumber{SpeedMode,1,ActiveMode}),25,100);
MissProj=zeros(max(MissTrialNumber{SpeedMode,1,ActiveMode}),25,100);
FAProj=zeros(max(FATrialNumber{SpeedMode,1}),25,100);

HitMeanScores=zeros(max(HitTrialNumber{SpeedMode,1}),25);
CRMeanScores=zeros(max(CRTrialNumber{SpeedMode,1,ActiveMode}),25);
MissMeanScores=zeros(max(MissTrialNumber{SpeedMode,1,ActiveMode}),25);
FAMeanScores=zeros(max(FATrialNumber{SpeedMode,1}),25);

NoiseDprime=zeros(2,25);
NoiseDprime_sh=zeros(2,25,1000);

for nr=1:100
    for td=1:25
        [~,maxInd]=max(mean(Info_Val_Fisher(area,td,1:20,:),4));
        W=MaxDecoders{area,td,maxInd,nr};
        
        index=HitDataset{SpeedMode,1};
        TNum=HitTrialNumber{SpeedMode,1};
        for i=ValTrials1{td,nr}
            if sum(TNum==i)>=td
                HitProj(i,td,nr)= W' * cellEvents(:,min(index(TNum==i))+td-1);
            end
        end
        
        index=CRDataset{SpeedMode,1,ActiveMode};
        TNum=CRTrialNumber{SpeedMode,1,ActiveMode};
        for i=ValTrials0{td,nr}
            if sum(TNum==i)>=td
                CRProj(i,td,nr)= W' * cellEvents(:,min(index(TNum==i))+td-1);
            end
        end
        
        index=FADataset{SpeedMode,1};
        TNum=FATrialNumber{SpeedMode,1};
        for i=unique(TNum)
            if sum(TNum==i)>=td
                FAProj(i,td,nr)= W' * cellEvents(:,min(index(TNum==i))+td-1);
            end
        end
        
        index=MissDataset{SpeedMode,1,ActiveMode};
        TNum=MissTrialNumber{SpeedMode,1,ActiveMode};
        for i=unique(TNum)
            if sum(TNum==i)>=td
                MissProj(i,td,nr)= W' * cellEvents(:,min(index(TNum==i))+td-1);
            end
        end
        
        
        
    end
end

for td=1:25
    ns=1;
    for i=1:size(HitProj,1)
        if sum(HitProj(i,td,:)~=0)>1
        HitMeanScores(ns,td)=mean(HitProj(i,td,HitProj(i,td,:)~=0));
        ns=ns+1;
        end
    end
    HitMeanScores=HitMeanScores(1:ns-1,:);
    ns=1;
    for i=1:size(CRProj,1)
        if sum(CRProj(i,td,:)~=0)>1
        CRMeanScores(i,td)=mean(CRProj(i,td,CRProj(i,td,:)~=0));
        ns=ns+1;
        end
    end
    CRMeanScores=CRMeanScores(1:ns-1,:);
    
    FAMeanScores(:,td)=mean(FAProj(:,td,:),3);
    MissMeanScores(:,td)=mean(MissProj(:,td,:),3);
    NoiseDprime(1,td)=abs(mean(CRMeanScores(:,td))-mean(FAMeanScores(:,td)))/...
        sqrt(0.5*(var(CRMeanScores(:,td))+var(FAMeanScores(:,td))));
    NoiseDprime(2,td)=abs(mean(HitMeanScores(:,td))-mean(MissMeanScores(:,td)))/...
        sqrt(0.5*(var(HitMeanScores(:,td))+var(MissMeanScores(:,td))));
    
    for sh=1:1000
        Group=[zeros(length(MissMeanScores(:,td)),1);ones(length(HitMeanScores(:,td)),1)];
        randGroup=Group(randperm(length(Group)));
        scores=[MissMeanScores(:,td);HitMeanScores(:,td)];
                   NoiseDprime_sh(2,td,sh)=abs(mean(scores(randGroup==0))-mean(scores(randGroup==1)))/...
                            sqrt(0.5*(var(scores(randGroup==0))+var(scores(randGroup==1))));
                        
        Group=[zeros(length(CRMeanScores(:,td)),1);ones(length(FAMeanScores(:,td)),1)];
        randGroup=Group(randperm(length(Group)));
        scores=[CRMeanScores(:,td);FAMeanScores(:,td)];
                   NoiseDprime_sh(1,td,sh)=abs(mean(scores(randGroup==0))-mean(scores(randGroup==1)))/...
                            sqrt(0.5*(var(scores(randGroup==0))+var(scores(randGroup==1))));
                      
        
    end
end

save(strcat('E:\Reports2\2018_09_15_NoiseBehaviorRelation\InstantHC\Session',num2str(Mouse),'_mode1'),'NoiseDprime','NoiseDprime_sh','HitMeanScores','CRMeanScores','FAMeanScores','MissMeanScores'...
    ,'HitProj','CRProj','FAProj','MissProj','area','ActiveMode','SpeedMode','-v7.3');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DPMtx=zeros(2,25,6);
DPMtx_sh=zeros(2,25,6);
nM=1;

for Mouse=[61,63:67]
    load(strcat('E:\Reports2\2018_09_15_NoiseBehaviorRelation\InstantHC\Session',num2str(Mouse),'_mode1.mat'));
   
    for CH=1:2
        DPMtx(CH,:,nM)=NoiseDprime(CH,:);
       DPMtx_sh(CH,:,nM)=max(NoiseDprime_sh(CH,:,:),[],3); 
    end
    figure();hold on
    plot(NoiseDprime(1,:),'r');
    plot(NoiseDprime(2,:),'b');
    plot(DPMtx_sh(1,:,nM),'r--');
    plot(DPMtx_sh(2,:,nM),'b--');
    
    nM=nM+1;
end
figure();hold on
shadedErrorBar(-0.4:0.1:2,mean(DPMtx(1,:,:),3),sqrt(var(squeeze(DPMtx(1,:,:))')),'r',1);
shadedErrorBar(-0.4:0.1:2,mean(DPMtx(2,:,:),3),sqrt(var(squeeze(DPMtx(2,:,:))')),'b',1);
shadedErrorBar(-0.4:0.1:2,mean(DPMtx_sh(1,:,:),3),sqrt(var(squeeze(DPMtx_sh(1,:,:))')),'r--',1);
shadedErrorBar(-0.4:0.1:2,mean(DPMtx_sh(2,:,:),3),sqrt(var(squeeze(DPMtx_sh(2,:,:))')),'b--',1);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Projection of HMCF on CCA Modes
Areas={'V1','LV','MV','PPC','A','S','M','RSC'};
ActiveMode=1;
SpeedMode=2;

for Mouse=[61,63:67]
    Mouse
    ProjectionsHit={};
    ProjectionsMiss={};
    ProjectionsCR={};
    ProjectionsFA={};
    BehaviorCorrelations=zeros(8,8,25,2,20);
    BehaviorCorrelations_sh=zeros(8,8,25,2,20,1000);
    BehaviorDprime=zeros(8,8,25,2,20);
    BehaviorDprime_sh=zeros(8,8,25,2,20,1000);
    loadpath=LoadPath{Mouse-10};
    load(strcat('E:\Reports2\2018_03_14_CannonicalCorrelations\ROI_Bin\Datasets--5to0\HitCR_L2_MixedStim\Session',num2str(Mouse),'_mode1'));
    load(strcat(loadpath,'/cellData.mat'));
    load(strcat(loadpath,'/cellData_ZS.mat'));
    load(strcat(loadpath,'/SessionLength.mat'));
    load(strcat(loadpath,'/Datasets/Datasets--5to0.mat'));
    cellEvents=cellData_ROI_bin;
    MW1={};MW2={};
    for a1=1:8%x
        a1
        for a2=a1+1:8%y
            a2
            if sum(CortexArea==a1)>5 && sum(CortexArea==a2)>5

                W1=zeros(size(Wx0{a1,a2,1},1),20,100);
                W2=zeros(size(Wy0{a1,a2,1},1),20,100);
                for nr=1:100
                    W1(:,:,nr)= (Wx0{a1,a2,nr}.* repmat(sign(mean(Wx0{a1,a2,nr})),[size(Wx0{a1,a2,nr},1),1]));
                    W2(:,:,nr)= (Wy0{a1,a2,nr}.* repmat(sign(mean(Wy0{a1,a2,nr})),[size(Wy0{a1,a2,nr},1),1]));
                end
                
                MW1{a1,a2}=squeeze(mean(W1,3));
                MW2{a1,a2}=squeeze(mean(W2,3));
                
                
                HitTrialMtx=zeros(1,length(CortexArea),25);
                n=1;
                index=HitDataset{SpeedMode,1};
                TNum=HitTrialNumber{SpeedMode,1};
                
                for i=unique(TNum)
                    if sum(TNum==i)==25
                        HitTrialMtx(n,:,1:25)= cellEvents(:,index(TNum==i));
                        n=n+1;
                    end
                end
                
                
                MissTrialMtx=zeros(1,length(CortexArea),25);
                n=1;
                index=MissDataset{SpeedMode,1,ActiveMode};
                TNum=MissTrialNumber{SpeedMode,1,ActiveMode};
                
                for i=unique(TNum)
                    if sum(TNum==i)==25
                        MissTrialMtx(n,:,1:25)= cellEvents(:,index(TNum==i));
                        n=n+1;
                    end
                end
                
                CRTrialMtx=zeros(1,length(CortexArea),25);
                n=1;
                index=CRDataset{SpeedMode,1,ActiveMode};
                TNum=CRTrialNumber{SpeedMode,1,ActiveMode};
                
                for i=unique(TNum)
                    if sum(TNum==i)==25
                        CRTrialMtx(n,:,1:25)= cellEvents(:,index(TNum==i));
                        n=n+1;
                    end
                end
                
                FATrialMtx=zeros(1,length(CortexArea),25);
                n=1;
                index=FADataset{SpeedMode,1};
                TNum=FATrialNumber{SpeedMode,1};
                
                for i=unique(TNum)
                    if sum(TNum==i)==25
                        FATrialMtx(n,:,1:25)= cellEvents(:,index(TNum==i));
                        n=n+1;
                    end
                end
                
                
                
                
                for ccnum=1:20
                    ccnum;

                    D1=squeeze(MW1{a1,a2}(:,ccnum)/norm(MW1{a1,a2}(:,ccnum)));
                    D2=squeeze(MW2{a1,a2}(:,ccnum)/norm(MW2{a1,a2}(:,ccnum)));
                    
                    
                    HitProj=zeros(size(HitTrialMtx,1),25,2);
                    CRProj=zeros(size(CRTrialMtx,1),25,2);
                    MissProj=zeros(size(MissTrialMtx,1),25,2);
                    FAProj=zeros(size(FATrialMtx,1),25,2);
                    for td=1:25
                        

                        HitProj(:,td,1)=D1'*squeeze(HitTrialMtx(:,CortexArea==a1,td))';
                        HitProj(:,td,2)=D2'*squeeze(HitTrialMtx(:,CortexArea==a2,td))';
                        
                        
                        MissProj(:,td,1)=D1'*squeeze(MissTrialMtx(:,CortexArea==a1,td))';
                        MissProj(:,td,2)=D2'*squeeze(MissTrialMtx(:,CortexArea==a2,td))';
                        
                        CRProj(:,td,1)=D1'*squeeze(CRTrialMtx(:,CortexArea==a1,td))';
                        CRProj(:,td,2)=D2'*squeeze(CRTrialMtx(:,CortexArea==a2,td))';

                        FAProj(:,td,1)=D1'*squeeze(FATrialMtx(:,CortexArea==a1,td))';
                        FAProj(:,td,2)=D2'*squeeze(FATrialMtx(:,CortexArea==a2,td))';
                       
                        
                        
                        
                
                    end
                    ProjectionsHit{a1,a2,ccnum}=HitProj;
                    ProjectionsMiss{a1,a2,ccnum}=MissProj;
                    ProjectionsCR{a1,a2,ccnum}=CRProj;
                    ProjectionsFA{a1,a2,ccnum}=FAProj;
                    for td=1:25
                        Group0=[zeros(length(ProjectionsCR{a1,a2,ccnum}(:,td,1)),1);ones(length(ProjectionsFA{a1,a2,ccnum}(:,td,1)),1)];
                        BehaviorCorrelations(a1,a2,td,1,ccnum)=similarity([[ProjectionsCR{a1,a2,ccnum}(:,td,1);ProjectionsFA{a1,a2,ccnum}(:,td,1)],Group0]);
                        BehaviorDprime(a1,a2,td,1,ccnum)=(mean(ProjectionsCR{a1,a2,ccnum}(:,td,1))-mean(ProjectionsFA{a1,a2,ccnum}(:,td,1)))/...
                            sqrt(0.5*(var(ProjectionsCR{a1,a2,ccnum}(:,td,1))+var(ProjectionsFA{a1,a2,ccnum}(:,td,1))));
                        BehaviorCorrelations(a2,a1,td,1,ccnum)=similarity([[ProjectionsCR{a1,a2,ccnum}(:,td,2);ProjectionsFA{a1,a2,ccnum}(:,td,2)],Group0]);
                        BehaviorDprime(a2,a1,td,1,ccnum)=(mean(ProjectionsCR{a1,a2,ccnum}(:,td,2))-mean(ProjectionsFA{a1,a2,ccnum}(:,td,2)))/...
                            sqrt(0.5*(var(ProjectionsCR{a1,a2,ccnum}(:,td,2))+var(ProjectionsFA{a1,a2,ccnum}(:,td,2))));
                      
                        Group1=[zeros(length(ProjectionsMiss{a1,a2,ccnum}(:,td,1)),1);ones(length(ProjectionsHit{a1,a2,ccnum}(:,td,1)),1)];
                        BehaviorCorrelations(a1,a2,td,2,ccnum)=similarity([[ProjectionsMiss{a1,a2,ccnum}(:,td,1);ProjectionsHit{a1,a2,ccnum}(:,td,1)],Group1]);
                        BehaviorDprime(a1,a2,td,2,ccnum)=(mean(ProjectionsMiss{a1,a2,ccnum}(:,td,1))-mean(ProjectionsHit{a1,a2,ccnum}(:,td,1)))/...
                            sqrt(0.5*(var(ProjectionsMiss{a1,a2,ccnum}(:,td,1))+var(ProjectionsHit{a1,a2,ccnum}(:,td,1))));
                        BehaviorCorrelations(a2,a1,td,2,ccnum)=similarity([[ProjectionsMiss{a1,a2,ccnum}(:,td,2);ProjectionsHit{a1,a2,ccnum}(:,td,2)],Group1]);
                        BehaviorDprime(a2,a1,td,2,ccnum)=(mean(ProjectionsMiss{a1,a2,ccnum}(:,td,2))-mean(ProjectionsHit{a1,a2,ccnum}(:,td,2)))/...
                            sqrt(0.5*(var(ProjectionsMiss{a1,a2,ccnum}(:,td,2))+var(ProjectionsHit{a1,a2,ccnum}(:,td,2))));
                      for sh=1:1000
                        randGroup=Group0(randperm(length(ProjectionsCR{a1,a2,ccnum}(:,td,1))+length(ProjectionsFA{a1,a2,ccnum}(:,td,1))));
                        scores=[ProjectionsCR{a1,a2,ccnum}(:,td,1);ProjectionsFA{a1,a2,ccnum}(:,td,1)];
                        BehaviorCorrelations_sh(a1,a2,td,1,ccnum,sh)=similarity([[ProjectionsCR{a1,a2,ccnum}(:,td,1);ProjectionsFA{a1,a2,ccnum}(:,td,1)],randGroup]);
                        BehaviorDprime_sh(a1,a2,td,1,ccnum,sh)=(mean(scores(randGroup==0))-mean(scores(randGroup==1)))/...
                            sqrt(0.5*(var(scores(randGroup==0))+var(scores(randGroup==1))));
                        
                        scores=[ProjectionsCR{a1,a2,ccnum}(:,td,2);ProjectionsFA{a1,a2,ccnum}(:,td,2)];
                        BehaviorCorrelations_sh(a2,a1,td,1,ccnum,sh)=similarity([[ProjectionsCR{a1,a2,ccnum}(:,td,2);ProjectionsFA{a1,a2,ccnum}(:,td,2)],randGroup]);
                        BehaviorDprime_sh(a2,a1,td,1,ccnum,sh)=(mean(scores(randGroup==0))-mean(scores(randGroup==1)))/...
                            sqrt(0.5*(var(scores(randGroup==0))+var(scores(randGroup==1))));
                        
                    
                        randGroup=Group1(randperm(length(ProjectionsMiss{a1,a2,ccnum}(:,td,1))+length(ProjectionsHit{a1,a2,ccnum}(:,td,1))));
                        scores=[ProjectionsMiss{a1,a2,ccnum}(:,td,1);ProjectionsHit{a1,a2,ccnum}(:,td,1)];
                        BehaviorCorrelations_sh(a1,a2,td,2,ccnum,sh)=similarity([[ProjectionsMiss{a1,a2,ccnum}(:,td,1);ProjectionsHit{a1,a2,ccnum}(:,td,1)],randGroup]);
                        BehaviorDprime_sh(a1,a2,td,2,ccnum,sh)=(mean(scores(randGroup==0))-mean(scores(randGroup==1)))/...
                            sqrt(0.5*(var(scores(randGroup==0))+var(scores(randGroup==1))));
                        
                        scores=[ProjectionsMiss{a1,a2,ccnum}(:,td,2);ProjectionsHit{a1,a2,ccnum}(:,td,2)];
                        BehaviorCorrelations_sh(a2,a1,td,2,ccnum,sh)=similarity([[ProjectionsMiss{a1,a2,ccnum}(:,td,2);ProjectionsHit{a1,a2,ccnum}(:,td,2)],randGroup]);
                        BehaviorDprime_sh(a2,a1,td,2,ccnum,sh)=(mean(scores(randGroup==0))-mean(scores(randGroup==1)))/...
                            sqrt(0.5*(var(scores(randGroup==0))+var(scores(randGroup==1))));
                      end
                    end
                    
                   
                end
               
            end

        end
    end
    save(strcat('E:\Reports2\2018_09_15_NoiseBehaviorRelation\CCAModes\Session',num2str(Mouse),'_mode1'),'BehaviorCorrelations','BehaviorCorrelations_sh',...
        'BehaviorDprime_sh','BehaviorDprime','ProjectionsCR','ProjectionsFA','ProjectionsHit','ProjectionsMiss','-v7.3');
end

%%%%%%%%%%%%%%%%%%%%

CorrMtx=zeros(8*25,8*20,2,6);
CorrMtx_sh=zeros(8*25,8*20,2,6);
nM=1;
for Mouse=[61,63:67]
    load(strcat('E:\Reports2\2018_09_15_NoiseBehaviorRelation\CCAModes\Session',num2str(Mouse),'_mode1.mat'));
   
    for CH=1:2
        test=abs(BehaviorDprime(:,:,:,CH,:));
        maxSh=max(abs(BehaviorDprime_sh(:,:,:,CH,:,:)),[],6);
        test(test<maxSh)=0;
        test=abs(squeeze(test));
        maxSh=abs(squeeze(maxSh));
        test=permute(test,[1,3,2,4]);
        maxSh=permute(maxSh,[1,3,2,4]);
        CorrMtx(:,:,CH,nM)=reshape(test,[8*25,8*20]);
        CorrMtx_sh(:,:,CH,nM)=reshape(maxSh,[8*25,8*20]);
    end
    
%     figure();
% realMean=[squeeze(CorrMtx(:,:,1,nM)),squeeze(CorrMtx(:,:,2,nM))];
% ShuffMean=[squeeze(CorrMtx_sh(:,:,1,nM)),squeeze(CorrMtx_sh(:,:,2,nM))];
% realMean(realMean<ShuffMean)=0;
% imagesc(realMean)
    
    
    nM=nM+1;
end
figure();
set(0,'DefaultAxesFontSize',14);
realMean=[squeeze(mean(CorrMtx(:,:,1,:),4)),zeros(200,10),squeeze(mean(CorrMtx(:,:,2,:),4))];
 %ShuffMean=[squeeze(mean(CorrMtx_sh(:,:,1,:),4)),squeeze(mean(CorrMtx_sh(:,:,2,:),4))];
 %realMean(realMean<ShuffMean)=0;
imagesc(realMean)


















