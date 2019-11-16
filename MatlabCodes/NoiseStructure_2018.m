AddressSetup;
for Mouse=[61,63:67];
    Mouse
S=load(strcat('E:\Reports2\2018_03_02_FisherInstant\ROI_Bin\Datasets--5to0\HitCR\Session',num2str(Mouse),'_mode1.mat'));
D=load(strcat('E:\Reports2\2018_03_02_FisherInstant\ROI_Bin\Datasets--5to0\HitCR\Session',num2str(Mouse),'_mode2.mat'));
R=load(strcat('E:\Reports2\2018_03_02_FisherInstant\ROI_Bin\Datasets--5to0\HitCR\Session',num2str(Mouse),'_mode3.mat'));
SC=load(strcat('E:\Reports2\2018_03_01_FisherInfo\ROI_Bin\Datasets--5to0\HitCR\Session',num2str(Mouse),'_mode1.mat'));
DC=load(strcat('E:\Reports2\2018_03_01_FisherInfo\ROI_Bin\Datasets--5to0\HitCR\Session',num2str(Mouse),'_mode2.mat'));
RC=load(strcat('E:\Reports2\2018_03_01_FisherInfo\ROI_Bin\Datasets--5to0\HitCR\Session',num2str(Mouse),'_mode3.mat'));
loadpath=LoadPath{Mouse-10};
SET0Name='C';
SET1Name='H';
CellT='IB';

SpeedMode=1;
Balanced=0;
LowDimVec=1:50;
NRepeate=1;
ActiveTrialNumber=20;
SpeedBalance=1;
ActiveMode=1;


warning('off','all');
load(strcat(loadpath,'/All_Sessions.mat'));
load(strcat(loadpath,'/cellData.mat'));
load(strcat(loadpath,'/cellData_ZS.mat'));
load(strcat(loadpath,'/SessionLength.mat'));
load(strcat(loadpath,'/Datasets/Datasets--5to0.mat'));
Speed=Speed-min(Speed);
progress='Data is loaded!'


clear cellData_BaseLine
clear cellData_Calcium
clear cellData_Noise
clear cellData_ROI_PerSession
clear cellData_Raw_PerSession
clear cellData_Raw
clear cellData_ROI
clear cellData_Baseline



switch CellT
    case 'R'
        cellEvents=cellData_Raw;
    case 'I'
        cellEvents=cellData_ROI;
    case 'O'
        load(strcat(loadpath,'/cellData_oopsi.mat'));
        cellEvents=cellData_oopsi;
    case 'RB'
        cellEvents=cellData_Raw_bin;
    case 'IB'
        cellEvents=cellData_ROI_bin;
    case 'RRT'
        cellEvents=cellData_Raw;
        cellEvents(cellData_Raw_bin==0)=0;
    case 'IRT'
        cellEvents=cellData_ROI;
        cellEvents(cellData_ROI_bin==0)=0;
        
end

cellCount=size(cellEvents,1);
CMtx=zeros(cellCount,cellCount,60);
CMtx_Sh=zeros(cellCount,cellCount);
MMtx=zeros(cellCount,2,60);


meanCorrG_Sh=zeros(9,9,60,100);
meanCorrNG_Sh=zeros(9,9,60,100);
meanCorrNGG_Sh=zeros(9,9,60,100);
meanCorrNN_Sh=zeros(9,9,60,100);
clear W
clear WC

W{9}=[];
WC{9}=[];
for mode=1:3
    %%% creating data sets
    
    TrialType0=zeros(1,max(CRTrialNumber{SpeedMode,mode,ActiveMode}));
    DataIndex0=CRDataset{SpeedMode,mode,ActiveMode};
    DataNumber0=CRTrialNumber{SpeedMode,mode,ActiveMode};
    
    TrialType1=zeros(1,max(HitTrialNumber{SpeedMode,mode}));
    DataIndex1=HitDataset{SpeedMode,mode};
    DataNumber1=HitTrialNumber{SpeedMode,mode};
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Speed Balance
    NBin=50;
    ActiveAnimal=ones(length(Lick),1);
    for i=1:(length(Lick)-75*ActiveTrialNumber)
        if max(Lick(i:(i+75*ActiveTrialNumber)))==0
            ActiveAnimal(i:(i+75*ActiveTrialNumber))=0;
        end
    end
    
    Speed=double(round(Speed*10))/10;
    [SQ,BQ]=SpeedQuantEqual(Speed,(Speed>0 & ActiveAnimal==1),NBin);
    
    Buf_DataNumber0=DataNumber0;
    Buf_DataNumber1=DataNumber1;
    
    
    shuff0=zeros(max(DataNumber0),NRepeate);
    shuff1=zeros(max(DataNumber1),NRepeate);
    
    nr=1;
    DataNumber0=Buf_DataNumber0;
    DataNumber1=Buf_DataNumber1;
    
    Sorder0=unique(DataNumber0);
    if Sorder0(1)==0 Sorder0=Sorder0(2:end);end
    Sorder0=Sorder0(randperm(length(Sorder0)));
    Sorder1=unique(DataNumber1);
    if Sorder1(1)==0 Sorder1=Sorder1(2:end);end
    Sorder1=Sorder1(randperm(length(Sorder1)));
    
    UMask0=zeros(size(DataIndex0));
    UMask1=zeros(size(DataIndex1));
    BalSign0=0;
    BalSign1=0;
    NotFinished=1;
    
    if SpeedBalance==1
        while (NotFinished==1)
            NotFinished=0;
            for k0=Sorder0
                if (TrialType0(k0) * BalSign0) <=0 && max(UMask0(DataNumber0==k0))==0
                    matchable=0;
                    Sbin0=max(SQ(DataIndex0(DataNumber0==k0)));
                    for k1=Sorder1
                        if max(UMask1(DataNumber1==k1))==0
                            Sbin1=max(SQ(DataIndex1(DataNumber1==k1)));
                            if Sbin0==Sbin1
                                matchable=1;
                                if (TrialType1(k1) * BalSign1) <=0
                                    UMask0(DataNumber0==k0)=1;
                                    UMask1(DataNumber1==k1)=1;
                                    BalSign0=BalSign0+TrialType0(k0);
                                    BalSign1=BalSign1+TrialType1(k1);
                                    NotFinished=1;
                                    break;
                                end
                            end
                        end
                    end
                    
                    if matchable==0
                        UMask0(DataNumber0==k0)=2;
                    end
                end
            end
            
        end
        DataNumber0(UMask0~=1)=0;
        DataNumber1(UMask1~=1)=0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%Data Ready
    %%%%%%%%%%Prepare Decoders
    
    rep=100;
    
    
    for area=1:9
        if mode==1
            TimeVec=1:25;
            OS=0;
            for itd=TimeVec
                SD=S.MaxDecoders{area,itd,1}/rep;
                SDC=SC.MaxDecoders{area,1}/rep;
                for r=2:rep
                    if length(S.MaxDecoders{area,itd,r})>0
                        SD=SD+(S.MaxDecoders{area,itd,r}/rep);
                        SDC=SDC+(SC.MaxDecoders{area,r}/rep);
                    end
                end
                W{area}=[W{area},(SD/norm(SD))];
                WC{area}=[WC{area},(SDC/norm(SDC))];
            end
        elseif mode==2
            TimeVec=1:5;
            OS=25;
            for itd=TimeVec
                SD=D.MaxDecoders{area,itd,1}/rep;
                SDC=DC.MaxDecoders{area,1}/rep;
                for r=2:rep
                    if length(D.MaxDecoders{area,itd,r})>0
                        SD=SD+(D.MaxDecoders{area,itd,r}/rep);
                        SDC=SDC+(DC.MaxDecoders{area,r}/rep);
                    end
                 
                end
                W{area}=[W{area},(SD/norm(SD))];
                WC{area}=[WC{area},(SDC/norm(SDC))];
            end
        elseif mode==3
             TimeVec=1:30;
             OS=30;
            for itd=TimeVec
                SD=R.MaxDecoders{area,itd,1}/rep;
                SDC=RC.MaxDecoders{area,1}/rep;
                for r=2:rep
                    if length(R.MaxDecoders{area,itd,r})>0
                        SD=SD+(R.MaxDecoders{area,itd,r}/rep);
                        SDC=SDC+(RC.MaxDecoders{area,r}/rep);
                    end
                end
                W{area}=[W{area},(SD/norm(SD))];
                WC{area}=[WC{area},(SDC/norm(SDC))];
            end
            
        end
    end
    
   % W=WC;
    
    
    
    Sorder0=unique(DataNumber0);
    if Sorder0(1)==0 Sorder0=Sorder0(2:end);end
    Sorder0=Sorder0(randperm(length(Sorder0)));
    Sorder1=unique(DataNumber1);
    if Sorder1(1)==0 Sorder1=Sorder1(2:end);end
    Sorder1=Sorder1(randperm(length(Sorder1)));
    
    
    for tind=TimeVec
        
        
        c0=1;
        c1=1;
        dd0=zeros(size(cellEvents,1),length(DataIndex0));
        for si=Sorder0
            trialLength=sum(DataNumber0==si);
            if(trialLength>=tind)
                dd0(:,c0)=cellEvents(:,min(DataIndex0(DataNumber0==si))+tind-1);
                c0=c0+1;
            end
        end
        c0=c0-1;
        dd0=dd0(:,1:c0);
        
        dd1=zeros(size(cellEvents,1),length(DataIndex1));
        for si=Sorder1
            trialLength=sum(DataNumber1==si);
            if(trialLength>=tind)
                dd1(:,c1)=cellEvents(:,min(DataIndex1(DataNumber1==si))+tind-1);
                c1=c1+1;
            end
        end
        c1=c1-1;
        dd1=dd1(:,1:c1);
        
        MMtx(:,1,tind+OS)=mean(dd0,2);
        MMtx(:,2,tind+OS)=mean(dd1,2);
        
        
        
        dd0_Sh=dd0;
        dd1_Sh=dd1;
        
%         dd0=dd0-repmat(mean(dd0,2),[1,c0]);
%         dd1=dd1-repmat(mean(dd1,2),[1,c1]);
%         dd=[dd0,dd1];
%         CMtx(:,:,tind+OS)=cov(dd');
        
        CMtx(:,:,tind+OS)=0.5*(cov(dd0')+cov(dd1'));
        
        
        %%%%%%%%%%

        
     
        
%         
%         for rsh=1%:100
%             dd0_Sh=dd0_Sh-repmat(mean(dd0_Sh,2),[1,c0]);
%             dd1_Sh=dd1_Sh-repmat(mean(dd1_Sh,2),[1,c1]);
%             dd_Sh=[dd0_Sh,dd1_Sh];
%             for cn=1:cellCount
%                 dd_Sh(cn,:)=dd_Sh(cn,randperm(c0+c1));
%             end
%             CMtx_Sh(:,:)=cov(dd_Sh');
%             M=abs(CMtx_Sh(:,:));
%             M(eye(cellCount)==1)=0;
%             
%             
%             for a1=1:8
%                 for a2=1:8
%                     B1=zeros(cellCount,1);
%                     B2=zeros(cellCount,1);
%                     if length(W{a1}) >0  && length(W{a2}) >0 
%                     B1(CortexArea==a1)=(W{a1}(:,tind+OS));
%                     B2(CortexArea==a2)=(W{a2}(:,tind+OS));
% 
%                     th1=2*sqrt(var(W{a1}(:,tind+OS)));
%                     th2=2*sqrt(var(W{a2}(:,tind+OS)));
%                     meanCorrG_Sh(a1,a2,tind+OS,rsh)=mean(mean(M(B1>th1,B2>th2),2));
%                     meanCorrNG_Sh(a1,a2,tind+OS,rsh)=mean(mean(M(B1<(-1*th1),B2<(-1*th2)),2));
%                     meanCorrNGG_Sh(a1,a2,tind+OS,rsh)=mean(mean(M(B1>th1,B2<(-1*th2)),2));
%                     meanCorrNN_Sh(a1,a2,tind+OS,rsh)=mean(mean(M(abs(B1)<th1 & abs(B1)>0 ,abs(B2)<(th2) & abs(B2)>0),2));
%                     end
%                 end
%             end
%             a1=9;a2=9;
%             th1=2*sqrt(var(W{a1}(:,tind+OS)));
%             th2=2*sqrt(var(W{a2}(:,tind+OS)));
%             meanCorrG_Sh(a1,a2,tind+OS,rsh)=mean(mean(M(W{a1}(:,tind+OS)>th1,W{a2}(:,tind+OS)>th2),2));
%             meanCorrNG_Sh(a1,a2,tind+OS,rsh)=mean(mean(M(W{a1}(:,tind+OS)<(-1*th1),W{a2}(:,tind+OS)<(-1*th2)),2));
%             meanCorrNGG_Sh(a1,a2,tind+OS,rsh)=mean(mean(M(W{a1}(:,tind+OS)>th1,W{a2}(:,tind+OS)<(-1*th2)),2));
%             meanCorrNN_Sh(a1,a2,tind+OS,rsh)=mean(mean(M(abs(W{a1}(:,tind+OS))<th1,abs(W{a2}(:,tind+OS))<(th2)),2));
% 
%         end
%         

        
        
        
    end%%%TimeVec
    
end

save(strcat('E:\Reports2\2018_01_26_NoiseCorrelation\ROI_Bin\NoiseStructure\HitCRSeparate\Mouse',num2str(Mouse)),'MMtx','CMtx','W','meanCorrG_Sh','meanCorrNG_Sh','meanCorrNGG_Sh','meanCorrNN_Sh','-v7.3')


end

Colors=linspecer(60);
meanVarT=zeros(9,60);
meanVarN=zeros(9,60);
meanVar=zeros(9,60);
meanCorrTN=zeros(9,9,60);
meanCorrTT=zeros(9,9,60);
meanCorrNN=zeros(9,9,60);

meanCorrG=zeros(9,9,60);
meanVarG=zeros(9,60);
meanCorrNG=zeros(9,9,60);
meanVarNG=zeros(9,60);
meanCorrNGG=zeros(9,9,60);



for i=1:60
    i
%     for a1=9
%         for a2=9
%             th1=1*sqrt(var(W{a1}(:,i)));
%             th2=1*sqrt(var(W{a2}(:,i)));
%             B1=zeros(cellCount,1);
%             B2=zeros(cellCount,1);
%             B1(CortexArea==a1)=(W{a1}(:,i));
%             B2(CortexArea==a2)=(W{a2}(:,i));
%             M=abs(CMtx(:,:,i));
%             V=diag(M);
%             meanVar(a1,i)=mean(V(CortexArea==a1));
%             meanVarT(a1,i)=mean(V(abs(B1)>th1));
%             meanVarN(a1,i)=mean(V(abs(B1)<th1 & CortexArea==a1));
%             meanVarG(a1,i)=mean(V(B1>th1));
%             meanVarNG(a1,i)=mean(V(B1<(-1*th1)));
%             M(eye(cellCount)==1)=0;
%             
%             meanCorrTT(a1,a2,i)=mean(mean(M(abs(B1)>th1,abs(B2)>th2),2));
%             meanCorrTN(a1,a2,i)=mean(mean(M(abs(B1)>th1,abs(B2)<th2),2));
%             meanCorrNN(a1,a2,i)=mean(mean(M(abs(B1)<th1,abs(B2)<th2),2));
%             meanCorrG(a1,a2,i)=mean(mean(M(B1>th1,B2>th2),2));
%             meanCorrNG(a1,a2,i)=mean(mean(M(B1<(-1*th1),B2<(-1*th2)),2));
%             meanCorrNGG(a1,a2,i)=mean(mean(M(B1>th1,B2<(-1*th2)),2));
%         end
%     end
    a1=9;a2=9;
    th=2*sqrt(var(W{a1}(:,i)));
    B1=abs(W{a1}(:,i));
    B2=abs(W{a2}(:,i));
    M=abs(CMtx(:,:,i));
    V=diag(M);
    meanVar(a1,i)=mean(V);
    meanVarT(a1,i)=mean(V(B1>th));
    meanVarN(a1,i)=mean(V(B1<th));
    meanVarG(a1,i)=mean(V(W{a1}(:,i)>th));
    meanVarNG(a1,i)=mean(V(W{a1}(:,i)<(-1*th)));
    M(eye(cellCount)==1)=0;
    
    meanCorrTT(a1,a2,i)=mean(mean(M(B1>th,B2>th),2));
    meanCorrTN(a1,a2,i)=mean(mean(M(B1>th,B2<th),2));
    meanCorrNN(a1,a2,i)=mean(mean(M(B1<th,B2<th),2));
    meanCorrG(a1,a2,i)=mean(mean(M(W{a1}(:,i)>th,W{a2}(:,i)>th),2));
    meanCorrNG(a1,a2,i)=mean(mean(M(W{a1}(:,i)<(-1*th),W{a2}(:,i)<(-1*th)),2));
    meanCorrNGG(a1,a2,i)=mean(mean(M(W{a1}(:,i)>th,W{a2}(:,i)<(-1*th)),2));
end





set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',14);
set(0,'defaultAxesFontName','Calibri');

a1=9;a2=9;
figure();hold on
xlabel('Time(s)');ylabel('Mean Absolute Covariance');
title('mean Correlations');
plot(-0.4:0.1:5.5,squeeze(meanCorrG(a1,a2,:)));
plot(-0.4:0.1:5.5,squeeze(meanCorrNG(a1,a2,:)),'r');
plot(-0.4:0.1:5.5,squeeze(meanCorrNGG(a1,a2,:)),'m');
% plot(-0.4:0.1:5.5,squeeze(meanCorrTT(a1,a2,:)),'g');
% plot(-0.4:0.1:5.5,squeeze(meanCorrTN(a1,a2,:)),'r');
% plot(-0.4:0.1:5.5,squeeze(meanCorrNN(a1,a2,:)),'k');
legend('Go','Nogo','Go and Nogo');

% figure();hold on
% %xlabel('Time(s)');ylabel('Mean Absolute Covariance');
% %title('mean Shuffled Correlations');
% %shadedErrorBar(dim,mean(V),sqrt(var(V)),'-b.',1);
% shadedErrorBar(-0.4:0.1:5.5,mean(meanCorrG_Sh(a1,a2,:,:),4),sqrt(var(squeeze(meanCorrG_Sh(a1,a2,:,:))')'),'b');
% shadedErrorBar(-0.4:0.1:5.5,mean(meanCorrNG_Sh(a1,a2,:,:),4),sqrt(var(squeeze(meanCorrNG_Sh(a1,a2,:,:))')'),'r');
% shadedErrorBar(-0.4:0.1:5.5,mean(meanCorrNGG_Sh(a1,a2,:,:),4),sqrt(var(squeeze(meanCorrNGG_Sh(a1,a2,:,:))')'),'m');
% 
% 
figure();hold on
title('Ind/corr');
xlabel('Time(s)');ylabel('Num Cell Estimation');
  plot(-0.4:0.1:5.5,squeeze(meanVarG(a1,:))'./squeeze(meanCorrG(a1,a2,:)),'b');
  plot(-0.4:0.1:5.5,squeeze(meanVarNG(a1,:))'./squeeze(meanCorrNG(a1,a2,:)),'r');
%plot(-0.4:0.1:2,meanVarT./meanCorrTT,'k');

% figure();hold on;
% title('fano factor');
% xlabel('Time(s)');ylabel('Variance / mean');
% %plot(meanVarT'./squeeze(mean(mean(MMtx(B>th,:,:),2),1)),'b');
% plot(-0.4:0.1:2,meanVar'./squeeze(mean(mean(MMtx(:,:,:),2),1)),'k');
% 



