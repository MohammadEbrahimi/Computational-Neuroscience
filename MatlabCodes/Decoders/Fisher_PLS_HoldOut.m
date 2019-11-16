function Fisher_PLS_HoldOut(loadpath,savePath,mode,train,TSin)
load(strcat(loadpath,'\All_Sessions.mat'));
load(strcat(loadpath,'\areas.mat'));
load(strcat(loadpath,'\cellImage.mat'));
%%% creating data sets
%%%%%%%%% Balance Trials
nse0=min(size(CRSE,1),size(FASE,1));
nse1=min(size(HitSE,1),size(MissSE,1));
BSE0=zeros(nse0*2,2);
BSE1=zeros(nse1*2,2);
for i=1:nse0
    BSE0(2*i -1,:)=CRSE(i,:);
    BSE0(2*i   ,:)=FASE(i,:);
end
for i=1:nse1
    BSE1(2*i -1,:)=HitSE(i,:);
    BSE1(2*i   ,:)=MissSE(i,:);
end
%%%%

SE0=CRSE(2:end,:);
SE1=FASE(2:end,:);
Balanced=0;
LowDimVec=1:50;
SpeedTh=1;
ActiveTrialNumber=20;
AreaVec=1:8;

if train==0
    R=open(strcat(savePath,'.mat'));
    BMAT=R.BMAT;
    InterceptMAT=R.InterceptMAT;
    IFisher_pre=R.IFisher;
    PLSRot=R.PLSRot;
end


cellEvents=cellData;

set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',12);


ActiveAnimal=ones(length(Lick),1);
for i=1:(length(Lick)-75*ActiveTrialNumber)
    if max(Lick(i:(i+75*ActiveTrialNumber)))==0
        ActiveAnimal(i:(i+75*ActiveTrialNumber))=0;
    end
end


%%%%Find cell Coordinates and assign area
cellIJ=zeros(cellCount,2);
CortexArea=zeros(cellCount,1);

for i=1:cellCount
    [row,ind]=max(cellImage{i});
    [~,indj]=max(row);
    indi=ind(indj);
    cellIJ(i,:)=TileTopLeft(i,:)+[indi-1,indj-1];
    for a=1:8
        
        if Area(cellIJ(i,1),cellIJ(i,2),a)
            CortexArea(i)=a;
            
        end
        
        
    end
    
end



%%%% creating data sets
%%%%%%%%%% Balance Trials
% nse0=min(size(CRSE,1),size(FASE,1));
% nse1=min(size(HitSE,1),size(MissSE,1));
% BSE0=zeros(nse0*2,2);
% BSE1=zeros(nse1*2,2);
% for i=1:nse0
%     BSE0(2*i -1,:)=CRSE(i,:);
%     BSE0(2*i   ,:)=FASE(i,:);
% end
% for i=1:nse1
%     BSE1(2*i -1,:)=HitSE(i,:);
%     BSE1(2*i   ,:)=MissSE(i,:);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mode==1
    group0=NogoTrials+GoTrials;
    group1=NogoTrials+GoTrials;
    LMax=20;
    tdVec=-5:19;
    trainTime=0:19;
elseif mode==2
    group0=Delay;
    group1=Delay;
    LMax=5;
    tdVec=-5:4;
    trainTime=0:4;
elseif mode==3
    group0=RewardWindow;
    group1=RewardWindow;
    LMax=30;
    tdVec=-5:29;
    trainTime=0:29;
end


%%%%%%%%%%%%%cut in time
%maxTime=106972;
% group0=group0(1:maxTime);
% group1=group1(1:maxTime);
% cellEvents=cellEvents(:,1:maxTime);
% SE0=SE0(SE0(:,2)>maxTime,:);
% SE1=SE1(SE1(:,2)>maxTime,:);



%%%%%%%%%%%%%%%%%%%%%%%%Initialize Variables
if train ==1
    Info_Val_Time_LR=zeros(8,length(tdVec),length(LowDimVec));
    Info_Val_Time_Fisher=zeros(8,length(tdVec),length(LowDimVec));
    Info_Train_Fisher=zeros(8,length(LowDimVec));
    Info_Val_Fisher=zeros(8,length(LowDimVec));
    Info_Train_LR=zeros(8,length(LowDimVec));
    Info_Val_LR=zeros(8,length(LowDimVec));
    Error_Val_LR=zeros(8,length(LowDimVec));
    Error_Train_LR=zeros(8,length(LowDimVec));
    
    Info_Train_Fisher_Sh=zeros(8,length(LowDimVec));
    Info_Val_Fisher_Sh=zeros(8,length(LowDimVec));
    Info_Train_LR_Sh=zeros(8,length(LowDimVec));
    Info_Val_LR_Sh=zeros(8,length(LowDimVec));
    Error_Val_LR_Sh=zeros(8,length(LowDimVec));
    Error_Train_LR_Sh=zeros(8,length(LowDimVec));
    
    BMAT=zeros(cellCount,8,length(LowDimVec));
    %PBMAT=zeros(max(LowDimVec),8,length(LowDimVec));
    InterceptMAT=zeros(8,1,length(LowDimVec));
    Signal=zeros(8,length(tdVec),length(LowDimVec));
    Noise=zeros(8,length(tdVec),length(LowDimVec));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%Start Learning

shuff0=randperm(size(SE0,1));
shuff1=randperm(size(SE1,1));

% X=cellEvents(CortexArea==area,:);
c0=1;
c1=1;
dd0=zeros(size(cellEvents,1),sum(group0));
for si=1:length(shuff0)
    
    k=shuff0(si);
    if(max(Speed(SE0(k,1):(SE0(k,1)+55)))<=SpeedTh && max(ActiveAnimal(SE0(k,1):SE0(k,2)))==1)
        j=0;
        while group0(SE0(k,1)+j)==0 && ((SE0(k,1)+j) < SE0(k,2))
            j=j+1;
        end
        stimL=0;
        for i=0:LMax
            if group0(SE0(k,1)+j+i)==0
                stimL=i;
                break;
            end
        end
        stimL=stimL-1;
        if ((SE0(k,1)+j+i) <= SE0(k,2)+1  && stimL>0)
            start=(SE0(k,1)+j+min(trainTime));
            finish=(SE0(k,1)+j+min(stimL,max(trainTime)));
            dd0(:,c0:c0+(finish-start))=cellEvents(:, start : finish );
            shuffInd0(c0:c0+(finish-start))=si;
            c0=c0+(finish-start)+1;
        end
    end
end
c0=c0-1;
dd0=dd0(:,1:c0);


dd1=zeros(size(cellEvents,1),sum(group1));
for si=1:length(shuff1)
    
    k=shuff1(si);
    if(max(Speed(SE1(k,1):(SE1(k,1)+55)))<=SpeedTh && max(ActiveAnimal(SE1(k,1):SE1(k,2)))==1)
        j=0;
        while group1(SE1(k,1)+j)==0 && ((SE1(k,1)+j) < SE1(k,2))
            j=j+1;
        end
        stimL=0;
        for i=0:LMax
            if group1(SE1(k,1)+j+i)==0
                stimL=i;
                break;
            end
        end
        stimL=stimL-1;
        if ((SE1(k,1)+j+i) <= SE1(k,2)+1  && stimL>0)
            start=(SE1(k,1)+j+min(trainTime));
            finish=(SE1(k,1)+j+min(stimL,max(trainTime)));
            dd1(:,c1:c1+(finish-start))=cellEvents(:, start : finish );
            shuffInd1(c1:c1+(finish-start))=si;
            c1=c1+(finish-start)+1;
        end
    end
end
c1=c1-1;
dd1=dd1(:,1:c1);


c=floor(min(c0,c1)/2);
if Balanced==0
    ResampleMask=ones(1,c);
elseif Balanced==1
    ResampleMask=zeros(1,c);
    ResampleMask(1:TSin)=1;
end
RM=[ResampleMask(randperm(c)),ResampleMask(randperm(c))];


TrainingSize=sum(RM)/2;
return; 

if TrainingSize<100
    return;
end

shuffInd = [shuffInd0(1:c),(shuffInd1(1:c)+shuffInd0(c))];

%%%%%%%%%%%%% Decoder
TrainSet=[dd0(:,1:c) dd1(:,1:c)];

PLSDataSet=TrainSet;
ValSet=[dd0(:,c+1:2*c) dd1(:,c+1:2*c)];

% PLSDataSet=[dd0(:,c+1:2*c) dd1(:,c+1:2*c)];
% ValSet=[dd0(:,2*c+1:3*c) dd1(:,2*c+1:3*c)];

PLSDataSetSh=PLSDataSet(:,randperm(2*c));
TrainSetSh=TrainSet(:,randperm(2*c));
ValSetSh=ValSet(:,randperm(2*c));

PLSDataSet=PLSDataSet(:,RM==1);
PLSDataSetSh=PLSDataSetSh(:,RM==1);
TrainSet=TrainSet(:,RM==1);
TrainSetSh=TrainSetSh(:,RM==1);
hc=sum(RM)/2;


group=[zeros(c,1);ones(c,1)];
group=group(RM==1);

for area=AreaVec
    TrainingArea=area
    Cells=(CortexArea==area);
    if sum(Cells)<5
        continue;end
    for indld=1:length(LowDimVec)
        if train==1
            
            
            [PLSRot{area,indld},~] = plsregress(PLSDataSet(Cells,:)',group,min(LowDimVec(indld),sum(Cells)));
            [PLSRotSh{area,indld},~] = plsregress(PLSDataSetSh(Cells,:)',group,min(LowDimVec(indld),sum(Cells)));
            
            
            
            %[B,intercept,ErrCurve,ErrVar,Lambda]=lassoglmcv(DataSet(:,RM==1),group(RM==1),shuffInd(RM==1),5,10,1);
            %[B,intercept,ErrCurve,ErrVar,Lambda]=elasticnetcv(DataSet(:,RM==1),group(RM==1),shuffInd(RM==1),5,10,1,0.5);
            [B,intercept,ErrCurve,~,~]=glmnoreg(PLSRot{area,indld}' * TrainSet(Cells,:),group,shuffInd,5,10,1);
            [BSh,interceptSh,ErrCurveSh,~,~]=glmnoreg(PLSRotSh{area,indld}' * TrainSetSh(Cells,:),group,shuffInd,5,10,1);
            
            %PBMAT(1:length(B),area,indld)=B;
            smallB=B;
            B=PLSRot{area,indld} *B;
            fullB=zeros(cellCount,1);
            fullB(CortexArea==area)=B;
            BMAT(:,area,indld)=fullB;
            InterceptMAT(area,indld)=intercept;
            
            CovMAT{area,indld}=CovSets(PLSRot{area,indld}' *TrainSet(Cells,1:hc),PLSRot{area,indld}' *TrainSet(Cells,(hc+1):(2*hc)));
            Error_Train_LR(area,indld)=ErrCurve;
            Info_Train_Fisher(area,indld)=FisherInformation(PLSRot{area,indld}' *TrainSet(Cells,1:hc),PLSRot{area,indld}' *TrainSet(Cells,(hc+1):(2*hc)));
            Info_Train_LR(area,indld)=FisherInformationDecoder(PLSRot{area,indld}' *TrainSet(Cells,1:hc),PLSRot{area,indld}' *TrainSet(Cells,(hc+1):(2*hc)),smallB);
            Info_Val_Fisher(area,indld)=FisherInformationVal(PLSRot{area,indld}' *TrainSet(Cells,1:hc),PLSRot{area,indld}' *TrainSet(Cells,(hc+1):(2*hc)),PLSRot{area,indld}' *ValSet(Cells,1:c),PLSRot{area,indld}' *ValSet(Cells,(c+1):(2*c)));
            Info_Val_LR(area,indld)=FisherInformationDecoder(PLSRot{area,indld}' *ValSet(Cells,1:c),PLSRot{area,indld}' *ValSet(Cells,(c+1):(2*c)),smallB);
            Error_Val_LR(area,indld)= mean( LRClassify([ValSet(Cells,:)'] , B , intercept) ~= [zeros(c,1);ones(c,1)]);
            
            Error_Train_LR_Sh(area,indld)=ErrCurveSh;
            Info_Train_Fisher_Sh(area,indld)=FisherInformation(PLSRotSh{area,indld}' *TrainSetSh(Cells,1:hc),PLSRotSh{area,indld}' *TrainSetSh(Cells,(hc+1):(2*hc)));
            Info_Train_LR_Sh(area,indld)=FisherInformationDecoder(PLSRotSh{area,indld}' *TrainSetSh(Cells,1:hc),PLSRotSh{area,indld}' *TrainSetSh(Cells,(hc+1):(2*hc)),BSh);
            Info_Val_Fisher_Sh(area,indld)=FisherInformationVal(PLSRotSh{area,indld}' *TrainSetSh(Cells,1:hc),PLSRotSh{area,indld}' *TrainSetSh(Cells,(hc+1):(2*hc)),PLSRotSh{area,indld}' *ValSetSh(Cells,1:c),PLSRotSh{area,indld}' *ValSetSh(Cells,(c+1):(2*c)));
            Info_Val_LR_Sh(area,indld)=FisherInformationDecoder(PLSRotSh{area,indld}' *ValSetSh(Cells,1:c),PLSRotSh{area,indld}' *ValSetSh(Cells,(c+1):(2*c)),BSh);
            Error_Val_LR_Sh(area,indld)= mean( LRClassify([ValSetSh(Cells,:)'] , PLSRotSh{area,indld} * BSh , interceptSh) ~= [zeros(c,1);ones(c,1)]);
            
            
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

VSE0=SE0(shuff0(shuffInd0(c)+1:end),:);
VSE1=SE1(shuff1(shuffInd1(c)+1:end),:);
Valsize=min(size(VSE0,1),size(VSE1,1));

for indtd=1:length(tdVec)
    
    
    td=tdVec(indtd);
    vc0=0;
    vc1=0;
    
    dv0=zeros(size(cellEvents,1),Valsize);
    for k=1:Valsize
        if(max(Speed(VSE0(k,1):(VSE0(k,1)+55)))<=SpeedTh && max(ActiveAnimal(VSE0(k,1):VSE0(k,2)))==1)
            j=0;
            while group0(VSE0(k,1)+j)==0 && ((VSE0(k,1)+j) < VSE0(k,2))
                j=j+1;
                
            end
            stimL=0;
            for i=0:LMax
                if group0(VSE0(k,1)+j+i)==0
                    stimL=i;
                    break;
                end
            end
            
            if (stimL>0 && td<stimL )
                vc0=vc0+1;
                dv0(:,vc0)=cellEvents(:,VSE0(k,1)+j+td);
            end
        end
    end
    vc0=vc0-1;
    dv0=dv0(:,1:vc0);
    
    dv1=zeros(size(cellEvents,1),Valsize);
    for k=1:Valsize
        if( max(Speed(VSE1(k,1):(VSE1(k,1)+55)))<=SpeedTh && max(ActiveAnimal(VSE1(k,1):VSE1(k,2)))==1)
            j=0;
            while group1(VSE1(k,1)+j)==0 && ((VSE1(k,1)+j) < VSE1(k,2))
                j=j+1;
            end
            
            stimL=0;
            for i=0:LMax
                if group1(VSE1(k,1)+j+i)==0
                    stimL=i;
                    break;
                end
            end
            
            if (stimL>0 && td<stimL)
                vc1=vc1+1;
                dv1(:,vc1)=cellEvents(:,VSE1(k,1)+j+td);
            end
        end
    end
    
    vc1=vc1-1;
    dv1=dv1(:,1:vc0);
    vc=min(vc0,vc1);
    
    for area=AreaVec
        
        Cells=(CortexArea==area);
        if sum(Cells)<5
            continue;end
        ddv0=dv0(Cells,1:vc);
        ddv1=dv1(Cells,1:vc);
        
        for indld=1:length(LowDimVec)
            B=BMAT(Cells,area,indld);
            intercept=InterceptMAT(area,indld);
            PLSRotation=PLSRot{area,indld};
            
            
            Scores0{area,indtd,indld}=(dv0(Cells,:)'*B) +InterceptMAT(area);
            Scores1{area,indtd,indld}= (dv1(Cells,:)'*B) +InterceptMAT(area);
            
            
            
            
            
            
            
            ValidationSize=vc;
            
            
            
            groupV=[zeros(vc,1);ones(vc,1)];
            
            
            md0=sum(ddv0,2)/vc;
            md1=sum(ddv1,2)/vc;
            
            S0=((ddv0-md0*ones(1,vc))*(ddv0-md0*ones(1,vc))')/(vc-1);
            S1=((ddv1-md1*ones(1,vc))*(ddv1-md1*ones(1,vc))')/(vc-1);
            
            S=(S0+S1)/2;
            dm=md1-md0;
            
            
            Info_Val_Time_LR(area,indtd,indld)=(B' * dm )^2 /(B'* S *B);
            
            Signal(area,indtd,indld)= (B' * dm )^2;
            Noise(area,indtd,indld)= (B'*S*B);
            
            %%%%%%%%%%PLS
            RS=CovMAT{area,indld};
            Rdv0=PLSRotation' * ddv0 ;
            Rdv1=PLSRotation' * ddv1 ;
            
            Rmd0=sum(Rdv0,2)/vc;
            Rmd1=sum(Rdv1,2)/vc;
            Rdm=Rmd1-Rmd0;
            Info_Val_Time_Fisher(area,indtd,indld)=  Rdm' * inv(RS) * Rdm;
            %%%%%%%%%%%%%%%%%
            
            Error_Val_Time_LR(area,indtd,indld)= mean( LRClassify([ddv0' ; ddv1'] , B , intercept) ~= groupV);
            
        end
    end
    
    
end

% 
% if(train==1 && TrainingSize>100)
%     save(strcat(savePath,'_mode',num2str(mode)),'Signal','Noise','BMAT','InterceptMAT','CovMAT',...
%         'SE0','SE1','Info_Train_LR','Info_Val_Fisher','Info_Val_LR','Error_Val_LR',...
%         'PLSRot','PLSRotSh','shuff0','shuff1','LowDimVec','TrainingSize','ValidationSize','hc'...
%         ,'Scores0','Scores1','Error_Val_Time_LR','Info_Val_Time_Fisher','Info_Val_Time_LR','Error_Train_LR','Info_Train_Fisher',...
%         'Error_Train_LR_Sh','Info_Train_Fisher_Sh','Info_Train_LR_Sh','Info_Val_Fisher_Sh','Info_Val_LR_Sh','Error_Val_LR_Sh');
% end
end




