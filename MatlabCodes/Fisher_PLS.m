%savePath='C:\Users\Sadegh\Documents\VLMReborn\Reports\2016_27_10_Fisher_PLS_TI_TB_ManualMap\DuringStim\HitCR_alldata\L354_12contrast_active_tr';
path='C:\Users\Sadegh\Documents\VLMReborn\L362bis\Data\allDays';
mode=1;
train=1;
%function Fisher_PLS(path,savePath,mode,train,TSin)

A=open(strcat(path,'\All_Sessions.mat'));
M=open(strcat(path,'\areas.mat'));
I=open(strcat(path,'\cellImage.mat'));


AirPuff=A.AirPuff;
Angle=A.Angle;
Area=M.Area;
CR=A.CR;
CRSE=A.CRSE;
Delay=A.Delay;
FA=A.FA;
FASE=A.FASE;
GoSE=A.GoSE;
GoTrials=A.GoTrials;
Hit=A.Hit;
HitSE=A.HitSE;
IC_time=A.IC_time;
Lick=A.Lick;
Miss=A.Miss;
MissSE=A.MissSE;
NoGoSE=A.NoGoSE;
NogoTrials=A.NogoTrials;
Punish=A.Punish;
RewardWindow=A.RewardWindow;
Speed=A.Speed;
TileTopLeft=I.TileTopLeft;
TimeOut=A.TimeOut;
WaterReward=A.WaterReward;
XSpeed=A.XSpeed;
%X_tile=A.X_tile;
YSpeed=A.YSpeed;
%Y_tile=A.Y_tile;
cellCount=A.cellCount;
cellImage=I.cellImage;
%eventBin=A.eventBin;
cellData=A.cellData;
%cellData=A.normCellData;


clear PLSRot

if train==0
    R=open(strcat(savePath,'.mat'));
    BMAT=R.BMAT;
    InterceptMAT=R.InterceptMAT;
    IFisher_pre=R.IFisher;
    PLSRot=R.PLSRot;
end

%%%%%%%%%%%%%%%%% Non Negative Matrix Deconvolution
% NNEvents=zeros(size(cellData));
% V.dt=0.1;
% V.Ncells=1;
% for i=1:cellCount
%     NNEvents(i,:)=(fast_oopsi(cellData(i,:),V))';
% end

%%%%%%%%%%%%%%%%%%


cellEvents=cellData;%padding(eventBin,0,2);

set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',12);
SpeedTh=50;
ActiveTrialNumber=20;

ActiveAnimal=ones(length(Lick),1);
% ActiveAnimal(1:204000)=0;
% ActiveAnimal(216000:end)=0;
for i=1:(length(Lick)-75*ActiveTrialNumber)
    if max(Lick(i:(i+75*ActiveTrialNumber)))==0
        ActiveAnimal(i:(i+75*ActiveTrialNumber))=0;
    end
end

%%%%%%% Erode area for conservative test
% erodst=strel('disk',100,8);
% for i=6:6
%     Area(:,:,i)=imerode(Area(:,:,i),erodst);
% end

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

clear HitScores
clear MissScores
clear FAScores
clear CRScores


%%%% creating data sets
%%%%%%%%%% Balance Trials
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SE0=CRSE(2:end,:);
SE1=HitSE(2:end,:);
TSE0=FASE(2:end,:);
TSE1=MissSE(2:end,:);
Balanced=0;
LowDimVec=1:20;
if mode==1
    group0=NogoTrials+GoTrials;
    group1=NogoTrials+GoTrials;
    LMax=20;
    tdVec=-5:19;
elseif mode==2
    group0=Delay;
    group1=Delay;
    LMax=5;
    tdVec=-5:4;
elseif mode==3
    group0=RewardWindow;
    group1=RewardWindow;
    LMax=30;
    tdVec=-5:29;
end


%%%%%%%%%%%%%cut in time
% maxTime=106972;
% % group0=group0(1:maxTime);
% % group1=group1(1:maxTime);
% % cellEvents=cellEvents(:,1:maxTime);
% SE0=SE0(SE0(:,2)>maxTime,:);
% SE1=SE1(SE1(:,2)>maxTime,:);
% TSE0=TSE0(TSE0(:,2)>maxTime,:);
% TSE1=TSE1(TSE1(:,2)>maxTime,:);


%%%%%%%%%%%%%%%%%%%%%%%%Initialize Variables
if train ==1
    IFisher=zeros(8,length(tdVec),length(LowDimVec));
    IFisherPLS=zeros(8,length(tdVec),length(LowDimVec));
    IFisherPLS_TrainNT=zeros(8,length(LowDimVec));
    IFisherPLS_ValNT=zeros(8,length(LowDimVec));
    IFisherPLSC_ValNT=zeros(8,length(LowDimVec));
    IFisher_TrainNT=zeros(8,length(LowDimVec));
    IFisher_ValNT=zeros(8,length(LowDimVec));
    errorValNT=zeros(8,length(LowDimVec));
    IFisherPLS_TrainNTSh=zeros(8,length(LowDimVec));
    IFisherPLS_ValNTSh=zeros(8,length(LowDimVec));
    IFisherPLSC_ValNTSh=zeros(8,length(LowDimVec));
    IFisher_TrainNTSh=zeros(8,length(LowDimVec));
    IFisher_ValNTSh=zeros(8,length(LowDimVec));
    errorValNTSh=zeros(8,length(LowDimVec));
    
    errorTrain=zeros(8,length(LowDimVec));
    errorVal=zeros(8,length(tdVec),length(LowDimVec));
    BMAT=zeros(cellCount,8,length(LowDimVec));
    PBMAT=zeros(max(LowDimVec),8,length(LowDimVec));
    InterceptMAT=zeros(8,1,length(LowDimVec));
    Signal=zeros(8,length(tdVec),length(LowDimVec));
    Noise=zeros(8,length(tdVec),length(LowDimVec));
    FAonCR=zeros(8,length(tdVec),length(LowDimVec));
    FAonFA=zeros(8,length(tdVec),length(LowDimVec));
    MissonHit=zeros(8,length(tdVec),length(LowDimVec));
    MissonMiss=zeros(8,length(tdVec),length(LowDimVec));
    FATrialCount=zeros(8,length(tdVec));
    MissTrialCount=zeros(8,length(tdVec));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%Start Learning

shuff0=randperm(size(SE0,1));
shuff1=randperm(size(SE1,1));

for area=1:3
    
    X=cellEvents(CortexArea==area,:);
    
    
    if size(X,1)<5
        continue;
    end
    
    c0=1;
    c1=1;
    
    
    dd0=zeros(size(X,1),length(group0));
    for si=1:length(shuff0)
        
        k=shuff0(si);
        if(max(Speed(SE0(k,1):(SE0(k,1)+55)))<=SpeedTh && max(ActiveAnimal(SE0(k,1):SE0(k,2)))==1 ...
                )%% && max(TimeOut(SE0(k,1):(SE0(k,1)+55))) ==0 )
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
                dd0(:,c0:c0+stimL)=X(:,SE0(k,1)+j:SE0(k,1)+j+stimL);
                shuffInd0(c0:c0+stimL)=si;
                c0=c0+stimL+1;
            end
        end
    end
    c0=c0-1;
    dd0=dd0(:,1:c0);
    
    
    dd1=zeros(size(X,1),length(group1));
    for si=1:length(shuff1)
        k=shuff1(si);
        if(max(Speed(SE1(k,1):(SE1(k,1)+55)))<=SpeedTh && max(ActiveAnimal(SE1(k,1):SE1(k,2)))==1 ...
                )%%&& max(TimeOut(SE1(k,1):(SE1(k,1)+55))) ==0 )
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
            if ((SE1(k,1)+j+i) <= SE1(k,2)+1 && stimL>0)
                dd1(:,c1:c1+stimL)=X(:,SE1(k,1)+j:SE1(k,1)+j+stimL);
                shuffInd1(c1:c1+stimL)=si;
                c1=c1+stimL+1;
            end
        end
    end
    c1=c1-1;
    dd1=dd1(:,1:c1);
    
    c=round(min(c0,c1)/2);
    if Balanced==0
        ResampleMask=ones(1,c);
    elseif Balanced==1
        ResampleMask=zeros(1,c);
        ResampleMask(1:TSin)=1;
    end
    RM=[ResampleMask(randperm(c)),ResampleMask(randperm(c))];
    
    
    
    
    TrainingSize=c;
    
    TSout=c;
    if c<100
        break;
    end
    
    shuffInd = [shuffInd0(1:c),(shuffInd1(1:c)+shuffInd0(c))];
    
    %%%%%%%%%%%%% Decoder
    xc=min(c1-c,c0-c);
    DataSet=[dd0(:,1:c) dd1(:,1:c)];
    DataSetSh=DataSet(:,randperm(2*c));
    DataSet=DataSet(:,RM==1);
    DataSetSh=DataSetSh(:,RM==1);
    hc=sum(RM)/2;
    
    ValDataSet=[dd0(:,c+1:c+xc) dd1(:,c+1:c+xc)];
    ValDataSetSh=ValDataSet(:,randperm(2*xc));
    
    group=[zeros(c,1);ones(c,1)];
    group=group(RM==1);
    
    area
    
    for indld=1:length(LowDimVec)
        if train==1
            
            LowDim=LowDimVec(indld);
            [PLSRotation,YL] = plsregress(DataSet',group,min(LowDim,size(DataSet,1)));
            [PLSRotationSh,YLSh] = plsregress(DataSetSh',group,min(LowDim,size(DataSet,1)));
            
            PLSRot{area,indld}=PLSRotation;
            
            %[B,intercept,ErrCurve,ErrVar,Lambda]=lassoglmcv(DataSet(:,RM==1),group(RM==1),shuffInd(RM==1),5,10,1);
            %[B,intercept,ErrCurve,ErrVar,Lambda]=elasticnetcv(DataSet(:,RM==1),group(RM==1),shuffInd(RM==1),5,10,1,0.5);
            [B,intercept,ErrCurve,ErrVar,Lambda]=glmnoreg(PLSRotation' * DataSet,group,shuffInd,5,10,1);
            [BSh,interceptSh,~,~,~]=glmnoreg(PLSRotationSh' * DataSetSh,group,shuffInd,5,10,1);

            PBMAT(1:length(B),area,indld)=B;
            smallB=B;
            B=PLSRotation *B;
            fullB=zeros(cellCount,1);
            fullB(CortexArea==area)=B;
            BMAT(:,area,indld)=fullB;
            InterceptMAT(area,indld)=intercept;
            
            
            errorTrain(area,indld)=ErrCurve;
            IFisherPLS_TrainNT(area,indld)=FisherInformation(PLSRotation' *DataSet(:,1:hc),PLSRotation' *DataSet(:,(hc+1):(2*hc)));
            IFisher_TrainNT(area,indld)=FisherInformationDecoder(PLSRotation' *DataSet(:,1:hc),PLSRotation' *DataSet(:,(hc+1):(2*hc)),smallB);
            IFisherPLS_ValNT(area,indld)=FisherInformation(PLSRotation' *ValDataSet(:,1:xc),PLSRotation' *ValDataSet(:,(xc+1):(2*xc)));
            IFisherPLSC_ValNT(area,indld)=FisherInformationVal(PLSRotation' *DataSet(:,1:hc),PLSRotation' *DataSet(:,(hc+1):(2*hc)),PLSRotation' *ValDataSet(:,1:xc),PLSRotation' *ValDataSet(:,(xc+1):(2*xc)));
            IFisher_ValNT(area,indld)=FisherInformationDecoder(PLSRotation' *ValDataSet(:,1:xc),PLSRotation' *ValDataSet(:,(xc+1):(2*xc)),smallB);
            errorValNT(area,indld)= mean( LRClassify([ValDataSet'] , B , intercept) ~= [zeros(xc,1);ones(xc,1)]);

            IFisherPLS_TrainNTSh(area,indld)=FisherInformation(PLSRotationSh' *DataSetSh(:,1:hc),PLSRotationSh' *DataSetSh(:,(hc+1):(2*hc)));
            IFisher_TrainNTSh(area,indld)=FisherInformationDecoder(PLSRotationSh' *DataSetSh(:,1:hc),PLSRotationSh' *DataSetSh(:,(hc+1):(2*hc)),BSh);
            IFisherPLS_ValNTSh(area,indld)=FisherInformation(PLSRotationSh' *ValDataSetSh(:,1:xc),PLSRotationSh' *ValDataSetSh(:,(xc+1):(2*xc)));
            IFisherPLSC_ValNTSh(area,indld)=FisherInformationVal(PLSRotationSh' *DataSetSh(:,1:hc),PLSRotationSh' *DataSetSh(:,(hc+1):(2*hc)),PLSRotationSh' *ValDataSetSh(:,1:xc),PLSRotationSh' *ValDataSetSh(:,(xc+1):(2*xc)));
            IFisher_ValNTSh(area,indld)=FisherInformationDecoder(PLSRotationSh' *ValDataSetSh(:,1:xc),PLSRotationSh' *ValDataSetSh(:,(xc+1):(2*xc)),BSh);
            errorValNTSh(area,indld)= mean( LRClassify([ValDataSetSh'] , PLSRotationSh *BSh , interceptSh) ~= [zeros(xc,1);ones(xc,1)]);

            

            
            
        end
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    VSE0=SE0(shuff0(shuffInd0(c)+1:end),:);
    VSE1=SE1(shuff1(shuffInd1(c)+1:end),:);
    Valsize=min(size(VSE0,1),size(VSE1,1));
    
    for indtd=1:length(tdVec)
        
        
        td=tdVec(indtd);
        vc0=0;
        vc1=0;
        
        dv0=zeros(size(X,1),Valsize);
        for k=1:Valsize
            if(max(Speed(VSE0(k,1):(VSE0(k,1)+55)))<=SpeedTh && max(ActiveAnimal(VSE0(k,1):VSE0(k,2)))==1 ...
                    )%% && max(TimeOut(VSE0(k,1):(VSE0(k,1)+55))) ==0 )
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
                    dv0(:,vc0)=X(:,VSE0(k,1)+j+td);
                end
            end
        end
        
        dv1=zeros(size(X,1),Valsize);
        for k=1:Valsize
            if( max(Speed(VSE1(k,1):(VSE1(k,1)+55)))<=SpeedTh && max(ActiveAnimal(VSE1(k,1):VSE1(k,2)))==1 ...
                    )%% && max(TimeOut(VSE1(k,1):(VSE1(k,1)+55))) ==0 )
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
                    dv1(:,vc1)=X(:,VSE1(k,1)+j+td);
                end
            end
        end
        

        tc0=0;
        tc1=0;
        
        tv0=zeros(size(X,1),size(TSE0,1));
        for k=1:size(TSE0,1)
            if(max(Speed(TSE0(k,1):(TSE0(k,1)+55)))<=SpeedTh && max(ActiveAnimal(TSE0(k,1):TSE0(k,2)))==1 ...
                    )%%&& max(TimeOut(TSE0(k,1):(TSE0(k,1)+55))) ==0 )
                
                j=0;
                while group0(TSE0(k,1)+j)==0 && ((TSE0(k,1)+j) < TSE0(k,2))
                    j=j+1;
                    
                end
                stimL=0;
                for i=0:LMax
                    if group0(TSE0(k,1)+j+i)==0
                        stimL=i;
                        break;
                    end
                end
                
                if (stimL>0 && td<stimL)
                    tc0=tc0+1;
                    tv0(:,tc0)=X(:,TSE0(k,1)+j+td);
                end
            end
        end
        tv0=tv0(:,1:tc0);
        
        tv1=zeros(size(X,1),size(TSE1,1));
        for k=1:size(TSE1,1)
            if(max(Speed(TSE1(k,1):(TSE1(k,1)+55)))<=SpeedTh && max(ActiveAnimal(TSE1(k,1):TSE1(k,2)))==1 ...
                    )%% && max(TimeOut(TSE1(k,1):(TSE1(k,1)+55))) ==0 )
                j=0;
                while group1(TSE1(k,1)+j)==0 && ((TSE1(k,1)+j) < TSE1(k,2))
                    j=j+1;
                end
                
                stimL=0;
                for i=0:LMax
                    if group1(TSE1(k,1)+j+i)==0
                        stimL=i;
                        break;
                    end
                end
                
                if (stimL>0 && td<stimL)
                    tc1=tc1+1;
                    tv1(:,tc1)=X(:,TSE1(k,1)+j+td);
                end
            end
        end
        tv1=tv1(:,1:tc1);
        
        
        for indld=1:length(LowDimVec)
            B=BMAT(CortexArea==area,area,indld);
            intercept=InterceptMAT(area,indld);
            PLSRotation=PLSRot{area,indld};
            
            
            HitScores{area,indtd,indld}=(dv1'*B) +InterceptMAT(area);
            MissScores{area,indtd,indld}= (tv1'*B) +InterceptMAT(area);
            FAScores{area,indtd,indld}= ((tv0'*B) +InterceptMAT(area));
            CRScores{area,indtd,indld}= ((dv0'*B) +InterceptMAT(area));
            
            
            
            
            FAonCR_buf=0;
            FAonFA_buf= 0;
            MissonMiss_buf=0;
            MissonHit_buf=0;
            
            FATrialCount(area,indtd)=min(vc0,tc0);
            MissTrialCount(area,indtd)=min(vc1,tc1);
            
            
            for rep=1:50
                CRmask=randperm(vc0,min(vc0,tc0));
                FAmask=randperm(tc0,min(vc0,tc0));
                Missmask=randperm(tc1,min(vc1,tc1));
                Hitmask=randperm(vc1,min(vc1,tc1));
                
                FAonCR_buf=FAonCR_buf+(mean( LRClassify([dv0(:,CRmask)'] , B, InterceptMAT(area)) ~= 0)/50);
                FAonFA_buf=FAonFA_buf+ ( mean( LRClassify([tv0(:,FAmask)'] , B , InterceptMAT(area)) ~= 0)/50);
                MissonMiss_buf=MissonMiss_buf+( mean( LRClassify([tv1(:,Missmask)'] , B , InterceptMAT(area)) ~= 1)/50);
                MissonHit_buf=MissonHit_buf+ ( mean( LRClassify([dv1(:,Hitmask)'] , B , InterceptMAT(area)) ~= 1)/50);
                
            end
            
            FAonCR(area,indtd,indld)=FAonCR_buf;
            FAonFA(area,indtd,indld)=  FAonFA_buf;
            MissonMiss(area,indtd,indld)= MissonMiss_buf;
            MissonHit(area,indtd,indld)=  MissonHit_buf;
            
            
            vc=min(vc0,vc1);
            
            ValidationSize=vc;
            
            ddv0=dv0(:,1:vc);
            ddv1=dv1(:,1:vc);
            
            groupV=[zeros(vc,1);ones(vc,1)];
            
            
            md0=sum(ddv0,2)/vc;
            md1=sum(ddv1,2)/vc;
            
            S0=((ddv0-md0*ones(1,vc))*(ddv0-md0*ones(1,vc))')/(vc-1);
            S1=((ddv1-md1*ones(1,vc))*(ddv1-md1*ones(1,vc))')/(vc-1);
            
            S=(S0+S1)/2;
            dm=md1-md0;
            
            
            IFisher(area,indtd,indld)=(B' * dm )^2 /(B'* S *B);
            
            Signal(area,indtd,indld)= (B' * dm )^2;
            Noise(area,indtd,indld)= (B'*S*B);
            
            %%%%%%%%%%PLS
            Rdv0=PLSRotation' * ddv0 ;
            Rdv1=PLSRotation' * ddv1 ;
            
            Rmd0=sum(Rdv0,2)/vc;
            Rmd1=sum(Rdv1,2)/vc;
            
            RS0=((Rdv0-Rmd0*ones(1,vc))*(Rdv0-Rmd0*ones(1,vc))')/(vc-1);
            RS1=((Rdv1-Rmd1*ones(1,vc))*(Rdv1-Rmd1*ones(1,vc))')/(vc-1);
            
            RS=(RS0+RS1)/2;
            Rdm=Rmd1-Rmd0;
            GausW{area,indtd,indld}= inv(RS)* Rdm;
            
            IFisherPLS(area,indtd,indld)=  Rdm' * inv(RS) * Rdm;
            %%%%%%%%%%%%%%%%%
            
            errorVal(area,indtd,indld)= mean( LRClassify([ddv0' ; ddv1'] , B , intercept) ~= groupV);
            
        end
    end
    
    
end


% if(train==1 && TrainingSize>100)
%     save(savePath,'IFisher','Signal','Noise','errorVal','errorTrain','BMAT','InterceptMAT'...
%         ,'FAonFA','FAonCR','MissonMiss','MissonHit','HitSE','MissSE','CRSE','FASE','IFisherPLS'...
%         ,'PLSRot','shuff0','shuff1','LowDimVec','TrainingSize','ValidationSize','xc','FATrialCount','MissTrialCount'...
%         ,'HitScores','MissScores','FAScores','CRScores','GausW','errorValNT','IFisherPLS_ValNT','IFisherPLSC_ValNT','IFisherPLS_TrainNT','IFisher_ValNT','IFisher_TrainNT'...
%         ,'errorValNTSh','IFisherPLS_ValNTSh','IFisherPLSC_ValNTSh','IFisherPLS_TrainNTSh','IFisher_ValNTSh','IFisher_TrainNTSh');
% end

%end



