% savePath='C:\Users\Sadegh\Documents\VLMReborn\Reports\2016_28_9_FisherTI_TB_ManualMap\DuringStim\HitMiss_alldata\L347_8_4_2015_active_tr';
% path='C:\Users\Sadegh\Documents\VLMReborn\L347\Data\allDays';
% mode=1;
% train=1;
function Fisher_PLS_TV(path,savePath,mode,train,TSin)

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
eventBin=A.eventBin;
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
LowDimVec=1:50;
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



%%%%%%%%%%%%%%%%%%%%%%%%Initialize Variables
if train ==1
    IFisher=zeros(8,length(tdVec));
    IFisherPLS=zeros(8,length(tdVec));
    BMAT_inst=zeros(cellCount,8,length(tdVec));
    InterceptMAT_inst=zeros(8,length(tdVec));
    NumPLSDim=zeros(8,length(tdVec));
    
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
    errorVal=zeros(8,length(tdVec));
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

for indtd=1:length(tdVec)
    indtd
     td=tdVec(indtd);
    for area=1:8
        
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
                if ( stimL>0 && td<=stimL)
                    dd0(:,c0)=X(:,SE0(k,1)+j+td);
                    shuffInd0(c0)=si;
                    c0=c0+1;
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
                if ( stimL>0 && td <= stimL)
                    dd1(:,c1)=X(:,SE1(k,1)+j+td);
                    shuffInd1(c1)=si;
                    c1=c1+1;
                end
            end
        end
        c1=c1-1;
        dd1=dd1(:,1:c1);
        
        c=round(min(c0,c1)/1.25);
        if Balanced==0
            ResampleMask=ones(1,c);
        elseif Balanced==1
            ResampleMask=zeros(1,c);
            ResampleMask(1:TSin)=1;
        end
        RM=[ResampleMask(randperm(c)),ResampleMask(randperm(c))];
        
        
        
        
        TrainingSize=c;
        
        TSout=c;
        if c<10
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
        
        area;
        
        for indld=1:length(LowDimVec)
            if train==1
                
                LowDim=LowDimVec(indld);
                [PLSRotation,YL] = plsregress(DataSet',group,min(LowDim,size(DataSet,1)));
                %[PLSRotationSh,YLSh] = plsregress(DataSetSh',group,min(LowDim,size(DataSet,1)));
                
                PLSRot{area,indld,indtd}=PLSRotation;
                
                %[B,intercept,ErrCurve,ErrVar,Lambda]=lassoglmcv(DataSet(:,RM==1),group(RM==1),shuffInd(RM==1),5,10,1);
                %[B,intercept,ErrCurve,ErrVar,Lambda]=elasticnetcv(DataSet(:,RM==1),group(RM==1),shuffInd(RM==1),5,10,1,0.5);
                [B,intercept,ErrCurve,ErrVar,Lambda]=glmnoreg(PLSRotation' * DataSet,group,shuffInd,5,10,1);
                %[BSh,interceptSh,~,~,~]=glmnoreg(PLSRotationSh' * DataSetSh,group,shuffInd,5,10,1);
                
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
                
%                 IFisherPLS_TrainNTSh(area,indld)=FisherInformation(PLSRotationSh' *DataSetSh(:,1:hc),PLSRotationSh' *DataSetSh(:,(hc+1):(2*hc)));
%                 IFisher_TrainNTSh(area,indld)=FisherInformationDecoder(PLSRotationSh' *DataSetSh(:,1:hc),PLSRotationSh' *DataSetSh(:,(hc+1):(2*hc)),BSh);
%                 IFisherPLS_ValNTSh(area,indld)=FisherInformation(PLSRotationSh' *ValDataSetSh(:,1:xc),PLSRotationSh' *ValDataSetSh(:,(xc+1):(2*xc)));
%                 IFisherPLSC_ValNTSh(area,indld)=FisherInformationVal(PLSRotationSh' *DataSetSh(:,1:hc),PLSRotationSh' *DataSetSh(:,(hc+1):(2*hc)),PLSRotationSh' *ValDataSetSh(:,1:xc),PLSRotationSh' *ValDataSetSh(:,(xc+1):(2*xc)));
%                 IFisher_ValNTSh(area,indld)=FisherInformationDecoder(PLSRotationSh' *ValDataSetSh(:,1:xc),PLSRotationSh' *ValDataSetSh(:,(xc+1):(2*xc)),BSh);
%                 errorValNTSh(area,indld)= mean( LRClassify([ValDataSetSh'] , PLSRotationSh *BSh , interceptSh) ~= [zeros(xc,1);ones(xc,1)]);
%                 
                
                
                
                
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%Pick the the best decoder 
        [~,maxind]=max(IFisherPLSC_ValNT(area,:));
        errorVal(area,indtd)=errorValNT(area,maxind);
        IFisherPLS(area,indtd)=IFisherPLSC_ValNT(area,maxind);
        IFisher(area,indtd)=IFisher_ValNT(area,maxind);
        
        NumPLSDim(area,indtd)=maxind;
        BMAT_inst(:,area,indtd)=BMAT(:,area,maxind);
        InterceptMAT_inst(area,indtd)=InterceptMAT(area,maxind);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
        
 
    end
    
    
end


 if(train==1 && TrainingSize>10)
     save(savePath,'errorVal','IFisherPLS','IFisher','NumPLSDim','BMAT_inst','InterceptMAT_inst','PLSRot');
%     save(savePath,'IFisher','Signal','Noise','errorVal','errorTrain','BMAT','InterceptMAT'...
%         ,'FAonFA','FAonCR','MissonMiss','MissonHit','HitSE','MissSE','CRSE','FASE','IFisherPLS'...
%         ,'PLSRot','shuff0','shuff1','LowDimVec','TrainingSize','ValidationSize','xc','FATrialCount','MissTrialCount'...
%         ,'HitScores','MissScores','FAScores','CRScores','GausW','errorValNT','IFisherPLS_ValNT','IFisherPLSC_ValNT','IFisherPLS_TrainNT','IFisher_ValNT','IFisher_TrainNT'...
%         ,'errorValNTSh','IFisherPLS_ValNTSh','IFisherPLSC_ValNTSh','IFisherPLS_TrainNTSh','IFisher_ValNTSh','IFisher_TrainNTSh');
 end

end



