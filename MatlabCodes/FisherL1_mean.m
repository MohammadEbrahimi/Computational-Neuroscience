% savePath='C:\Users\Sadegh\Documents\VLMReborn\Reports\2016_22_8_FisherTI_balanced\DuringStim\MissCR\L347_8_4_2015_active_tr.mat';
% path='C:\Users\Sadegh\Documents\VLMReborn\L347\Data\8_4_2015';
% mode=3;
% train=1;
function [TSout]=FisherL1_mean(path,savePath,mode,train,TSin)

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


cellEvents=padding(eventBin,0,2);

set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',12);
SpeedTh=50;
ActiveTrialNumber=20;

ActiveAnimal=ones(length(Lick),1);
for i=1:(length(Lick)-75*ActiveTrialNumber)
    if max(Lick(i:(i+75*ActiveTrialNumber)))==0
        ActiveAnimal(i:(i+75*ActiveTrialNumber))=0;
    end
end

%%%%%%% Erode area for conservative test
% erodst=strel('disk',70,8);
% for i=1:8
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

%%%% creating data sets
if mode==1
    group0=NogoTrials;
    group1=GoTrials;
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

if train ==1
    IFisher=zeros(8,length(tdVec));
    IFisherPLS=zeros(8,length(tdVec));
    errorD=zeros(8,1);
    errorVal=zeros(8,length(tdVec));
    errorDVar=zeros(8,length(tdVec));
    BMAT=zeros(cellCount,8);
    InterceptMAT=zeros(8,1);
    Signal=zeros(8,length(tdVec));
    Noise=zeros(8,length(tdVec));
    FAonCR=zeros(8,length(tdVec));
    FAonFA=zeros(8,length(tdVec));
    MissonHit=zeros(8,length(tdVec));
    MissonMiss=zeros(8,length(tdVec));
    FARatio=zeros(8,2);
    MissRatio=zeros(8,2);  
    FATrialCount=zeros(8,length(tdVec));
    MissTrialCount=zeros(8,length(tdVec));
end

%%%%%%%%%% Balance Trials
%     nse0=min(size(CRSE,1),size(FASE,1));
%     nse1=min(size(HitSE,1),size(MissSE,1));
%     BSE0=zeros(nse0*2,2);
%     BSE1=zeros(nse1*2,2);
%     for i=1:nse0
%         BSE0(2*i -1,:)=CRSE(i,:);
%         BSE0(2*i   ,:)=FASE(i,:);
%     end
%     for i=1:nse1
%         BSE1(2*i -1,:)=HitSE(i,:);
%         BSE1(2*i   ,:)=MissSE(i,:);
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






    SE0=CRSE(2:end,:);
    SE1=HitSE(2:end,:);
       
    shuff0=randperm(size(SE0,1));
    shuff1=randperm(size(SE1,1)); 

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
            if ((SE0(k,1)+j+i) <= SE0(k,2)+1  && stimL>0)
                dd0(:,c0)=mean(X(:,SE0(k,1)+j:SE0(k,1)+j+stimL),2);
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
            if ((SE1(k,1)+j+i) <= SE1(k,2)+1 && stimL>0)
                dd1(:,c1)=mean(X(:,SE1(k,1)+j:SE1(k,1)+j+stimL),2);
                shuffInd1(c1)=si;
                c1=c1+1;
            end
        end
    end
    c1=c1-1;
    dd1=dd1(:,1:c1);

        c=round(min(c0,c1)/2);
        ResampleMask=ones(1,c);
%        ResampleMask=zeros(1,c);
%        ResampleMask(1:TSin)=1;
        RM=[ResampleMask(randperm(c)),ResampleMask(randperm(c))];

    TSout=c;
    
    
    TrainingSize=c;
    
    
    
    shuffInd = [shuffInd0(1:c),(shuffInd1(1:c)+shuffInd0(c))];
    
    %%%%%%%%%%%%% Decoder
    DataSet=[dd0(:,1:c) dd1(:,1:c)];
    group=[zeros(c,1);ones(c,1)];
    area
    if train==1
        [B,intercept,ErrCurve,ErrVar,Lambda]=lassoglmcv(DataSet(:,RM==1),group(RM==1),shuffInd(RM==1),5,10,1);
       % [B,intercept,ErrCurve,ErrVar,Lambda]=glmnoreg(DataSet(:,RM==1),group(RM==1),shuffInd(RM==1),5,10,1);

        ErrorCurve{area}=ErrCurve;
        fullB=zeros(cellCount,1);
        fullB(CortexArea==area)=B;
        BMAT(:,area)=fullB;
        InterceptMAT(area)=intercept;
        errorD(area)=min(ErrCurve);
        
        LowDim=5;
        [PLSRotation,YL] = plsregress(DataSet',group,min(LowDim,size(DataSet,1)));
        PLSRot{area}=PLSRotation;
        
        
    else
        
        B=BMAT(CortexArea==area,area);
        intercept=InterceptMAT(area);
        PLSRotation=PLSRot{area};
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    foc=0;CRCount=0;
    fof=0;FACount=0;
    moh=0;HitCount=0;
    mom=0;MissCount=0;
    
    
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
                
                if (stimL>0 && td<stimL)
                    vc0=vc0+1;
                    if td < 0
                       dv0(:,vc0)=mean(X(:,VSE0(k,1)+j-5:VSE0(k,1)+j-1),2);
                    else
                        dv0(:,vc0)=mean(X(:,VSE0(k,1)+j:VSE0(k,1)+j+stimL),2);
                    end
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
                    if td <0
                        dv1(:,vc1)=mean(X(:,VSE1(k,1)+j-5:VSE1(k,1)+j-1),2);
                    else
                        dv1(:,vc1)=mean(X(:,VSE1(k,1)+j:VSE1(k,1)+j+stimL),2);
                    end
                end
            end
        end
        
        TSE0=FASE(2:end,:);
        TSE1=MissSE(2:end,:);
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
                    if td <0
                        tv0(:,tc0)=mean(X(:,TSE0(k,1)+j-5:TSE0(k,1)+j-1),2);
                    else 
                        tv0(:,tc0)=mean(X(:,TSE0(k,1)+j:TSE0(k,1)+j+stimL),2);
                    end
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
                    if td<0
                        tv1(:,tc1)=mean(X(:,TSE1(k,1)+j-5:TSE1(k,1)+j-1),2);
                    else
                        tv1(:,tc1)=mean(X(:,TSE1(k,1)+j:TSE1(k,1)+j+stimL),2);
                    end
                end
            end
        end
        tv1=tv1(:,1:tc1);
        
        if (td > 0)
            HitCount=HitCount+vc1;
            MissCount=MissCount+tc1;
            CRCount=CRCount+vc0;
            FACount=FACount+tc0;
            mom=mom+sum( LRClassify([tv1'] , B , InterceptMAT(area)) ~= 1);
            moh=moh+sum( LRClassify([dv1'] , B , InterceptMAT(area)) ~= 1);
            fof=fof+sum( LRClassify([tv0'] , B , InterceptMAT(area)) ~= 0);
            foc=foc+sum( LRClassify([dv0'] , B, InterceptMAT(area)) ~= 0);
        end
        
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
        
        FAonCR(area,indtd)=FAonCR_buf;
        FAonFA(area,indtd)=  FAonFA_buf;
        MissonMiss(area,indtd)= MissonMiss_buf;
        MissonHit(area,indtd)=  MissonHit_buf;
        
        
        vc=min(vc0,vc1);
        
        ValidationSize=vc;
        
        dv0=dv0(:,1:vc);
        dv1=dv1(:,1:vc);
        
        groupV=[zeros(vc,1);ones(vc,1)];
        
        
        md0=sum(dv0,2)/vc;
        md1=sum(dv1,2)/vc;
        
        S0=((dv0-md0*ones(1,vc))*(dv0-md0*ones(1,vc))')/(vc-1);
        S1=((dv1-md1*ones(1,vc))*(dv1-md1*ones(1,vc))')/(vc-1);
        
        S=(S0+S1)/2;
        dm=md1-md0;
        
        
        IFisher(area,indtd)=(B' * dm )^2 /(B'* S *B);
        
        Signal(area,indtd)= (B' * dm )^2;
        Noise(area,indtd)= (B'*S*B);
        
        %%%%%%%%%%PLS
        Rdv0=PLSRotation' * dv0 ;
        Rdv1=PLSRotation' * dv1 ;
        
        Rmd0=sum(Rdv0,2)/vc;
        Rmd1=sum(Rdv1,2)/vc;
        
        RS0=((Rdv0-Rmd0*ones(1,vc))*(Rdv0-Rmd0*ones(1,vc))')/(vc-1);
        RS1=((Rdv1-Rmd1*ones(1,vc))*(Rdv1-Rmd1*ones(1,vc))')/(vc-1);
        
        RS=(RS0+RS1)/2;
        Rdm=Rmd1-Rmd0;
        
        IFisherPLS(area,indtd)=  Rdm' * inv(RS) * Rdm;
        %%%%%%%%%%%%%%%%%
        
        errorVal(area,indtd)= mean( LRClassify([dv0' ; dv1'] , B , intercept) ~= groupV);
        
    end
    
    FARatio(area,:)=[fof/FACount,foc/CRCount];
    MissRatio(area,:)=[mom/MissCount,moh/HitCount];
    
end


if(train==1)
    save(savePath,'IFisher','Signal','Noise','errorVal','ErrorCurve','BMAT','InterceptMAT'...
        ,'FAonFA','FAonCR','MissonMiss','MissonHit','HitSE','MissSE','CRSE','FASE','IFisherPLS'...
        ,'PLSRot','shuff0','shuff1','LowDim','TrainingSize','ValidationSize','FATrialCount','MissTrialCount');
end

end


