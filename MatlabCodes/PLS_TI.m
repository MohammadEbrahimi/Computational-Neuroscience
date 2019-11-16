% clear all
% path='C:\Users\Sadegh\Documents\VLMReborn\L347\Data\8_4_2015'
% A=open(strcat(path,'\All_Sessions.mat'));
% M=open(strcat(path,'\areas.mat'));
% I=open(strcat(path,'\cellImage.mat'));
% mode=1;
% LowDim=10;

function [IFisherPLS, PLSRot]=FisherPLS(path,LowDim,mode) 


AirPuff=A.AirPuff;
Angle=A.Angle
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
X_tile=A.X_tile;
YSpeed=A.YSpeed;
Y_tile=A.Y_tile;
cellCount=A.cellCount;
cellData=A.cellData;
cellImage=I.cellImage;
eventBin=A.eventBin;


cellEvents=cellData;%padding(eventBin,0,2);

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
    

IFisherPLS=zeros(8,length(tdVec));





for area=1:8
    
    
    
        SE0=CRSE(2:end,:);
        SE1=HitSE(2:end,:);

    X=cellEvents(CortexArea==area,:);

    
    shuff0=randperm(size(SE0,1));
    shuff1=randperm(size(SE1,1));
    
    c0=1;
    c1=1;
    
    clear shuffInd0
    dd0=zeros(size(X,1),length(group0));
    for si=1:length(shuff0)
        
        k=shuff0(si);
        if(max(Speed(SE0(k,1):SE0(k,2)))<=SpeedTh && max(ActiveAnimal(SE0(k,1):SE0(k,2)))==1  )
            j=0;
            while group0(SE0(k,1)+j)==0 && ((SE0(k,1)+j) < SE0(k,2))
                j=j+1;
            end
            stimL=0;
            for i=0:LMax
                if group0(SE0(k,1)+j+i)==0 && ((SE0(k,1)+j+i) < SE0(k,2))
                    stimL=i;
                    break;
                end
            end
            stimL=stimL-1;
            if ((SE0(k,1)+j+i) < SE0(k,2))
                dd0(:,c0:c0+stimL)=X(:,SE0(k,1)+j:SE0(k,1)+j+stimL);
                shuffInd0(c0:c0+stimL)=si;
                c0=c0+stimL+1;
            end
        end
    end
    c0=c0-1;
    dd0=dd0(:,1:c0);
    
    clear shuffInd1
    dd1=zeros(size(X,1),length(group1));
    for si=1:length(shuff1)
        k=shuff1(si);
        if(max(Speed(SE1(k,1):SE1(k,2)))<=SpeedTh && max(ActiveAnimal(SE1(k,1):SE1(k,2)))==1 )
            j=0;
            while group1(SE1(k,1)+j)==0 && ((SE1(k,1)+j) < SE1(k,2))
                j=j+1;
            end
            
            stimL=0;
            for i=0:LMax
                if group1(SE1(k,1)+j+i)==0 && ((SE1(k,1)+j+i) < SE1(k,2))
                    stimL=i;
                    break;
                end
            end
            
            stimL=stimL-1;
            if ((SE1(k,1)+j+i) < SE1(k,2))
                dd1(:,c1:c1+stimL)=X(:,SE1(k,1)+j:SE1(k,1)+j+stimL);
                shuffInd1(c1:c1+stimL)=si;
                c1=c1+stimL+1;
            end
        end
    end
    c1=c1-1;
    dd1=dd1(:,1:c1);
    

    v=round(min(c0,c1)/2);
    c=min(c0,c1)-v;
    

    shuffInd = [shuffInd0(1:c),(shuffInd1(1:c)+shuffInd0(c))];
    
%     % %%%%%%%%%%%%% Decoder
    DataSet=[dd0(:,1:c) dd1(:,1:c)];
    group=[zeros(c,1);ones(c,1)];
    area

    [PLSRotation,YL] = plsregress(DataSet',group,LowDim); 
    PLSRot{area}=PLSRotation;

%     % %%%% Information should be calculated on a seperate validation set

    
    
    
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
            if(max(Speed(VSE0(k,1):VSE0(k,2)))<=SpeedTh && max(ActiveAnimal(VSE0(k,1):VSE0(k,2)))==1 )
                j=0;
                while group0(VSE0(k,1)+j)==0 && ((VSE0(k,1)+j) < VSE0(k,2))
                    j=j+1;
                    
                end
                stimL=0;
                for i=0:LMax
                    if group0(VSE0(k,1)+j+i)==0 && ((VSE0(k,1)+j+i) < VSE0(k,2))
                        stimL=i;
                        break;
                    end
                end
                
                if (td<stimL && stimL>0)
                    vc0=vc0+1;
                    dv0(:,vc0)=X(:,VSE0(k,1)+j+td);
                end
            end
        end
        
        dv1=zeros(size(X,1),Valsize);
        for k=1:Valsize
            if(max(Speed(VSE1(k,1):VSE1(k,2)))<=SpeedTh && max(ActiveAnimal(VSE1(k,1):VSE1(k,2)))==1 )
                j=0;
                while group1(VSE1(k,1)+j)==0 && ((VSE1(k,1)+j) < VSE1(k,2))
                    j=j+1;
                end
                
                stimL=0;
                for i=0:LMax
                    if group1(VSE1(k,1)+j+i)==0 && ((VSE1(k,1)+j+i) < VSE1(k,2))
                        stimL=i;
                        break;
                    end
                end
                
                if (td<stimL && stimL>0)
                    vc1=vc1+1;
                    dv1(:,vc1)=X(:,VSE1(k,1)+j+td);
                end
            end
        end
        
        TSE0=FASE(2:end,:);
        TSE1=MissSE(2:end,:);
        tc0=0;
        tc1=0;
        
        tv0=zeros(size(X,1),size(TSE0,1));
        for k=1:size(TSE0,1)
            if(max(Speed(TSE0(k,1):TSE0(k,2)))<=SpeedTh && max(ActiveAnimal(TSE0(k,1):TSE0(k,2)))==1 )
                j=0;
                while group0(TSE0(k,1)+j)==0 && ((TSE0(k,1)+j) < TSE0(k,2))
                    j=j+1;
                    
                end
                stimL=0;
                for i=0:LMax
                    if group0(TSE0(k,1)+j+i)==0 && ((TSE0(k,1)+j+i) < TSE0(k,2))
                        stimL=i;
                        break;
                    end
                end
                
                if (stimL>0)
                    tc0=tc0+1;
                    tv0(:,tc0)=X(:,TSE0(k,1)+j+td);
                end
            end
        end
        tv0=tv0(:,1:tc0);
        
        tv1=zeros(size(X,1),size(TSE1,1));
        for k=1:size(TSE1,1)
            if(max(Speed(TSE1(k,1):TSE1(k,2)))<=SpeedTh && max(ActiveAnimal(TSE1(k,1):TSE1(k,2)))==1 )
                j=0;
                while group1(TSE1(k,1)+j)==0 && ((TSE1(k,1)+j) < TSE1(k,2))
                    j=j+1;
                end
                
                stimL=0;
                for i=0:LMax
                    if group1(TSE1(k,1)+j+i)==0 && ((TSE1(k,1)+j+i) < TSE1(k,2))
                        stimL=i;
                        break;
                    end
                end
                
                if (stimL>0)
                    tc1=tc1+1;
                    tv1(:,tc1)=X(:,TSE1(k,1)+j+td);
                end
            end
        end
        tv1=tv1(:,1:tc1);
        
            
        vc=min(vc0,vc1);
        
        dv0=dv0(:,1:vc);
        dv1=dv1(:,1:vc);

        groupV=[zeros(vc,1);ones(vc,1)];
        %%%%%%%%%%PLS
        dv0=PLSRotation' * dv0 ;
        dv1=PLSRotation' * dv1 ;
        %%%%%%%%%%%%%%%%%
        
        md0=sum(dv0,2)/vc;
        md1=sum(dv1,2)/vc;
        
        S0=((dv0-md0*ones(1,vc))*(dv0-md0*ones(1,vc))')/(vc-1);
        S1=((dv1-md1*ones(1,vc))*(dv1-md1*ones(1,vc))')/(vc-1);
        
        S=(S0+S1)/2;
        dm=md1-md0;
        
        
        IFisherPLS(area,indtd)= dm' * inv(S) * dm;
 
        
        
        
    end

    
end

end


