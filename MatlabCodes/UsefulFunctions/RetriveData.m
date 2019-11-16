%function [X,CortexArea,HitC,MissC,CRC,FAC,Delay,RewardWindow,Speed]=RetriveData(path,mode)
    function [X,CortexArea,HitC,MissC,CRC,FAC,HitSE,MissSE,CRSE,FASE]=RetriveData(path,mode)
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


X=cellData;%padding(eventBin,0,2);

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


HitC=zeros(size(Hit));
MissC=zeros(size(Hit));
CRC=zeros(size(Hit));
FAC=zeros(size(Hit));

%%%%%%%%%%%%%%%%%%%%%%%%%%%Start Learning



    

    
    c0=1;
    c1=1;
    
    for k=1:size(SE0,1)
        
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
                CRC(SE0(k,1)+j:SE0(k,1)+j+stimL)=1;
                c0=c0+stimL+1;
            end
        end
    end
    c0=c0-1;

    
    

    for k=1:size(SE1,1)

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
                HitC(SE1(k,1)+j:SE1(k,1)+j+stimL)=1;
                c1=c1+stimL+1;
            end
        end
    end
    c1=c1-1;
    
    c=round(min(c0,c1)/2);

 
    tc0=1;
    tc1=1;
    
    for k=1:size(TSE0,1)
        
        if(max(Speed(TSE0(k,1):(TSE0(k,1)+55)))<=SpeedTh && max(ActiveAnimal(TSE0(k,1):TSE0(k,2)))==1 ...
                )%% && max(TimeOut(SE0(k,1):(SE0(k,1)+55))) ==0 )
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
            stimL=stimL-1;
            if ((TSE0(k,1)+j+i) <= TSE0(k,2)+1  && stimL>0)
                FAC(TSE0(k,1)+j:TSE0(k,1)+j+stimL)=1;
                tc0=tc0+stimL+1;
            end
        end
    end
    tc0=tc0-1;

    
    

    for k=1:size(TSE1,1)

        if(max(Speed(TSE1(k,1):(TSE1(k,1)+55)))<=SpeedTh && max(ActiveAnimal(TSE1(k,1):TSE1(k,2)))==1 ...
                )%%&& max(TimeOut(SE1(k,1):(SE1(k,1)+55))) ==0 )
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
            
            stimL=stimL-1;
            if ((TSE1(k,1)+j+i) <= TSE1(k,2)+1 && stimL>0)
                MissC(TSE1(k,1)+j:TSE1(k,1)+j+stimL)=1;
                tc1=tc1+stimL+1;
            end
        end
    end
    tc1=tc1-1;
    
    
    
 
end