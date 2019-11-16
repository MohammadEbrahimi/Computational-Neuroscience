cellEvents=cellData;%padding(eventBin,0,2);
set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',12);

Hit=zeros(size(GoTrials));
Miss=zeros(size(GoTrials));
CR=zeros(size(GoTrials));
FA=zeros(size(GoTrials));
nH=0;nM=0;nC=0;nF=0;
clear HitSE MissSE FASE CRSE
for i=1:size(GoSE,1)
    if(max(WaterReward(GoSE(i,1):GoSE(i,2)))==1)
        Hit(GoSE(i,1):GoSE(i,2))=GoTrials(GoSE(i,1):GoSE(i,2));
        nH=nH+1;
        HitSE(nH,:)=GoSE(i,:);
    else
        Miss(GoSE(i,1):GoSE(i,2))=GoTrials(GoSE(i,1):GoSE(i,2));
        nM=nM+1;
        MissSE(nM,:)=GoSE(i,:);
    end
end
Punish=RewardWindow & AirPuff;
for i=1:size(NoGoSE,1)
    if(max(Punish(NoGoSE(i,1):NoGoSE(i,2)+1))==1)
        FA(NoGoSE(i,1):NoGoSE(i,2))=NogoTrials(NoGoSE(i,1):NoGoSE(i,2));
        nF=nF+1;
        FASE(nF,:)=NoGoSE(i,:);
    else
        CR(NoGoSE(i,1):NoGoSE(i,2))=NogoTrials(NoGoSE(i,1):NoGoSE(i,2));
        nC=nC+1;
        CRSE(nC,:)=NoGoSE(i,:);
    end
end

%%%%Find cell Coordinates and assign area
cellIJ=zeros(cellCount,2);
CortexArea=zeros(cellCount,1);

for i=1:cellCount
    [row,ind]=max(cellImage{i});
    [~,indj]=max(row);
    indi=ind(indj);
    cellIJ(i,:)=TileTopLeft(i,:)+[indi-1,indj-1]+[5,5];
end

for b=1:8
    a=9-b;
    vertex=cortexMap{a};
    [in,on] = inpolygon(cellIJ(:,1),cellIJ(:,2),vertex(:,1),vertex(:,2));
    CortexArea(in==1)=a;
end

%%%% creating data sets
group0=CR;
group1=Hit;
LMax=20;
avgWin=0;
tdVec=-5:19;
IFisher=zeros(8,length(tdVec));
errorD=zeros(8,length(tdVec));
errorVal=zeros(8,length(tdVec));
errorDVar=zeros(8,length(tdVec));
BMAT=zeros(cellCount,length(tdVec),8);
InterceptMAT=zeros(8,length(tdVec));
Signal=zeros(8,length(tdVec));
Noise=zeros(8,length(tdVec));

for area=1:3
for indtd=1:length(tdVec)
    
td=tdVec(indtd);
c0=0;
c1=0;

%%%%Active sesseions
% SE0=CRSE(1:83,:);
% SE1=HitSE(1:114,:);
% SE0=NoGoSE(1:110,:);
% SE1=GoSE(1:119,:);
%%%%% Inactive session
SE0=CRSE;
SE1=HitSE;

X=cellEvents(CortexArea==area,:);
%X=cellEvents;

dd0=zeros(size(X,1),size(SE0,1));
for k=2:length(SE0)
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
    
    if (td<stimL && stimL>0)
        c0=c0+1;
        t0(c0)=SE0(k,1)+j+td;
        dd0(:,c0)=sum(X(:,SE0(k,1)+j+td-avgWin:SE0(k,1)+j+td),2);
    end
end

dd1=zeros(size(X,1),size(SE1,1));
for k=2:length(SE1)
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
    
    if (td<stimL && stimL>0)
        c1=c1+1;
        t1(c1)=SE1(k,1)+j+td;
        dd1(:,c1)=sum(X(:,SE1(k,1)+j+td-avgWin:SE1(k,1)+j+td),2);
    end
end

v=round(min(c0,c1)/3);
c=min(c0,c1)-v;
shuff=randperm(c0);
trials0=[ones(c,1);zeros(c0-c,1)];
trials0=(trials0(shuff)==1);
shuff=randperm(c1);
trials1=[ones(c,1);zeros(c1-c,1)];
trials1=(trials1(shuff)==1);


%%%%%%%%%%%%% Decoder
DataSet=[dd0(:,trials0) dd1(:,trials1)];
group=[zeros(c,1);ones(c,1)];
area
[B,intercept,ErrCurve,ErrVar,Lambda]=lassoglmcv(DataSet(:,:),group,5);

fullB=zeros(cellCount,1);
fullB(CortexArea==area)=B;
BMAT(:,indtd,area)=fullB;
InterceptMAT(area,indtd)=intercept;
%%%% Information should be calculated on a seperate validation set
% v=min(c0,c1);
% dv0=dd0(:,1:v);
% dv1=dd1(:,1:v);
% B=BMAT(:,indtd,area);

dv0=dd0(:,(trials0==0));
dv0=dv0(:,1:v);
dv1=dd1(:,(trials1==0));
dv1=dv1(:,1:v);
groupV=[zeros(v,1);ones(v,1)];

md0=sum(dv0,2)/v;
md1=sum(dv1,2)/v;

S0=(dv0-md0*ones(1,v))*(dv0-md0*ones(1,v))';
S1=(dv1-md1*ones(1,v))*(dv1-md1*ones(1,v))';

S=(S0+S1)/(2*(v-1));
dm=md1-md0;


IFisher(area,indtd)= 2*(B' * dm )^2 /(B'*(S1+S0)*B); 
Signal(area,indtd)= 2*(B' * dm )^2;
Noise(area,indtd)= (B'*(S1+S0)*B);

errorD(area,indtd)=min(ErrCurve);
errorDVar(area,indtd)=ErrVar;

errorVal(area,indtd)= mean( LRClassify([dv0' ; dv1'] , B , intercept) ~= groupV) 
indtd
end

end


%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%Visualization
IFisher(isnan(IFisher))=0;
%IFisher_i(isnan(IFisher_i))=0;

A1=eye(25)/3;
for i=1:24
A1(i,i+1)=1/3;
A1(i+1,i)=1/3;
end
A1(1:2,1)=0.5;
A1(24:25,25)=0.5;

td=-0.4:0.1:2;
color={'b','k','m','r','r','m','k','b'};
  TitleV={'V1','LV','MV','PTLP','A','S','M','RSC'};
  figure();
    hold on;
for area=1:4
 
    vect=(100*IFisher(area,:)*A1/max(IFisher(area,:)*A1));
    %vect=IFisher(area,:)*A1;
    %vect2=IFisher_i(area,:)*A1;
    %vect=errorD(area,:);
    plot(td,vect,color{area});
    %plot(td,vect2,'k--');
    %title(TitleV(area));
    xlabel('time (s)');
    ylabel('Fisher Information');
end

legend('V1','LV','MV','PTLP');
%legend('S','M','RSC');
grid on
















