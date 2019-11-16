cellEvents=padding(eventBin,0,2);
set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',12);
%%%%Find cell Coordinates and assign area
cellIJ=zeros(cellCount,2);
CortexArea=zeros(cellCount,1);

for i=1:cellCount
    [row,ind]=max(cellImage{i});
    [~,indj]=max(row);
    indi=ind(indj);
    cellIJ(i,:)=TileTopLeft(i,:)+[indi-1,indj-1];
end

for a=1:8
    vertex=cortexMap{a};
    [in,on] = inpolygon(cellIJ(:,1),cellIJ(:,2),vertex(:,1),vertex(:,2));
    CortexArea(in==1)=a;
end

%%%% creating data sets
group0=NogoTrials;
group1=GoTrials;
LMax=20;
tdVec=-5:20;
IFisher=zeros(8,length(tdVec));

area=1;
for indtd=1:length(tdVec)
    
td=tdVec(indtd);
c0=0;
c1=0;

SE0=NoGoSE;
SE1=GoSE;

X=cellEvents(CortexArea==area,:);

dd0=zeros(size(X,1),size(SE0,1));
for k=2:length(SE0)
    j=0;
    while group0(SE0(k,1)+j)==0
        j=j+1;
    end
    
    stimL=0;
    for i=0:LMax
        if group0(SE0(k,1)+j+i)==0
            stimL=i;
            break;
        end
    end
    
    if (td<stimL)
        c0=c0+1;
        dd0(:,c0)=sum(X(:,SE0(k,1)+j+td-1:SE0(k,1)+j+td+1),2);
    end
end

dd1=zeros(size(X,1),size(SE1,1));
for k=2:length(SE1)
    j=0;
    while group1(SE1(k,1)+j)==0
        j=j+1;
    end
    
    stimL=0;
    for i=0:LMax
        if group1(SE1(k,1)+j+i)==0
            stimL=i;
            break;
        end
    end
    
    if (td<stimL)
        c1=c1+1;
        dd1(:,c1)=sum(X(:,SE1(k,1)+j+td-1:SE1(k,1)+j+td+1),2);
    end
end

c=min(c0,c1);
shuff=randperm(c0);
trials0=[ones(c,1);zeros(c0-c,1)];
trials0=(trials0(shuff)==1);
shuff=randperm(c1);
trials1=[ones(c,1);zeros(c1-c,1)];
trials1=(trials1(shuff)==1);

md0=sum(dd0(:,trials0),2)/c;
md1=sum(dd1(:,trials1),2)/c;

S0=(dd0(:,trials0)-md0*ones(1,c))*(dd0(:,trials0)-md0*ones(1,c))';
S1=(dd1(:,trials1)-md1*ones(1,c))*(dd1(:,trials1)-md1*ones(1,c))';

S=(S0+S1)/(2*(c-1));
dm=md1-md0;


DataSet=[dd0(:,trials0) dd1(:,trials1)];
group=[zeros(c,1);ones(c,1)];
[B,intercept,MissCurve,FalseAlarmCurve,Lambda]=lassoglmcv(DataSet,group,5);

%%%% Information should be calculated on a seperate validation set

IFisher(area,indtd)= 2*(B' * dm )^2 /(B'*(S1+S0)*B) 

end



