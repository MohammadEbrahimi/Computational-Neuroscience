cellEvents=padding(eventBin,0,2);
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
group0=NogoTrials;
group1=GoTrials;
LMax=20;
avgWin=0;
tdVec=-5:19;
IFisher=zeros(8,length(tdVec));
NDim=zeros(8,1);
Ncells=zeros(8,1);
eigMul=zeros(8,length(tdVec));
eigSum=zeros(8,length(tdVec));
BMAT=zeros(cellCount,8);
InterceptMAT=zeros(8,1);
Signal=zeros(8,length(tdVec));
Noise=zeros(8,length(tdVec));

for area=1:8
    
    
    
    
    %%%%Active sesseions
%     SE0=CRSE(1:83,:);
%     SE1=HitSE(1:114,:);
%     SE0=NoGoSE(111:end,:);
%     SE1=GoSE(120:end,:);
    %%%%% Inactive session
    SE0=CRSE(2:end,:);
    SE1=HitSE(2:end,:);
    
    X=cellEvents(CortexArea==area,:);
    %X=cellEvents;
    
    shuff0=randperm(size(SE0,1));
    shuff1=randperm(size(SE1,1));
    
    c0=1;
    c1=1;
    
    clear shuffInd0
    dd0=zeros(size(X,1),length(group0));
    for si=1:length(shuff0)
        
        k=shuff0(si);
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
            dd0(:,c0:c0+stimL-5)=X(:,SE0(k,1)+j+5:SE0(k,1)+j+stimL);
            shuffInd0(c0:c0+stimL-5)=si;
            c0=c0+stimL+1-5;
        end
    end
    c0=c0-1;

    clear shuffInd1
    dd1=zeros(size(X,1),length(group1));
    for si=1:length(shuff1)
        k=shuff1(si);
        j=0;
        while group1(SE1(k,1)+j)==0 && ((SE1(k,1)+j) < SE1(k,2))
            j=j+1;
        end
        si
        stimL=0;
        for i=0:LMax
            if group1(SE1(k,1)+j+i)==0 && ((SE1(k,1)+j+i) < SE1(k,2))
                stimL=i;
                break;
            end
        end
        
        stimL=stimL-1
        if ((SE1(k,1)+j+i) < SE1(k,2) )
            dd1(:,c1:c1+stimL-5)=X(:,SE1(k,1)+j+5:SE1(k,1)+j+stimL);
            shuffInd1(c1:c1+stimL-5)=si;
            c1=c1+stimL+1-5;
        end
    end
    c1=c1-1;
    
    
    c=min(c0,c1);
    dd0=dd0(:,1:c);
    
    dd1=dd1(:,1:c);
    % %%%%%%%%%%%%% Decoder
    md0=sum(dd0,2)/c;
    md1=sum(dd1,2)/c;
        
S0=(dd0-md0*ones(1,c))*(dd0-md0*ones(1,c))';
S1=(dd1-md1*ones(1,c))*(dd1-md1*ones(1,c))';

S=(S0+S1)/(2*(c-1));
dm=md1-md0;

[U,D,V]=svd(S);
lamb= diag(D);

threshold= max(2.858 * median(lamb),1e-3);
for k=1:size(S,1)
    if (lamb(k) < threshold && sum ((lamb(1:k))/sum(lamb)) >0.99)
        NDim(area) = k-1;
        break
    end
end

Ncells(area)=sum(CortexArea==area);
% %%%% Information should be calculated on a seperate validation set
    
    VSE0=SE0;
    VSE1=SE1;
    Valsize=min(size(VSE0,1),size(VSE1,1));
    
    for indtd=1:length(tdVec)
        
        td=tdVec(indtd);
        vc0=0;
        vc1=0;
        
        dv0=zeros(size(X,1),Valsize);
        for k=1:Valsize
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
        
        dv1=zeros(size(X,1),Valsize);
        for k=1:Valsize
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
        
        
        
        vc=min(vc0,vc1)
        
        dv0=dv0(:,1:vc);
        dv1=dv1(:,1:vc);
        
        md0=sum(dv0,2)/vc;
        md1=sum(dv1,2)/vc;
        
        dm=md1-md0;
        
        gm= V'*dm;
        
        IFX=(gm .^2)./lamb;
        
       eigMul(area,indtd)=prod(lamb);
       eigSum(area,indtd)=sum(lamb); 
       IFisher(area,indtd)= sum(IFX(1:NDim(area))); 
        
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


IFisher(isnan(IFisher))=0;

L=25
A=eye(L)/3;
for i=1:L-1
A(i,i+1)=1/3;
A(i+1,i)=1/3;
end
A(1:2,1)=0.5;
A(L-1:L,L)=0.5;


td=-0.4:0.1:2;
color={'b','b--','k--','r','g','m','k','y'};
  TitleV={'V1','LV','MV','PTLP','A','S','M','RSC'};
  figure();
    hold on;
for area=1:8

    vect=(100*IFisher(area,:)*A/max(IFisher(area,:)*A));
    %vect=IFisher(area,:)*A;
    %vect2=IFisher_i(area,:)*A1;
    %vect=errorD(area,:);
    plot(td,vect,color{area});
    %plot(td,vect2,'k--');
    title('L362');
    xlabel('time (s)');
    ylabel('validation');
end

legend('V1','LV','MV','PTLP','A','S','M','RSC');
grid on







