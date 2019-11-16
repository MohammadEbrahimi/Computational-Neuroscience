 path='C:\Users\Sadegh\Documents\VLMReborn\L354\Data\allDays';
 mode=1;
 train=1;
%function Fisher_PLS(path,savePath,mode,train,TSin)

load(strcat(path,'\All_Sessions.mat'));
load(strcat(path,'\areas.mat'));
load(strcat(path,'\cellImage.mat'));

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

%%%%%%%%%%%
X=cellEvents(CortexArea==6,:);
N=size(X,1);

TrialSE=[HitSE;MissSE;CRSE;FASE];
trial=NogoTrials+GoTrials;
L=20;
K=size(TrialSE,1);

T=zeros(N,L+6,K);
kp=1;

for k=1:K
    os=0;
    while trial(TrialSE(k,1)+os)==0 && TrialSE(k,1)+os<=TrialSE(k,2)
        os=os+1;
    end
    if trial(TrialSE(k,1)+os)==0 || TrialSE(k,1)<5 continue; end
    
    T(:,:,kp)=X(:,TrialSE(k,1)+os-5:TrialSE(k,1)+os+L);
    
    if k<size(HitSE,1) C(kp)=1;
    elseif k<size([HitSE;MissSE],1) C(kp)=2;
    elseif k<size([HitSE;MissSE;FASE],1) C(kp)=3;
    else C(kp)=4;
    end  
    kp=1+kp;
end

kp=kp-1;
T=T(:,:,1:kp);
C=C(1:kp);
%R=rankest(T);
R=10;
Uest=cpd(T,R);

[SCA,NAreaIndex]=sort(CortexArea);

TimeFactor=Uest{2};
NeuronFactor=Uest{1};
TrialFactor=Uest{3};
tk=1:kp;

color={'b','c','c','r','g','m','k','y'};
colord={'b.','m.','r.','k.'};
set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',16);
set(0,'defaultAxesFontName','Calibri');
for i=1:R
    figure();
    title(strcat('L347_',int2str(i)));
    subplot(1,3,1);
    plot(-0.5:0.1:2,TimeFactor(:,i));
    grid on
    subplot(1,3,2);
    hold on
%     for a=1:8
%         plot((SCA==a).*NeuronFactor(NAreaIndex,i),color{a});
%     end
plot(NeuronFactor(:,i),color{1});
    subplot(1,3,3);
    hold on
    for c=1:4
        plot(tk(C==c),TrialFactor(C==c,i),colord{c});
    end
end
    
    
    








