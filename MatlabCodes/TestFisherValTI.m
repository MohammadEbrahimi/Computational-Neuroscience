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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
group0=RewardWindow;
group1=RewardWindow;
LMax=30;
avgWin=0;
tdVec=-5:29;
IFisher_test=zeros(8,length(tdVec));
errorVal_test=zeros(8,length(tdVec));
FA_test=zeros(8,length(tdVec));
Miss_test=zeros(8,length(tdVec));
Signal_test=zeros(8,length(tdVec));
Noise_test=zeros(8,length(tdVec));

 VSE0=FASE(2:end,:);
 VSE1=HitSE(2:end,:);
 Valsize=min(size(VSE0,1),size(VSE1,1));
 
for area=1:8   
    X=cellEvents(CortexArea==area,:);
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
                if group0(VSE0(k,1)+j+i)==0 && ((VSE0(k,1)+j+i) < VSE0(k,2)+3)
                    stimL=i+1;
                    break;
                end
            end
            
            if (td<stimL)
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
                if group1(VSE1(k,1)+j+i)==0 && ((VSE1(k,1)+j+i) < VSE1(k,2)+3)
                    stimL=i+1;
                    break;
                end
            end
            
            if (td<stimL)
                vc1=vc1+1;
                dv1(:,vc1)=X(:,VSE1(k,1)+j+td);
            end
        end
        
        
        B= BMAT((CortexArea==area),area);
        intercept= InterceptMAT(area); 
        vc=min(vc0,vc1);
        
        dv0=dv0(:,1:vc);
        dv1=dv1(:,1:vc);
        groupV=[zeros(vc,1);ones(vc,1)];
        
        md0=sum(dv0,2)/vc;
        md1=sum(dv1,2)/vc;
        
        S0=(dv0-md0*ones(1,vc))*(dv0-md0*ones(1,vc))';
        S1=(dv1-md1*ones(1,vc))*(dv1-md1*ones(1,vc))';
        
        S=(S0+S1)/(2*(vc-1));
        dm=md1-md0;
        
        
        IFisher_test(area,indtd)= 2*(B' * dm )^2 /(B'*(S1+S0)*B);
        Signal_test(area,indtd)= 2*(B' * dm )^2;
        Noise_test(area,indtd)= (B'*(S1+S0)*B);
        
        
        errorVal_test(area,indtd)= mean( LRClassify([dv0' ; dv1'] , B , intercept) ~= groupV);
        FA_test(area,indtd)= mean( LRClassify([dv0'] , B , intercept));
        Miss_test(area,indtd)= 1-mean( LRClassify([dv1'] , B , intercept));
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L=25
% A=eye(L)/3;
% for i=1:L-1
% A(i,i+1)=1/3;
% A(i+1,i)=1/3;
% end
% A(1:2,1)=0.5;
% A(L-1:L,L)=0.5;
% 
% 
% td=-0.4:0.1:2;
% color={'b','k','m','r','r','m','k','b'};
%   TitleV={'V1','LV','MV','PTLP','A','S','M','RSC'};
%   figure();
%     hold on;
% for area=1:4
%     %vect=(errorVal_test(area,:)*A);
%     %vect=(100*IFisher(area,:)*A/max(IFisher(area,:)*A));
%     vect=IFisher_test(area,:)*A;
%     %vect2=IFisher_i(area,:)*A1;
%     %vect=errorD(area,:);
%     plot(td,vect,color{area});
%     %plot(td,vect2,'k--');
%     title('L347 inactive sessions during stim');
%     xlabel('time (s)');
%     ylabel('Validation Error');
% end
% 
% legend('V1','LV','MV','PTLP');
% %legend('A','S','M','RSC');
% grid on



















