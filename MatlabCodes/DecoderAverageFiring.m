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

Delay=zeros(size(GoTrials));
for i=2:length(Delay)
    if (RewardWindow(i)==1 && RewardWindow(i-1)==0)
        Delay(i-5:i-1)=1;
    end
end


%%%% creating data sets
group0=Delay;
group1=Delay;
LMax=5;
avgWin=0;


%%%initialization
IFisher=zeros(8,1);
errorD=zeros(8,1);
errorVal=zeros(8,1);
errorDVar=zeros(8,1);
BMAT=zeros(cellCount,8);
InterceptMAT=zeros(8,1);
Signal=zeros(8,1);
Noise=zeros(8,1);
% FARatio=zeros(8,3);
% MissRatio=zeros(8,3);

for area=1:8
    
    
    
    
    %     %%%%Active sesseions
        SE0=CRSE(1:83,:);
        SE1=HitSE(1:114,:);
    %     SE0=NoGoSE(111:end,:);
    %     SE1=GoSE(120:end,:);
    %%%%% Inactive session
%     SE0=CRSE(2:end,:);
%     SE1=HitSE(2:end,:);
%     
    X=cellEvents(CortexArea==area,:);
    %X=cellEvents;
    
    shuff0=randperm(size(SE0,1));
    shuff1=randperm(size(SE1,1));
    
    c0=1;
    c1=1;
    
    dd0=zeros(size(X,1),size(SE0,1));
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
            dd0(:,c0)=mean(X(:,SE0(k,1)+j:SE0(k,1)+j+stimL),2);
            c0=c0+1;
        end
    end
    c0=c0-1;
    dd0=dd0(:,1:c0);
    
    dd1=zeros(size(X,1),size(SE1,1));
    for si=1:length(shuff1)
        k=shuff1(si)
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
            dd1(:,c1)=mean(X(:,SE1(k,1)+j:SE1(k,1)+j+stimL),2);
            c1=c1+1;
        end
    end
    c1=c1-1;
    dd1=dd1(:,1:c1);
    
    
    v=round(min(c0,c1)/2);
    c=min(c0,c1)-v;
    
    % %%%%%%%%%%%%% Decoder
    DataSet=[dd0(:,1:c) dd1(:,1:c)];
    group=[zeros(c,1);ones(c,1)];
    area
    [B,intercept,ErrCurve,ErrVar,Lambda]=lassoglmcv(DataSet,group,5,5,1,1);
    
    fullB=zeros(cellCount,1);
    fullB(CortexArea==area)=B;
    BMAT(:,area)=fullB;
    InterceptMAT(area)=intercept;
    
    errorD(area)=min(ErrCurve);
    %%%% Information should be calculated on a seperate validation set
    
    VSE0=SE0(shuff0(c+1:end),:);
    VSE1=SE1(shuff1(c+1:end),:);
    Valsize=min(size(VSE0,1),size(VSE1,1));
    
    
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
        
        if (stimL>0)
            vc0=vc0+1;
            dv0(:,vc0)=mean(X(:,VSE0(k,1)+j:VSE0(k,1)+j+stimL),2);
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
        
        if (stimL>0)
            vc1=vc1+1;
            dv1(:,vc1)=mean(X(:,VSE1(k,1)+j:VSE1(k,1)+j+stimL),2);
        end
    end
        
        
   TSE0=FASE;
   TSE1=MissSE;
    tc0=0;
    tc1=0;
    
    tv0=zeros(size(X,1),size(TSE0,1));
    for k=1:size(TSE0,1)
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
            tv0(:,tc0)=mean(X(:,TSE0(k,1)+j:TSE0(k,1)+j+stimL),2);
        end
    end
    
    tv1=zeros(size(X,1),size(TSE1,1));
    for k=1:size(TSE1,1)
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
            tv1(:,tc1)=mean(X(:,TSE1(k,1)+j:TSE1(k,1)+j+stimL),2);
        end
    end
        
 B=BMAT(CortexArea==area,area);
FARatio(area,2)=  mean( LRClassify([tv0'] , B , InterceptMAT(area)) ~= 0)/...
    mean( LRClassify([dv0'] , B, InterceptMAT(area)) ~= 0); 

MissRatio(area,2)= mean( LRClassify([tv1'] , B , InterceptMAT(area)) ~= 1)/...
    mean( LRClassify([dv1'] , B , InterceptMAT(area)) ~= 1);
                    

        
        vc=min(vc0,vc1)
        
        dv0=dv0(:,1:vc);
        dv1=dv1(:,1:vc);
        groupV=[zeros(vc,1);ones(vc,1)];
        
        md0=sum(dv0,2)/vc;
        md1=sum(dv1,2)/vc;
        
        S0=(dv0-md0*ones(1,vc))*(dv0-md0*ones(1,vc))';
        S1=(dv1-md1*ones(1,vc))*(dv1-md1*ones(1,vc))';
        
        S=(S0+S1)/(2*(vc-1));
        dm=md1-md0;
        
        
        IFisher(area)= 2*(B' * dm )^2 /(B'*(S1+S0)*B);
        Signal(area)= 2*(B' * dm )^2;
        Noise(area)= (B'*(S1+S0)*B);
        
        
        errorVal(area)= mean( LRClassify([dv0' ; dv1'] , B , intercept) ~= groupV);
        
end






% [sortedW cellIndex]=sort(B,'descend');
% 
% figure();hold on
% plot(GoTrials);
% plot(NogoTrials,'r');
% plot(X(cellIndex(1),:),'k')
%
% save('C:\Users\Sadegh\Documents\VLMReborn\Reports\AvergeFiringDecoder\during Stim\L354_8_3_2015_all_ev',...
%    'IFisher','Signal','Noise','errorVal','BMAT','InterceptMAT');
% % % %%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%
% % % % %%%%%%%%%%%%%%%%%%Visualization
% cellEvents=padding(eventBin,0,2);
% area=2;
% W=BMAT(:,area);
% group=GoTrials;
% Lmax=20;
% 
% Raster=zeros(sum(W>0),Lmax+5);
% Counter=zeros(sum(W>0),Lmax+5);
% 
% X=cellEvents(W>0,:);
% [s,ind]=sort(W(W>0),'descend');
% X=X(ind,:);
% 
% for i=1:length(GoTrials)-1
%     if (group(i)==0 && group(i+1)==1)
%         Raster(:,1:5)=Raster(:,1:5)+X(:,i-4:i);
%         Counter(:,1:5) = Counter(:,1:5) + 1;
%         
%         for j=1:Lmax
%             Raster(:,j+5)=Raster(:,j+5)+X(:,i+j);
%             Counter(:,j+5) = Counter(:,j+5) + 1;
%             if (group(i+j+1)==0)
%                 break;
%             end
%         end
%     end
%     
% end
% 
% Raster = Raster ./ Counter; 
% figure();imagesc(Raster);
% 
% 
% 

