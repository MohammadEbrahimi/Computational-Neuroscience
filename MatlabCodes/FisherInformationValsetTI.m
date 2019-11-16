
% NNEvents=zeros(size(cellData));
% V.dt=0.1;
% V.Ncells=1;
% for i=1:cellCount
%     NNEvents(i,:)=(fast_oopsi(cellData(i,:),V))';
% end

cellEvents=risingEvents(cellData,eventBin,0.7);%padding(eventBin,0,2);
cellEvents(cellEvents<0)=0;
set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',12);
Speed = sqrt(XSpeed.^2 + YSpeed.^2);

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

Delay=zeros(size(GoTrials));
for i=2:length(Delay)
    if (RewardWindow(i)==1 && RewardWindow(i-1)==0)
        Delay(i-5:i-1)=1;
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
group0=Delay;
group1=Delay;
LMax=5;
avgWin=0;
tdVec=-5:4;
IFisher=zeros(8,length(tdVec));
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

for area=1:8
    
    
    
    
    %%%%Active sesseions
        SE0=CRSE(2:83,:);
        SE1=HitSE(2:114,:);
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
    
    clear shuffInd0
    dd0=zeros(size(X,1),length(group0));
    for si=1:length(shuff0)
        
        k=shuff0(si);
        if(max(Speed(SE0(k,1):SE0(k,2)))<=10)
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
        k=shuff1(si)
        if(max(Speed(SE1(k,1):SE1(k,2)))<=10)
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
    [B,intercept,ErrCurve,ErrVar,Lambda]=lassoglmcv(DataSet,group,shuffInd,5,10,1);
    
    fullB=zeros(cellCount,1);
    fullB(CortexArea==area)=B;
    BMAT(:,area)=fullB;
    InterceptMAT(area)=intercept;
    errorD(area)=min(ErrCurve);
%     % %%%% Information should be calculated on a seperate validation set
     B=BMAT(CortexArea==area,area);
     intercept=InterceptMAT(area);
    
    
    
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
            if(max(Speed(VSE0(k,1):VSE0(k,2)))<=10)
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
            if(max(Speed(VSE1(k,1):VSE1(k,2)))<=10)
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
            if(max(Speed(TSE0(k,1):TSE0(k,2)))<=10)
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
            if(max(Speed(TSE1(k,1):TSE1(k,2)))<=10)
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
        
        
        IFisher(area,indtd)= 2*(B' * dm )^2 /(B'*(S1+S0)*B);
        Signal(area,indtd)= 2*(B' * dm )^2;
        Noise(area,indtd)= (B'*(S1+S0)*B);
        
        
        errorVal(area,indtd)= mean( LRClassify([dv0' ; dv1'] , B , intercept) ~= groupV);
        
    end
    
    FARatio(area,:)=[fof/FACount,foc/CRCount];
    MissRatio(area,:)=[mom/MissCount,moh/HitCount];
    
end
% 
% save('C:\Users\Sadegh\Documents\VLMReborn\Reports\2016_3_3_FisherInformationTI\DuringStim\HitCR\L362_8_3_2015_all_rev',...
%   'IFisher','Signal','Noise','errorVal','BMAT','InterceptMAT','FAonFA','FAonCR','MissonMiss','MissonHit');
% % % % %%%%%%%%%%%%%%%%%%
% % % %%%%%%%%%%%%%%%%%%
% % % % %%%%%%%%%%%%%%%%%%Visualization
% IFisher(isnan(IFisher))=0;
% %IFisher_i(isnan(IFisher_i))=0;
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

% 
%  
% for area=1:2
%      figure();
%     hold on;
%     %vect1=(errorVal(area,:));
%     %vect=(100*IFisher(area,:)*A/max(IFisher(area,:)*A));
%      vect1=errorVal(area,:);
%      vect2=EV62(area,:);
%     %vect=errorD(area,:);
%     plot(td,vect1,'r');
%      plot(td,vect2,'b');
%     %plot(td,vect2,'k--');
%     title(TitleV{area});
%     xlabel('time (s)');
%     ylabel('validation');
% end

% legend('V1','LV','MV','PTLP','A','S','M','RSC');
% grid on
%IFisher=[IFisher_S,IFisher_D(:,6:10),IFisher_R(:,6:35)];


%
%
% IFisher_62=[IFisher_S,IFisher_D(:,6:10),IFisher_R(:,6:35)];
% errorVal_62=[errorVal_S,errorVal_D(:,6:10),errorVal_R(:,6:35)];
%
% D1=FF{1};
% D2=FF{2};
% D3=FF{3};
% 
% M1=FC{1};
% M2=FC{2};
% M3=FC{3};
% 
% 
% for a=1:8
%  figure();
%     hold on
% 
%     v1=(D1(a,:));% / max(D1(a,:));
%     v2=(D2(a,:));% / max(D2(a,:));
%     v3=(D3(a,:));%/ max(D3(a,:));
%     x=td;
%     y=(v1+v2+v3)/3;
%     errBar=zeros(2,length(x));
%     errBar(1,:)=max([v1;v2;v3]) - y;
%     errBar(2,:)=y-min([v1;v2;v3]);
%     shadedErrorBar(x,y,errBar,'r',1)
%     
% %     plot(x,y,color{a});
% %     plot(zeros(1,40),0.025:0.025:1,'b.')
% %     plot(2*ones(1,40),0.025:0.025:1,'r.')
% %     plot(2.5*ones(1,40),0.025:0.025:1,'g.')
% 
%     v1=M1(a,:);% / max(IF47(a,:));
%     v2=M2(a,:);% / max(IF54(a,:));
%     v3=M3(a,:);% / max(IF62(a,:));
%     x=td;
%     y=(v1+v2+v3)/3;
%     errBar=zeros(2,length(x));
%     errBar(1,:)=max([v1;v2;v3]) - y;
%     errBar(2,:)=y-min([v1;v2;v3]);
%     shadedErrorBar(x,y,errBar,'b',1)
% 
%     xlabel('time (s)');
%     ylabel('Validation error');
%     title(TitleV{a});
%     grid on
% end
% 
% 
% % 
% % 
% % 
% i=2;
% EV{i}=errorVal;
% IF{i}=IFisher;
% FF{i}=FAonFA;
% FC{i}=FAonCR;
% MM{i}=MissonMiss;
% MH{i}=MissonHit;
% 
% 




