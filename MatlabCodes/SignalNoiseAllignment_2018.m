
AddressSetup;

for  Mouse=[61,63:67];
    loadpath=LoadPath{Mouse-10};
    
    
    
    SET0Name='C';
    SET1Name='H';
    CellT='IB';
    mode=1;
    TSin=200;
    SpeedMode=2;
    
    Balanced=0;
    NPC=10;
    NRepeate=100;
    ActiveTrialNumber=20;
    SpeedBalance=0;
    ActiveMode=1;
    
    
    % mycluster=parcluster('local');
    % mycluster.NumWorkers=8;
    % saveProfile(mycluster);
    %
    % poolobj=parpool(mycluster,8)
    
    warning('off','all');
    load(strcat(loadpath,'/All_Sessions.mat'));
    load(strcat(loadpath,'/areas.mat'));
    load(strcat(loadpath,'/cellData.mat'));
    load(strcat(loadpath,'/cellData_ZS.mat'));
    load(strcat(loadpath,'/SessionLength.mat'));
    load(strcat(loadpath,'/Datasets/Datasets--5to0.mat'));
    Speed=Speed-min(Speed);
    progress='Data is loaded!'
    
    if mode==1
        TimeVec=1:25;
    elseif mode==2
        TimeVec=1:5;
    elseif mode==3
        TimeVec=1:30;
        
    end
    
    
    switch CellT
        case 'R'
            cellEvents=cellData_Raw;
        case 'I'
            cellEvents=cellData_ROI;
        case 'O'
            load(strcat(loadpath,'/cellData_oopsi.mat'));
            cellEvents=cellData_oopsi;
        case 'RB'
            cellEvents=cellData_Raw_bin;
        case 'IB'
            cellEvents=cellData_ROI_bin;
        case 'RRT'
            cellEvents=cellData_Raw;
            cellEvents(cellData_Raw_bin==0)=0;
        case 'IRT'
            cellEvents=cellData_ROI;
            cellEvents(cellData_ROI_bin==0)=0;
            
    end
    
    cellCount=size(cellEvents,1);
    
    %%% creating data sets
    switch SET0Name
        case 'H'
            TrialType0=zeros(1,max(HitTrialNumber{SpeedMode,mode}));
            DataIndex0=HitDataset{SpeedMode,mode};
            DataNumber0=HitTrialNumber{SpeedMode,mode};
        case 'M'
            TrialType0=zeros(1,max(MissTrialNumber{SpeedMode,mode,ActiveMode}));
            DataIndex0=MissDataset{SpeedMode,mode,ActiveMode};
            DataNumber0=MissTrialNumber{SpeedMode,mode,ActiveMode};
        case 'C'
            TrialType0=zeros(1,max(CRTrialNumber{SpeedMode,mode,ActiveMode}));
            DataIndex0=CRDataset{SpeedMode,mode,ActiveMode};
            DataNumber0=CRTrialNumber{SpeedMode,mode,ActiveMode};
        case 'F'
            TrialType0=zeros(1,max(FATrialNumber{SpeedMode,mode}));
            DataIndex0=FADataset{SpeedMode,mode};
            DataNumber0=FATrialNumber{SpeedMode,mode};
            
            
            
    end
    
    switch SET1Name
        case 'H'
            TrialType1=zeros(1,max(HitTrialNumber{SpeedMode,mode}));
            DataIndex1=HitDataset{SpeedMode,mode};
            DataNumber1=HitTrialNumber{SpeedMode,mode};
        case 'M'
            TrialType1=zeros(1,max(MissTrialNumber{SpeedMode,mode,ActiveMode}));
            DataIndex1=MissDataset{SpeedMode,mode,ActiveMode};
            DataNumber1=MissTrialNumber{SpeedMode,mode,ActiveMode};
        case 'C'
            TrialType1=zeros(1,max(CRTrialNumber{SpeedMode,mode,ActiveMode}));
            DataIndex1=CRDataset{SpeedMode,mode,ActiveMode};
            DataNumber1=CRTrialNumber{SpeedMode,mode,ActiveMode};
        case 'F'
            TrialType1=zeros(1,max(FATrialNumber{SpeedMode,mode}));
            DataIndex1=FADataset{SpeedMode,mode};
            DataNumber1=FATrialNumber{SpeedMode,mode};
            
            
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Speed Balance
    NBin=50;
    ActiveAnimal=ones(length(Lick),1);
    for i=1:(length(Lick)-75*ActiveTrialNumber)
        if max(Lick(i:(i+75*ActiveTrialNumber)))==0
            ActiveAnimal(i:(i+75*ActiveTrialNumber))=0;
        end
    end
    
    Speed=double(round(Speed*10))/10;
    [SQ,BQ]=SpeedQuantEqual(Speed,(Speed>0 & ActiveAnimal==1),NBin);
    
    Buf_DataNumber0=DataNumber0;
    Buf_DataNumber1=DataNumber1;
    
    
    shuff0=zeros(max(DataNumber0),NRepeate);
    shuff1=zeros(max(DataNumber1),NRepeate);
    
    HitFR=zeros(length(TimeVec),NRepeate);
    CRFR=zeros(length(TimeVec),NRepeate);
    MFR=zeros(length(TimeVec),NRepeate);
    
    for  nr=1:NRepeate;
        Rpeate=nr
        
        DataNumber0=Buf_DataNumber0;
        DataNumber1=Buf_DataNumber1;
        
        
        Sorder0=unique(DataNumber0);
        if Sorder0(1)==0 Sorder0=Sorder0(2:end);end
        Sorder0=Sorder0(randperm(length(Sorder0)));
        Sorder1=unique(DataNumber1);
        if Sorder1(1)==0 Sorder1=Sorder1(2:end);end
        Sorder1=Sorder1(randperm(length(Sorder1)));
        
        UMask0=zeros(size(DataIndex0));
        UMask1=zeros(size(DataIndex1));
        BalSign0=0;
        BalSign1=0;
        NotFinished=1;
        
        if SpeedBalance==1
            while (NotFinished==1)
                NotFinished=0;
                for k0=Sorder0
                    if (TrialType0(k0) * BalSign0) <=0 && max(UMask0(DataNumber0==k0))==0
                        matchable=0;
                        Sbin0=max(SQ(DataIndex0(DataNumber0==k0)));
                        for k1=Sorder1
                            if max(UMask1(DataNumber1==k1))==0
                                Sbin1=max(SQ(DataIndex1(DataNumber1==k1)));
                                if Sbin0==Sbin1
                                    matchable=1;
                                    if (TrialType1(k1) * BalSign1) <=0
                                        UMask0(DataNumber0==k0)=1;
                                        UMask1(DataNumber1==k1)=1;
                                        BalSign0=BalSign0+TrialType0(k0);
                                        BalSign1=BalSign1+TrialType1(k1);
                                        NotFinished=1;
                                        break;
                                    end
                                end
                            end
                        end
                        
                        if matchable==0
                            UMask0(DataNumber0==k0)=2;
                        end
                    end
                end
                
            end
            DataNumber0(UMask0~=1)=0;
            DataNumber1(UMask1~=1)=0;
        end
        
        
        
        Sorder0=unique(DataNumber0);
        if Sorder0(1)==0 Sorder0=Sorder0(2:end);end
        Sorder0=Sorder0(randperm(length(Sorder0)));
        Sorder1=unique(DataNumber1);
        if Sorder1(1)==0 Sorder1=Sorder1(2:end);end
        Sorder1=Sorder1(randperm(length(Sorder1)));
        
        
        
        for tind=TimeVec
            
            
            c0=1;
            c1=1;
            dd0=zeros(size(cellEvents,1),length(DataIndex0));
            for si=Sorder0
                trialLength=sum(DataNumber0==si);
                if(trialLength>=tind)
                    dd0(:,c0)=cellEvents(:,min(DataIndex0(DataNumber0==si))+tind-1);
                    c0=c0+1;
                end
            end
            c0=c0-1;
            dd0=dd0(:,1:c0);
            
            dd1=zeros(size(cellEvents,1),length(DataIndex1));
            for si=Sorder1
                trialLength=sum(DataNumber1==si);
                if(trialLength>=tind)
                    dd1(:,c1)=cellEvents(:,min(DataIndex1(DataNumber1==si))+tind-1);
                    c1=c1+1;
                end
            end
            c1=c1-1;
            dd1=dd1(:,1:c1);
            
            
            
            
            c=min(size(dd0,2),size(dd1,2));
            
            X0=dd0(:,1:c);
            X1=dd1(:,1:c);
            
            
            %                 HitFR(tind,nr)=mean(mean(dd1));
            %                 CRFR(tind,nr)=mean(mean(dd0));
            %                 MFR(tind,nr)=mean(mean([dd0,dd1]));
            
            % [coeff,~,~] = pca([(X0-repmat(mean(X0,2),[1,c])),(X1-repmat(mean(X1,2),[1,c]))]');
            %
            % M0=mean(X0'*coeff(:,1:NPC),1);
            % M1=mean(X1'*coeff(:,1:NPC),1);
            %
            % S0{Mouse-60,tind,nr}=cov(X0'*coeff(:,1:NPC));
            % S1{Mouse-60,tind,nr}=cov(X1'*coeff(:,1:NPC));
            %
            % S=0.5*(cov(X0'*coeff(:,1:NPC))+cov(X1'*coeff(:,1:NPC)));
            %
            % [V,D]=eig(S);
            %
            %
            % dM=(M1-M0)/norm(M1-M0);
            % Res=dM*V';
            %
            %
            % Lamb{Mouse-60,tind,nr}=diag(D);
            % dMe{Mouse-60,tind,nr}=Res;
            
            
            %%%%%%%%%
            hc=floor(c/2);
            [coeff,~] = plsregress([X0(:,1:hc),X1(:,1:hc)]',[zeros(hc,1);ones(hc,1)],NPC);
            %[coeff,~,~] = pca([(X0(:,1:hc)-repmat(mean(X0(:,1:hc),2),[1,hc])),(X1(:,1:hc)-repmat(mean(X1(:,1:hc),2),[1,hc]))]');
            M0=mean(X0(:,hc+1:end)'*coeff(:,1:NPC),1);
            M1=mean(X1(:,hc+1:end)'*coeff(:,1:NPC),1);
            
            
            S0{Mouse-60,tind,nr}=cov(X0(:,hc+1:end)'*coeff(:,1:NPC));
            S1{Mouse-60,tind,nr}=cov(X1(:,hc+1:end)'*coeff(:,1:NPC));
            
            S=0.5*(cov(X0(:,hc+1:end)'*coeff(:,1:NPC))+cov(X1(:,hc+1:end)'*coeff(:,1:NPC)));
            
            [V,D]=eig(S);
            
            dM=(M1-M0)/norm(M1-M0);
            Res=dM*V';
            
            
            Lamb{Mouse-60,tind,nr}=diag(D);
            dMe{Mouse-60,tind,nr}=Res;
            
            
        end
    end
    %    figure();hold on
    %    xlabel('Time (s)');
    %    ylabel('Mean Firing Rate');
    %    shadedErrorBar(-0.4:0.1:2,squeeze(mean(HitFR,2))',sqrt(var(HitFR')),'b',1);
    %    shadedErrorBar(-0.4:0.1:2,squeeze(mean(CRFR,2))',sqrt(var(CRFR')),'r',1);
    %    shadedErrorBar(-0.4:0.1:2,squeeze(mean(MFR,2))',sqrt(var(MFR')),{'Color',[0.6,0.6,0.6]},1);
    
    
    
end

save('E:\Reports2\2018_08_07_SignalNoiseAllignment\Instant_IB_Speed2_Mode1_PLS2_Val_10PC','Lamb','dMe','S0','S1','-v7.3');
%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Colors=linspecer(10);
 Colorst=linspecer(10);
 MList={'o','*','s','^','+','x','.','.','.','.'};
 Mouse=[1,3,4,5,6,7];
 dMeMtx=zeros(6,25,100,10);
 LambMtx=zeros(6,25,100,10);
 TimeV=14:25;

 for M=1:6
     for tind=1:25
         for nr=1:100
            [LambMtx(Mouse(M),tind,nr,:),ord]=sort(Lamb{Mouse(M),tind,nr},'descend');
            dMeMtx(Mouse(M),tind,nr,:)=abs(dMe{Mouse(M),tind,nr}(ord));
         end
     end

     figure();hold on
     for pc=1:10

         plot(squeeze(mean(mean(dMeMtx(Mouse(M),TimeV,:,pc),1),2)),squeeze(mean(mean(LambMtx(Mouse(M),TimeV,:,pc),1),2)),'Color',Colorst(pc,:),'Marker',MList{pc},'LineStyle','none');

     end
 end

 for M=1:6
figure();hold on
     for pc=1:10

         shadedErrorBar(-0.4:0.1:2,squeeze(mean(dMeMtx(Mouse(M),:,:,pc),3)),sqrt(var(squeeze(dMeMtx(Mouse(M),:,:,pc))')),{'Color',Colorst(pc,:)},1);

     end
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% KL Divergence of two noise distributions
nM=1;
KLDiv=zeros(25,6);
CovCorr=zeros(25,6);
for M=[1,3:7]
    
    for tind=1:25
        cov0=0;
        cov1=0;
        for nr=1:100
            cov0=cov0+(S0{M,tind,nr}/100);
            cov1=cov1+(S1{M,tind,nr}/100);
            
        end
        NPC=size(cov0,1);
        
        KLDiv(tind,nM)=0.5*( log(det(cov0)/det(cov1)) +sum(diag(inv(cov1)*cov0)) - NPC);
        CovCorr(tind,nM)=similarity([cov0(eye(NPC)==0),cov1(eye(NPC)==0)]);
    end
    nM=nM+1;
end
figure();
plot(mean(CovCorr,2));

%%%%%%%%%%%%%%%%%%%%%%%% Correlation of Covariance Matrices elements

nM=1;

cov0=zeros(NPC,NPC,25,100);
cov1=zeros(NPC,NPC,25,100);
StimCorr=zeros(3*NPC*(NPC-1),25,2);
TimeCorr=zeros(25,1);
RVal=zeros(25,1);

figure();
hold on

for M=[1,3:7]
    
    
    for tind=1:25
   
        for nr=1:100
            cov0(:,:,tind,nr)=S0{M,tind,nr}./sqrt(diag(S0{M,tind,nr})*diag(S0{M,tind,nr})');
            cov1(:,:,tind,nr)=S1{M,tind,nr}./sqrt(diag(S1{M,tind,nr})*diag(S1{M,tind,nr})');
        end
        
        ns=1;
        for i=1:NPC
            for j=(i+1):NPC
                
                    if tind==20
                        ex=sqrt(var(cov0(i,j,tind,:)));
                        ey=sqrt(var(cov1(i,j,tind,:)));
                       plot(mean(cov0(i,j,tind,:),4),mean(cov1(i,j,tind,:),4),'bo');
                        plot([mean(cov0(i,j,tind,:),4)-ex,mean(cov0(i,j,tind,:),4)+ex],[mean(cov1(i,j,tind,:),4),mean(cov1(i,j,tind,:),4)],'k');
                       plot([mean(cov0(i,j,tind,:),4),mean(cov0(i,j,tind,:),4)],[mean(cov1(i,j,tind,:),4)-ey,mean(cov1(i,j,tind,:),4)+ey],'k');
                       

                    end
                        StimCorr(0.5*NPC*(NPC-1)*(nM-1)+ns,tind,1)=mean(cov0(i,j,tind,:),4);
                        StimCorr(0.5*NPC*(NPC-1)*(nM-1)+ns,tind,2)=mean(cov1(i,j,tind,:),4);
                        ns=ns+1;
                
            end
        end
        TimeCorr(tind)=similarity(squeeze(StimCorr(:,tind,:)));
    end
    nM=nM+1;
end
figure();
plot(TimeCorr);
















