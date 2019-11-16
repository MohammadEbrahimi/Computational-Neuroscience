AddressSetup;

for  Mouse=[66:67];
    loadpath=LoadPath{Mouse-10};
    
    
    
    SET0Name='C';
    SET1Name='H';
    CellT='IB';
    mode=1;
    TSin=200;
    SpeedMode=2;
    
    Balanced=0;
    NPC=20;
    NRepeate=25;
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
    
    
    PCACoefficients={};
    PCACoefficients{9,NRepeate}=[];
    PCAVariances=zeros(9,NPC,NRepeate);
    
    
    
    
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
        
   
            c0=1;
            c1=1;
            dd0=zeros(size(cellEvents,1),length(DataIndex0));
            for si=Sorder0
                trialLength=sum(DataNumber0==si)-10;
                if(trialLength>=10)
                    dd0(:,c0:c0+trialLength-1)=cellEvents(:,min(DataIndex0(DataNumber0==si)):min(DataIndex0(DataNumber0==si))+trialLength-1);
                    c0=c0+trialLength;
                end
            end
            c0=c0-1;
            dd0=dd0(:,1:c0);
            
            dd1=zeros(size(cellEvents,1),length(DataIndex1));
            for si=Sorder1
                trialLength=sum(DataNumber1==si)-10;
                if(trialLength>=10)
                    dd1(:,c1:c1+trialLength-1)=cellEvents(:,min(DataIndex1(DataNumber1==si)):min(DataIndex1(DataNumber1==si))+trialLength-1);
                    c1=c1+1;
                end
            end
            c1=c1-1;
            dd1=dd1(:,1:c1);
            
            
            
            
            c=min(size(dd0,2),size(dd1,2));
            
            X0=dd0(:,1:c);
            X1=dd1(:,1:c);  
            
            %%%%%%%%%
            for area=1:9
                if area<9
                    Cells=(CortexArea==area);
                else
                    Cells=(CortexArea<area);
                end
                
                if (sum(Cells)>20)
                    [coeff,~,latent] = pca([(X0(Cells,:)-repmat(mean(X0(Cells,:),2),[1,c])),(X1(Cells,:)-repmat(mean(X1(Cells,:),2),[1,c]))]');
                
                    PCAVariances(area,:,nr)=latent(1:NPC);
                    PCACoefficients{area,nr}=coeff(:,1:NPC);
                end
                
                
                
            end
            
            
            
       
    end
   
  save(strcat('E:\Reports2\2018_11_17_PCANoiseAnalysis\ConsensusPC5to0\Mouse',num2str(Mouse-60)),...
      'PCAVariances','PCACoefficients','NPC','-v7.3');
  
    
end
