%function SparseCanonicalCorrelation_2018(loadpath,savepath,SET0Name,SET1Name)
for Mouse=61;
    AddressSetup;
    CVpath='E:\Reports2\2018_03_14_CannonicalCorrelations\ROI_Bin\Datasets--5to0_par\HitCR';
    savepath='E:\Reports2\2018_03_14_CannonicalCorrelations\ROI_Bin\Datasets--5to0\HitCR_L2';
    loadpath=LoadPath{Mouse-10};
    GAP=0;
    ccNum=20;
    reg=2;
    Shuffle=0;
    
    SET0Name='C';
    SET1Name='H';
    CellT='IB';
    mode=1;
    TSin=200;
    SpeedMode=1;
    cellSelect=0;
    
    
    Balanced=0;
    LowDimVec=1:50;
    NRepeate=25;
    ActiveTrialNumber=20;
    SpeedBalance=1;
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
    par=load(strcat(CVpath,'/Session',num2str(Mouse),'_mode',num2str(mode)));
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
    NearCortexArea=zeros(cellCount,8);
if GAP==1

Area=Area>0;
        se=strel('disk',15,8);
        Mask=zeros(size(Area));
        for a=1:8
            Mask(:,:,a)=imdilate(Area(:,:,a),se);
            Mask(:,:,a)=Mask(:,:,a)>0;
            for cnum=1:cellCount
                if Mask(cellIJ(cnum,1),cellIJ(cnum,2),a)>0

                        NearCortexArea(cnum,a)=1;

                end

            end
        end

end

    
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
        case 'G'
            TrialType0=[(-1*ones(1,max(HitTrialNumber{SpeedMode,mode}))),ones(1,max(MissTrialNumber{SpeedMode,mode,ActiveMode}))];
            DataIndex0=[HitDataset{SpeedMode,mode},MissDataset{SpeedMode,mode,ActiveMode}];
            DataNumber0=[HitTrialNumber{SpeedMode,mode},(MissTrialNumber{SpeedMode,mode,ActiveMode}+ max(HitTrialNumber{SpeedMode,mode}))];
            
        case 'N'
            TrialType0=[(-1*ones(1,max(FATrialNumber{SpeedMode,mode}))),ones(1,max(CRTrialNumber{SpeedMode,mode,ActiveMode}))];
            DataIndex0=[FADataset{SpeedMode,mode},CRDataset{SpeedMode,mode,ActiveMode}];
            DataNumber0=[FATrialNumber{SpeedMode,mode},(CRTrialNumber{SpeedMode,mode,ActiveMode}+max(FATrialNumber{SpeedMode,mode}))];
        case 'L'
            TrialType0=[(-1*ones(1,max(FATrialNumber{SpeedMode,mode}))),ones(1,max(HitTrialNumber{SpeedMode,mode}))];
            DataIndex0=[FADataset{SpeedMode,mode},HitDataset{SpeedMode,mode}];
            DataNumber0=[FATrialNumber{SpeedMode,mode},(HitTrialNumber{SpeedMode,mode}+max(FATrialNumber{SpeedMode,mode}))];
        case 'NL'
            TrialType0=[(-1*ones(1,max(MissTrialNumber{SpeedMode,mode,ActiveMode}))),ones(1,max(CRTrialNumber{SpeedMode,mode,ActiveMode}))];
            DataIndex0=[MissDataset{SpeedMode,mode,ActiveMode},CRDataset{SpeedMode,mode,ActiveMode}];
            DataNumber0=[MissTrialNumber{SpeedMode,mode,ActiveMode},(CRTrialNumber{SpeedMode,mode,ActiveMode}+max(MissTrialNumber{SpeedMode,mode,ActiveMode}))];
        case 'Cor'
            TrialType0=[(-1*ones(1,max(CRTrialNumber{SpeedMode,mode,ActiveMode}))),ones(1,max(HitTrialNumber{SpeedMode,mode}))];
            DataIndex0=[CRDataset{SpeedMode,mode,ActiveMode},HitDataset{SpeedMode,mode}];
            DataNumber0=[CRTrialNumber{SpeedMode,mode,ActiveMode},(HitTrialNumber{SpeedMode,mode}+max(CRTrialNumber{SpeedMode,mode,ActiveMode}))];
        case 'Err'
            TrialType0=[(-1*ones(1,max(MissTrialNumber{SpeedMode,mode,ActiveMode}))),ones(1,max(FATrialNumber{SpeedMode,mode}))];
            DataIndex0=[MissDataset{SpeedMode,mode,ActiveMode},FADataset{SpeedMode,mode}];
            DataNumber0=[MissTrialNumber{SpeedMode,mode,ActiveMode},(FATrialNumber{SpeedMode,mode}+max(MissTrialNumber{SpeedMode,mode,ActiveMode}))];
            
            
            
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
        case 'G'
            TrialType1=[(-1*ones(1,max(HitTrialNumber{SpeedMode,mode}))),ones(1,max(MissTrialNumber{SpeedMode,mode,ActiveMode}))];
            DataIndex1=[HitDataset{SpeedMode,mode},MissDataset{SpeedMode,mode,ActiveMode}];
            DataNumber1=[HitTrialNumber{SpeedMode,mode},(MissTrialNumber{SpeedMode,mode,ActiveMode}+ max(HitTrialNumber{SpeedMode,mode}))];
            
        case 'N'
            TrialType1=[(-1*ones(1,max(FATrialNumber{SpeedMode,mode}))),ones(1,max(CRTrialNumber{SpeedMode,mode,ActiveMode}))];
            DataIndex1=[FADataset{SpeedMode,mode},CRDataset{SpeedMode,mode,ActiveMode}];
            DataNumber1=[FATrialNumber{SpeedMode,mode},(CRTrialNumber{SpeedMode,mode,ActiveMode}+max(FATrialNumber{SpeedMode,mode}))];
        case 'L'
            TrialType1=[(-1*ones(1,max(FATrialNumber{SpeedMode,mode}))),ones(1,max(HitTrialNumber{SpeedMode,mode}))];
            DataIndex1=[FADataset{SpeedMode,mode},HitDataset{SpeedMode,mode}];
            DataNumber1=[FATrialNumber{SpeedMode,mode},(HitTrialNumber{SpeedMode,mode}+max(FATrialNumber{SpeedMode,mode}))];
        case 'NL'
            TrialType1=[(-1*ones(1,max(MissTrialNumber{SpeedMode,mode,ActiveMode}))),ones(1,max(CRTrialNumber{SpeedMode,mode,ActiveMode}))];
            DataIndex1=[MissDataset{SpeedMode,mode,ActiveMode},CRDataset{SpeedMode,mode,ActiveMode}];
            DataNumber1=[MissTrialNumber{SpeedMode,mode,ActiveMode},(CRTrialNumber{SpeedMode,mode,ActiveMode}+max(MissTrialNumber{SpeedMode,mode,ActiveMode}))];
        case 'Cor'
            TrialType1=[(-1*ones(1,max(CRTrialNumber{SpeedMode,mode,ActiveMode}))),ones(1,max(HitTrialNumber{SpeedMode,mode}))];
            DataIndex1=[CRDataset{SpeedMode,mode,ActiveMode},HitDataset{SpeedMode,mode}];
            DataNumber1=[CRTrialNumber{SpeedMode,mode,ActiveMode},(HitTrialNumber{SpeedMode,mode}+max(CRTrialNumber{SpeedMode,mode,ActiveMode}))];
        case 'Err'
            TrialType1=[(-1*ones(1,max(MissTrialNumber{SpeedMode,mode,ActiveMode}))),ones(1,max(FATrialNumber{SpeedMode,mode}))];
            DataIndex1=[MissDataset{SpeedMode,mode,ActiveMode},FADataset{SpeedMode,mode}];
            DataNumber1=[MissTrialNumber{SpeedMode,mode,ActiveMode},(FATrialNumber{SpeedMode,mode}+max(MissTrialNumber{SpeedMode,mode,ActiveMode}))];
            
            
            
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
    
    TimeVec=0;
    rVal0={};rVal0{8,8,NRepeate}=[];
    rVal1={};rVal1{8,8,NRepeate}=[];
    rTrain0={};rTrain0{8,8,NRepeate}=[];
    rTrain1={};rTrain1{8,8,NRepeate}=[];
    Wx0={}; Wx1={};
    Wy0={}; Wy1={};
    Wx0{8,8,NRepeate}=[]; Wx1{8,8,NRepeate}=[];
    Wy0{8,8,NRepeate}=[]; Wy1{8,8,NRepeate}=[];
    CanCorr=zeros(8,8,NRepeate,ccNum);
    CanCorr0_train=zeros(8,8,NRepeate,ccNum);
    CanCorr1_train=zeros(8,8,NRepeate,ccNum);
    
    
    CanCorr0_val=zeros(8,8,NRepeate,ccNum);
    CanCorr1_val=zeros(8,8,NRepeate,ccNum);
    
    
    
    for nr=1:NRepeate
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
        
        
        offset=0;
        if mode==1
            offset=10;
        end
        
        c0=1;
        c1=1;
        dd0=zeros(size(cellEvents,1),length(DataIndex0));
        for si=Sorder0
            trialLength=sum(DataNumber0==si);
            if(trialLength>0)
                trialLength=trialLength-offset;
                trialIndex=DataIndex0(DataNumber0==si);
                
                dd0(:,c0:c0+trialLength-1)=cellEvents(:,trialIndex(offset+1:end));
                c0=c0+trialLength;
            end
            
        end
        c0=c0-1;
        dd0=dd0(:,1:c0);
        
        dd1=zeros(size(cellEvents,1),length(DataIndex1));
        for si=Sorder1
            trialLength=sum(DataNumber1==si);
            if(trialLength>0)
                trialLength=trialLength-offset;
                trialIndex=DataIndex1(DataNumber1==si);
                
                dd1(:,c1:c1+trialLength-1)=cellEvents(:,trialIndex(offset+1:end));
                c1=c1+trialLength;
            end
        end
        c1=c1-1;
        dd1=dd1(:,1:c1);
        
        
        
        
        
        c=floor(min(c0,c1)/2)
        dd0v=dd0(:,(c+1):(2*c));
        dd1v=dd1(:,(c+1):(2*c));
        dd0=dd0(:,1:c);
        dd1=dd1(:,1:c);
        
        if Shuffle==1
            for cnum=1:size(dd0,1)
                dd0v(cnum,:)=dd0v(cnum,randperm(c));
                dd1v(cnum,:)=dd1v(cnum,randperm(c));
                
                dd0(cnum,:)=dd0(cnum,randperm(c));
                dd1(cnum,:)=dd1(cnum,randperm(c));
                
            end
        end
        
        
        
        if reg==2
            cxpar=1;
            cypar=1;
        end
        
        
 for af=1:8
            af

            if(sum(CortexArea==af)<5) continue;end
            parfor am=af+1:8
                if(sum(CortexArea==am)>5)

                    [Wx0{af,am,nr}, Wy0{af,am,nr},rTrain0{af,am,nr}]=SparseCCA(dd0(CortexArea==af & NearCortexArea(:,am)==0,:)',dd0(CortexArea==am & NearCortexArea(:,af)==0 ,:)',1,1,reg,ccNum);
                end
            end
            parfor am=af+1:8

                if(sum(CortexArea==am)>5)
                    [Wx1{af,am,nr}, Wy1{af,am,nr},rTrain1{af,am,nr}]=SparseCCA(dd1(CortexArea==af & NearCortexArea(:,am)==0 ,:)',dd1(CortexArea==am & NearCortexArea(:,af)==0,:)',1,1,reg,ccNum);
                end
            end
            for am=af+1:8

                if(sum(CortexArea==am)>5)
                    rVal0{af,am,nr}=corrcoef([dd0v(CortexArea==af & NearCortexArea(:,am)==0,:)'*Wx0{af,am,nr},dd0v(CortexArea==am & NearCortexArea(:,af)==0,:)'*Wy0{af,am,nr}]);
                    CanCorr0_train(af,am,nr,:)=rTrain0{af,am,nr};
                    CanCorr0_val(af,am,nr,:)=diag(rVal0{af,am,nr}(ccNum+1:end,1:ccNum));


                    rVal1{af,am,nr}=corrcoef([dd1v(CortexArea==af & NearCortexArea(:,am)==0,:)'*Wx1{af,am,nr},dd1v(CortexArea==am & NearCortexArea(:,af)==0,:)'*Wy1{af,am,nr}]);
                    CanCorr1_train(af,am,nr,:)=rTrain1{af,am,nr};
                    CanCorr1_val(af,am,nr,:)=diag(rVal1{af,am,nr}(ccNum+1:end,1:ccNum));

                end


            end
end

        
    end
    
    
    sh='';
    if Shuffle==1
        sh='Shuffled';
    end
    
    save(strcat(savepath,'/Session',num2str(Mouse),'_mode',num2str(mode),sh),'CanCorr0_train','CanCorr1_train','NRepeate','NearCortexArea'...
        ,'CanCorr0_val','CanCorr1_val','TimeVec','CellT','c','rVal0','rVal1','rTrain0','rTrain1','Wx0','Wy0','Wx1','Wy1','-v7.3');
end

