MouseName={'L347','L362bis','L364','L365','L367','L368'};

for Mouse=[1]
    reg=2;
    loadpath=strcat('C:\Users\Sadegh\Documents\VLMReborn\',MouseName{Mouse},'\Data\ConcatDays');
    SET0Name='C';
    SET1Name='H';
    CellT='IB';
    mode=1;
    TSin=200;
    SpeedMode=1;
    cellSelect=0;
    
    
    Balanced=0;
    LowDimVec=1:50;
    NRepeate=10;
    ActiveTrialNumber=20;
    SpeedBalance=1;
    ActiveMode=1;
    
    warning('off','all');
    load(strcat(loadpath,'/All_Sessions.mat'));
    load(strcat(loadpath,'/cellData.mat'));
    load(strcat(loadpath,'/cellData_ZS.mat'));
    load(strcat(loadpath,'/SessionLength.mat'));
    load(strcat(loadpath,'/Datasets/Datasets--5to0.mat'));
    Speed=Speed-min(Speed);
    progress='Data is loaded!'
    
    if mode==1
        TimeVector=1:25;
    elseif mode==2
        TimeVector=1:5;
    elseif mode==3
        TimeVector=1:30;
        
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
    
    

    rVal0={};
    rVal1={};
    rTrain0={};
    rTrain1={};
    CanCorr=zeros(8,8,length(TimeVector),NRepeate);
    CanCorr0_train=zeros(8,8,length(TimeVector),NRepeate);
    CanCorr1_train=zeros(8,8,length(TimeVector),NRepeate);

    
    CanCorr0_val=zeros(8,8,length(TimeVector),NRepeate);
    CanCorr1_val=zeros(8,8,length(TimeVector),NRepeate);


    
    for nr=1:NRepeate
        
        Repeate=nr
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


        for indtd=TimeVector
        
        c0=1;
        c1=1;
        dd0=zeros(size(cellEvents,1),length(DataIndex0));
        for si=Sorder0
            trialLength=sum(DataNumber0==si);
            if(trialLength>0)
                        dd0(:,c0)=cellEvents(:,min(DataIndex0(DataNumber0==si))+indtd-1);
                        c0=c0+1;
            end
        end
        c0=c0-1;
        
        
        dd1=zeros(size(cellEvents,1),length(DataIndex1));
        for si=Sorder1
            trialLength=sum(DataNumber1==si);
            if(trialLength>0)
                    dd1(:,c1)=cellEvents(:,min(DataIndex1(DataNumber1==si))+indtd-1);
                    c1=c1+1;
            end
        end
        c1=c1-1;
        
        
        
        
        
        
        c=floor(min(c0,c1)/2)
        dd0v=dd0(:,(c+1):(2*c));
        dd1v=dd1(:,(c+1):(2*c));
        dd0=dd0(:,1:c);
        dd1=dd1(:,1:c);
        
        Sh1=randperm(c);
        Sh2=randperm(c);
        cxVec={};
        cyVec={};
        

            for af=1:8
                af
                CF=CortexArea==af;
                if (sum(CF)<5) continue;end 
                    for am=af+1:8
                        am
                        CM=CortexArea==am;
                        if (sum(CM)<5) continue;end 

                            
                            [cxVec{af,am,nr},cyVec{af,am,nr},rVal0{af,am,nr},rTrain0{af,am,nr}]=SparseCCA_Par(dd0(CF,:)',dd0(CM,:)',dd0v(CF,:)',dd0v(CM,:)',reg);
                            
                            CanCorr0_train(af,am,indtd,nr)=max(max(rTrain0{af,am,nr}));
                            CanCorr0_val(af,am,indtd,nr)=max(max(rVal0{af,am,nr}));
                            
                            [cxVec{af,am,nr},cyVec{af,am,nr},rVal1{af,am,nr},rTrain1{af,am,nr}]=SparseCCA_Par(dd1(CF,:)',dd1(CM,:)',dd1v(CF,:)',dd1v(CM,:)',reg);
                            
                            CanCorr1_train(af,am,indtd,nr)=max(max(rTrain1{af,am,nr}));
                            CanCorr1_val(af,am,indtd,nr)=max(max(rVal1{af,am,nr}));
                            
                            
                            
                            
                        end
                    end
                end
                
            end

    
    
    
    save(strcat('E:\Reports2\2018_03_14_CannonicalCorrelations\HitCR\SparseInstant\ParOpt\',MouseName{Mouse}),'CanCorr0_train','CanCorr1_train','CanCorr0_Sh_train','CanCorr1_Sh_train'...
        ,'CanCorr0_val','CanCorr1_val','CanCorr0_Sh_val','CanCorr1_Sh_val','TimeVec','CellT','c','cellSelect','rVal0','rVal1','rTrain0','rTrain1','cxVec','cyVec','-v7.3');
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MouseVec=[61,63:67];
% p='E:\Reports2\2018_03_14_CannonicalCorrelations\HitCR\Sparse\ParOpt\Session'
% 
% mean0=zeros(8,8,6);
% mean1=zeros(8,8,6);
% for M=1:6
%     load(strcat(p,num2str(MouseVec(M)),'.mat'));
% sumrVal0=zeros(8,8,20,20);
% sumrTrain0=zeros(8,8,20,20);
% nrVal0=zeros(8,8);
% sumrVal1=zeros(8,8,20,20);
% sumrTrain1=zeros(8,8,20,20);
% nrVal1=zeros(8,8);
% 
% for a1=1:size(rVal0,1)
%     for a2=a1+1:size(rVal0,2)
%         for r=1:size(rVal0,3)
%             if max(max(rVal0{a1,a2,r}))>0
%                 sumrVal0(a1,a2,:,:)=squeeze(sumrVal0(a1,a2,:,:))+rVal0{a1,a2,r};
%                 sumrTrain0(a1,a2,:,:)=squeeze(sumrTrain0(a1,a2,:,:))+rTrain0{a1,a2,r};
%                 nrVal0(a1,a2)=nrVal0(a1,a2)+1;
%             end
%             if max(max(rVal1{a1,a2,r}))>0
%                 sumrVal1(a1,a2,:,:)=squeeze(sumrVal1(a1,a2,:,:))+rVal1{a1,a2,r};
%                 sumrTrain1(a1,a2,:,:)=squeeze(sumrTrain1(a1,a2,:,:))+rTrain1{a1,a2,r};
%                 nrVal1(a1,a2)=nrVal1(a1,a2)+1;
%             end
%         end
%         sumrVal0(a1,a2,:,:)=squeeze(sumrVal0(a1,a2,:,:))/ nrVal0(a1,a2);
%         sumrTrain0(a1,a2,:,:)=squeeze(sumrTrain0(a1,a2,:,:))/ nrVal0(a1,a2);
%         mean0(a1,a2,M)=max(max(sumrVal0(a1,a2,:,:)));
%         mean0(a2,a1,M)=max(max(sumrVal0(a1,a2,:,:)));
% 
%         sumrVal1(a1,a2,:,:)=squeeze(sumrVal1(a1,a2,:,:))/ nrVal1(a1,a2);
%         sumrTrain1(a1,a2,:,:)=squeeze(sumrTrain1(a1,a2,:,:))/ nrVal1(a1,a2);
%         mean1(a1,a2,M)=max(max(sumrVal1(a1,a2,:,:)));
%         mean1(a2,a1,M)=max(max(sumrVal1(a1,a2,:,:)));
%     end
% end
% 
% end
% 
% 
% CCAMtx0=zeros(8,8);
% CCAMtx1=zeros(8,8);
% for a1=1:8
%     for a2=1:8
%         x=squeeze(mean0(a1,a2,:));
%         mask= x~=0 & isnan(x)==0;
%         CCAMtx0(a1,a2)=mean(x(mask));
% 
%         x=squeeze(mean1(a1,a2,:));
%         mask= x~=0 & isnan(x)==0;
%         CCAMtx1(a1,a2)=mean(x(mask));
%     end
% end
% 









