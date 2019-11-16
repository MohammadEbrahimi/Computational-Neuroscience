MouseName={'L347','L362bis','L364','L365','L367','L368'};

for Mouse=[2]

loadpath=strcat('C:\Users\Sadegh\Documents\VLMReborn\',MouseName{Mouse},'\Data\ConcatDays');
SET0Name='M';
SET1Name='H';
CellT='RB';
mode=1;
TSin=200;
SpeedMode=1;



Balanced=0;
LowDimVec=1:50;
NRepeate=100;
ActiveTrialNumber=20;
SpeedBalance=1;
ActiveMode=3;

warning('off','all');
load(strcat(loadpath,'/All_Sessions.mat'));
load(strcat(loadpath,'/cellData.mat'));
load(strcat(loadpath,'/cellData_ZS.mat'));
load(strcat(loadpath,'/SessionLength.mat'));
load(strcat(loadpath,'/Datasets/Datasets--5to0.mat'));
Speed=Speed-min(Speed);
progress='Data is loaded!'


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


shuff0=zeros(max(DataNumber0),NRepeate);
shuff1=zeros(max(DataNumber1),NRepeate);

 TimeVec=-5:5;
cellVector=20:50:1000;
    CanCorr=zeros(8,8,length(TimeVec),NRepeate,length(cellVector));
    CanCorr0_train=zeros(8,8,length(TimeVec),NRepeate,length(cellVector));
    CanCorr1_train=zeros(8,8,length(TimeVec),NRepeate,length(cellVector));
    CanCorr0_Sh_train=zeros(8,8,length(TimeVec),NRepeate,length(cellVector));
    CanCorr1_Sh_train=zeros(8,8,length(TimeVec),NRepeate,length(cellVector));
    
    CanCorr0_val=zeros(8,8,length(TimeVec),NRepeate,length(cellVector));
    CanCorr1_val=zeros(8,8,length(TimeVec),NRepeate,length(cellVector));
    CanCorr0_Sh_val=zeros(8,8,length(TimeVec),NRepeate,length(cellVector));
    CanCorr1_Sh_val=zeros(8,8,length(TimeVec),NRepeate,length(cellVector));
for CellNumIdx=1:length(cellVector)
    CellNum=cellVector(CellNumIdx)
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
    
    
    
    c0=1;
    c1=1;
    dd0=zeros(size(cellEvents,1),length(DataIndex0));
    for si=Sorder0
        trialLength=sum(DataNumber0==si);
        if(trialLength>0)
            dd0(:,c0:c0+trialLength-1)=cellEvents(:,DataIndex0(DataNumber0==si));
            c0=c0+trialLength;
        end
    end
    c0=c0-1;
    
    
    dd1=zeros(size(cellEvents,1),length(DataIndex1));
    for si=Sorder1
        trialLength=sum(DataNumber1==si);
        if(trialLength>0)
            dd1(:,c1:c1+trialLength-1)=cellEvents(:,DataIndex1(DataNumber1==si));
            c1=c1+trialLength;
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
    
    for tind=1:length(TimeVec)
        td=TimeVec(tind)
        for af=1:8
            CF=CortexArea==af;
            if sum(CF)>CellNum
                idx=find(CF);
                idx=idx(randperm(length(idx)));
                CF(idx(CellNum:end))=false;
                for am=af+1:8
                    
                    CM=CortexArea==am;
                    if sum(CM)>CellNum
                        idx=find(CM);
                        idx=idx(randperm(length(idx)));
                        CM(idx(CellNum:end))=false;
                        
                        sm=2*max(TimeVec)+1-td;
                        em=c-2*max(TimeVec)-td;
                        sf=2*max(TimeVec)+1;
                        ef=c-2*max(TimeVec);
                        
                        
                        
                        
                        
                        
                        [A,B,r] = canoncorr(dd0(CF,sf:ef)',dd0(CM,sm:em)');
                        rv=corrcoef([dd0v(CF,sf:ef)'*A(:,1),dd0v(CM,sm:em)'*B(:,1)]);
                        CanCorr0_train(af,am,tind,nr,CellNumIdx)=r(1);
                        CanCorr0_val(af,am,tind,nr,CellNumIdx)=rv(1,2);
                        [A,B,r] = canoncorr(dd1(CF,sf:ef)',dd1(CM,sm:em)');
                        rv=corrcoef([dd1v(CF,sf:ef)'*A(:,1),dd1v(CM,sm:em)'*B(:,1)]);
                        CanCorr1_train(af,am,tind,nr,CellNumIdx)=r(1);
                        CanCorr1_val(af,am,tind,nr,CellNumIdx)=rv(1,2);
%                         
%                         [~,~,r] = canoncorr(dd0(CF,Sh1(sf:ef))',dd0(CM,Sh2(sm:em))');
%                         rv=corrcoef([dd0v(CF,Sh1(sf:ef))'*A(:,1),dd0v(CM,Sh2(sm:em))'*B(:,1)]);
%                         CanCorr0_Sh_train(af,am,tind,nr,CellNumIdx)=r(1);
%                         CanCorr0_Sh_val(af,am,tind,nr,CellNumIdx)=rv(1,2);
%                         [~,~,r] = canoncorr(dd1(CF,Sh1(sf:ef))',dd1(CM,Sh2(sm:em))');
%                         rv=corrcoef([dd1v(CF,Sh1(sf:ef))'*A(:,1),dd1v(CM,Sh2(sm:em))'*B(:,1)]);
%                         CanCorr1_Sh_train(af,am,tind,nr,CellNumIdx)=r(1);
%                         CanCorr1_Sh_val(af,am,tind,nr,CellNumIdx)=rv(1,2);
%                         
                        
                        
                        
                    end
                end
            end
            
        end
    end
    
    
    
end %%%%% end repeate
end


save(strcat('E:\Reports2\2018_03_14_CannonicalCorrelations\',MouseName{Mouse}),'CanCorr0_train','CanCorr1_train','CanCorr0_Sh_train','CanCorr1_Sh_train'...
    ,'CanCorr0_val','CanCorr1_val','CanCorr0_Sh_val','CanCorr1_Sh_val','TimeVec','CellT','cellVector','c','-v7.3');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  p='E:\Reports2\2018_03_14_CannonicalCorrelations\HitMiss';
% mean0=zeros(8,8,6);
% meanTime0=zeros(8,8,6);
% n0=zeros(8,8,6);
% mean1=zeros(8,8,6);
% meanTime1=zeros(8,8,6);
% n1=zeros(8,8,6);
% T=1;
% C=1:20;
% for M=[1:6]
%     load(strcat(p,'\',MouseName{M},'.mat'));
%     for a1=1:8
%         if max(max(mean(CanCorr0_val(a1,:,T,:,C),4),[],5),[],2)==0 && max(max(mean(CanCorr0_val(:,a1,T,:,C),4),[],5),[],1)==0
%             mean0(a1,:,M)=0.03;
%             mean1(a1,:,M)=0.03;
%             mean0(:,a1,M)=0.03;
%             mean1(:,a1,M)=0.03;
%             continue;
%         end
%         for a2=a1+1:8
%            % if max(mean(CanCorr0_val(a1,a2,T,:,C),4),[],5)>=max(max(CanCorr0_Sh_val(a1,a2,T,:,C),[],5),[],4)
%                 n0(a1,a2,M)=n0(a1,a2,M)+1;
%                 mean0(a1,a2,M)=max(mean(CanCorr0_val(a1,a2,T,:,C),4),[],5);
%             %end
%             %if max(mean(CanCorr1_val(a1,a2,T,:,C),4),[],5)>=max(max(CanCorr1_Sh_val(a1,a2,T,:,C),[],5),[],4)
%                 n1(a1,a2,M)=n1(a1,a2,M)+1;
%                 mean1(a1,a2,M)=max(mean(CanCorr1_val(a1,a2,T,:,C),4),[],5);
%             %end
% 
%         end
%     end
%     mean0(:,:,M)=mean0(:,:,M)+mean0(:,:,M)';
%     mean1(:,:,M)=mean1(:,:,M)+mean1(:,:,M)';
%     figure();
%     
%     imagesc([mean0(:,:,M),0.05*ones(8,2),mean1(:,:,M)]);
% end
% 
% Average0=zeros(8,8);
% Average1=zeros(8,8);
% Time0=zeros(8,8);
% Time1=zeros(8,8);
% sigC=zeros(8,8);
% sigT=zeros(8,8);
% 
% for a1=1:8
%     for a2=1:8
%         mask0=mean0(a1,a2,:)~=0.06;
%         Average0(a1,a2)=mean(mean0(a1,a2,mask0));
%         mask1=mean1(a1,a2,:)~=0.06;
%         Average1(a1,a2)=mean(mean1(a1,a2,mask1));
% %         Time0(a1,a2)=mean(meanTime0(a1,a2,mask));
% %         Time1(a1,a2)=mean(meanTime1(a1,a2,mask));
%         mask=mask0 | mask1;
%         sigC(a1,a2)=max(sum(mean0(a1,a2,mask)> mean1(a1,a2,mask)) ,...
%             sum(mean0(a1,a2,mask)< mean1(a1,a2,mask)))/sum(mask);
%        % sigT(a1,a2)=sum((meanTime0(a1,a2,mask)./n0(a1,a2,mask))>(meanTime1(a1,a2,mask)./n1(a1,a2,mask)))==sum(mask) ||...
%         %    sum((meanTime0(a1,a2,mask)./n0(a1,a2,mask))<(meanTime1(a1,a2,mask)./n1(a1,a2,mask)))==sum(mask);
%     end
% end

% figure();
% imagesc([(meanTime0./n0),(meanTime1./n1)])
% figure();
% imagesc([(mean0./n0),(mean1./n1)])
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%train /val
AreaColors=double([0,0,128;0,191,255;100,149,237;138,43,226;34,139,34;255,69,0;112,128,144;244,164,96])/256;
figure();hold on;
title('During Miss Training- Validation--')
xlabel('Number of Neurons');
ylabel('Canonical Correlation with V1'); 
set(0,'defaultlinelinewidth',1.5);
set(0,'DefaultAxesFontSize',16);
set(0,'defaultAxesFontName','Calibri');
for a=[3,4,6,8]
    IT=squeeze(mean(CanCorr1_train(1,a,1,:,:),4));
    IT=IT(IT>0);
    IT=[0;IT];
    IV=squeeze(mean(CanCorr1_val(1,a,1,:,:),4));
    IV=IV(IV>0);
    IV=[0;IV];
    
    %plot([0,cellVector(1:length(IT)-1)],IT,'Color',AreaColors(a,:));
    plot([0,cellVector(1:length(IV)-1)],IV,'Color',AreaColors(a,:),'LineStyle','--');
    
end
legend('MV','PPC','S','RSC'); 
legend('LV','MV','PPC','A','S','M','RSC'); 
    
    
figure();hold on;
for a=1:8
a1=min(4,a);
a2=max(4,a);
%plot(-0.5:0.1:0.5,squeeze(max(mean(CanCorr1_val(a1,a2,:,:,:),4),[],5)),'Color',AreaColors(a,:))
plot(-0.5:0.1:0.5,squeeze(max(max(CanCorr1_Sh_val(a1,a2,:,:,:),[],4),[],5)),'Color',AreaColors(a,:),'LineStyle','--')
end
legend('V1','LV','MV','PPC','A','S','M','RSC');








