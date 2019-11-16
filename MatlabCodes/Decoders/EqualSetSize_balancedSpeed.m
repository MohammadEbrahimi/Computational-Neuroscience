%function [DataSize]=EqualSetSize(loadpath,SET0Name,SET1Name,SET2Name,SET3Name,mode,SpeedMode,ActiveMode,TPin)
ActiveTrialNumber=20;
SET0Name='H';
SET1Name='M';
SET2Name='F';
SET3Name='C';
mode=1;
SpeedMode=1;
ActiveMode=1;
TPin='-0to0';

loadpath='C:\Users\Sadegh\Documents\VLMReborn\L364\Data\allDays\'
load(strcat(loadpath,'\All_Sessions.mat'));
load(strcat(loadpath,'\Datasets',TPin,'.mat'));
%%% creating data sets
%%% creating data sets
switch SET0Name
        case 'H'
                DataIndex0=HitDataset{SpeedMode,mode};
        DataNumber0=HitTrialNumber{SpeedMode,mode};
        case 'M'
                DataIndex0=MissDataset{SpeedMode,mode,ActiveMode};
        DataNumber0=MissTrialNumber{SpeedMode,mode,ActiveMode};
        case 'C'
                DataIndex0=CRDataset{SpeedMode,mode,ActiveMode};
        DataNumber0=CRTrialNumber{SpeedMode,mode,ActiveMode};
        case 'F'
                DataIndex0=FADataset{SpeedMode,mode};
        DataNumber0=FATrialNumber{SpeedMode,mode};
        case 'G'
        MinNumTrial=min( max(MissTrialNumber{SpeedMode,mode,ActiveMode}) , max(HitTrialNumber{SpeedMode,mode}));
                DataIndex0=[HitDataset{SpeedMode,mode}(HitTrialNumber{SpeedMode,mode}<=MinNumTrial),...
            MissDataset{SpeedMode,mode,ActiveMode}(MissTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)];
        DataNumber0=[HitTrialNumber{SpeedMode,mode}(HitTrialNumber{SpeedMode,mode}<=MinNumTrial),...
            (MissTrialNumber{SpeedMode,mode,ActiveMode}(MissTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)+MinNumTrial)];
        case 'N'
                MinNumTrial=min( max(CRTrialNumber{SpeedMode,mode,ActiveMode}) , max(FATrialNumber{SpeedMode,mode}));
                DataIndex0=[FADataset{SpeedMode,mode}(FATrialNumber{SpeedMode,mode}<=MinNumTrial),...
            CRDataset{SpeedMode,mode,ActiveMode}(CRTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)];
        DataNumber0=[FATrialNumber{SpeedMode,mode}(FATrialNumber{SpeedMode,mode}<=MinNumTrial),...
            (CRTrialNumber{SpeedMode,mode,ActiveMode}(CRTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)+MinNumTrial)];
        case 'L'
                MinNumTrial=min( max(HitTrialNumber{SpeedMode,mode}) , max(FATrialNumber{SpeedMode,mode}));
                DataIndex0=[FADataset{SpeedMode,mode}(FATrialNumber{SpeedMode,mode}<=MinNumTrial),...
            HitDataset{SpeedMode,mode}(HitTrialNumber{SpeedMode,mode}<=MinNumTrial)];
        DataNumber0=[FATrialNumber{SpeedMode,mode}(FATrialNumber{SpeedMode,mode}<=MinNumTrial),...
            (HitTrialNumber{SpeedMode,mode}(HitTrialNumber{SpeedMode,mode}<=MinNumTrial)+MinNumTrial)];
        case 'NL'
                MinNumTrial=min( max(CRTrialNumber{SpeedMode,mode,ActiveMode}) , max(MissTrialNumber{SpeedMode,mode,ActiveMode}));
                DataIndex0=[MissDataset{SpeedMode,mode,ActiveMode}(MissTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial),...
            CRDataset{SpeedMode,mode,ActiveMode}(CRTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)];
        DataNumber0=[MissTrialNumber{SpeedMode,mode,ActiveMode}(MissTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial),...
            (CRTrialNumber{SpeedMode,mode,ActiveMode}(CRTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)+MinNumTrial)];
        case 'MAI'
                MinNumTrial=max(MissTrialNumber{SpeedMode,mode,1});
                DataIndex0=[MissDataset{SpeedMode,mode,1},MissDataset{SpeedMode,mode,2}];
        DataNumber0=[MissTrialNumber{SpeedMode,mode,1},(MissTrialNumber{SpeedMode,mode,2}+MinNumTrial)];


end

switch SET1Name
        case 'H'
                DataIndex1=HitDataset{SpeedMode,mode};
        DataNumber1=HitTrialNumber{SpeedMode,mode};
        case 'M'
                DataIndex1=MissDataset{SpeedMode,mode,ActiveMode};
        DataNumber1=MissTrialNumber{SpeedMode,mode,ActiveMode};
        case 'C'
                DataIndex1=CRDataset{SpeedMode,mode,ActiveMode};
        DataNumber1=CRTrialNumber{SpeedMode,mode,ActiveMode};
        case 'F'
                DataIndex1=FADataset{SpeedMode,mode};
        DataNumber1=FATrialNumber{SpeedMode,mode};
        case 'G'
        MinNumTrial=min( max(MissTrialNumber{SpeedMode,mode,ActiveMode}) , max(HitTrialNumber{SpeedMode,mode}));
                DataIndex1=[HitDataset{SpeedMode,mode}(HitTrialNumber{SpeedMode,mode}<=MinNumTrial),...
            MissDataset{SpeedMode,mode,ActiveMode}(MissTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)];
        DataNumber1=[HitTrialNumber{SpeedMode,mode}(HitTrialNumber{SpeedMode,mode}<=MinNumTrial),...
            (MissTrialNumber{SpeedMode,mode,ActiveMode}(MissTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)+MinNumTrial)];
        case 'N'
                MinNumTrial=min( max(CRTrialNumber{SpeedMode,mode,ActiveMode}) , max(FATrialNumber{SpeedMode,mode}));
                DataIndex1=[FADataset{SpeedMode,mode}(FATrialNumber{SpeedMode,mode}<=MinNumTrial),...
            CRDataset{SpeedMode,mode,ActiveMode}(CRTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)];
        DataNumber1=[FATrialNumber{SpeedMode,mode}(FATrialNumber{SpeedMode,mode}<=MinNumTrial),...
            (CRTrialNumber{SpeedMode,mode,ActiveMode}(CRTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)+MinNumTrial)];
        case 'L'
                MinNumTrial=min( max(HitTrialNumber{SpeedMode,mode}) , max(FATrialNumber{SpeedMode,mode}));
                DataIndex1=[FADataset{SpeedMode,mode}(FATrialNumber{SpeedMode,mode}<=MinNumTrial),...
            HitDataset{SpeedMode,mode}(HitTrialNumber{SpeedMode,mode}<=MinNumTrial)];
        DataNumber1=[FATrialNumber{SpeedMode,mode}(FATrialNumber{SpeedMode,mode}<=MinNumTrial),...
            (HitTrialNumber{SpeedMode,mode}(HitTrialNumber{SpeedMode,mode}<=MinNumTrial)+MinNumTrial)];
        case 'NL'
                MinNumTrial=min( max(CRTrialNumber{SpeedMode,mode,ActiveMode}) , max(MissTrialNumber{SpeedMode,mode,ActiveMode}));
                DataIndex1=[MissDataset{SpeedMode,mode,ActiveMode}(MissTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial),...
            CRDataset{SpeedMode,mode,ActiveMode}(CRTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)];
        DataNumber1=[MissTrialNumber{SpeedMode,mode,ActiveMode}(MissTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial),...
            (CRTrialNumber{SpeedMode,mode,ActiveMode}(CRTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)+MinNumTrial)];
        case 'MAI'
                MinNumTrial=max(MissTrialNumber{SpeedMode,mode,1});
                DataIndex1=[MissDataset{SpeedMode,mode,1},MissDataset{SpeedMode,mode,2}];
        DataNumber1=[MissTrialNumber{SpeedMode,mode,1},(MissTrialNumber{SpeedMode,mode,2}+MinNumTrial)];


end

%%%%%%%%%%%%%%%%%%%%%%%%% Speed Balance
NBin=50;

ActiveAnimal=ones(length(Lick),1);
for i=1:(length(Lick)-75*ActiveTrialNumber)
    if max(Lick(i:(i+75*ActiveTrialNumber)))==0
        ActiveAnimal(i:(i+75*ActiveTrialNumber))=0;
    end
end

[SQ,BQ]=SpeedQuantEqual(Speed,(Speed>0 & ActiveAnimal==1),NBin);
UMask=zeros(size(DataIndex1));
Sorder0=randperm(max(DataNumber0));
Sorder1=randperm(max(DataNumber1));
for kt=1:max(DataNumber0)
        st=Sorder0(kt);
        Sbin0=max(SQ(DataIndex0(DataNumber0==st)));
        matched=0;
        for st1=Sorder1
                Sbin1=max(SQ(DataIndex1(DataNumber1==st1)));
                if Sbin1==Sbin0 && max(UMask(DataNumber1==st1))==0 && matched==0
                        UMask(DataNumber1==st1)=1;
                        matched=1;
                        break;
                end
        end
        if matched==0
                DataNumber0(DataNumber0==st)=0;
        end
end

DataNumber1(UMask==0)=0;

DataSize0=min([sum(DataNumber0>0),sum(DataNumber1>0)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch SET2Name
        case 'H'
                DataIndex0=HitDataset{SpeedMode,mode};
        DataNumber0=HitTrialNumber{SpeedMode,mode};
        case 'M'
                DataIndex0=MissDataset{SpeedMode,mode,ActiveMode};
        DataNumber0=MissTrialNumber{SpeedMode,mode,ActiveMode};
        case 'C'
                DataIndex0=CRDataset{SpeedMode,mode,ActiveMode};
        DataNumber0=CRTrialNumber{SpeedMode,mode,ActiveMode};
        case 'F'
                DataIndex0=FADataset{SpeedMode,mode};
        DataNumber0=FATrialNumber{SpeedMode,mode};
        case 'G'
        MinNumTrial=min( max(MissTrialNumber{SpeedMode,mode,ActiveMode}) , max(HitTrialNumber{SpeedMode,mode}));
                DataIndex0=[HitDataset{SpeedMode,mode}(HitTrialNumber{SpeedMode,mode}<=MinNumTrial),...
            MissDataset{SpeedMode,mode,ActiveMode}(MissTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)];
        DataNumber0=[HitTrialNumber{SpeedMode,mode}(HitTrialNumber{SpeedMode,mode}<=MinNumTrial),...
            (MissTrialNumber{SpeedMode,mode,ActiveMode}(MissTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)+MinNumTrial)];
        case 'N'
                MinNumTrial=min( max(CRTrialNumber{SpeedMode,mode,ActiveMode}) , max(FATrialNumber{SpeedMode,mode}));
                DataIndex0=[FADataset{SpeedMode,mode}(FATrialNumber{SpeedMode,mode}<=MinNumTrial),...
            CRDataset{SpeedMode,mode,ActiveMode}(CRTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)];
        DataNumber0=[FATrialNumber{SpeedMode,mode}(FATrialNumber{SpeedMode,mode}<=MinNumTrial),...
            (CRTrialNumber{SpeedMode,mode,ActiveMode}(CRTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)+MinNumTrial)];
        case 'L'
                MinNumTrial=min( max(HitTrialNumber{SpeedMode,mode}) , max(FATrialNumber{SpeedMode,mode}));
                DataIndex0=[FADataset{SpeedMode,mode}(FATrialNumber{SpeedMode,mode}<=MinNumTrial),...
            HitDataset{SpeedMode,mode}(HitTrialNumber{SpeedMode,mode}<=MinNumTrial)];
        DataNumber0=[FATrialNumber{SpeedMode,mode}(FATrialNumber{SpeedMode,mode}<=MinNumTrial),...
            (HitTrialNumber{SpeedMode,mode}(HitTrialNumber{SpeedMode,mode}<=MinNumTrial)+MinNumTrial)];
        case 'NL'
                MinNumTrial=min( max(CRTrialNumber{SpeedMode,mode,ActiveMode}) , max(MissTrialNumber{SpeedMode,mode,ActiveMode}));
                DataIndex0=[MissDataset{SpeedMode,mode,ActiveMode}(MissTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial),...
            CRDataset{SpeedMode,mode,ActiveMode}(CRTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)];
        DataNumber0=[MissTrialNumber{SpeedMode,mode,ActiveMode}(MissTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial),...
            (CRTrialNumber{SpeedMode,mode,ActiveMode}(CRTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)+MinNumTrial)];
        case 'MAI'
                MinNumTrial=max(MissTrialNumber{SpeedMode,mode,1});
                DataIndex0=[MissDataset{SpeedMode,mode,1},MissDataset{SpeedMode,mode,2}];
        DataNumber0=[MissTrialNumber{SpeedMode,mode,1},(MissTrialNumber{SpeedMode,mode,2}+MinNumTrial)];


end

switch SET3Name
        case 'H'
                DataIndex1=HitDataset{SpeedMode,mode};
        DataNumber1=HitTrialNumber{SpeedMode,mode};
        case 'M'
                DataIndex1=MissDataset{SpeedMode,mode,ActiveMode};
        DataNumber1=MissTrialNumber{SpeedMode,mode,ActiveMode};
        case 'C'
                DataIndex1=CRDataset{SpeedMode,mode,ActiveMode};
        DataNumber1=CRTrialNumber{SpeedMode,mode,ActiveMode};
        case 'F'
                DataIndex1=FADataset{SpeedMode,mode};
        DataNumber1=FATrialNumber{SpeedMode,mode};
        case 'G'
        MinNumTrial=min( max(MissTrialNumber{SpeedMode,mode,ActiveMode}) , max(HitTrialNumber{SpeedMode,mode}));
                DataIndex1=[HitDataset{SpeedMode,mode}(HitTrialNumber{SpeedMode,mode}<=MinNumTrial),...
            MissDataset{SpeedMode,mode,ActiveMode}(MissTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)];
        DataNumber1=[HitTrialNumber{SpeedMode,mode}(HitTrialNumber{SpeedMode,mode}<=MinNumTrial),...
            (MissTrialNumber{SpeedMode,mode,ActiveMode}(MissTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)+MinNumTrial)];
        case 'N'
                MinNumTrial=min( max(CRTrialNumber{SpeedMode,mode,ActiveMode}) , max(FATrialNumber{SpeedMode,mode}));
                DataIndex1=[FADataset{SpeedMode,mode}(FATrialNumber{SpeedMode,mode}<=MinNumTrial),...
            CRDataset{SpeedMode,mode,ActiveMode}(CRTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)];
        DataNumber1=[FATrialNumber{SpeedMode,mode}(FATrialNumber{SpeedMode,mode}<=MinNumTrial),...
            (CRTrialNumber{SpeedMode,mode,ActiveMode}(CRTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)+MinNumTrial)];
        case 'L'
                MinNumTrial=min( max(HitTrialNumber{SpeedMode,mode}) , max(FATrialNumber{SpeedMode,mode}));
                DataIndex1=[FADataset{SpeedMode,mode}(FATrialNumber{SpeedMode,mode}<=MinNumTrial),...
            HitDataset{SpeedMode,mode}(HitTrialNumber{SpeedMode,mode}<=MinNumTrial)];
        DataNumber1=[FATrialNumber{SpeedMode,mode}(FATrialNumber{SpeedMode,mode}<=MinNumTrial),...
            (HitTrialNumber{SpeedMode,mode}(HitTrialNumber{SpeedMode,mode}<=MinNumTrial)+MinNumTrial)];
        case 'NL'
                MinNumTrial=min( max(CRTrialNumber{SpeedMode,mode,ActiveMode}) , max(MissTrialNumber{SpeedMode,mode,ActiveMode}));
                DataIndex1=[MissDataset{SpeedMode,mode,ActiveMode}(MissTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial),...
            CRDataset{SpeedMode,mode,ActiveMode}(CRTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)];
        DataNumber1=[MissTrialNumber{SpeedMode,mode,ActiveMode}(MissTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial),...
            (CRTrialNumber{SpeedMode,mode,ActiveMode}(CRTrialNumber{SpeedMode,mode,ActiveMode}<=MinNumTrial)+MinNumTrial)];
        case 'MAI'
                MinNumTrial=max(MissTrialNumber{SpeedMode,mode,1});
                DataIndex1=[MissDataset{SpeedMode,mode,1},MissDataset{SpeedMode,mode,2}];
        DataNumber1=[MissTrialNumber{SpeedMode,mode,1},(MissTrialNumber{SpeedMode,mode,2}+MinNumTrial)];


end

%%%%%%%%%%%%%%%%%%%%%%%%% Speed Balance
UMask=zeros(size(DataIndex1));
Sorder0=randperm(max(DataNumber0));
Sorder1=randperm(max(DataNumber1));
for kt=1:max(DataNumber0)
        st=Sorder0(kt);
        Sbin0=max(SQ(DataIndex0(DataNumber0==st)));
        matched=0;
        for st1=Sorder1
                Sbin1=max(SQ(DataIndex1(DataNumber1==st1)));
                if Sbin1==Sbin0 && max(UMask(DataNumber1==st1))==0 && matched==0
                        UMask(DataNumber1==st1)=1;
                        matched=1;
                        break;
                end
        end
        if matched==0
                DataNumber0(DataNumber0==st)=0;
        end
end

DataNumber1(UMask==0)=0;

DataSize1=min([sum(DataNumber0>0),sum(DataNumber1>0)]);





DataSize=min([DataSize0,DataSize1])

%end


