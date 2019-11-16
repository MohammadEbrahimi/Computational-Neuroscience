function [DataSize]=EqualSetSize(loadpath,SET0Name,SET1Name,SET2Name,SET3Name,mode,SpeedMode,ActiveMode)
load(strcat(loadpath,'/Datasets.mat'));
Sets={SET0Name,SET1Name,SET2Name,SET3Name};
setsize=zeros(4,1);
%%% creating data sets
for i=1:4
    SETName=Sets{i}
switch SETName
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

    setsize(i)=length(DataIndex0);
end

DataSize=min(setsize)
end


