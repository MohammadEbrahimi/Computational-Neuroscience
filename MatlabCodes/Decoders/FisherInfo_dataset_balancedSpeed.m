function FisherInfo_dataset_balancedSpeed(loadpath,savePath,SET0Name,SET1Name,mode,SpeedMode,ActiveMode,TPin,TSin)
Balanced=1;
LowDimVec=1:50;
AreaVec=1:8;
NRepeate=100;
ActiveTrialNumber=20;


mycluster=parcluster('local');
mycluster.NumWorkers=32;
saveProfile(mycluster);
warning('off','all');
load(strcat(loadpath,'/All_Sessions.mat'));
load(strcat(loadpath,'/Calcium.mat'));
load(strcat(loadpath,'/Datasets',TPin,'.mat'));
progress='Data is loaded!'

cellEvents=cellData_Calcium;
cellCount=size(cellEvents,1);
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

%%%%%%%%%%%%%%%%%%%%%%%%Initialize Variables

    Info_Train_Fisher=zeros(8,length(LowDimVec),NRepeate);
    Info_Val_Fisher=zeros(8,length(LowDimVec),NRepeate);
    Info_Train_LR=zeros(8,length(LowDimVec),NRepeate);
    Info_Val_LR=zeros(8,length(LowDimVec),NRepeate);
    Error_Val_LR=zeros(8,length(LowDimVec),NRepeate);
    Error_Train_LR=zeros(8,length(LowDimVec),NRepeate);
    
    Info_Train_Fisher_Sh=zeros(8,length(LowDimVec),NRepeate);
    Info_Val_Fisher_Sh=zeros(8,length(LowDimVec),NRepeate);
    Info_Train_LR_Sh=zeros(8,length(LowDimVec),NRepeate);
    Info_Val_LR_Sh=zeros(8,length(LowDimVec),NRepeate);
    Error_Val_LR_Sh=zeros(8,length(LowDimVec),NRepeate);
    Error_Train_LR_Sh=zeros(8,length(LowDimVec),NRepeate);
    
    BMAT={}
    InterceptMAT=zeros(8,length(LowDimVec),NRepeate);




%%%%%%%%%%%%%%%%%%%%%%%%%%%Start Learning
%%%%DataIndex0,1  DataNumber0,1
poolobj=parpool(mycluster,32)

shuff0=zeros(max(DataNumber0),NRepeate);
shuff1=zeros(max(DataNumber1),NRepeate);

for nr=1:NRepeate
Repeate=nr

shuff0(:,nr)=randperm(max(DataNumber0));
shuff1(:,nr)=randperm(max(DataNumber1));

c0=1;
c1=1;
dd0=zeros(size(cellEvents,1),length(DataIndex0));
for si=shuff0(:,nr)'
    trialLength=sum(DataNumber0==si);
    if(trialLength>0)
    dd0(:,c0:c0+trialLength-1)=cellEvents(:,DataIndex0(DataNumber0==si));
    c0=c0+trialLength; 
    end
end
c0=c0-1;

dd1=zeros(size(cellEvents,1),length(DataIndex1));
for si=shuff1(:,nr)'
    trialLength=sum(DataNumber1==si);
    if(trialLength>0)
    dd1(:,c1:c1+trialLength-1)=cellEvents(:,DataIndex1(DataNumber1==si));
    c1=c1+trialLength;  
    end
end
c1=c1-1;


c=floor(min(c0,c1)/2);
if Balanced==0
    ResampleMask=ones(1,c);
elseif Balanced==1
    ResampleMask=zeros(1,c);
    TSinHalf=floor(TSin/2);
    ResampleMask(1:TSinHalf)=1;
    if(TSinHalf>c)
	strcat('Error 1: TSin=',num2str(TSinHalf),' > c=',num2str(c))
	%return;
    end 
end
RM=[ResampleMask(randperm(c)),ResampleMask(randperm(c))];


TrainingSize=sum(RM);


if TrainingSize<100
strcat('error: ',loadpath,': Small Data Size',num2str(TrainingSize),',',num2str(TSin))
delete(poolobj);
    return;
end


%%%%%%%%%%%%% Decoder
TrainSet=[dd0(:,1:c) dd1(:,1:c)];

PLSDataSet=TrainSet;
ValSet=[dd0(:,c+1:2*c) dd1(:,c+1:2*c)];

% PLSDataSet=[dd0(:,c+1:2*c) dd1(:,c+1:2*c)];
% ValSet=[dd0(:,2*c+1:3*c) dd1(:,2*c+1:3*c)];

PLSDataSetSh=PLSDataSet(:,randperm(2*c));
TrainSetSh=TrainSet(:,randperm(2*c));
ValSetSh=ValSet(:,randperm(2*c));

PLSDataSet=PLSDataSet(:,RM==1);
PLSDataSetSh=PLSDataSetSh(:,RM==1);
TrainSet=TrainSet(:,RM==1);
TrainSetSh=TrainSetSh(:,RM==1);
ValSet=ValSet(:,RM==1);
ValSetSh=ValSetSh(:,RM==1);

hc=sum(RM)/2;


group=[zeros(c,1);ones(c,1)];
group=group(RM==1);

%%%%%%%%% Predefine because matlab parfor is dumb
PLSRot={};PLSRotSh={};
CovMAT={};FisherTuning={};


for area=AreaVec
    TrainingArea=area
    Cells=(CortexArea==area);
    B=zeros(sum(Cells),length(LowDimVec));BSh={};
    intercept=zeros(length(LowDimVec),1);interceptSh=zeros(length(LowDimVec),1);
    ErrCurve=zeros(length(LowDimVec),1);ErrCurveSh=zeros(length(LowDimVec),1);
    smallB={};fullB=zeros(cellCount,length(LowDimVec));

    if sum(Cells)<5
        continue;end
	
       maxDim=min(max(LowDimVec),sum(max(PLSDataSet(Cells,:)')~=min(PLSDataSet(Cells,:)')));
       parfor indld=1:maxDim
    	    [PLSRot{area,indld},~] = plsregress(PLSDataSet(Cells,:)',group,min(LowDimVec(indld),sum(Cells)));
            [PLSRotSh{area,indld},~] = plsregress(PLSDataSetSh(Cells,:)',group,min(LowDimVec(indld),sum(Cells)));
            
            
            
%            [smallB{indld},intercept(indld),ErrCurve(indld),~,~]=glmnoreg(PLSRot{area,indld}' * TrainSet(Cells,:),group,0,5,10,1);
 %           [BSh{indld},interceptSh(indld),ErrCurveSh(indld),~,~]=glmnoreg(PLSRotSh{area,indld}' * TrainSetSh(Cells,:),group,0,5,10,1);
            
            
%            B(:,indld)=PLSRot{area,indld} *smallB{indld};
%            BMAT{area,indld,nr}=B(:,indld);
%            InterceptMAT(area,indld,nr)=intercept(indld);
            
            CovMAT{area,indld,nr}=CovSets(PLSRot{area,indld}' *TrainSet(Cells,1:hc),PLSRot{area,indld}' *TrainSet(Cells,(hc+1):(2*hc)));

            
%            Error_Train_LR(area,indld,nr)=ErrCurve(indld);
            Info_Train_Fisher(area,indld,nr)=FisherInformation(PLSRot{area,indld}' *TrainSet(Cells,1:hc),PLSRot{area,indld}' *TrainSet(Cells,(hc+1):(2*hc)));
%            Info_Train_LR(area,indld,nr)=FisherInformationDecoder(PLSRot{area,indld}' *TrainSet(Cells,1:hc),PLSRot{area,indld}' *TrainSet(Cells,(hc+1):(2*hc)),smallB{indld});
            [Info_Val_Fisher(area,indld,nr),FisherTuning{area,indld,nr}]=FisherInformationVal(PLSRot{area,indld}' *TrainSet(Cells,1:hc),PLSRot{area,indld}' *TrainSet(Cells,(hc+1):(2*hc)),PLSRot{area,indld}' *ValSet(Cells,1:hc),PLSRot{area,indld}' *ValSet(Cells,(hc+1):(2*hc)));
%            Info_Val_LR(area,indld,nr)=FisherInformationDecoder(PLSRot{area,indld}' *ValSet(Cells,1:hc),PLSRot{area,indld}' *ValSet(Cells,(hc+1):(2*hc)),smallB{indld});
%            Error_Val_LR(area,indld,nr)= mean( LRClassify([ValSet(Cells,:)'] , B(:,indld) , intercept(indld)) ~= [zeros(hc,1);ones(hc,1)]);
            
%           Error_Train_LR_Sh(area,indld,nr)=ErrCurveSh(indld);
            Info_Train_Fisher_Sh(area,indld,nr)=FisherInformation(PLSRotSh{area,indld}' *TrainSetSh(Cells,1:hc),PLSRotSh{area,indld}' *TrainSetSh(Cells,(hc+1):(2*hc)));
%            Info_Train_LR_Sh(area,indld,nr)=FisherInformationDecoder(PLSRotSh{area,indld}' *TrainSetSh(Cells,1:hc),PLSRotSh{area,indld}' *TrainSetSh(Cells,(hc+1):(2*hc)),BSh{indld});
            Info_Val_Fisher_Sh(area,indld,nr)=FisherInformationVal(PLSRotSh{area,indld}' *TrainSetSh(Cells,1:hc),PLSRotSh{area,indld}' *TrainSetSh(Cells,(hc+1):(2*hc)),PLSRotSh{area,indld}' *ValSetSh(Cells,1:hc),PLSRotSh{area,indld}' *ValSetSh(Cells,(hc+1):(2*hc)));
%            Info_Val_LR_Sh(area,indld,nr)=FisherInformationDecoder(PLSRotSh{area,indld}' *ValSetSh(Cells,1:hc),PLSRotSh{area,indld}' *ValSetSh(Cells,(hc+1):(2*hc)),BSh{indld});
%            Error_Val_LR_Sh(area,indld,nr)= mean( LRClassify([ValSetSh(Cells,:)'] , PLSRotSh{area,indld} * BSh{indld} , interceptSh(indld)) ~= [zeros(hc,1);ones(hc,1)]);

            %{ 
            %}
            
        end
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end



progress='Done!'




end %%%%% end repeate 


delete(poolobj)

if(TrainingSize>=100)
%{
    save(strcat(savePath,'_mode',num2str(mode)),'BMAT','CovMAT','InterceptMAT','CovMAT',...
        'Info_Train_LR','Info_Val_Fisher','FisherTuning','Info_Val_LR','Error_Val_LR',...
        'PLSRot','PLSRotSh','shuff0','shuff1','LowDimVec','TrainingSize','hc'...
        ,'Error_Train_LR','Info_Train_Fisher',...
        'Error_Train_LR_Sh','Info_Train_Fisher_Sh','Info_Train_LR_Sh','Info_Val_Fisher_Sh','Info_Val_LR_Sh','Error_Val_LR_Sh');
%}
    save(strcat(savePath,'_mode',num2str(mode)),'CovMAT','CovMAT','Info_Val_Fisher','FisherTuning',...
        'PLSRot','PLSRotSh','shuff0','shuff1','LowDimVec','TrainingSize','hc','Info_Train_Fisher','Info_Train_Fisher_Sh','Info_Val_Fisher_Sh','-v7.3');
    progress='Results are saved successfully!'
else
Progress='Insufficient data'
end
end



