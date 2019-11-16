function FisherInfo_Instant(loadpath,savePath,SET0Name,SET1Name,mode,SpeedMode,TSin,CellT)

Balanced=1;
LowDimVec=1:50;
AreaVec=1:9;
NRepeate=100;
ActiveTrialNumber=20;
SpeedBalance=1;
ActiveMode=1;

mycluster=parcluster('local');
mycluster.NumWorkers=32;
saveProfile(mycluster);
warning('off','all');
load(strcat(loadpath,'/All_Sessions.mat'));
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

%%%%%%%%%%%%%%%%%%%%%%%%Initialize Variables

    Info_Train_Fisher=zeros(9,length(TimeVec),length(LowDimVec),NRepeate);
    Info_Val_Fisher=zeros(9,length(TimeVec),length(LowDimVec),NRepeate);
    Error_Val_Fisher=zeros(9,length(TimeVec),length(LowDimVec),NRepeate);
    Noise_Val_Fisher=zeros(9,length(TimeVec),length(LowDimVec),NRepeate);
    Signal_Val_Fisher=zeros(9,length(TimeVec),length(LowDimVec),NRepeate);
    Info_Train_LR=zeros(9,length(LowDimVec),NRepeate);
    Info_Val_LR=zeros(9,length(LowDimVec),NRepeate);
    Error_Val_LR=zeros(9,length(LowDimVec),NRepeate);
    Error_Train_LR=zeros(9,length(LowDimVec),NRepeate);
    
    Info_Train_Fisher_Sh=zeros(9,length(TimeVec),length(LowDimVec),NRepeate);
    Info_Val_Fisher_Sh=zeros(9,length(TimeVec),length(LowDimVec),NRepeate);
    Noise_Val_Fisher_Sh=zeros(9,length(TimeVec),length(LowDimVec),NRepeate);
    Signal_Val_Fisher_Sh=zeros(9,length(TimeVec),length(LowDimVec),NRepeate);
    Info_Train_LR_Sh=zeros(9,length(LowDimVec),NRepeate);
    Info_Val_LR_Sh=zeros(9,length(LowDimVec),NRepeate);
    Error_Val_LR_Sh=zeros(9,length(LowDimVec),NRepeate);
    Error_Train_LR_Sh=zeros(9,length(LowDimVec),NRepeate);

    MaxDecoders={};

    
    BMAT={}
    InterceptMAT=zeros(8,length(LowDimVec),NRepeate);
CovMAT={};FisherTuning={};





%%%%%%%%%%%%%%%%%%%%%%%%%%%Start Learning
%%%%DataIndex0,1  DataNumber0,1
poolobj=parpool(mycluster,32)

shuff0=zeros(max(DataNumber0),NRepeate);
shuff1=zeros(max(DataNumber1),NRepeate);

for nr=1:NRepeate

	Repeate=nr
	DataNumber0=Buf_DataNumber0;
	DataNumber1=Buf_DataNumber1;
	shuff0(:,nr)=randperm(max(DataNumber0));
	shuff1(:,nr)=randperm(max(DataNumber1));

	Sorder0=shuff0(:,nr)';
	Sorder1=shuff1(:,nr)';
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

	
	for tind=TimeVec


		c0=1;
		c1=1;
		dd0=zeros(size(cellEvents,1),length(DataIndex0));
		for si=shuff0(:,nr)'
		    trialLength=sum(DataNumber0==si);
		    if(trialLength>=tind)
		    	dd0(:,c0)=cellEvents(:,min(DataIndex0(DataNumber0==si))+tind-1);
		    	c0=c0+1; 
		    end
		end
		c0=c0-1;
		dd0=dd0(:,1:c0);

		dd1=zeros(size(cellEvents,1),length(DataIndex1));
		for si=shuff1(:,nr)'
		    trialLength=sum(DataNumber1==si);
		    if(trialLength>=tind)
		    dd1(:,c1)=cellEvents(:,min(DataIndex1(DataNumber1==si))+tind-1);
		    c1=c1+1;  
		    end
		end
		c1=c1-1;
		dd1=dd1(:,1:c1);

		devideSets=3;

		c=floor(min(c0,c1)/devideSets);
		if Balanced==0
		    ResampleMask=ones(1,c);
		elseif Balanced==1
		    ResampleMask=zeros(1,c);
		    TSinHalf=floor(TSin/devideSets);
		    ResampleMask(1:TSinHalf)=1;
		    if(TSinHalf>c)
			strcat('Error 1: TSin=',num2str(TSinHalf),' > c=',num2str(c))
			%return;
		    end 
		end
		RM=[ResampleMask(randperm(c)),ResampleMask(randperm(c))];


		TrainingSize=sum(RM);


		if TrainingSize<max(LowDimVec)+2
		strcat('error: ',loadpath,': Small Data Size',num2str(TrainingSize),',',num2str(TSin))
		delete(poolobj);
		    return;
		end


		%%%%%%%%%%%%% Decoder
		TrainSet=[dd0(:,1:c) dd1(:,1:c)];

		%PLSDataSet=TrainSet;
		%ValSet=[dd0(:,c+1:2*c) dd1(:,c+1:2*c)];

		PLSDataSet=[dd0(:,c+1:2*c) dd1(:,c+1:2*c)];
		ValSet=[dd0(:,2*c+1:3*c) dd1(:,2*c+1:3*c)];

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


		for area=AreaVec
		    TrainingArea=area
		    if area<9
			Cells=(CortexArea==area);
		    else
			Cells=(CortexArea<=area);
		    end
		    B=zeros(sum(Cells),length(LowDimVec));BSh={};
		    intercept=zeros(length(LowDimVec),1);interceptSh=zeros(length(LowDimVec),1);
		    ErrCurve=zeros(length(LowDimVec),1);ErrCurveSh=zeros(length(LowDimVec),1);
		    smallB={};fullB=zeros(cellCount,length(LowDimVec));

		    if sum(Cells)<=20
			continue;end
			
		       maxDim=min(max(LowDimVec),sum(max(PLSDataSet(Cells,:)')~=min(PLSDataSet(Cells,:)')));
		       parfor indld=1:maxDim
			    [PLSRot{area,indld},~] = plsregress(PLSDataSet(Cells,:)',group,min(LowDimVec(indld),sum(Cells)));
			    [PLSRotSh{area,indld},~] = plsregress(PLSDataSetSh(Cells,:)',group,min(LowDimVec(indld),sum(Cells)));
			    
			   Info_Train_Fisher(area,tind,indld,nr)=FisherInformation(PLSRot{area,indld}' *TrainSet(Cells,1:hc),PLSRot{area,indld}' *TrainSet(Cells,(hc+1):(2*hc)));
			   
		           [Info_Val_Fisher(area,tind,indld,nr),FisherTuning{area,indld},Signal_Val_Fisher(area,tind,indld,nr),Noise_Val_Fisher(area,tind,indld,nr)]=FisherInformationVal(PLSRot{area,indld}' *TrainSet(Cells,1:hc),PLSRot{area,indld}' *TrainSet(Cells,(hc+1):(2*hc)),PLSRot{area,indld}' *ValSet(Cells,1:hc),PLSRot{area,indld}' *ValSet(Cells,(hc+1):(2*hc)));
			    
		           Info_Val_Fisher_Sh(area,tind,indld,nr)=FisherInformationVal(PLSRotSh{area,indld}' *TrainSetSh(Cells,1:hc),PLSRotSh{area,indld}' *TrainSetSh(Cells,(hc+1):(2*hc)),PLSRotSh{area,indld}' *ValSetSh(Cells,1:hc),PLSRotSh{area,indld}' *ValSetSh(Cells,(hc+1):(2*hc)));
		   

			  Error_Val_Fisher(area,tind,indld,nr)=sum((( FisherTuning{area,indld}' * [PLSRot{area,indld}' *ValSet(Cells,1:hc),PLSRot{area,indld}' *ValSet(Cells,(hc+1):(2*hc))]) > mean( FisherTuning{area,indld}' * [PLSRot{area,indld}' *ValSet(Cells,1:hc),PLSRot{area,indld}' *ValSet(Cells,(hc+1):(2*hc))])) ~= [zeros(1,hc),ones(1,hc)] ) / (2*hc);
 
		    end
		

		[~,maxIdx]=max(Info_Val_Fisher(area,:,nr));	   
	        if maxDim>5  && min(size(PLSRot{area,maxIdx}))>=1 && size(PLSRot{area,maxIdx} ,2) == size(FisherTuning{area,maxIdx},1)  
		%	size( FisherTuning{area,maxIdx})
		%	size(PLSRot{area,maxIdx}) 
			MaxDecoders{area,tind,nr}= PLSRot{area,maxIdx} *  FisherTuning{area,maxIdx};
		end

 
	       end

	progress='Done!'

	end%%%TimeVec




end %%%%% end repeate 


delete(poolobj)

if(TrainingSize> max(LowDimVec)+2)

    save(strcat(savePath,'_mode',num2str(mode)),'Info_Val_Fisher','Info_Val_Fisher_Sh','Signal_Val_Fisher','Noise_Val_Fisher',...
        'Error_Val_Fisher','MaxDecoders','LowDimVec','TrainingSize','hc','Info_Train_Fisher','-v7.3');
    progress='Results are saved successfully!'
else
Progress='Insufficient data'
end
end



