AddressSetup;
clear HDatasets;
clear BMtx;
HLength=10;
NRepeate=100;
HDatasets{4,HLength}=[];
BMtx{HLength,NRepeate,2}=[];
PredictionError=zeros(HLength,NRepeate,2);
PredictionError_Sh=zeros(HLength,NRepeate,2);

HistoryCount=zeros(4,HLength,4);
TrialCount=zeros(4,1);


for Mouse=61:67;
    load(strcat(LoadPath{Mouse-10},'\All_Sessions.mat'));
    
    ActiveTrialNumber=20;
    ActiveAnimal=ones(length(Lick),1);
    for i=1:(length(Lick)-75*ActiveTrialNumber)
        if max(Lick(i:(i+75*ActiveTrialNumber)))==0
            ActiveAnimal(i:(i+75*ActiveTrialNumber))=0;
        end
    end
    
    %%%%% Hit=1 Miss=2 CR=3 FA=4 In-Active=0 %%%%%
    TrialSE=[HitSE;MissSE;CRSE;FASE];
    TrialType=[ones(size(HitSE,1),1);2*ones(size(MissSE,1),1);3*ones(size(CRSE,1),1);4*ones(size(FASE,1),1)];
    StimuliType=[ones(size(HitSE,1),1);ones(size(MissSE,1),1);zeros(size(CRSE,1),1);zeros(size(FASE,1),1)];
    
    ResponseType=[ones(size(HitSE,1),1);zeros(size(MissSE,1),1);zeros(size(CRSE,1),1);ones(size(FASE,1),1)];
    
    for i=1:size(TrialSE,1)
        if max(ActiveAnimal(TrialSE(i,1):TrialSE(i,2)))==0
            TrialType(i)=0;
        end
    end
    
    [~,sortedIndex]=sort(TrialSE(:,1),'ascend');
    
    TrialSE=TrialSE(sortedIndex,:);
    TrialType=TrialType(sortedIndex);
    StimuliType=StimuliType(sortedIndex);
    ResponseType=ResponseType(sortedIndex);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for HL=1:HLength
        for i=HL+1:length(TrialType)
            if min(TrialType((i-HL):i))>0
                HDatasets{TrialType(i),HL}=[HDatasets{TrialType(i),HL};[StimuliType(i-HL:i-1)'-ResponseType(i-HL:i-1)']];
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for HL=1:HLength
        for nr=1:NRepeate
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%No Go
            Nt=min(size(HDatasets{3,HL},1),size(HDatasets{4,HL},1));
            GoData=[HDatasets{3,HL}(randperm(size(HDatasets{3,HL},1),Nt),:);HDatasets{4,HL}(randperm(size(HDatasets{4,HL},1),Nt),:)];
            Group=[ones(Nt,1);2*ones(Nt,1)];
            TrainSet=GoData([1:round(Nt/2),Nt+1:Nt+round(Nt/2)],:);
            TrainGroup=Group([1:round(Nt/2),Nt+1:Nt+round(Nt/2)]);
            B=mnrfit(TrainSet,TrainGroup);
            B_Sh=mnrfit(TrainSet,TrainGroup(randperm(length(TrainGroup))));
            
            
            ValSet=GoData([round(Nt/2)+1:Nt,Nt+1+round(Nt/2):2*Nt],:);
            ValGroup=Group([round(Nt/2)+1:Nt,Nt+1+round(Nt/2):2*Nt]);
            pihat=mnrval(B,ValSet);
            [~,Prediction]=max(pihat,[],2);
            pihat=mnrval(B_Sh,ValSet);
            [~,Prediction_Sh]=max(pihat,[],2);
            BMtx{HL,nr,1}=B;
            PredictionError(HL,nr,1)=sum(ValGroup~=Prediction)/length(ValGroup);
            PredictionError_Sh(HL,nr,1)=sum(ValGroup~=Prediction_Sh)/length(ValGroup);
            
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Go
            Nt=min(size(HDatasets{2,HL},1),size(HDatasets{2,HL},1));
            GoData=[HDatasets{2,HL}(randperm(size(HDatasets{2,HL},1),Nt),:);HDatasets{1,HL}(randperm(size(HDatasets{1,HL},1),Nt),:)];
            Group=[ones(Nt,1);2*ones(Nt,1)];
            TrainSet=GoData([1:round(Nt/2),Nt+1:Nt+round(Nt/2)],:);
            TrainGroup=Group([1:round(Nt/2),Nt+1:Nt+round(Nt/2)]);
            B=mnrfit(TrainSet,TrainGroup);
            B_Sh=mnrfit(TrainSet,TrainGroup(randperm(length(TrainGroup))));
            
            
            ValSet=GoData([round(Nt/2)+1:Nt,Nt+1+round(Nt/2):2*Nt],:);
            ValGroup=Group([round(Nt/2)+1:Nt,Nt+1+round(Nt/2):2*Nt]);
            pihat=mnrval(B,ValSet);
            [~,Prediction]=max(pihat,[],2);
            pihat=mnrval(B_Sh,ValSet);
            [~,Prediction_Sh]=max(pihat,[],2);
            
            BMtx{HL,nr,2}=B;
            PredictionError(HL,nr,2)=sum(ValGroup~=Prediction)/length(ValGroup);
            PredictionError_Sh(HL,nr,2)=sum(ValGroup~=Prediction_Sh)/length(ValGroup);
            
            
            
            
            
        end
    end
    save(strcat('C:\Users\Sadegh\Documents\VLMReborn\MatlabCodes\SensoryPrior\2018_11_14_BehaviorPriors\ErrorMouse',num2str(Mouse-60)), ...
        'PredictionError','PredictionError_Sh','BMtx','-v7.3');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MeanPredictionError=zeros(HLength,2,7);
MeanShufflePredictionError=zeros(HLength,2,7);
SDPredictionError=zeros(HLength,2,7);
SDShufflePredictionError=zeros(HLength,2,7);
HPriorMtx{HLength,2,7}=[];
figure()
for Mouse=61:67
    load(strcat('C:\Users\Sadegh\Documents\VLMReborn\MatlabCodes\SensoryPrior\2018_11_14_BehaviorPriors\ErrorMouse',num2str(Mouse-60),'.mat'))
    for HL=1:HLength
            HPriorMtx{HL,1,Mouse-60}=(BMtx{HL,1,1}/NRepeate);
            HPriorMtx{HL,2,Mouse-60}=(BMtx{HL,1,2}/NRepeate);
        
        
%         for nr=2:NRepeate
%             HPriorMtx{HL,1,Mouse-60}=HPriorMtx{HL,1,Mouse-60}+(BMtx{HL,nr,1}/NRepeate);
%             HPriorMtx{HL,2,Mouse-60}=HPriorMtx{HL,2,Mouse-60}+(BMtx{HL,nr,2}/NRepeate);
%         end
    end
    MeanPredictionError(:,:,Mouse-60)=mean(PredictionError,2);
    MeanShufflePredictionError(:,:,Mouse-60)=mean(PredictionError_Sh(:,:,:),2);
    
    SDPredictionError(:,1,Mouse-60)=sqrt(var(squeeze(PredictionError(:,1,:))')')/10;
    SDShufflePredictionError(:,1,Mouse-60)=sqrt(var(squeeze(PredictionError_Sh(:,1,:))')')/10;
    
    SDPredictionError(:,2,Mouse-60)=sqrt(var(squeeze(PredictionError(:,2,:))')')/10;
    SDShufflePredictionError(:,2,Mouse-60)=sqrt(var(squeeze(PredictionError_Sh(:,2,:))')')/10;
    
    
    subplot(1,7,Mouse-60);hold on
    shadedErrorBar(1:10,MeanPredictionError(:,1,Mouse-60),SDPredictionError(:,1,Mouse-60),{'Color','k'},1);
    shadedErrorBar(1:10,MeanPredictionError(:,2,Mouse-60),SDPredictionError(:,2,Mouse-60),{'Color','b'},1);
    shadedErrorBar(1:10,MeanShufflePredictionError(:,1,Mouse-60),SDShufflePredictionError(:,1,Mouse-60),{'Color','k','LineStyle','--'},1);
    shadedErrorBar(1:10,MeanShufflePredictionError(:,2,Mouse-60),SDShufflePredictionError(:,2,Mouse-60),{'Color','b','LineStyle','--'},1);
     title(strcat('Mouse ',num2str(Mouse-60)));
     if Mouse==61
    xlabel('History length');
    ylabel('prediction error for next response');
     end
    ylim([0.2,0.55]);
    %     plot(MeanPredictionError(:,2,Mouse-60),'b');
    %     plot(MinShufflePredictionError(:,1,Mouse-60),'k--');
    %     plot(MinShufflePredictionError(:,2,Mouse-60),'b--');
end



figure();
HL=10;
for M=1:7
    subplot(1,7,M);hold on
    
    plot(HPriorMtx{HL,1,M}(2:1+HL),'k')

    
    plot(HPriorMtx{HL,2,M}(2:1+HL),'b')

    
end




