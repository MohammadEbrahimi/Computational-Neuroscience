
AddressSetup;

HLength=10;
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
    
    for i=1:size(TrialSE,1)
        if max(ActiveAnimal(TrialSE(i,1):TrialSE(i,2)))==0
            TrialType(i)=0;
        end
    end
    
    [~,sortedIndex]=sort(TrialSE(:,1),'ascend');
    
    TrialSE=TrialSE(sortedIndex,:);
    TrialType=TrialType(sortedIndex);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    start=1;
    for i=1:length(TrialType)
        if  TrialType(i)==0
            start=i+1;
            continue;
        end
        TrialCount(TrialType(i))=TrialCount(TrialType(i))+1;
        if i>start
            HistoryCount(TrialType(start),i-start,TrialType(i))=HistoryCount(TrialType(start),i-start,TrialType(i))+1;
        end
        if (TrialType(i)~=TrialType(start) || (i-start)==HLength)
            start=i;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HistoryProb=HistoryCount;
TrialProb=TrialCount;
TrialProb(1:2)=TrialProb(1:2)/sum(TrialProb(1:2));
TrialProb(3:4)=TrialProb(3:4)/sum(TrialProb(3:4));
HistoryProb(:,:,1:2)=HistoryProb(:,:,1:2) ./ repmat(sum(HistoryCount(:,:,1:2),3),[1,1,2]);
HistoryProb(:,:,3:4)=HistoryProb(:,:,3:4) ./ repmat(sum(HistoryCount(:,:,3:4),3),[1,1,2]);

HistoryProbNorm(:,:,1)=HistoryProb(:,:,1)/TrialProb(1);
HistoryProbNorm(:,:,2)=HistoryProb(:,:,2)/TrialProb(2);
HistoryProbNorm(:,:,3)=HistoryProb(:,:,3)/TrialProb(3);
HistoryProbNorm(:,:,4)=HistoryProb(:,:,4)/TrialProb(4);









