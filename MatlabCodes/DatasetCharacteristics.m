StimDelay=-20;
StimCut=15;
ActiveTrialNumber=40;
%%%% Speed * 0.4 is CM/S
SpeedTh=2.5;

AddressSetup;
Days=[51:57];

% NumTrials=zeros(length(Days),6);
% Note='NumTrials: (Sessions,H/M/C/F/Inactive Go/Inactive Nogo)';


for Day=Days;

Day

path=LoadPath{Day};
load(strcat(path,'\All_Sessions.mat'));

Speed=Speed-min(Speed);

% cutOff=173976;
% CRSE=CRSE(CRSE(:,1)<cutOff,:);
% HitSE=HitSE(HitSE(:,1)<cutOff,:);
% FASE=FASE(FASE(:,1)<cutOff,:);
% MissSE=MissSE(MissSE(:,1)<cutOff,:);


CRSE=CRSE(2:end,:);
HitSE=HitSE(2:end,:);
FASE=FASE(2:end,:);
MissSE=MissSE(2:end,:);

ActiveAnimal=zeros(length(Lick),1);
for i=1:(length(Lick)-75*ActiveTrialNumber)
    if max(Lick(i:(i+75*ActiveTrialNumber)))==1
        ActiveAnimal(i:(i+75*ActiveTrialNumber))=1;
    end
end
% 
% ActiveAnimal=ones(length(Lick),1);
% for i=1:(length(Lick)-75*ActiveTrialNumber)
%     if max(Lick(i:(i+75*ActiveTrialNumber)))==0
%         ActiveAnimal(i:(i+75*ActiveTrialNumber))=0;
%     end
% end

clear HitDataset
clear HitTrialNumber  
clear MissDataset
clear MissTrialNumber
clear CRDataset
clear CRTrialNumber
clear FADataset
clear FATrialNumber



DatasetStats=zeros(5,6,3);%(S,H/M/C/F/IG/IN,S/D/R)
HitDataset{2,3}=[];%{Speed,SDR}
HitTrialNumber{2,3}=[];%{Speed,SDR}
MissDataset{2,3,3}=[];%{Speed,SDR,A/I}
MissTrialNumber{2,3,3}=[];%{Speed,SDR,A/I}
CRDataset{2,3,3}=[];%{Speed,SDR,A/I}
CRTrialNumber{2,3,3}=[];%{Speed,SDR,A/I}
FADataset{2,3}=[];%{Speed,SDR}
FATrialNumber{2,3}=[];%{Speed,SDR}

n11=1;
n21=1;
n12=1;
n22=1;
n13=1;
n23=1;

for k=1:size(HitSE,1)
    TrialTime=HitSE(k,1):HitSE(k,2);


    StimStart=TrialTime(1)+StimDelay+(1-Hit(TrialTime(1)));
    StimEnd=StimStart-StimDelay+19-StimCut;
    %StimEnd=TrialTime(diff(Hit(TrialTime))==-1)-StimCut;
    
    RewardStart=TrialTime(diff(RewardWindow(TrialTime))==1)+1;
    RewardEnd=RewardStart+29;
    %RewardEnd=TrialTime(diff(RewardWindow(TrialTime))==-1);
    
    if length(StimEnd)==1 && min(Hit((StimStart-StimDelay):StimEnd+1))==1
        DatasetStats(1,1,1)= DatasetStats(1,1,1)+1;
        HitDataset{1,1}=[HitDataset{1,1},StimStart:StimEnd];
        HitTrialNumber{1,1}=[HitTrialNumber{1,1},n11*ones(1,1+StimEnd-StimStart)];
        n11=n11+1;
        if max(Speed(StimStart:StimEnd))>50
            DatasetStats(2,1,1)= DatasetStats(2,1,1)+1;
        elseif max(Speed(StimStart:StimEnd))>10
            DatasetStats(3,1,1)= DatasetStats(3,1,1)+1;
        elseif max(Speed(StimStart:StimEnd))>SpeedTh
            DatasetStats(4,1,1)= DatasetStats(4,1,1)+1;
        else
            DatasetStats(5,1,1)= DatasetStats(5,1,1)+1;
            HitDataset{2,1}=[HitDataset{2,1},StimStart:StimEnd];
            HitTrialNumber{2,1}=[HitTrialNumber{2,1},n21*ones(1,1+StimEnd-StimStart)];
            n21=n21+1;
        end
    end
    
    StimEnd=TrialTime(diff(Hit(TrialTime))==-1);
    DelayTime=StimEnd(end)+1:StimEnd(end)+5;
    if max(Lick(DelayTime))==0
       DatasetStats(1,1,2)= DatasetStats(1,1,2)+1;
       HitDataset{1,2}=[HitDataset{1,2},DelayTime];
       HitTrialNumber{1,2}=[HitTrialNumber{1,2},n12*ones(1,length(DelayTime))];
       n12=n12+1;
       if max(Speed(DelayTime))>50
           DatasetStats(2,1,2)= DatasetStats(2,1,2)+1;
       elseif max(Speed(DelayTime))>10
           DatasetStats(3,1,2)= DatasetStats(3,1,2)+1;
       elseif max(Speed(DelayTime))>SpeedTh
           DatasetStats(4,1,2)= DatasetStats(4,1,2)+1;
       else
           DatasetStats(5,1,2)= DatasetStats(5,1,2)+1;
           HitDataset{2,2}=[HitDataset{2,2},DelayTime];
           HitTrialNumber{2,2}=[HitTrialNumber{2,2},n22*ones(1,length(DelayTime))];
           n22=n22+1;
       end       
    end
    
    if length(RewardEnd)==1 
       DatasetStats(1,1,3)= DatasetStats(1,1,3)+1;
       HitDataset{1,3}=[HitDataset{1,3},RewardStart:RewardEnd];
       HitTrialNumber{1,3}=[HitTrialNumber{1,3},n13*ones(1,1+RewardEnd-RewardStart)];
       n13=n13+1;
       if max(Speed(RewardStart:RewardEnd))>50
           DatasetStats(2,1,3)= DatasetStats(2,1,3)+1;
       elseif max(Speed(RewardStart:RewardEnd))>10
           DatasetStats(3,1,3)= DatasetStats(3,1,3)+1;
       elseif max(Speed(RewardStart:RewardEnd))>SpeedTh
           DatasetStats(4,1,3)= DatasetStats(4,1,3)+1;
       else
           DatasetStats(5,1,3)= DatasetStats(5,1,3)+1;
           HitDataset{2,3}=[HitDataset{2,3},RewardStart:RewardEnd];
           HitTrialNumber{2,3}=[HitTrialNumber{2,3},n23*ones(1,1+RewardEnd-RewardStart)];
           n23=n23+1;
       end   
    
    end
end
    
    
n11=[1,1,1];
n21=[1,1];
n12=[1,1,1];
n22=[1,1];
n13=[1,1,1];
n23=[1,1];


for k=1:size(CRSE,1)
    TrialTime=CRSE(k,1):CRSE(k,2);
   
    StimStart=TrialTime(1)+StimDelay+(1-CR(TrialTime(1)));
    StimEnd=StimStart-StimDelay+19-StimCut;
    %StimEnd=TrialTime(diff(CR(TrialTime))==-1)-StimCut;
    
    RewardStart=TrialTime(diff(RewardWindow(TrialTime))==1)+1;
    RewardEnd=RewardStart+29;
    %RewardEnd=TrialTime(diff(RewardWindow(TrialTime))==-1);
    Index=3;
    Act=1;
    if max(ActiveAnimal(StimStart:RewardEnd))==0
        Index=[3,6];
        Act=2;
    end
    
    if length(StimEnd)==1 && min(CR((StimStart-StimDelay):StimEnd+1))==1
       DatasetStats(1,Index,1)= DatasetStats(1,Index,1)+1;
       CRDataset{1,1,Act}=[CRDataset{1,1,Act},StimStart:StimEnd];
       CRTrialNumber{1,1,Act}=[CRTrialNumber{1,1,Act},n11(Act)*ones(1,1+StimEnd-StimStart)];
       n11(Act)=n11(Act)+1;
       CRDataset{1,1,3}=[CRDataset{1,1,3},StimStart:StimEnd];
       CRTrialNumber{1,1,3}=[CRTrialNumber{1,1,3},n11(3)*ones(1,1+StimEnd-StimStart)];
       n11(3)=n11(3)+1;
       
       if max(Speed(StimStart:StimEnd))>50
           DatasetStats(2,Index,1)= DatasetStats(2,Index,1)+1;
       elseif max(Speed(StimStart:StimEnd))>10
           DatasetStats(3,Index,1)= DatasetStats(3,Index,1)+1;
       elseif max(Speed(StimStart:StimEnd))>SpeedTh
           DatasetStats(4,Index,1)= DatasetStats(4,Index,1)+1;
       else
           DatasetStats(5,Index,1)= DatasetStats(5,Index,1)+1;
           CRDataset{2,1,Act}=[CRDataset{2,1,Act},StimStart:StimEnd];
           CRTrialNumber{2,1,Act}=[CRTrialNumber{2,1,Act},n21(Act)*ones(1,1+StimEnd-StimStart)];
           n21(Act)=n21(Act)+1;
       end       
    end
    StimEnd=TrialTime(diff(CR(TrialTime))==-1);
    DelayTime=StimEnd(end)+1:StimEnd(end)+5;
    if max(Lick(DelayTime))==0
       DatasetStats(1,Index,2)= DatasetStats(1,Index,2)+1;
       CRDataset{1,2,Act}=[CRDataset{1,2,Act},DelayTime];
       CRTrialNumber{1,2,Act}=[CRTrialNumber{1,2,Act},n12(Act)*ones(1,length(DelayTime))];
       n12(Act)=n12(Act)+1;
       
       CRDataset{1,2,3}=[CRDataset{1,2,3},DelayTime];
       CRTrialNumber{1,2,3}=[CRTrialNumber{1,2,3},n12(3)*ones(1,length(DelayTime))];
       n12(3)=n12(3)+1;
       
       if max(Speed(DelayTime))>50
           DatasetStats(2,Index,2)= DatasetStats(2,Index,2)+1;
       elseif max(Speed(DelayTime))>10
           DatasetStats(3,Index,2)= DatasetStats(3,Index,2)+1;
       elseif max(Speed(DelayTime))>SpeedTh
           DatasetStats(4,Index,2)= DatasetStats(4,Index,2)+1;
       else
           DatasetStats(5,Index,2)= DatasetStats(5,Index,2)+1;
           CRDataset{2,2,Act}=[CRDataset{2,2,Act},DelayTime];
           CRTrialNumber{2,2,Act}=[CRTrialNumber{2,2,Act},n22(Act)*ones(1,length(DelayTime))];
           n22(Act)=n22(Act)+1;
       end       
    end
    
    if length(RewardEnd)==1 
       DatasetStats(1,Index,3)= DatasetStats(1,Index,3)+1;
       CRDataset{1,3,Act}=[CRDataset{1,3,Act},RewardStart:RewardEnd];
       CRTrialNumber{1,3,Act}=[CRTrialNumber{1,3,Act},n13(Act)*ones(1,1+RewardEnd-RewardStart)];
       n13(Act)=n13(Act)+1;
       
       CRDataset{1,3,3}=[CRDataset{1,3,3},RewardStart:RewardEnd];
       CRTrialNumber{1,3,3}=[CRTrialNumber{1,3,3},n13(3)*ones(1,1+RewardEnd-RewardStart)];
       n13(3)=n13(3)+1;
       
       if max(Speed(RewardStart:RewardEnd))>50
           DatasetStats(2,Index,3)= DatasetStats(2,Index,3)+1;
       elseif max(Speed(RewardStart:RewardEnd))>10
           DatasetStats(3,Index,3)= DatasetStats(3,Index,3)+1;
       elseif max(Speed(RewardStart:RewardEnd))>SpeedTh
           DatasetStats(4,Index,3)= DatasetStats(4,Index,3)+1;
       else
           DatasetStats(5,Index,3)= DatasetStats(5,Index,3)+1;
           CRDataset{2,3,Act}=[CRDataset{2,3,Act},RewardStart:RewardEnd];
           CRTrialNumber{2,3,Act}=[CRTrialNumber{2,3,Act},n23(Act)*ones(1,1+RewardEnd-RewardStart)];
           n23(Act)=n23(Act)+1;
       end   
    
    end
end
    
n11=[1,1,1];
n21=[1,1];
n12=[1,1,1];
n22=[1,1];
n13=[1,1,1];
n23=[1,1];

for k=1:size(MissSE,1)
    TrialTime=MissSE(k,1):MissSE(k,2);
   
    StimStart=TrialTime(1)+StimDelay+(1-Miss(TrialTime(1)));
    StimEnd=StimStart-StimDelay+19-StimCut;
    %StimEnd=TrialTime(diff(Miss(TrialTime))==-1)-StimCut;
    
    RewardStart=TrialTime(diff(RewardWindow(TrialTime))==1)+1;
    RewardEnd=RewardStart+29;
    %RewardEnd=TrialTime(diff(RewardWindow(TrialTime))==-1);
    Index=2;
    Act=1;
    if max(ActiveAnimal(StimStart:RewardEnd))==0
        Index=[2,5];
        Act=2;
    end
    
    if length(StimEnd)==1 && min(Miss((StimStart-StimDelay):StimEnd+1))==1
       DatasetStats(1,Index,1)= DatasetStats(1,Index,1)+1;
       MissDataset{1,1,Act}=[MissDataset{1,1,Act},StimStart:StimEnd];
       MissTrialNumber{1,1,Act}=[MissTrialNumber{1,1,Act},n11(Act)*ones(1,1+StimEnd-StimStart)];
       n11(Act)=n11(Act)+1;
       
       MissDataset{1,1,3}=[MissDataset{1,1,3},StimStart:StimEnd];
       MissTrialNumber{1,1,3}=[MissTrialNumber{1,1,3},n11(3)*ones(1,1+StimEnd-StimStart)];
       n11(3)=n11(3)+1;
       
       if max(Speed(StimStart:StimEnd))>50
           DatasetStats(2,Index,1)= DatasetStats(2,Index,1)+1;
       elseif max(Speed(StimStart:StimEnd))>10
           DatasetStats(3,Index,1)= DatasetStats(3,Index,1)+1;
       elseif max(Speed(StimStart:StimEnd))>SpeedTh
           DatasetStats(4,Index,1)= DatasetStats(4,Index,1)+1;
       else
           DatasetStats(5,Index,1)= DatasetStats(5,Index,1)+1;
           MissDataset{2,1,Act}=[MissDataset{2,1,Act},StimStart:StimEnd];
           MissTrialNumber{2,1,Act}=[MissTrialNumber{2,1,Act},n21(Act)*ones(1,1+StimEnd-StimStart)];
           n21(Act)=n21(Act)+1;
       end       
    end
    StimEnd=TrialTime(diff(Miss(TrialTime))==-1);
    DelayTime=StimEnd(end)+1:StimEnd(end)+5;
    if max(Lick(DelayTime))==0
       DatasetStats(1,Index,2)= DatasetStats(1,Index,2)+1;
       MissDataset{1,2,Act}=[MissDataset{1,2,Act},DelayTime];
       MissTrialNumber{1,2,Act}=[MissTrialNumber{1,2,Act},n12(Act)*ones(1,length(DelayTime))];
       n12(Act)=n12(Act)+1;
       
       MissDataset{1,2,3}=[MissDataset{1,2,3},DelayTime];
       MissTrialNumber{1,2,3}=[MissTrialNumber{1,2,3},n12(3)*ones(1,length(DelayTime))];
       n12(3)=n12(3)+1;
       
       if max(Speed(DelayTime))>50
           DatasetStats(2,Index,2)= DatasetStats(2,Index,2)+1;
       elseif max(Speed(DelayTime))>10
           DatasetStats(3,Index,2)= DatasetStats(3,Index,2)+1;
       elseif max(Speed(DelayTime))>SpeedTh
           DatasetStats(4,Index,2)= DatasetStats(4,Index,2)+1;
       else
           DatasetStats(5,Index,2)= DatasetStats(5,Index,2)+1;
           MissDataset{2,2,Act}=[MissDataset{2,2,Act},DelayTime];
           MissTrialNumber{2,2,Act}=[MissTrialNumber{2,2,Act},n22(Act)*ones(1,length(DelayTime))];
           n22(Act)=n22(Act)+1;
       end       
    end
    
    if length(RewardEnd)==1 
       DatasetStats(1,Index,3)= DatasetStats(1,Index,3)+1;
       MissDataset{1,3,Act}=[MissDataset{1,3,Act},RewardStart:RewardEnd];
       MissTrialNumber{1,3,Act}=[MissTrialNumber{1,3,Act},n13(Act)*ones(1,1+RewardEnd-RewardStart)];
       n13(Act)=n13(Act)+1;
       
       MissDataset{1,3,3}=[MissDataset{1,3,3},RewardStart:RewardEnd];
       MissTrialNumber{1,3,3}=[MissTrialNumber{1,3,3},n13(3)*ones(1,1+RewardEnd-RewardStart)];
       n13(3)=n13(3)+1;
       
       if max(Speed(RewardStart:RewardEnd))>50
           DatasetStats(2,Index,3)= DatasetStats(2,Index,3)+1;
       elseif max(Speed(RewardStart:RewardEnd))>10
           DatasetStats(3,Index,3)= DatasetStats(3,Index,3)+1;
       elseif max(Speed(RewardStart:RewardEnd))>SpeedTh
           DatasetStats(4,Index,3)= DatasetStats(4,Index,3)+1;
       else
           DatasetStats(5,Index,3)= DatasetStats(5,Index,3)+1;
           MissDataset{2,3,Act}=[MissDataset{2,3,Act},RewardStart:RewardEnd];
           MissTrialNumber{2,3,Act}=[MissTrialNumber{2,3,Act},n23(Act)*ones(1,1+RewardEnd-RewardStart)];
           n23(Act)=n23(Act)+1;
       end   
    
    end
end


n11=1;
n21=1;
n12=1;
n22=1;
n13=1;
n23=1;


for k=1:size(FASE,1)
    TrialTime=FASE(k,1):FASE(k,2);
   
    StimStart=TrialTime(1)+StimDelay+(1-FA(TrialTime(1)));
    StimEnd=StimStart-StimDelay+19-StimCut;
    %StimEnd=TrialTime(diff(FA(TrialTime))==-1)-StimCut;
    
    RewardStart=TrialTime(diff(RewardWindow(TrialTime))==1)+1;
    RewardEnd=RewardStart+29;
    %RewardEnd=TrialTime(diff(RewardWindow(TrialTime))==-1);
    Index=4;
    
    if length(StimEnd)==1 && min(FA((StimStart-StimDelay):StimEnd+1))==1
       DatasetStats(1,Index,1)= DatasetStats(1,Index,1)+1;
       FADataset{1,1}=[FADataset{1,1},StimStart:StimEnd];
       FATrialNumber{1,1}=[FATrialNumber{1,1},n11*ones(1,1+StimEnd-StimStart)];
       n11=n11+1;
       if max(Speed(StimStart:StimEnd))>50
           DatasetStats(2,Index,1)= DatasetStats(2,Index,1)+1;
       elseif max(Speed(StimStart:StimEnd))>10
           DatasetStats(3,Index,1)= DatasetStats(3,Index,1)+1;
       elseif max(Speed(StimStart:StimEnd))>SpeedTh
           DatasetStats(4,Index,1)= DatasetStats(4,Index,1)+1;
       else
           DatasetStats(5,Index,1)= DatasetStats(5,Index,1)+1;
           FADataset{2,1}=[FADataset{2,1},StimStart:StimEnd];
           FATrialNumber{2,1}=[FATrialNumber{2,1},n21*ones(1,1+StimEnd-StimStart)];
           n21=n21+1;
       end       
    end
    StimEnd=TrialTime(diff(FA(TrialTime))==-1);
    DelayTime=StimEnd(end)+1:StimEnd(end)+5;
    if max(Lick(DelayTime))==0
       DatasetStats(1,Index,2)= DatasetStats(1,Index,2)+1;
       FADataset{1,2}=[FADataset{1,2},DelayTime];
       FATrialNumber{1,2}=[FATrialNumber{1,2},n12*ones(1,length(DelayTime))];
       n12=n12+1;
       if max(Speed(DelayTime))>50
           DatasetStats(2,Index,2)= DatasetStats(2,Index,2)+1;
       elseif max(Speed(DelayTime))>10
           DatasetStats(3,Index,2)= DatasetStats(3,Index,2)+1;
       elseif max(Speed(DelayTime))>SpeedTh
           DatasetStats(4,Index,2)= DatasetStats(4,Index,2)+1;
       else
           DatasetStats(5,Index,2)= DatasetStats(5,Index,2)+1;
           FADataset{2,2}=[FADataset{2,2},DelayTime];
           FATrialNumber{2,2}=[FATrialNumber{2,2},n22*ones(1,length(DelayTime))];
           n22=n22+1;
       end       
    end
    
    if length(RewardEnd)==1 
       DatasetStats(1,Index,3)= DatasetStats(1,Index,3)+1;
       FADataset{1,3}=[FADataset{1,3},RewardStart:RewardEnd];
       FATrialNumber{1,3}=[FATrialNumber{1,3},n13*ones(1,1+RewardEnd-RewardStart)];
       n13=n13+1;
       if max(Speed(RewardStart:RewardEnd))>50
           DatasetStats(2,Index,3)= DatasetStats(2,Index,3)+1;
       elseif max(Speed(RewardStart:RewardEnd))>10
           DatasetStats(3,Index,3)= DatasetStats(3,Index,3)+1;
       elseif max(Speed(RewardStart:RewardEnd))>SpeedTh
           DatasetStats(4,Index,3)= DatasetStats(4,Index,3)+1;
       else
           DatasetStats(5,Index,3)= DatasetStats(5,Index,3)+1;
           FADataset{2,3}=[FADataset{2,3},RewardStart:RewardEnd];
           FATrialNumber{2,3}=[FATrialNumber{2,3},n23*ones(1,1+RewardEnd-RewardStart)];
           n23=n23+1;
       end   
    
    end
end
Note='{Speed,S/D/R,A/I}';
NumTrials(Day-Days(1)+1,:)=DatasetStats(5,:,1);
  save(strcat(path,'/Datasets/Datasets-',num2str(StimDelay),'to',num2str(StimCut)),'SpeedTh','HitDataset','HitTrialNumber','CRDataset','CRTrialNumber',...
    'MissDataset','MissTrialNumber','FADataset','FATrialNumber','Note','-v7.3');
%save('C:\Users\Sadegh\Documents\VLMReborn\Reports2\2017_10_26_UsefulData\Mouse7_S2','NumTrials','Note');

end