clear all
StimDelay=-10;
ActiveTrialNumber=20;
SpeedTh=2;
Mouse=1;
AddressSetup;
ShFlag=0;
num_shuff=1000;

Days=[17,18,19,47,48,49,50];
areaNames={'V1','LV','MV','PPC','A','S','M','RSC','Border'};
Day=Days(Mouse);

path=LoadPath{Day};
load(strcat(path,'\All_Sessions.mat'));
load(strcat(path,'\registerationUnion.mat'));
load(strcat(path,'\Calcium.mat'));
%  load(strcat(path,'\areas.mat'));

cellCount=size(cellData_Raw,1);
cellIJ=zeros(cellCount,2);
cellIm={};
for i=1:cellCount
    cnum=cellNumbers(i);
    cellIm{i}=cellImage{cnum};
    
    [row,ind]=max(cellImage{cnum});
    [~,indj]=max(row);
    indi=ind(indj);
    cellIJ(i,:)=TileTopLeft(i,:)+[indi-1,indj-1];
end



X=cellData_Calcium;
% X=[Speed';Speed'];
% cellCount=2;
ind_sh={};
for r=1:num_shuff
    ind_sh{r}=randperm(size(X,2));
end



CRSE=CRSE(2:end,:);
HitSE=HitSE(2:end,:);
FASE=FASE(2:end,:);
MissSE=MissSE(2:end,:);


ActiveAnimal=ones(length(Lick),1);
for i=1:(length(Lick)-75*ActiveTrialNumber)
    if max(Lick(i:(i+75*ActiveTrialNumber)))==0
        ActiveAnimal(i:(i+75*ActiveTrialNumber))=0;
    end
end



RespH=zeros(cellCount,65);
RespM=zeros(cellCount,65);
RespC=zeros(cellCount,65);
RespF=zeros(cellCount,65);
qRespH=zeros(cellCount,65);
qRespM=zeros(cellCount,65);
qRespC=zeros(cellCount,65);
qRespF=zeros(cellCount,65);
cH=zeros(cellCount,65);
cM=zeros(cellCount,65);
cC=zeros(cellCount,65);
cF=zeros(cellCount,65);

RespH_binned=zeros(cellCount,13);
RespM_binned=zeros(cellCount,13);
RespC_binned=zeros(cellCount,13);
RespF_binned=zeros(cellCount,13);
qRespH_binned=zeros(cellCount,13);
qRespM_binned=zeros(cellCount,13);
qRespC_binned=zeros(cellCount,13);
qRespF_binned=zeros(cellCount,13);

cH_binned=zeros(cellCount,13);
cM_binned=zeros(cellCount,13);
cC_binned=zeros(cellCount,13);
cF_binned=zeros(cellCount,13);

RespH_binned_sh=zeros(cellCount,13,num_shuff);
RespM_binned_sh=zeros(cellCount,13,num_shuff);
RespC_binned_sh=zeros(cellCount,13,num_shuff);
RespF_binned_sh=zeros(cellCount,13,num_shuff);
qRespH_binned_sh=zeros(cellCount,13,num_shuff);
qRespM_binned_sh=zeros(cellCount,13,num_shuff);
qRespC_binned_sh=zeros(cellCount,13,num_shuff);
qRespF_binned_sh=zeros(cellCount,13,num_shuff);



semH=zeros(cellCount,65);
semM=zeros(cellCount,65);
semC=zeros(cellCount,65);
semF=zeros(cellCount,65);


progress='data loaded'


for k=1:size(HitSE,1)
    floor(100*k/size(HitSE,1))
    TrialTime=HitSE(k,1):HitSE(k,2);

   
    StimStart=TrialTime(1)+1+StimDelay;
    %StimEnd=TrialTime(diff(Hit(TrialTime))==-1);
    StimEnd=StimStart+29;
    
    RewardStart=TrialTime(diff(RewardWindow(TrialTime))==1)+1;
    RewardEnd=TrialTime(diff(RewardWindow(TrialTime))==-1);
    
    StimTime=StimStart:StimEnd;
    if length(StimEnd)==1 && Hit(TrialTime(1)+1)==1 && max(Lick(StimTime))==0
       if max(Speed(StimTime))<SpeedTh && max(ActiveAnimal(StimTime))==1
                RespH(:,1:length(StimTime))=RespH(:,1:length(StimTime))+X(:,StimTime);
                qRespH(:,1:length(StimTime))=qRespH(:,1:length(StimTime))+(X(:,StimTime).^2);
                cH(:,1:length(StimTime))=cH(:,1:length(StimTime))+1;
                for bin=1:6
                    RespH_binned(:,bin)=RespH_binned(:,bin)+mean(X(:,StimStart+(bin-1)*5:min(StimStart+bin*5,StimEnd)),2);
                    qRespH_binned(:,bin)=qRespH_binned(:,bin)+(mean(X(:,StimStart+(bin-1)*5:min(StimStart+bin*5,StimEnd)),2).^2);
                    cH_binned(:,bin)=cH_binned(:,bin)+1;
                    if ShFlag==1
                    for r=1:num_shuff
                        time_sh=ind_sh{r}(StimStart+(bin-1)*5:min(StimStart+bin*5,StimEnd));
                        RespH_binned_sh(:,bin,r)=RespH_binned_sh(:,bin,r)+mean(X(:,time_sh),2);
                        qRespH_binned_sh(:,bin,r)=qRespH_binned_sh(:,bin,r)+(mean(X(:,time_sh),2).^2);
                    end    
                    end
                end
       end       
    end
    StimEnd=TrialTime(diff(Hit(TrialTime))==-1);
    DelayTime=StimEnd(end)+1:RewardStart-1;
    if max(Lick(DelayTime))==0
        if max(Speed(DelayTime))<SpeedTh && max(ActiveAnimal(DelayTime))==1
                RespH(:,31:30+length(DelayTime))=RespH(:,31:30+length(DelayTime))+X(:,DelayTime);
                qRespH(:,31:30+length(DelayTime))=qRespH(:,31:30+length(DelayTime))+(X(:,DelayTime).^2);
                cH(:,31:30+length(DelayTime))=cH(:,31:30+length(DelayTime))+1;
                bin=7;
                    RespH_binned(:,bin)=RespH_binned(:,bin)+mean(X(:,DelayTime),2);
                    qRespH_binned(:,bin)=qRespH_binned(:,bin)+(mean(X(:,DelayTime),2).^2);
                    cH_binned(:,bin)=cH_binned(:,bin)+1;
                    if ShFlag==1
                    for r=1:num_shuff
                        time_sh=ind_sh{r}(DelayTime);
                        RespH_binned_sh(:,bin,r)=RespH_binned_sh(:,bin,r)+mean(X(:,time_sh),2);
                        qRespH_binned_sh(:,bin,r)=qRespH_binned_sh(:,bin,r)+(mean(X(:,time_sh),2).^2);
                    end 
                    end
       end       
    end
    
    RewardTime=RewardStart:RewardEnd;
    RewardTime=RewardTime(1:min(length(RewardTime),31));
    if length(RewardEnd)==1 
        if max(Speed(RewardTime))<SpeedTh && max(ActiveAnimal(RewardTime))==1
                RespH(:,35:34+length(RewardTime))=RespH(:,35:34+length(RewardTime))+X(:,RewardTime);
                qRespH(:,35:34+length(RewardTime))=qRespH(:,35:34+length(RewardTime))+(X(:,RewardTime).^2);
                cH(:,35:34+length(RewardTime))=cH(:,35:34+length(RewardTime))+1;
                for bin=8:13
                    RespH_binned(:,bin)=RespH_binned(:,bin)+mean(X(:,RewardStart+(bin-8)*5:min(RewardStart+(bin-7)*5,RewardEnd)),2);
                    qRespH_binned(:,bin)=qRespH_binned(:,bin)+(mean(X(:,RewardStart+(bin-8)*5:min(RewardStart+(bin-7)*5,RewardEnd)),2).^2);
                    cH_binned(:,bin)=cH_binned(:,bin)+1;
                    if ShFlag==1
                    for r=1:num_shuff
                        time_sh=ind_sh{r}(RewardTime);
                        RespH_binned_sh(:,bin,r)=RespH_binned_sh(:,bin,r)+mean(X(:,time_sh),2);
                        qRespH_binned_sh(:,bin,r)=qRespH_binned_sh(:,bin,r)+(mean(X(:,time_sh),2).^2);
                    end   
                    end
                end
       end    
    
    end
end
progress='Hit done!'
for k=1:size(MissSE,1)
    floor(100*k/size(MissSE,1))
    TrialTime=MissSE(k,1):MissSE(k,2);

   
    StimStart=TrialTime(1)+1+StimDelay;
    %StimEnd=TrialTime(diff(Miss(TrialTime))==-1);
    StimEnd=StimStart+29;
    
    RewardStart=TrialTime(diff(RewardWindow(TrialTime))==1)+1;
    RewardEnd=TrialTime(diff(RewardWindow(TrialTime))==-1);
    
    StimTime=StimStart:StimEnd;
    if length(StimEnd)==1 && Miss(TrialTime(1)+1)==1 && max(Lick(StimTime))==0
       if max(Speed(StimTime))<SpeedTh && max(ActiveAnimal(StimTime))==1
                RespM(:,1:length(StimTime))=RespM(:,1:length(StimTime))+X(:,StimTime);
                qRespM(:,1:length(StimTime))=qRespM(:,1:length(StimTime))+(X(:,StimTime).^2);
                cM(:,1:length(StimTime))=cM(:,1:length(StimTime))+1;
                for bin=1:6
                    RespM_binned(:,bin)=RespM_binned(:,bin)+mean(X(:,StimStart+(bin-1)*5:min(StimStart+bin*5,StimEnd)),2);
                    qRespM_binned(:,bin)=qRespM_binned(:,bin)+(mean(X(:,StimStart+(bin-1)*5:min(StimStart+bin*5,StimEnd)),2).^2);
                    cM_binned(:,bin)=cM_binned(:,bin)+1;
                    if ShFlag==1
                    for r=1:num_shuff
                        time_sh=ind_sh{r}(StimStart+(bin-1)*5:min(StimStart+bin*5,StimEnd));
                        RespM_binned_sh(:,bin,r)=RespM_binned_sh(:,bin,r)+mean(X(:,time_sh),2);
                        qRespM_binned_sh(:,bin,r)=qRespM_binned_sh(:,bin,r)+(mean(X(:,time_sh),2).^2);
                    end
                    end
                end
       end       
    end
    StimEnd=TrialTime(diff(Miss(TrialTime))==-1);
    DelayTime=StimEnd(end)+1:RewardStart-1;
    if max(Lick(DelayTime))==0
        if max(Speed(DelayTime))<SpeedTh && max(ActiveAnimal(DelayTime))==1
                RespM(:,31:30+length(DelayTime))=RespM(:,31:30+length(DelayTime))+X(:,DelayTime);
                qRespM(:,31:30+length(DelayTime))=qRespM(:,31:30+length(DelayTime))+(X(:,DelayTime).^2);
                cM(:,31:30+length(DelayTime))=cM(:,31:30+length(DelayTime))+1;
                bin=7;
                    RespM_binned(:,bin)=RespM_binned(:,bin)+mean(X(:,DelayTime),2);
                    qRespM_binned(:,bin)=qRespM_binned(:,bin)+(mean(X(:,DelayTime),2).^2);
                    cM_binned(:,bin)=cM_binned(:,bin)+1;
                    if ShFlag==1
                    for r=1:num_shuff
                        time_sh=ind_sh{r}(DelayTime);
                        RespM_binned_sh(:,bin,r)=RespM_binned_sh(:,bin,r)+mean(X(:,time_sh),2);
                        qRespM_binned_sh(:,bin,r)=qRespM_binned_sh(:,bin,r)+(mean(X(:,time_sh),2).^2);
                    end
                    end
       end       
    end
    
    RewardTime=RewardStart:RewardEnd;
    RewardTime=RewardTime(1:min(length(RewardTime),31));
    if length(RewardEnd)==1 
        if max(Speed(RewardTime))<SpeedTh && max(ActiveAnimal(RewardTime))==1
                RespM(:,35:34+length(RewardTime))=RespM(:,35:34+length(RewardTime))+X(:,RewardTime);
                qRespM(:,35:34+length(RewardTime))=qRespM(:,35:34+length(RewardTime))+(X(:,RewardTime).^2);
                cM(:,35:34+length(RewardTime))=cM(:,35:34+length(RewardTime))+1;
                for bin=8:13
                    RespM_binned(:,bin)=RespM_binned(:,bin)+mean(X(:,RewardStart+(bin-8)*5:min(RewardStart+(bin-7)*5,RewardEnd)),2);
                    qRespM_binned(:,bin)=qRespM_binned(:,bin)+(mean(X(:,RewardStart+(bin-8)*5:min(RewardStart+(bin-7)*5,RewardEnd)),2).^2);
                    cM_binned(:,bin)=cM_binned(:,bin)+1;
                    if ShFlag==1
                    for r=1:num_shuff
                        time_sh=ind_sh{r}(RewardTime);
                        RespM_binned_sh(:,bin,r)=RespM_binned_sh(:,bin,r)+mean(X(:,time_sh),2);
                        qRespM_binned_sh(:,bin,r)=qRespM_binned_sh(:,bin,r)+(mean(X(:,time_sh),2).^2);
                    end
                    end
                end
       end    
    
    end
end
progress='Miss done!'
for k=1:size(CRSE,1)
    floor(100*k/size(CRSE,1))
    TrialTime=CRSE(k,1):CRSE(k,2);

   
    StimStart=TrialTime(1)+1+StimDelay;
    %StimEnd=TrialTime(diff(CR(TrialTime))==-1);
    StimEnd=StimStart+29;
    
    RewardStart=TrialTime(diff(RewardWindow(TrialTime))==1)+1;
    RewardEnd=TrialTime(diff(RewardWindow(TrialTime))==-1);
    
    StimTime=StimStart:StimEnd;
    if length(StimEnd)==1 && CR(TrialTime(1)+1)==1 && max(Lick(StimTime))==0
       if max(Speed(StimTime))<SpeedTh && max(ActiveAnimal(StimTime))==1
                RespC(:,1:length(StimTime))=RespC(:,1:length(StimTime))+X(:,StimTime);
                qRespC(:,1:length(StimTime))=qRespC(:,1:length(StimTime))+(X(:,StimTime).^2);
                cC(:,1:length(StimTime))=cC(:,1:length(StimTime))+1;
                for bin=1:6
                    RespC_binned(:,bin)=RespC_binned(:,bin)+mean(X(:,StimStart+(bin-1)*5:min(StimStart+bin*5,StimEnd)),2);
                    qRespC_binned(:,bin)=qRespC_binned(:,bin)+(mean(X(:,StimStart+(bin-1)*5:min(StimStart+bin*5,StimEnd)),2).^2);
                    cC_binned(:,bin)=cC_binned(:,bin)+1;
                    if ShFlag==1
                    for r=1:num_shuff
                        time_sh=ind_sh{r}(StimStart+(bin-1)*5:min(StimStart+bin*5,StimEnd));
                        RespC_binned_sh(:,bin,r)=RespC_binned_sh(:,bin,r)+mean(X(:,time_sh),2);
                        qRespC_binned_sh(:,bin,r)=qRespC_binned_sh(:,bin,r)+(mean(X(:,time_sh),2).^2);
                    end
                    end
                end
       end       
    end
    
    StimEnd=TrialTime(diff(CR(TrialTime))==-1);
    DelayTime=StimEnd(end)+1:RewardStart-1;
    if max(Lick(DelayTime))==0
        if max(Speed(DelayTime))<SpeedTh && max(ActiveAnimal(DelayTime))==1
                RespC(:,31:30+length(DelayTime))=RespC(:,31:30+length(DelayTime))+X(:,DelayTime);
                qRespC(:,31:30+length(DelayTime))=qRespC(:,31:30+length(DelayTime))+(X(:,DelayTime).^2);
                cC(:,31:30+length(DelayTime))=cC(:,31:30+length(DelayTime))+1;
                bin=7;
                    RespC_binned(:,bin)=RespC_binned(:,bin)+mean(X(:,DelayTime),2);
                    qRespC_binned(:,bin)=qRespC_binned(:,bin)+(mean(X(:,DelayTime),2).^2);
                    cC_binned(:,bin)=cC_binned(:,bin)+1;
                    if ShFlag==1
                    for r=1:num_shuff
                        time_sh=ind_sh{r}(DelayTime);
                        RespC_binned_sh(:,bin,r)=RespC_binned_sh(:,bin,r)+mean(X(:,time_sh),2);
                        qRespC_binned_sh(:,bin,r)=qRespC_binned_sh(:,bin,r)+(mean(X(:,time_sh),2).^2);
                    end
                    end
                
       end       
    end
    
    RewardTime=RewardStart:RewardEnd;
    RewardTime=RewardTime(1:min(length(RewardTime),31));
    if length(RewardEnd)==1 
        if max(Speed(RewardTime))<SpeedTh && max(ActiveAnimal(RewardTime))==1
                RespC(:,35:34+length(RewardTime))=RespC(:,35:34+length(RewardTime))+X(:,RewardTime);
                qRespC(:,35:34+length(RewardTime))=qRespC(:,35:34+length(RewardTime))+(X(:,RewardTime).^2);
                cC(:,35:34+length(RewardTime))=cC(:,35:34+length(RewardTime))+1;
                for bin=8:13
                    RespC_binned(:,bin)=RespC_binned(:,bin)+mean(X(:,RewardStart+(bin-8)*5:min(RewardStart+(bin-7)*5,RewardEnd)),2);
                    qRespC_binned(:,bin)=qRespC_binned(:,bin)+(mean(X(:,RewardStart+(bin-8)*5:min(RewardStart+(bin-7)*5,RewardEnd)),2).^2);
                    cC_binned(:,bin)=cC_binned(:,bin)+1;
                    if ShFlag==1
                    for r=1:num_shuff
                        time_sh=ind_sh{r}(RewardTime);
                        RespC_binned_sh(:,bin,r)=RespC_binned_sh(:,bin,r)+mean(X(:,time_sh),2);
                        qRespC_binned_sh(:,bin,r)=qRespC_binned_sh(:,bin,r)+(mean(X(:,time_sh),2).^2);
                    end
                    end
                end
       end    
    
    end
end
progress='CR done!'
for k=1:size(FASE,1)
    floor(100*k/size(FASE,1))
    TrialTime=FASE(k,1):FASE(k,2);

   
    StimStart=TrialTime(1)+1+StimDelay;
    %StimEnd=TrialTime(diff(FA(TrialTime))==-1);
    StimEnd=StimStart+29;
    
    RewardStart=TrialTime(diff(RewardWindow(TrialTime))==1)+1;
    RewardEnd=TrialTime(diff(RewardWindow(TrialTime))==-1);
    
    StimTime=StimStart:StimEnd;
    if length(StimEnd)==1 && FA(TrialTime(1)+1)==1 && max(Lick(StimTime))==0
       if max(Speed(StimTime))<SpeedTh && max(ActiveAnimal(StimTime))==1
                RespF(:,1:length(StimTime))=RespF(:,1:length(StimTime))+X(:,StimTime);
                qRespF(:,1:length(StimTime))=qRespF(:,1:length(StimTime))+(X(:,StimTime).^2);
                cF(:,1:length(StimTime))=cF(:,1:length(StimTime))+1;
                for bin=1:6
                    RespF_binned(:,bin)=RespF_binned(:,bin)+mean(X(:,StimStart+(bin-1)*5:min(StimStart+bin*5,StimEnd)),2);
                    qRespF_binned(:,bin)=qRespF_binned(:,bin)+(mean(X(:,StimStart+(bin-1)*5:min(StimStart+bin*5,StimEnd)),2).^2);
                    cF_binned(:,bin)=cF_binned(:,bin)+1;
                    if ShFlag==1
                    for r=1:num_shuff
                        time_sh=ind_sh{r}(StimStart+(bin-1)*5:min(StimStart+bin*5,StimEnd));
                        RespF_binned_sh(:,bin,r)=RespF_binned_sh(:,bin,r)+mean(X(:,time_sh),2);
                        qRespF_binned_sh(:,bin,r)=qRespF_binned_sh(:,bin,r)+(mean(X(:,time_sh),2).^2);
                    end
                    end
                end
       end       
    end
    StimEnd=TrialTime(diff(FA(TrialTime))==-1);
    DelayTime=StimEnd(end)+1:RewardStart-1;
    if max(Lick(DelayTime))==0
        if max(Speed(DelayTime))<SpeedTh && max(ActiveAnimal(DelayTime))==1
                RespF(:,31:30+length(DelayTime))=RespF(:,31:30+length(DelayTime))+X(:,DelayTime);
                qRespF(:,31:30+length(DelayTime))=qRespF(:,31:30+length(DelayTime))+(X(:,DelayTime).^2);
                cF(:,31:30+length(DelayTime))=cF(:,31:30+length(DelayTime))+1;
                bin=7;
                    RespF_binned(:,bin)=RespF_binned(:,bin)+mean(X(:,DelayTime),2);
                    qRespF_binned(:,bin)=qRespF_binned(:,bin)+(mean(X(:,DelayTime),2).^2);
                    cF_binned(:,bin)=cF_binned(:,bin)+1;
                    if ShFlag==1
                    for r=1:num_shuff
                        time_sh=ind_sh{r}(DelayTime);
                        RespF_binned_sh(:,bin,r)=RespF_binned_sh(:,bin,r)+mean(X(:,time_sh),2);
                        qRespF_binned_sh(:,bin,r)=qRespF_binned_sh(:,bin,r)+(mean(X(:,time_sh),2).^2);
                    end
                    end
                
       end       
    end
    
    
    RewardTime=RewardStart:RewardEnd;
    RewardTime=RewardTime(1:min(length(RewardTime),31));
    if length(RewardEnd)==1 
        if max(Speed(RewardTime))<SpeedTh && max(ActiveAnimal(RewardTime))==1
                RespF(:,35:34+length(RewardTime))=RespF(:,35:34+length(RewardTime))+X(:,RewardTime);
                qRespF(:,35:34+length(RewardTime))=qRespF(:,35:34+length(RewardTime))+(X(:,RewardTime).^2);
                cF(:,35:34+length(RewardTime))=cF(:,35:34+length(RewardTime))+1;
                for bin=8:13
                    RespF_binned(:,bin)=RespF_binned(:,bin)+mean(X(:,RewardStart+(bin-8)*5:min(RewardStart+(bin-7)*5,RewardEnd)),2);
                    qRespF_binned(:,bin)=qRespF_binned(:,bin)+(mean(X(:,RewardStart+(bin-8)*5:min(RewardStart+(bin-7)*5,RewardEnd)),2).^2);
                    cF_binned(:,bin)=cF_binned(:,bin)+1;
                    if ShFlag==1
                    for r=1:num_shuff
                        time_sh=ind_sh{r}(RewardTime);
                        RespF_binned_sh(:,bin,r)=RespF_binned_sh(:,bin,r)+mean(X(:,time_sh),2);
                        qRespF_binned_sh(:,bin,r)=qRespF_binned_sh(:,bin,r)+(mean(X(:,time_sh),2).^2);
                    end
                    end
                end
       end    
    
    end
end
progress='FA done!'
AverageH = RespH ./ cH;
semH=((qRespH./(cH))- AverageH.^2)./sqrt(cH);

AverageM = RespM ./ cM;
semM=((qRespM./(cM))- AverageM.^2)./sqrt(cM);

AverageC = RespC ./ cC;
semC=((qRespC./(cC))- AverageC.^2)./sqrt(cC);

AverageF = RespF ./ cF;
semF=((qRespF./(cF))- AverageF.^2)./sqrt(cF);

AverageH_binned = RespH_binned ./ cH_binned;
semH_binned=((qRespH_binned./(cH_binned))- AverageH_binned.^2)./sqrt(cH_binned);
if ShFlag==1
for r=1:num_shuff
    AverageH_binned_sh(:,:,r) = RespH_binned_sh(:,:,r) ./ cH_binned;
    semH_binned_sh(:,:,r)=((qRespH_binned_sh(:,:,r)./(cH_binned))- AverageH_binned_sh(:,:,r).^2)./sqrt(cH_binned);
end
end
AverageM_binned = RespM_binned ./ cM_binned;
semM_binned=((qRespM_binned./(cM_binned))- AverageM_binned.^2)./sqrt(cM_binned);
if ShFlag==1
for r=1:num_shuff
    AverageM_binned_sh(:,:,r) = RespM_binned_sh(:,:,r) ./ cM_binned;
    semM_binned_sh(:,:,r)=((qRespM_binned_sh(:,:,r)./(cM_binned))- AverageM_binned_sh(:,:,r).^2)./sqrt(cM_binned);
end
end

AverageC_binned = RespC_binned ./ cC_binned;
semC_binned=((qRespC_binned./(cC_binned))- AverageC_binned.^2)./sqrt(cC_binned);
if ShFlag==1
for r=1:num_shuff
    AverageC_binned_sh(:,:,r) = RespC_binned_sh(:,:,r) ./ cC_binned;
    semC_binned_sh(:,:,r)=((qRespC_binned_sh(:,:,r)./(cC_binned))- AverageC_binned_sh(:,:,r).^2)./sqrt(cC_binned);
end
end

AverageF_binned = RespF_binned ./ cF_binned;
semF_binned=((qRespF_binned./(cF_binned))- AverageF_binned.^2)./sqrt(cF_binned);
if ShFlag==1
for r=1:num_shuff
    AverageF_binned_sh(:,:,r) = RespF_binned_sh(:,:,r) ./ cF_binned;
    semF_binned_sh(:,:,r)=((qRespF_binned_sh(:,:,r)./(cF_binned))- AverageF_binned_sh(:,:,r).^2)./sqrt(cF_binned);
end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Dprime 
%  Average_binned{1}=AverageH_binned;
% Average_binned_sh{1}=AverageH_binned_sh;
%  sem_binned{1}=semH_binned;
% sem_binned_sh{1}=semH_binned_sh;
%  c_binned{1}=cH_binned;
% 
%  Average_binned{2}=AverageM_binned;
% Average_binned_sh{2}=AverageM_binned_sh;
%  sem_binned{2}=semM_binned;
% sem_binned_sh{2}=semM_binned_sh;
%  c_binned{2}=cM_binned;
% 
%  Average_binned{3}=AverageC_binned;
% Average_binned_sh{3}=AverageC_binned_sh;
%  sem_binned{3}=semC_binned;
% sem_binned_sh{3}=semC_binned_sh;
%  c_binned{3}=cC_binned;
% 
%  Average_binned{4}=AverageF_binned;
% Average_binned_sh{4}=AverageF_binned_sh;
%  sem_binned{4}=semF_binned;
% sem_binned_sh{4}=semF_binned_sh;
%  c_binned{4}=cF_binned;
% 

% 
%  Note='HC,FM,HF,CM,HM,CF';
%  POA=[1,3;4,2;1,4;3,2;1,2;3,4];
%  cells_POA=zeros(cellCount,size(POA,1),3);
%  psh=zeros(cellCount,13,size(POA,1));
% timeBin{1}=[3:6];
% timeBin{2}=[7];
% timeBin{3}=[8:13];
%  pvalue=0.001;
% 
% 
%  Dprime=zeros(cellCount,13,size(POA,1));
%  for p=1:size(POA,1)
% 
% 
% Dprime_sh=zeros(cellCount,13,100);
% Fstat=zeros(cellCount,13);
% DegFstat=zeros(cellCount,13);
% cellType=zeros(cellCount,13);
% 
% for i= 1:cellCount
%     M0=Average_binned{POA(p,1)}(i,:);
%     C0=c_binned{POA(p,1)}(i,:);
%     V0=sem_binned{POA(p,1)}(i,:) .* sqrt(C0);
%     SS0=V0 .* C0;
%     
%     M1=Average_binned{POA(p,2)}(i,:);
%     C1=c_binned{POA(p,2)}(i,:);
%     V1=sem_binned{POA(p,2)}(i,:) .* sqrt(C1);
%     SS1=V1 .* C1;
%     
%     
%    for bin=1:13
%        Dprime(i,bin,p)=abs(M0(bin)-M1(bin)) / sqrt(0.5*(V0(bin)+V1(bin)));
%        for r=1:num_shuff
%             M0_sh=Average_binned_sh{POA(p,1)}(i,bin,r);
%             V0_sh=sem_binned_sh{POA(p,1)}(i,bin,r) .* sqrt(C0(bin));
%             
%             M1_sh=Average_binned_sh{POA(p,2)}(i,bin,r);
%             V1_sh=sem_binned_sh{POA(p,2)}(i,bin,r) .* sqrt(C1(bin));
%             Dprime_sh(i,bin,r)=abs(M0_sh-M1_sh) / sqrt(0.5*(V0_sh+V1_sh));
%        end
%        
%        
%        
% %        
% %        SSW=(SS0(bin)+SS1(bin)) / (C1(bin)+C0(bin)-2);
% %        MT=((C0(bin)*M0(bin))+(C1(bin)*M1(bin)))/(C0(bin)+C1(bin));
% %        SSB=(C0(bin)*((M0(bin)-MT)^2))  +  (C1(bin)*((M1(bin)-MT)^2));
% %     
% %        Fstat(i,bin)=SSB/SSW;
% %        DegFstat(i,bin)=(C1(bin)+C0(bin)-2); 
%        
%        if M0(bin)>M1(bin)
%            cellType(i,bin)=1;
%        else
%            cellType(i,bin)=2;
%        end
%        
%    end
% end
% 
% 
% for r=1:num_shuff
%     psh(:,:,p)=psh(:,:,p)+(Dprime_sh(:,:,r)>Dprime(:,:,p));
% end
% psh(:,:,p)=psh(:,:,p)/num_shuff;
% 
% for mode=1:3
% cells_POA(:,p,mode)=max(psh(:,timeBin{mode},p)<pvalue,[],2);
% 
% end
% 
%  end


%  save(strcat('C:\Users\Sadegh\Documents\VLMReborn\Reports2\2017_10_23_AverageActivity_CellType\Speed\Mouse',num2str(Mouse)),'POA','Note','cells_POA','psh','pvalue','Dprime','CortexArea')
%%%%%%%%%%%%%%%%%%%%Automatic Cell Detection
% CodingTable=zeros(8,8,8);
% clear DPAll
% DPAll{8,3,8}=[];
% meanDP=zeros(8,3,8);
% ErrorTable=zeros(8,8,8);
% for M=1:7
%     load(strcat('C:\Users\Sadegh\Documents\VLMReborn\Reports2\2017_10_23_AverageActivity_CellType\Coding\MouseDenoised',num2str(M)));
% 
%     S= cells_POA(:,1,1)==1;
%     D= cells_POA(:,1,2)==1;
%     R= cells_POA(:,1,3)==1;
%     
%     
%     for area=1:8
%         DPAll{area,1,M}= max(Dprime(S & CortexArea==area,3:6),[],2);
%         DPAll{area,2,M}= Dprime(D & CortexArea==area,7);
%         DPAll{area,3,M}= max(Dprime(R & CortexArea==area,8:13),[],2);
%         meanDP(area,1,M)=mean(DPAll{area,1,M});
%         meanDP(area,2,M)=mean(DPAll{area,2,M});
%         meanDP(area,3,M)=mean(DPAll{area,3,M});
%         DPAll{area,1,8}= [DPAll{area,1,8};max(Dprime(S & CortexArea==area,3:6),[],2)];
%         DPAll{area,2,8}= [DPAll{area,2,8};Dprime(D & CortexArea==area,7)];
%         DPAll{area,3,8}= [DPAll{area,3,8};max(Dprime(R & CortexArea==area,8:13),[],2)];
%         %%SDR
%         CodingTable(area,1,M)=sum(S & D & R & CortexArea==area);
%         %%SD
%         CodingTable(area,2,M)=sum(S & D & ~R & CortexArea==area);
%         %%SR
%         CodingTable(area,3,M)=sum(S & ~D & R & CortexArea==area);
%         %%S
%         CodingTable(area,4,M)=sum(S & ~D & ~R & CortexArea==area);
%         %%DR
%         CodingTable(area,5,M)=sum(~S & D & R & CortexArea==area);
%         %%D
%         CodingTable(area,6,M)=sum(~S & D & ~R & CortexArea==area);
%         %%R
%         CodingTable(area,7,M)=sum(~S & ~D & R & CortexArea==area);
%         %%cells
%         CodingTable(area,8,M)=sum( CortexArea==area);
%         
%         %%Stimuli HM
%         ErrorTable(area,1,M)=sum( cells_POA(:,5,1)==1 & cells_POA(:,6,1)==0 & CortexArea==area);
%         %%Stimuli FC
%         ErrorTable(area,2,M)=sum( cells_POA(:,5,1)==0 & cells_POA(:,6,1)==1 & CortexArea==area);
%         %%Stimuli HM-FC
%         ErrorTable(area,3,M)=sum( cells_POA(:,5,1)==1 & cells_POA(:,6,1)==1 & CortexArea==area);
%         %%Delay HM
%         ErrorTable(area,4,M)=sum( cells_POA(:,5,2)==1 & cells_POA(:,6,2)==0 & CortexArea==area);
%         %%Delay  FC
%         ErrorTable(area,5,M)=sum( cells_POA(:,5,2)==0 & cells_POA(:,6,2)==1 & CortexArea==area);
%         %%Delay  HM-FC
%         ErrorTable(area,6,M)=sum( cells_POA(:,5,2)==1 & cells_POA(:,6,2)==1 & CortexArea==area);
%         %%Reward  MC
%         ErrorTable(area,7,M)=sum( cells_POA(:,4,3)==1 & CortexArea==area);
%         %%Error Predicting
%         ErrorTable(area,8,M)=sum( (cells_POA(:,5,1)==1 | cells_POA(:,6,1)==1 | cells_POA(:,5,2)==1 | cells_POA(:,6,2)==1) &    CortexArea==area);
%         
% 
%     end
% 
% CodingTable(:,:,8)=CodingTable(:,:,8)+CodingTable(:,:,M);
% ErrorTable(:,:,8)=ErrorTable(:,:,8)+ErrorTable(:,:,M);    
% end
% for area=1:8
%     meanDP(area,1,8)=mean(DPAll{area,1,8});
%     meanDP(area,2,8)=mean(DPAll{area,2,8});
%     meanDP(area,3,8)=mean(DPAll{area,3,8});
% end
% 
% 
% dpTh=0.15;
% for i=1:1%cellCount
%     %if Dprime(i,3)>dpTh && Dprime(i,8)>dpTh && cellType(i,3)~=cellType(i,8)
%     %if cells_POA(i,6,1)==1
%         figure();
%             
%         %title(strcat('Mouse: ',num2str(Mouse),'  cell Number: ',num2str(i),' Area: ',areaNames(CortexArea(i))));
%         td=-0.9:0.1:5.5;
%         hold on;grid on
%         shadedErrorBar(td,AverageH(i,1:length(td)),semH(i,1:length(td)),'b',1);
%         shadedErrorBar(td,AverageC(i,1:length(td)),semC(i,1:length(td)),'k',1);
%         
%         %figure();imagesc(cellImage{i});
%          shadedErrorBar(td,AverageM(i,1:length(td)),semM(i,1:length(td)),'m',1);
%          shadedErrorBar(td,AverageF(i,1:length(td)),semF(i,1:length(td)),'r',1);
%     %end
% end


%%%%%%%%Manual Cell Detection

%%%right Hit , left CR , up both , down none 
%%right:MC>MH, down: not
load(strcat('C:\Users\Sadegh\Documents\VLMReborn\Reports2\2017_10_23_AverageActivity_CellType\SinosuidalCells\Session',num2str(Days(Mouse)),'_Raw.mat'));
%load(strcat('C:\Users\Sadegh\Documents\VLMReborn\Reports2\2017_10_23_AverageActivity_CellType\Coding\MouseDenoised',num2str(Mouse)));
%cellType=zeros(cellCount,1);
for cnum=1:cellCount
    cnum
    if cellType(cnum)==3
figHandle = figure(1);



clf(figHandle);
returnMap = containers.Map;
set(figHandle, 'KeyPressFcn', ...
    @(fig_obj , eventDat) readInput(fig_obj, eventDat,returnMap));

td=-0.9:0.1:5.5;
hold on;%grid on
title(strcat('Mouse: ',num2str(Mouse),'  cell Number: ',num2str(cnum),' Area: ',areaNames(CortexArea(cnum))));
shadedErrorBar(td,AverageH(cnum,1:length(td)),semH(cnum,1:length(td)),'b',1);
shadedErrorBar(td,AverageC(cnum,1:length(td)),semC(cnum,1:length(td)),'k',1);
% shadedErrorBar(td,AverageM(cnum,1:length(td)),semM(cnum,1:length(td)),'m',1);
% shadedErrorBar(td,AverageF(cnum,1:length(td)),semF(cnum,1:length(td)),'r',1);

waitfor(figHandle);

 newN= returnMap('newN')
cellType(cnum)=newN;
    end
end

% Note='1:Hit,2:CR,3:both,0:none';
% % Note='1: MC>MH';
%  save(strcat('C:\Users\Sadegh\Documents\VLMReborn\Reports2\2017_10_23_AverageActivity_CellType\SinosuidalCells\Session17_Raw'),'cellType','CortexArea','Note');
% 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Sinosuidal phase
% relativePhase=zeros(cellCount,2);
% ref1=1014;
% ref2=2363;
% sig=AverageH;
% ref=sig(ref1,10:30);
% ref=ref-mean(ref);
% figure();hold on;
% plot(ref,'r')
% for i=1:cellCount
% 
%     if cellType(i)==1 || cellType(i)==3
%         [CC,~,relativePhase(i,1)]=CrossCorrelation(sig(i,10:30)-mean(sig(i,10:30)),ref,3);
%         
%         plot(1-relativePhase(i,1):21-relativePhase(i,1),sig(i,10:30))
%     end   
% end
% 
% sig=AverageC;
% figure();hold on;
% ref=sig(ref2,10:30);
% ref=ref-mean(ref);
% plot(ref,'r')
% for i=1:cellCount
% 
%     if cellType(i)==2 || cellType(i)==3
%         [CC,~,relativePhase(i,2)]=CrossCorrelation(sig(i,10:30)-mean(sig(i,10:30)),ref,3);
%         
%         plot(1-relativePhase(i,2):21-relativePhase(i,2),sig(i,10:30))
%     end   
% end
% 
% save('C:\Users\Sadegh\Documents\VLMReborn\Reports2\2017_10_23_AverageActivity_CellType\SinosuidalCells\relPhase_Session17','relativePhase','ref1','ref2')
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Visualization
% cellVec=1:cellCount;
% 
% for cellNum=cellVec(cellType==2)
% 
% figure();
% title(cellNum)
% td=-0.9:0.1:5.5;
% hold on;grid on
% shadedErrorBar(td,AverageH(cellNum,1:length(td)),semH(cellNum,1:length(td)),'b',1);
% shadedErrorBar(td,AverageC(cellNum,1:length(td)),semC(cellNum,1:length(td)),'k',1);
% 
% 
% end
% 
% 
% 
cellMap=zeros(1017,1017)*5;
RGBarea=[0 0 3 ; 0 1 3 ; 0 1 3 ; 3 0 0 ; 0 3 0 ; 2 0 2 ;2 2 2; 3 2 0 ];
wh=5;
col=2;
for i=1:1800

    %if (cellType(i)==col || cellType(i)==3  ) 
    
       cellMap(max(1,cellIJ(i,1)-wh):min(1017,cellIJ(i,1)+wh),max(1,cellIJ(i,2)-wh):min(1017,cellIJ(i,2)+wh))=CortexArea(i);
        
    %end
    
end
figure()
imagesc(cellMap);







% 
