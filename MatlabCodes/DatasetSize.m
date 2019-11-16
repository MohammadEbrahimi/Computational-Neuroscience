%function [TrainingSize]=DatasetSize(loadpath,SET0Name,SET1Name,mode)
 loadpath='C:\Users\Sadegh\Documents\VLMReborn\L364\Data\allDays';
 mode=1;
 SET0Name='C';
 SET1Name='M';
warning off;
load(strcat(loadpath,'\All_Sessions.mat'));
progress='Data is loaded!'
%%% creating data sets
%%%%%%%%% Balance Trials
CRSE=CRSE(2:end,:);
FASE=FASE(2:end,:);
HitSE=HitSE(2:end,:);
MissSE=MissSE(2:end,:);
 nse0=min(size(CRSE,1),size(FASE,1));
 nse1=min(size(HitSE,1),size(MissSE,1));
 BSE0=zeros(nse0*2,2);
 BSE1=zeros(nse1*2,2);
 for i=1:nse0
     BSE0(2*i -1,:)=CRSE(i,:);
     BSE0(2*i   ,:)=FASE(i,:);
 end
 for i=1:nse1
     BSE1(2*i -1,:)=HitSE(i,:);
     BSE1(2*i   ,:)=MissSE(i,:);
 end

 nse0=min(size(CRSE,1),size(MissSE,1));
 nse1=min(size(HitSE,1),size(FASE,1));
 LSE0=zeros(nse0*2,2);
 LSE1=zeros(nse1*2,2);
 for i=1:nse0
     LSE0(2*i -1,:)=CRSE(i,:);
     LSE0(2*i   ,:)=MissSE(i,:);
 end
 for i=1:nse1
     LSE1(2*i -1,:)=HitSE(i,:);
     LSE1(2*i   ,:)=FASE(i,:);
 end
%%%%
switch SET0Name
        case 'H'
                SE0=HitSE;
        case 'M'
                SE0=MissSE;
        case 'C'
                SE0=CRSE;
        case 'F'
                SE0=FASE;
        case 'G'
                SE0=BSE1;
        case 'N'
                SE0=BSE0;
        case 'L'
                SE0=LSE1;
        case 'NL'
                SE0=LSE0;
        case 'SG'
                SE0=[HitSE;MissSE];
        case 'SN'
                SE0=[CRSE;FASE];

end

switch SET1Name
        case 'H'
                SE1=HitSE;
        case 'M'
                SE1=MissSE;
        case 'C'
                SE1=CRSE;
        case 'F'
                SE1=FASE;
        case 'G'
                SE1=BSE1;
        case 'N'
                SE1=BSE0;
        case 'L'
                SE1=LSE1;
        case 'NL'
                SE1=LSE0;
        case 'SG'
                SE1=[HitSE;MissSE];
        case 'SN'
                SE1=[CRSE;FASE];

end


SpeedTh=500;
ActiveTrialNumber=20;



cellEvents=cellData;


ActiveAnimal=ones(length(Lick),1);
for i=1:(length(Lick)-75*ActiveTrialNumber)
    if max(Lick(i:(i+75*ActiveTrialNumber)))==0
        ActiveAnimal(i:(i+75*ActiveTrialNumber))=0;
    end
end

ActiveAnimal=1-ActiveAnimal;
%%%%%%%%%%%%%%%%%%%%%%
Delay(Lick==1)=0;
if mode==1
    group0=NogoTrials+GoTrials;
    group1=NogoTrials+GoTrials;
    LMax=20;
    tdVec=-5:19;
    trainTime=0:19;
elseif mode==2
    group0=Delay;
    group1=Delay;
    LMax=5;
    tdVec=-5:4;
    trainTime=0:4;
elseif mode==3
    group0=RewardWindow;
    group1=RewardWindow;
    LMax=30;
    tdVec=-5:29;
    trainTime=0:29;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%Start Learning

shuff0=randperm(size(SE0,1));
shuff1=randperm(size(SE1,1));

c0=1;
c1=1;
dd0=zeros(size(cellEvents,1),sum(group0));
for si=1:length(shuff0)
    
    k=shuff0(si);
    if(max(Speed(SE0(k,1):(SE0(k,1)+25)))<=SpeedTh && max(ActiveAnimal(SE0(k,1):SE0(k,2)))==1 && max(Lick(SE0(k,1):(SE0(k,1)+25)))==0)
        j=0;
        while group0(SE0(k,1)+j)==0 && ((SE0(k,1)+j) < SE0(k,2))
            j=j+1;
        end
        stimL=0;
        for i=0:LMax
            if group0(SE0(k,1)+j+i)==0
                stimL=i;
                break;
            end
        end
        stimL=stimL-1;
        if ((SE0(k,1)+j+i) <= SE0(k,2)+1  && stimL>min(trainTime))
            start=(SE0(k,1)+j+min(trainTime));
            finish=(SE0(k,1)+j+min(stimL,max(trainTime)));
            dd0(:,c0:c0+(finish-start))=cellEvents(:, start : finish );
            shuffInd0(c0:c0+(finish-start))=si;
            c0=c0+(finish-start)+1;
        end
    end
end
c0=c0-1
dd0=dd0(:,1:c0);


dd1=zeros(size(cellEvents,1),sum(group1));
for si=1:length(shuff1)
    
    k=shuff1(si);
    if(max(Speed(SE1(k,1):(SE1(k,1)+25)))<=SpeedTh && max(ActiveAnimal(SE1(k,1):SE1(k,2)))==1 && max(Lick(SE1(k,1):(SE1(k,1)+25)))==0)
        j=0;
        while group1(SE1(k,1)+j)==0 && ((SE1(k,1)+j) < SE1(k,2))
            j=j+1;
        end
        stimL=0;
        for i=0:LMax
            if group1(SE1(k,1)+j+i)==0
                stimL=i;
                break;
            end
        end
        stimL=stimL-1;
        if ((SE1(k,1)+j+i) <= SE1(k,2)+1  && stimL>min(trainTime))
            start=(SE1(k,1)+j+min(trainTime));
            finish=(SE1(k,1)+j+min(stimL,max(trainTime)));
            dd1(:,c1:c1+(finish-start))=cellEvents(:, start : finish );
            shuffInd1(c1:c1+(finish-start))=si;
            c1=c1+(finish-start)+1;
        end
    end
end
c1=c1-1
dd1=dd1(:,1:c1);


c=floor(min(c0,c1)/2);




TrainingSize=c;


%end




