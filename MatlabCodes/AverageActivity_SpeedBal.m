clear all
SET0Name='C';
SET1Name='F';
StimDelay=-10;
ActiveTrialNumber=40;
MouseVec=1:7;
AddressSetup;
ShFlag=1;
num_repeat=1;
num_shuff=1000;
spr=num_shuff/num_repeat;

Days=[51:57];
areaNames={'V1','LV','MV','PPC','A','S','M','RSC','All'};
for Mouse=MouseVec
Day=Days(Mouse);

path=LoadPath{Day};
load(strcat(path,'\All_Sessions.mat'));
load(strcat(path,'\cellData.mat'));
load(strcat(path,'\cellData_ZS.mat'));
%load(strcat(path,'\areas.mat'));
Speed=Speed-min(Speed);
CRSE=CRSE(2:end,:);
HitSE=HitSE(2:end,:);
FASE=FASE(2:end,:);
MissSE=MissSE(2:end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch SET0Name
	case 'H'
        TrialType0=zeros(1,size(HitSE,1));
		SE0=HitSE;       
	case 'M'
        TrialType0=zeros(1,size(MissSE,1));
		SE0=MissSE;
	case 'C'
        TrialType0=zeros(1,size(CRSE,1));
		SE0=CRSE;
	case 'F'
        TrialType0=zeros(1,size(FASE,1));
		SE0=FASE;
	case 'G'
        TrialType0=[-1*ones(1,size(HitSE,1)),ones(1,size(MissSE,1))];
		SE0=[HitSE;MissSE];
	case 'N'
		TrialType0=[-1*ones(1,size(FASE,1)),ones(1,size(CRSE,1))];
		SE0=[FASE;CRSE];
	case 'L'
		TrialType0=[-1*ones(1,size(FASE,1)),ones(1,size(HitSE,1))];
		SE0=[FASE;HitSE];
	case 'NL'
		TrialType0=[-1*ones(1,size(MissSE,1)),ones(1,size(CRSE,1))];
		SE0=[MissSE;CRSE];
	case 'Cor'
		TrialType0=[-1*ones(1,size(CRSE,1)),ones(1,size(HitSE,1))];
		SE0=[CRSE;HitSE];
	case 'Err'
		TrialType0=[-1*ones(1,size(MissSE,1)),ones(1,size(FASE,1))];
		SE0=[MissSE;FASE];
end

switch SET1Name
	case 'H'
        TrialType1=zeros(1,size(HitSE,1));
		SE1=HitSE;
	case 'M'
        TrialType1=zeros(1,size(MissSE,1));
		SE1=MissSE;
	case 'C'
        TrialType1=zeros(1,size(CRSE,1));
		SE1=CRSE;
	case 'F'
        TrialType1=zeros(1,size(FASE,1));
		SE1=FASE;
	case 'G'
        TrialType1=[-1*ones(1,size(HitSE,1)),ones(1,size(MissSE,1))];
		SE1=[HitSE;MissSE];
	case 'N'
		TrialType1=[-1*ones(1,size(FASE,1)),ones(1,size(CRSE,1))];
		SE1=[FASE;CRSE];
	case 'L'
		TrialType1=[-1*ones(1,size(FASE,1)),ones(1,size(HitSE,1))];
		SE1=[FASE;HitSE];
	case 'NL'
		TrialType1=[-1*ones(1,size(MissSE,1)),ones(1,size(CRSE,1))];
		SE1=[MissSE;CRSE];
	case 'Cor'
		TrialType1=[-1*ones(1,size(CRSE,1)),ones(1,size(HitSE,1))];
		SE1=[CRSE;HitSE];
	case 'Err'
		TrialType1=[-1*ones(1,size(MissSE,1)),ones(1,size(FASE,1))];
		SE1=[MissSE;FASE];
end

Resp0=zeros(cellCount,65,num_repeat);
Resp1=zeros(cellCount,65,num_repeat);
qResp0=zeros(cellCount,65,num_repeat);
qResp1=zeros(cellCount,65,num_repeat);
c0=zeros(cellCount,65,num_repeat);
c1=zeros(cellCount,65,num_repeat);
Resp0_binned=zeros(cellCount,13,num_repeat);
Resp1_binned=zeros(cellCount,13,num_repeat);
qResp0_binned=zeros(cellCount,13,num_repeat);
qResp1_binned=zeros(cellCount,13,num_repeat);
c0_binned=zeros(cellCount,13,num_repeat);
c1_binned=zeros(cellCount,13,num_repeat);
Resp0_binned_sh=zeros(cellCount,13,num_shuff);
Resp1_binned_sh=zeros(cellCount,13,num_shuff);
qResp0_binned_sh=zeros(cellCount,13,num_shuff);
qResp1_binned_sh=zeros(cellCount,13,num_shuff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Stimuli=GoTrials+NogoTrials;
cellCount=size(cellData_Raw,1);
X=cellData_Raw_bin;

ind_sh={};
for r=1:num_shuff
    ind_sh{r}=randperm(size(X,2));
end


NBin=50;

ActiveAnimal=zeros(length(Lick),1);
for i=1:(length(Lick)-75*ActiveTrialNumber)
    if max(Lick(i:(i+75*ActiveTrialNumber)))==1
        ActiveAnimal(i:(i+75*ActiveTrialNumber))=1;
    end
end
Speed=double(round(Speed*10))/10;
[SQ,BQ]=SpeedQuantEqual(Speed,(Speed>0 & ActiveAnimal==1),NBin);

for rep=1:num_repeat

%%%%%%Select Clean trials for Stimuli window
nsse0=0;
clear SSE0
clear STrialType0
clear TW0
for i=1:size(SE0,1)     
    StimStart=SE0(i,1)+(1-Stimuli(SE0(i,1)));
    StimEnd=StimStart+19;
    
    if max(Lick(StimStart:StimEnd))==0 && Stimuli(StimStart)==1 &&  max(ActiveAnimal(SE0(i,1):SE0(i,2)))==1
        nsse0=nsse0+1;
        SSE0(nsse0,1:2)=SE0(i,:);
        STrialType0(1,nsse0)=TrialType0(i); 
        TW0{nsse0}=StimStart:StimEnd;
    end  
end

nsse1=0;
clear SSE1
clear STrialType1
clear TW1
for i=1:size(SE1,1)  
    StimStart=SE1(i,1)+(1-Stimuli(SE1(i,1)));
    StimEnd=StimStart+19;
    if max(Lick(StimStart:StimEnd))==0 && Stimuli(StimStart)==1 && max(ActiveAnimal(SE1(i,1):SE1(i,2)))==1
        nsse1=nsse1+1;
        SSE1(nsse1,1:2)=SE1(i,:);
        STrialType1(1,nsse1)=TrialType1(i);
        TW1{nsse1}=StimStart:StimEnd;
    end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%Balance Speed/type
Sorder0=randperm(nsse0);
Sorder1=randperm(nsse1);
UMask0=zeros(nsse0,1);
UMask1=zeros(nsse1,1);
BalSign0=0;
BalSign1=0;
NotFinished=1;

while (NotFinished==1)
   NotFinished=0; 
    for k0=Sorder0
        if (STrialType0(k0) * BalSign0) <=0 && UMask0(k0)==0
            matchable=0;
            Sbin0=max(SQ(TW0{k0}));
            for k1=Sorder1
                if UMask1(k1)==0
                    Sbin1=max(SQ(TW1{k1}));
                    if Sbin0==Sbin1
                       matchable=1;
                       if (STrialType1(k1) * BalSign1) <=0
                           UMask0(k0)=1;
                           UMask1(k1)=1;
                           BalSign0=BalSign0+STrialType0(k0);
                           BalSign1=BalSign1+STrialType1(k1);
                           NotFinished=1;
                           break;
                       end
                    end
                end
            end
            
            if matchable==0
                UMask0(k0)=2;
            end
        end
    end

end
SSE0=SSE0(UMask0==1,:);
SSE1=SSE1(UMask1==1,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeVector=[];
for i=1:size(SSE0,1)
    StimStart=SSE0(i,1)+(1-Stimuli(SSE0(i,1)))+StimDelay;
    StimEnd=StimStart+29;
    StimTime=StimStart:StimEnd;
    timeVector=[timeVector,StimTime];
end
for i=1:size(SSE1,1)
    StimStart=SSE1(i,1)+(1-Stimuli(SSE1(i,1)))+StimDelay;
    StimEnd=StimStart+29;
    StimTime=StimStart:StimEnd;
    timeVector=[timeVector,StimTime];
end
timeVector=timeVector(randperm(length(timeVector)));


for k=1:size(SSE0,1)
    floor(100*k/size(SSE0,1))
    TrialTime=SSE0(k,1):SSE0(k,2);
    
    StimStart=TrialTime(1)+(1-Stimuli(TrialTime(1)))+StimDelay;
    StimEnd=StimStart+29;
    
    StimTime=StimStart:StimEnd;
    
    Resp0(:,1:length(StimTime),rep)=Resp0(:,1:length(StimTime),rep)+X(:,StimTime);
    qResp0(:,1:length(StimTime),rep)=qResp0(:,1:length(StimTime),rep)+(X(:,StimTime).^2);
    c0(:,1:length(StimTime),rep)=c0(:,1:length(StimTime),rep)+1;
    for bin=1:6
        Resp0_binned(:,bin,rep)=Resp0_binned(:,bin,rep)+mean(X(:,StimStart+(bin-1)*5:min(StimStart+bin*5,StimEnd)),2);
        qResp0_binned(:,bin,rep)=qResp0_binned(:,bin,rep)+(mean(X(:,StimStart+(bin-1)*5:min(StimStart+bin*5,StimEnd)),2).^2);
        c0_binned(:,bin,rep)=c0_binned(:,bin,rep)+1;
        if ShFlag==1 
            for r=((rep-1)*spr+1):(rep*spr)
                %time_sh=ind_sh{r}(StimStart+(bin-1)*5:min(StimStart+bin*5,StimEnd));
                time_sh=timeVector(randperm(length(timeVector),5));
                Resp0_binned_sh(:,bin,r)=Resp0_binned_sh(:,bin,r)+mean(X(:,time_sh),2);
                qResp0_binned_sh(:,bin,r)=qResp0_binned_sh(:,bin,r)+(mean(X(:,time_sh),2).^2);
            end
        end
    end
    
end


for k=1:size(SSE1,1)
    floor(100*k/size(SSE1,1))
    TrialTime=SSE1(k,1):SSE1(k,2);
    
    StimStart=TrialTime(1)+(1-Stimuli(TrialTime(1)))+StimDelay;
    StimEnd=StimStart+29;
    
    StimTime=StimStart:StimEnd;
    
    Resp1(:,1:length(StimTime),rep)=Resp1(:,1:length(StimTime),rep)+X(:,StimTime);
    qResp1(:,1:length(StimTime),rep)=qResp1(:,1:length(StimTime),rep)+(X(:,StimTime).^2);
    c1(:,1:length(StimTime),rep)=c1(:,1:length(StimTime),rep)+1;
    for bin=1:6
        Resp1_binned(:,bin,rep)=Resp1_binned(:,bin,rep)+mean(X(:,StimStart+(bin-1)*5:min(StimStart+bin*5,StimEnd)),2);
        qResp1_binned(:,bin,rep)=qResp1_binned(:,bin,rep)+(mean(X(:,StimStart+(bin-1)*5:min(StimStart+bin*5,StimEnd)),2).^2);
        c1_binned(:,bin,rep)=c1_binned(:,bin,rep)+1;
        if ShFlag==1
            for r=((rep-1)*spr+1):(rep*spr)
                %time_sh=ind_sh{r}(StimStart+(bin-1)*5:min(StimStart+bin*5,StimEnd));
                time_sh=timeVector(randperm(length(timeVector),5));
                Resp1_binned_sh(:,bin,r)=Resp1_binned_sh(:,bin,r)+mean(X(:,time_sh),2);
                qResp1_binned_sh(:,bin,r)=qResp1_binned_sh(:,bin,r)+(mean(X(:,time_sh),2).^2);
            end
        end
    end
    
end

%%%%%%Select Clean trials for Delay window
nsse0=0;
clear SSE0
clear STrialType0
clear TW0
for i=1:size(SE0,1)
    TrialTime=SE0(i,1):SE0(i,2);
    StimEnd=max(TrialTime(Stimuli(TrialTime)==1));
    DelayTime=StimEnd(end)+1:StimEnd(end)+5;
    if  max(Lick(DelayTime))==0 &&  max(ActiveAnimal(SE0(i,1):SE0(i,2)))==1
        nsse0=nsse0+1;
        SSE0(nsse0,1:2)=SE0(i,:);
        STrialType0(1,nsse0)=TrialType0(i); 
        TW0{nsse0}=DelayTime;
    end  
end

nsse1=0;
clear SSE1
clear STrialType1
clear TW1
for i=1:size(SE1,1)     
    TrialTime=SE1(i,1):SE1(i,2);
    StimEnd=max(TrialTime(Stimuli(TrialTime)==1));
    DelayTime=StimEnd(end)+1:StimEnd(end)+5;
    if max(Lick(DelayTime))==0 && max(ActiveAnimal(SE1(i,1):SE1(i,2)))==1
        nsse1=nsse1+1;
        SSE1(nsse1,1:2)=SE1(i,:);
        STrialType1(1,nsse1)=TrialType1(i);
        TW1{nsse1}=DelayTime;
    end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%Balance Speed/type
Sorder0=randperm(nsse0);
Sorder1=randperm(nsse1);
UMask0=zeros(nsse0,1);
UMask1=zeros(nsse1,1);
BalSign0=0;
BalSign1=0;
NotFinished=1;

while (NotFinished==1)
   NotFinished=0; 
    for k0=Sorder0
        if (STrialType0(k0) * BalSign0) <=0 && UMask0(k0)==0
            matchable=0;
            Sbin0=max(SQ(TW0{k0}));
            for k1=Sorder1
                if UMask1(k1)==0
                    Sbin1=max(SQ(TW1{k1}));
                    if Sbin0==Sbin1
                       matchable=1;
                       if (STrialType1(k1) * BalSign1) <=0
                           UMask0(k0)=1;
                           UMask1(k1)=1;
                           BalSign0=BalSign0+STrialType0(k0);
                           BalSign1=BalSign1+STrialType1(k1);
                           NotFinished=1;
                           break;
                       end
                    end
                end
            end
            
            if matchable==0
                UMask0(k0)=2;
            end
        end
    end

end
SSE0=SSE0(UMask0==1,:);
SSE1=SSE1(UMask1==1,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeVector=[];
for i=1:size(SSE0,1)
    TrialTime=SSE0(i,1):SSE0(i,2);
    StimEnd=max(TrialTime(Stimuli(TrialTime)==1));
    DelayTime=StimEnd(end)+1:StimEnd(end)+5;
    timeVector=[timeVector,DelayTime];
end
for i=1:size(SSE1,1)
    TrialTime=SSE1(i,1):SSE1(i,2);
    StimEnd=max(TrialTime(Stimuli(TrialTime)==1));
    DelayTime=StimEnd(end)+1:StimEnd(end)+5;
    timeVector=[timeVector,DelayTime];
end
timeVector=timeVector(randperm(length(timeVector)));


for k=1:size(SSE0,1)
    floor(100*k/size(SSE0,1))
    
    TrialTime=SSE0(k,1):SSE0(k,2);
    StimEnd=max(TrialTime(Stimuli(TrialTime)==1));
    DelayTime=StimEnd(end)+1:StimEnd(end)+5;
    
    Resp0(:,31:30+length(DelayTime),rep)=Resp0(:,31:30+length(DelayTime),rep)+X(:,DelayTime);
    qResp0(:,31:30+length(DelayTime),rep)=qResp0(:,31:30+length(DelayTime),rep)+(X(:,DelayTime).^2);
    c0(:,31:30+length(DelayTime),rep)=c0(:,31:30+length(DelayTime),rep)+1;
    bin=7;
    Resp0_binned(:,bin,rep)=Resp0_binned(:,bin,rep)+mean(X(:,DelayTime),2);
    qResp0_binned(:,bin,rep)=qResp0_binned(:,bin,rep)+(mean(X(:,DelayTime),2).^2);
    c0_binned(:,bin,rep)=c0_binned(:,bin,rep)+1;
    if ShFlag==1 
        for r=((rep-1)*spr+1):(rep*spr)
            %time_sh=ind_sh{r}(DelayTime);
            time_sh=timeVector(randperm(length(timeVector),5));
            Resp0_binned_sh(:,bin,r)=Resp0_binned_sh(:,bin,r)+mean(X(:,time_sh),2);
            qResp0_binned_sh(:,bin,r)=qResp0_binned_sh(:,bin,r)+(mean(X(:,time_sh),2).^2);
        end
    end
    
end


for k=1:size(SSE1,1)
    floor(100*k/size(SSE1,1))
    
    TrialTime=SSE1(k,1):SSE1(k,2);
    StimEnd=max(TrialTime(Stimuli(TrialTime)==1));
    DelayTime=StimEnd(end)+1:StimEnd(end)+5;
    
    Resp1(:,31:30+length(DelayTime),rep)=Resp1(:,31:30+length(DelayTime),rep)+X(:,DelayTime);
    qResp1(:,31:30+length(DelayTime),rep)=qResp1(:,31:30+length(DelayTime),rep)+(X(:,DelayTime).^2);
    c1(:,31:30+length(DelayTime),rep)=c1(:,31:30+length(DelayTime),rep)+1;
    bin=7;
    Resp1_binned(:,bin,rep)=Resp1_binned(:,bin,rep)+mean(X(:,DelayTime),2);
    qResp1_binned(:,bin,rep)=qResp1_binned(:,bin,rep)+(mean(X(:,DelayTime),2).^2);
    c1_binned(:,bin,rep)=c1_binned(:,bin,rep)+1;
    if ShFlag==1 
        for r=((rep-1)*spr+1):(rep*spr)
            %time_sh=ind_sh{r}(DelayTime);
            time_sh=timeVector(randperm(length(timeVector),5));
            Resp1_binned_sh(:,bin,r)=Resp1_binned_sh(:,bin,r)+mean(X(:,time_sh),2);
            qResp1_binned_sh(:,bin,r)=qResp1_binned_sh(:,bin,r)+(mean(X(:,time_sh),2).^2);
        end
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Select Clean trials for Rward window
nsse0=0;
clear SSE0
clear STrialType0
clear TW0
for i=1:size(SE0,1)
    TrialTime=SE0(i,1):SE0(i,2);
    RewardStart=min(TrialTime(RewardWindow(TrialTime)==1));
    RewardTime=RewardStart:RewardStart+29;
    if length(RewardStart)==1 &&  max(ActiveAnimal(SE0(i,1):SE0(i,2)))==1
        nsse0=nsse0+1;
        SSE0(nsse0,1:2)=SE0(i,:);
        STrialType0(1,nsse0)=TrialType0(i); 
        TW0{nsse0}=RewardTime;
    end  
end

nsse1=0;
clear SSE1
clear STrialType1
clear TW1
for i=1:size(SE1,1)     
    TrialTime=SE1(i,1):SE1(i,2);
    RewardStart=min(TrialTime(RewardWindow(TrialTime)==1));
    RewardTime=RewardStart:RewardStart+29;
    if length(RewardStart)==1 && max(ActiveAnimal(SE1(i,1):SE1(i,2)))==1
        nsse1=nsse1+1;
        SSE1(nsse1,1:2)=SE1(i,:);
        STrialType1(1,nsse1)=TrialType1(i);
        TW1{nsse1}=RewardTime;
    end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%Balance Speed/type
Sorder0=randperm(nsse0);
Sorder1=randperm(nsse1);
UMask0=zeros(nsse0,1);
UMask1=zeros(nsse1,1);
BalSign0=0;
BalSign1=0;
NotFinished=1;

while (NotFinished==1)
   NotFinished=0; 
    for k0=Sorder0
        if (STrialType0(k0) * BalSign0) <=0 && UMask0(k0)==0
            matchable=0;
            Sbin0=max(SQ(TW0{k0}));
            for k1=Sorder1
                if UMask1(k1)==0
                    Sbin1=max(SQ(TW1{k1}));
                    if Sbin0==Sbin1
                       matchable=1;
                       if (STrialType1(k1) * BalSign1) <=0
                           UMask0(k0)=1;
                           UMask1(k1)=1;
                           BalSign0=BalSign0+STrialType0(k0);
                           BalSign1=BalSign1+STrialType1(k1);
                           NotFinished=1;
                           break;
                       end
                    end
                end
            end
            
            if matchable==0
                UMask0(k0)=2;
            end
        end
    end

end
SSE0=SSE0(UMask0==1,:);
SSE1=SSE1(UMask1==1,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeVector=[];
for i=1:size(SSE0,1)
    TrialTime=SSE0(i,1):SSE0(i,2);
    RewardStart=min(TrialTime(RewardWindow(TrialTime)==1));
    RewardTime=RewardStart:RewardStart+29;
    timeVector=[timeVector,RewardTime];
end
for i=1:size(SSE1,1)
    TrialTime=SSE1(i,1):SSE1(i,2);
    RewardStart=min(TrialTime(RewardWindow(TrialTime)==1));
    RewardTime=RewardStart:RewardStart+29;
    timeVector=[timeVector,RewardTime];
end
timeVector=timeVector(randperm(length(timeVector)));


for k=1:size(SSE0,1)
    floor(100*k/size(SSE0,1))
    TrialTime=SSE0(k,1):SSE0(k,2);
    RewardStart=min(TrialTime(RewardWindow(TrialTime)==1));
    RewardTime=RewardStart:RewardStart+29;
    
    Resp0(:,36:35+length(RewardTime),rep)=Resp0(:,36:35+length(RewardTime),rep)+X(:,RewardTime);
    qResp0(:,36:35+length(RewardTime),rep)=qResp0(:,36:35+length(RewardTime),rep)+(X(:,RewardTime).^2);
    c0(:,36:35+length(RewardTime),rep)=c0(:,36:35+length(RewardTime),rep)+1;
    for bin=8:13
        Resp0_binned(:,bin,rep)=Resp0_binned(:,bin,rep)+mean(X(:,RewardStart+(bin-8)*5:min(RewardStart+(bin-7)*5,RewardStart+30)),2);
        qResp0_binned(:,bin,rep)=qResp0_binned(:,bin,rep)+(mean(X(:,RewardStart+(bin-8)*5:min(RewardStart+(bin-7)*5,RewardStart+30)),2).^2);
        c0_binned(:,bin,rep)=c0_binned(:,bin,rep)+1;
        if ShFlag==1 
            for r=((rep-1)*spr+1):(rep*spr)
                %time_sh=ind_sh{r}(RewardTime);
                time_sh=timeVector(randperm(length(timeVector),5));
                Resp0_binned_sh(:,bin,r)=Resp0_binned_sh(:,bin,r)+mean(X(:,time_sh),2);
                qResp0_binned_sh(:,bin,r)=qResp0_binned_sh(:,bin,r)+(mean(X(:,time_sh),2).^2);
            end
        end
        
        
    end
end


for k=1:size(SSE1,1)
    floor(100*k/size(SSE1,1))
    TrialTime=SSE1(k,1):SSE1(k,2);
    RewardStart=min(TrialTime(RewardWindow(TrialTime)==1));
    RewardTime=RewardStart:RewardStart+29;
    
    Resp1(:,36:35+length(RewardTime),rep)=Resp1(:,36:35+length(RewardTime),rep)+X(:,RewardTime);
    qResp1(:,36:35+length(RewardTime),rep)=qResp1(:,36:35+length(RewardTime),rep)+(X(:,RewardTime).^2);
    c1(:,36:35+length(RewardTime),rep)=c1(:,36:35+length(RewardTime),rep)+1;
    for bin=8:13
        Resp1_binned(:,bin,rep)=Resp1_binned(:,bin,rep)+mean(X(:,RewardStart+(bin-8)*5:min(RewardStart+(bin-7)*5,RewardStart+30)),2);
        qResp1_binned(:,bin,rep)=qResp1_binned(:,bin,rep)+(mean(X(:,RewardStart+(bin-8)*5:min(RewardStart+(bin-7)*5,RewardStart+30)),2).^2);
        c1_binned(:,bin,rep)=c1_binned(:,bin,rep)+1;
        if ShFlag==1 
            for r=((rep-1)*spr+1):(rep*spr)
                %time_sh=ind_sh{r}(RewardTime);
                time_sh=timeVector(randperm(length(timeVector),5));
                Resp1_binned_sh(:,bin,r)=Resp1_binned_sh(:,bin,r)+mean(X(:,time_sh),2);
                qResp1_binned_sh(:,bin,r)=qResp1_binned_sh(:,bin,r)+(mean(X(:,time_sh),2).^2);
            end
        end
        
        
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Average0 = mean(Resp0 ./ c0,3);
sem0=(mean(qResp0./(c0),3)- Average0.^2)./sqrt(mean(c0,3));

Average1 = mean(Resp1 ./ c1,3);
sem1=(mean(qResp1./(c1),3)- Average1.^2)./sqrt(mean(c1,3));



Average0_binned = mean(Resp0_binned ./ c0_binned,3);
sem0_binned=(mean(qResp0_binned./(c0_binned),3)- Average0_binned.^2)./sqrt(mean(c0_binned,3));
Average0_binned_sh=zeros(size(Resp0_binned_sh));
sem0_binned_sh=zeros(size(Resp0_binned_sh));
if ShFlag==1
for r=1:num_shuff
    Average0_binned_sh(:,:,r) = Resp0_binned_sh(:,:,r) ./ c0_binned(:,:,1);
    sem0_binned_sh(:,:,r)=((qResp0_binned_sh(:,:,r)./(c0_binned(:,:,1)))- Average0_binned_sh(:,:,r).^2)./sqrt(c0_binned(:,:,1));
end
end


Average1_binned = mean(Resp1_binned ./ c1_binned,3);
sem1_binned=(mean(qResp1_binned./(c1_binned),3)- Average1_binned.^2)./sqrt(mean(c1_binned,3));
Average1_binned_sh=zeros(size(Resp0_binned_sh));
sem1_binned_sh=zeros(size(Resp0_binned_sh));
if ShFlag==1
for r=1:num_shuff
    Average1_binned_sh(:,:,r) = Resp1_binned_sh(:,:,r) ./ c1_binned(:,:,1);
    sem1_binned_sh(:,:,r)=((qResp1_binned_sh(:,:,r)./(c1_binned(:,:,1)))- Average1_binned_sh(:,:,r).^2)./sqrt(c1_binned(:,:,1));
end
end
save(strcat('E:\Reports2\2018_03_09_AverageActivity\RawBin\FACR\ConcatMouse',num2str(Mouse)),'cellCount','num_shuff',...
    'Average0','sem0','Resp0','qResp0','c0','Average0_binned','sem0_binned','Resp0_binned',...
    'qResp0_binned','c0_binned','Average0_binned_sh','sem0_binned_sh','Resp0_binned_sh','qResp0_binned_sh',...
    'Average1','sem1','Resp1','qResp1','c1','Average1_binned','sem1_binned','Resp1_binned',...
    'qResp1_binned','c1_binned','Average1_binned_sh','sem1_binned_sh','Resp1_binned_sh','qResp1_binned_sh')

end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Dprime 
% Mouse=1;
% load(strcat('E:\Reports2\2018_03_09_AverageActivity\RawBin\HitCR\ConcatMouse',num2str(M),'.mat'));
% load(strcat(LoadPath{Mouse+50},'\cellData.mat'));
% cells_POA=zeros(cellCount,3);
% psh=zeros(cellCount,13);
% timeBin{1}=[3:6];
% timeBin{2}=[7];
% timeBin{3}=[8:13];
% pvalue=0.001;
% 
% 
% Dprime=zeros(cellCount,13);
% 
% Dprime_sh=zeros(cellCount,13,100);
% Fstat=zeros(cellCount,13);
% DegFstat=zeros(cellCount,13);
% cellType=zeros(cellCount,13);
% 
% for i= 1:cellCount
%     M0=Average0_binned(i,:);
%     C0=mean(c0_binned(i,:,:),3);
%     V0=sem0_binned(i,:) .* sqrt(C0);
% 
%     
%     M1=Average1_binned(i,:);
%     C1=mean(c1_binned(i,:,:),3);
%     V1=sem1_binned(i,:) .* sqrt(C1);
%     
%     
%     for bin=1:13
%         Dprime(i,bin)=abs(M0(bin)-M1(bin)) / sqrt(0.5*(V0(bin)+V1(bin)));
%         for r=1:num_shuff
%             M0_sh=Average0_binned_sh(i,bin,r);
%             V0_sh=sem0_binned_sh(i,bin,r) .* sqrt(C0(bin));
%             
%             M1_sh=Average1_binned_sh(i,bin,r);
%             V1_sh=sem1_binned_sh(i,bin,r) .* sqrt(C1(bin));
%             Dprime_sh(i,bin,r)=abs(M0_sh-M1_sh) / sqrt(0.5*(V0_sh+V1_sh));
%         end
%         
%         
%         if M0(bin)>M1(bin)
%             cellType(i,bin)=1;
%         else
%             cellType(i,bin)=2;
%         end
%         
%     end
% end
% 
% 
% for r=1:num_shuff
%     psh(:,:)=psh(:,:)+(Dprime_sh(:,:,r)>Dprime(:,:));
% end
% psh(:,:)=psh(:,:)/num_shuff;
% 
% for mode=1:3
%     cells_POA(:,mode)=max(psh(:,timeBin{mode})<pvalue,[],2);
%     
% end
%     
% 
% save(strcat('E:\Reports2\2018_03_09_AverageActivity\RawBin\HitCR\DPMouse',num2str(Mouse)),'POA','Note','cells_POA','psh','pvalue','Dprime','CortexArea')
%%%%%%%%%%%%%%%%%%%%Automatic Cell Detection
% 
% % 
% % % dpTh=0.15;
% for i=1:cellCount
%     %if Dprime(i,3)>dpTh && Dprime(i,8)>dpTh && cellType(i,3)~=cellType(i,8)
%     if psh(i,3)<0.001
%         figure();
%             
%         %title(strcat('Mouse: ',num2str(Mouse),'  cell Number: ',num2str(i),' Area: ',areaNames(CortexArea(i))));
%         td=-0.9:0.1:5.5;
%         hold on;grid on
%         shadedErrorBar(td,Average1(i,1:length(td)),sem1(i,1:length(td)),'b',1);
%         shadedErrorBar(td,Average0(i,1:length(td)),sem0(i,1:length(td)),'k',1);
% 
%     end
% end


%%%%%%%%Manual Cell Detection

%%%right Hit , left CR , up both , down none 
%%right:MC>MH, down: not
% load(strcat('C:\Users\Sadegh\Documents\VLMReborn\Reports2\2017_10_23_AverageActivity_CellType\SinosuidalCells\Session',num2str(Days(Mouse)),'_Raw.mat'));
% %load(strcat('C:\Users\Sadegh\Documents\VLMReborn\Reports2\2017_10_23_AverageActivity_CellType\Coding\MouseDenoised',num2str(Mouse)));
% %cellType=zeros(cellCount,1);
% for cnum=1:cellCount
%     cnum
%     if cellType(cnum)==3
% figHandle = figure(1);
% 
% 
% 
% clf(figHandle);
% returnMap = containers.Map;
% set(figHandle, 'KeyPressFcn', ...
%     @(fig_obj , eventDat) readInput(fig_obj, eventDat,returnMap));
% 
% td=-0.9:0.1:5.5;
% hold on;%grid on
% title(strcat('Mouse: ',num2str(Mouse),'  cell Number: ',num2str(cnum),' Area: ',areaNames(CortexArea(cnum))));
% shadedErrorBar(td,AverageH(cnum,1:length(td)),semH(cnum,1:length(td)),'b',1);
% shadedErrorBar(td,AverageC(cnum,1:length(td)),semC(cnum,1:length(td)),'k',1);
% % shadedErrorBar(td,AverageM(cnum,1:length(td)),semM(cnum,1:length(td)),'m',1);
% % shadedErrorBar(td,AverageF(cnum,1:length(td)),semF(cnum,1:length(td)),'r',1);
% 
% waitfor(figHandle);
% 
%  newN= returnMap('newN')
% cellType(cnum)=newN;
%     end
% end

% Note='1:Hit,2:CR,3:both,0:none';
% % Note='1: MC>MH';
%  save(strcat('C:\Users\Sadegh\Documents\VLMReborn\Reports2\2017_10_23_AverageActivity_CellType\SinosuidalCells\Session17_Raw'),'cellType','CortexArea','Note');
% 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Sinosuidal phase


% 
