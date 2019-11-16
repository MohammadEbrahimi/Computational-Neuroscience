global StorageArray 
cellObjectNumber=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];
X_tile=[1 1 1 1 246 246 246 246 491 491 491 491 736 736 736 736];
Y_tile=[1 246 491 736 1 246 491 736 1 246 491 736 1 246 491 736];

savePath='C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_8_2015';

timeScale=0.5;



for tile=1:length(cellObjectNumber)
   tileCount=length(StorageArray{cellObjectNumber(tile)}.Children);
   tileMask{tile}=ones(1,tileCount);
end

 cellMask= tileMask{1};
 for i=2:16
         cellMask = [cellMask,tileMask{i}];
 end
cellCount= sum(cellMask);

IC_time=StorageArray{1}.Children{1}.Object.Trace.XVector * timeScale;

Ns=length(IC_time);
cellData=zeros(cellCount,Ns);
TileTopLeft=zeros(cellCount,2);

clear cellImage
n=1;
 for tile=1:length(cellObjectNumber)
     Mask=tileMask{tile};
     for i=1:length(StorageArray{cellObjectNumber(tile)}.Children)
         if Mask(i)==1
        cellData(n,:)=StorageArray{cellObjectNumber(tile)}.Children{i}.Object.Trace.Data;
        cellImage{n}=StorageArray{cellObjectNumber(tile)}.Children{i}.Object.Image.Data;
        TileTopLeft(n,:)=[X_tile(tile) Y_tile(tile)];
        n=n+1;
         end
     end
 end
 [eventTimes, peakTimes, eventAmps, eventBin, cellTraceSigmas] = detectEvents(cellData);


     

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%  Ns=11998+ds(2);
%  trialObjectNumber=18;
%  
%  
%     Angle_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 2}.Object.Data;
%     Lick_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 3}.Object.Data;
%     RW_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 4}.Object.Data;
%     GoTrials_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 5}.Object.Data;
%     NoGoTrials_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 6}.Object.Data;
%     AP_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 7}.Object.Data;
%     WR_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 8}.Object.Data;
%     XSpeed_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 9}.Object.Data;
%     YSpeed_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 10}.Object.Data;
%     TO_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 13}.Object.Data;
%     UA_time=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 2}.Object.XVector;
% 
% n=1;
% Angle=zeros(Ns,1);
% pre=1;
% for i=1:length(UA_time)
%     if(UA_time(i)>IC_time(n))
% 
%         Angle(n)=min(Angle_AU(pre:i-1));
%         pre=i;
%         n=n+1;
%         if(n>Ns)
%             break;
%         end
%     end
% 
% end
% 
% n=1;
% Lick=zeros(Ns,1);
% pre=1;
% for i=1:length(UA_time)
%     if(UA_time(i)>IC_time(n))
% 
%         Lick(n)=max(Lick_AU(max(pre-12,1):i-1));
%         pre=i;
%         n=n+1;
%         if(n>Ns)
%             break;
%         end
%     end
% 
% end
% 
% 
% n=1;
% RewardWindow=zeros(Ns,1);
% pre=1;
% for i=1:length(UA_time)
%     if(UA_time(i)>IC_time(n))
%         RewardWindow(n)=round(mean(RW_AU(pre:i-1)));
%         pre=i;
%         n=n+1;
%         if(n>Ns)
%             break;
%         end
%     end
% end
% 
% n=1;
% GoTrials=zeros(Ns,1);
% pre=1;
% for i=1:length(UA_time)
%     if(UA_time(i)>IC_time(n))
%         GoTrials(n)=round(mean(GoTrials_AU(pre:i-1)));
%         pre=i;
%         n=n+1;
%         if(n>Ns)
%             break;
%         end
%     end
% end
%  
% n=1;
% NogoTrials=zeros(Ns,1);
% pre=1;
% for i=1:length(UA_time)
%     if(UA_time(i)>IC_time(n))
%         NogoTrials(n)=round(mean(NoGoTrials_AU(pre:i-1)));
%         pre=i;
%         n=n+1;
%         if(n>Ns)
%             break;
%         end
%     end
% end
% 
% n=1;
% AirPuff=zeros(Ns,1);
% pre=1;
% for i=1:length(UA_time)
%     if(UA_time(i)>IC_time(n))
%         AirPuff(n)=max(AP_AU(pre:i-1));
%         pre=i;
%         n=n+1;
%         if(n>Ns)
%             break;
%         end
%     end
% end
%  
% n=1;
% WaterReward=zeros(Ns,1);
% pre=1;
% for i=1:length(UA_time)
%     if(UA_time(i)>IC_time(n))
%         WaterReward(n)=max(WR_AU(pre:i-1));
%         pre=i;
%         n=n+1;
%         if(n>Ns)
%             break;
%         end
%     end
% end
%  
% n=1;
% XSpeed=zeros(Ns,1);
% pre=1;
% for i=1:length(UA_time)
%     if(UA_time(i)>IC_time(n))
%         XSpeed(n)=mean(XSpeed_AU(pre:i-1));
%         pre=i;
%         n=n+1;
%         if(n>Ns)
%             break;
%         end
%     end
% end
% 
%  
% n=1;
% YSpeed=zeros(Ns,1);
% pre=1;
% for i=1:length(UA_time)
%     if(UA_time(i)>IC_time(n))
%         YSpeed(n)=mean(YSpeed_AU(pre:i-1));
%         pre=i;
%         n=n+1;
%         if(n>Ns)
%             break;
%         end
%     end
% end
% 
% n=1;
% TimeOut=zeros(Ns,1);
% pre=1;
% for i=1:length(UA_time)
%     if(UA_time(i)>IC_time(n))
%         TimeOut(n)=round(mean(TO_AU(pre:i-1)));
%         pre=i;
%         n=n+1;
%         if(n>Ns)
%             break;
%         end
%     end
% end
%  
%  
%  
%  T=1:length(GoTrials);
%     
%     Angle(1)=-1;
%     S3=sum(diff(Angle)== -1);%#go
%     S4=sum(diff(Angle)== -91);%#nogo
%     
%     
%     GoSE=zeros(S3,2);
%     NoGoSE=zeros(S4,2);
%     buf=T(diff(Angle)== 1);
%     GoSE=[buf(1:S3)',T(diff(Angle)== -1)'];
%     buf=T(diff(Angle)== 91);
%     NoGoSE=[buf(1:S4)',T(diff(Angle)== -91)'];
% 
% s=2;
% AirPuffo{s}=AirPuff;
% Angleo{s}=Angle;
% GoSEo{s}=GoSE;
% GoTrialso{s}=GoTrials;
% IC_timeo{s}=IC_time;
% Licko{s}=Lick;
% NoGoSEo{s}=NoGoSE;
% NogoTrialso{s}=NogoTrials;
% RewardWindowo{s}=RewardWindow;
% TileTopLefto{s}=TileTopLeft;
% TimeOuto{s}=TimeOut;
% WaterRewardo{s}=WaterReward;
% XSpeedo{s}=XSpeed;
% %X_tileo{s}=X_tile;
% YSpeedo{s}=YSpeed;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

 
 
 Ns=11998+ds(3);
 trialObjectNumber=19;
 
    Angle_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 2}.Object.Data;
    Lick_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 3}.Object.Data;
    RW_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 4}.Object.Data;
    GoTrials_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 5}.Object.Data;
    NoGoTrials_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 6}.Object.Data;
    AP_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 7}.Object.Data;
    WR_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 8}.Object.Data;
    XSpeed_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 9}.Object.Data;
    YSpeed_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 10}.Object.Data;
    TO_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 13}.Object.Data;
    UA_time=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 2}.Object.XVector;

n=1;
Angle=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))

        Angle(n)=min(Angle_AU(pre:i-1));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end

end

n=1;
Lick=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))

        Lick(n)=max(Lick_AU(max(pre-12,1):i-1));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end

end


n=1;
RewardWindow=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))
        RewardWindow(n)=round(mean(RW_AU(pre:i-1)));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end
end

n=1;
GoTrials=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))
        GoTrials(n)=round(mean(GoTrials_AU(pre:i-1)));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end
end
 
n=1;
NogoTrials=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))
        NogoTrials(n)=round(mean(NoGoTrials_AU(pre:i-1)));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end
end

n=1;
AirPuff=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))
        AirPuff(n)=max(AP_AU(pre:i-1));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end
end
 
n=1;
WaterReward=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))
        WaterReward(n)=max(WR_AU(pre:i-1));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end
end
 
n=1;
XSpeed=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))
        XSpeed(n)=mean(XSpeed_AU(pre:i-1));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end
end

 
n=1;
YSpeed=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))
        YSpeed(n)=mean(YSpeed_AU(pre:i-1));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end
end

n=1;
TimeOut=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))
        TimeOut(n)=round(mean(TO_AU(pre:i-1)));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end
end
 
 
 
 T=1:length(GoTrials);
    
    Angle(1)=-1;
    S3=sum(diff(Angle)== -1);%#go
    S4=sum(diff(Angle)== -91);%#nogo
    
    
    GoSE=zeros(S3,2);
    NoGoSE=zeros(S4,2);
    buf=T(diff(Angle)== 1);
    GoSE=[buf(1:S3)',T(diff(Angle)== -1)'];
    buf=T(diff(Angle)== 91);
    NoGoSE=[buf(1:S4)',T(diff(Angle)== -91)'];

s=3;
AirPuffo{s}=AirPuff;
Angleo{s}=Angle;
GoSEo{s}=GoSE;
GoTrialso{s}=GoTrials;
IC_timeo{s}=IC_time;
Licko{s}=Lick;
NoGoSEo{s}=NoGoSE;
NogoTrialso{s}=NogoTrials;
RewardWindowo{s}=RewardWindow;
TileTopLefto{s}=TileTopLeft;
TimeOuto{s}=TimeOut;
WaterRewardo{s}=WaterReward;
XSpeedo{s}=XSpeed;
%X_tileo{s}=X_tile;
YSpeedo{s}=YSpeed;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 Ns=11998+ds(1);
 trialObjectNumber=17;
 
    Angle_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 2}.Object.Data;
    Lick_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 3}.Object.Data;
    RW_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 4}.Object.Data;
    GoTrials_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 5}.Object.Data;
    NoGoTrials_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 6}.Object.Data;
    AP_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 7}.Object.Data;
    WR_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 8}.Object.Data;
    XSpeed_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 9}.Object.Data;
    YSpeed_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 10}.Object.Data;
    TO_AU=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 13}.Object.Data;
    UA_time=StorageArray.Children{1, trialObjectNumber}.Object.Children{1, 2}.Object.XVector;

n=1;
Angle=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))

        Angle(n)=min(Angle_AU(pre:i-1));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end

end

n=1;
Lick=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))

        Lick(n)=max(Lick_AU(max(pre-12,1):i-1));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end

end


n=1;
RewardWindow=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))
        RewardWindow(n)=round(mean(RW_AU(pre:i-1)));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end
end

n=1;
GoTrials=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))
        GoTrials(n)=round(mean(GoTrials_AU(pre:i-1)));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end
end
 
n=1;
NogoTrials=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))
        NogoTrials(n)=round(mean(NoGoTrials_AU(pre:i-1)));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end
end

n=1;
AirPuff=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))
        AirPuff(n)=max(AP_AU(pre:i-1));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end
end
 
n=1;
WaterReward=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))
        WaterReward(n)=max(WR_AU(pre:i-1));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end
end
 
n=1;
XSpeed=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))
        XSpeed(n)=mean(XSpeed_AU(pre:i-1));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end
end

 
n=1;
YSpeed=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))
        YSpeed(n)=mean(YSpeed_AU(pre:i-1));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end
end

n=1;
TimeOut=zeros(Ns,1);
pre=1;
for i=1:length(UA_time)
    if(UA_time(i)>IC_time(n))
        TimeOut(n)=round(mean(TO_AU(pre:i-1)));
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end
end
 
 
 
 T=1:length(GoTrials);
    
    Angle(1)=-1;
    S3=sum(diff(Angle)== -1);%#go
    S4=sum(diff(Angle)== -91);%#nogo
    
    
    GoSE=zeros(S3,2);
    NoGoSE=zeros(S4,2);
    buf=T(diff(Angle)== 1);
    GoSE=[buf(1:S3)',T(diff(Angle)== -1)'];
    buf=T(diff(Angle)== 91);
    NoGoSE=[buf(1:S4)',T(diff(Angle)== -91)'];



    
    
for i=3:3
    AirPuff=[AirPuff ; AirPuffo{i}];
    Angle=[Angle ; Angleo{i}];
    GoSE=[GoSE ; (GoSEo{i} + length(GoTrials))];
    GoTrials=[GoTrials ; GoTrialso{i}];
    IC_time=[IC_time , (IC_timeo{i}+max(IC_time))];
    Lick=[Lick;Licko{i}];
    NoGoSE=[NoGoSE ; (NoGoSEo{i} + length(NogoTrials))];
    NogoTrials=[NogoTrials ; NogoTrialso{i}];
    RewardWindow=[RewardWindow; RewardWindowo{i}];
    TimeOut=[TimeOut;TimeOuto{i}];
    WaterReward=[WaterReward; WaterRewardo{i}];
    XSpeed=[XSpeed;XSpeedo{i}];
    YSpeed=[YSpeed;YSpeedo{i}];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    Speed = sqrt(XSpeed.^2 + YSpeed.^2);

Hit=zeros(size(GoTrials));
Miss=zeros(size(GoTrials));
CR=zeros(size(GoTrials));
FA=zeros(size(GoTrials));
nH=0;nM=0;nC=0;nF=0;
clear HitSE MissSE FASE CRSE
for i=1:size(GoSE,1)
    if(max(WaterReward(GoSE(i,1):GoSE(i,2)))==1)
        Hit(GoSE(i,1):GoSE(i,2))=GoTrials(GoSE(i,1):GoSE(i,2));
        nH=nH+1;
        HitSE(nH,:)=GoSE(i,:);
    else
        Miss(GoSE(i,1):GoSE(i,2))=GoTrials(GoSE(i,1):GoSE(i,2));
        nM=nM+1;
        MissSE(nM,:)=GoSE(i,:);
    end
end
Punish=RewardWindow & AirPuff;
for i=1:size(NoGoSE,1)
    if(max(Punish(NoGoSE(i,1):NoGoSE(i,2)+1))==1)
        FA(NoGoSE(i,1):NoGoSE(i,2))=NogoTrials(NoGoSE(i,1):NoGoSE(i,2));
        nF=nF+1;
        FASE(nF,:)=NoGoSE(i,:);
    else
        CR(NoGoSE(i,1):NoGoSE(i,2))=NogoTrials(NoGoSE(i,1):NoGoSE(i,2));
        nC=nC+1;
        CRSE(nC,:)=NoGoSE(i,:);
    end
end

Delay=zeros(size(GoTrials));
for i=2:length(Delay)
    if (RewardWindow(i)==1 && RewardWindow(i-1)==0)
        Delay(i-5:i-1)=1;
    end
end


    
    
    
    save(strcat(savePath,'\All_Sessions'),'X_tile','Y_tile',...
    'cellCount','IC_time','cellData','eventBin','Angle','Lick','RewardWindow',...
    'GoTrials','NogoTrials','AirPuff','WaterReward','XSpeed','YSpeed','Speed','TimeOut'...
    ,'GoSE','NoGoSE','Hit','Miss','CR','FA','HitSE','MissSE','CRSE','FASE','Delay','Punish','-v7.3');
    
    save(strcat(savePath,'\cellImage'),'cellImage','TileTopLeft','-v7.3')




