AddressSetup;
Mouse=57;

cellObjectNumber=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16];
X_tile=[1 1 1 1 246 246 246 246 491 491 491 491 736 736 736 736];
Y_tile=[1 246 491 736 1 246 491 736 1 246 491 736 1 246 491 736];
%%%L347, L354, L362 : BehaviorTraceOld
MouseBehavior{1}=[1:14];
MouseBehavior{2}=[16:20,22:33];
MouseBehavior{3}=[34:48];
%%%L364, L365, L367, L368 :BehaviorTrace
MouseBehavior{4}=[1:5];
MouseBehavior{5}=[6:12];
MouseBehavior{6}=[13:19];
MouseBehavior{7}=[20:26];

%load(strcat(LoadPath{Mouse},'/CellSorting/AllCells_pcaica.mat'));
load(strcat(LoadPath{Mouse},'/CellSorting/cellSelect.mat'));
load(strcat(LoadPath{Mouse},'/areas.mat'));
load(strcat(LoadPath{Mouse},'/SessionLength.mat'));
cellCount=sum(cellSelect);
%TimeLength=size(TraceOut,2);
%cellData_Raw=zeros([cellCount,TimeLength],'single');
%cellData_Calcium=zeros([cellCount,TimeLength],'single');
%cellData_Noise=zeros([cellCount,TimeLength],'single');
%cellData_Baseline=zeros([cellCount,TimeLength],'single');
%TileTopLeft=zeros(cellCount,2);
%cellIJ=zeros(cellCount,2);
%CortexArea=zeros(cellCount,1);
%cellImage={};


% 
% n=1;
% for i=1:16
%     for k=sum(Nc(1:i-1))+1:sum(Nc(1:i))
%         if cellSelect(k)==1
%             n
%             cellImage{n}=ImageOut{k};
%             cellData_Raw(n,:)=TraceOut(k,:);
%             %if(n>3575)
%                [cellData_Calcium(n,:),cellData_Noise(n,:),cellData_Baseline(n,:)]=CaDecomposition(TraceOut(k,:));
% %             else
% %                cellData_Calcium(n,:)=cC(n,:); 
% %                cellData_Noise(n,:)=cN(n,:); 
% %                cellData_Baseline(n,:)=cB(n,:); 
% %             end
%             
%             TileTopLeft(n,:)=[X_tile(i),Y_tile(i)];
%             [row,ind]=max(cellImage{n});
%             [~,indj]=max(row);
%             indi=ind(indj);
%             cellIJ(n,:)=TileTopLeft(n,:)+[indi-1,indj-1];
%             n=n+1;
%         end
%         
%     end
% end
% 
% se=strel('disk',10,8);
% Mask=zeros(size(Area));
% for a=1:8
%     Mask(:,:,a)=Area(:,:,a)>0;
%     Mask(:,:,a)=imdilate(Mask(:,:,a),se);
% end
% for n=1:cellCount
%     for a=[8,7,5,1,6,2,3,4]
%         if Mask(cellIJ(n,1),cellIJ(n,2),a)>0
%             CortexArea(n)=a;
%         end
%     end
% end
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BehaviorTrace=BehaviorTraceOld;
s=1;
for Session=MouseBehavior{Mouse-50};
    Session

Ns=SessionLength(s);
Dim=[cellCount,Ns];

load(strcat(BehaviorTrace{Session},'\Obj_2 - Stimulation angle.mat'));
Angle_AU=Object.Data;
load(strcat(BehaviorTrace{Session},'\Obj_3 - Licking.mat'));
Lick_AU=Object.Data;
load(strcat(BehaviorTrace{Session},'\Obj_4 - Response Window.mat'));
RW_AU=Object.Data;
load(strcat(BehaviorTrace{Session},'\Obj_5 - Go trials.mat'));
GoTrials_AU=Object.Data;
load(strcat(BehaviorTrace{Session},'\Obj_6 - No Go trials.mat'));
NoGoTrials_AU=Object.Data;
load(strcat(BehaviorTrace{Session},'\Obj_7 - Air Puff.mat'));
AP_AU=Object.Data;
load(strcat(BehaviorTrace{Session},'\Obj_8 - Water Reward.mat'));
WR_AU=Object.Data;
load(strcat(BehaviorTrace{Session},'\Obj_9 - Mouse X Speed.mat'));
XSpeed_AU=Object.Data;
load(strcat(BehaviorTrace{Session},'\Obj_10 - Mouse Y Speed.mat'));
YSpeed_AU=Object.Data;
load(strcat(BehaviorTrace{Session},'\Obj_13 - Time out.mat'));
TO_AU=Object.Data;
load(strcat(BehaviorTrace{Session},'\Obj_12 - ImagingFrame.mat'));
ImgFrame_AU=Object.Data;
load(strcat(BehaviorTrace{Session},'\Obj_11 - Contrast level.mat'));
Contrast_AU=Object.Data;

UA_time=Object.XVector;


%IC_time=(1:Ns)*0.1;
IC_time=linspace(0.1,max(UA_time),Ns);


progress=0;

n=1;
Angle=zeros(Ns,1);
Lick=zeros(Ns,1);
RewardWindow=zeros(Ns,1);
GoTrials=zeros(Ns,1);
NogoTrials=zeros(Ns,1);
AirPuff=zeros(Ns,1);
WaterReward=zeros(Ns,1);
XSpeed=zeros(Ns,1);
YSpeed=zeros(Ns,1);
TimeOut=zeros(Ns,1);
ImageFrame=zeros(Ns,1);
Contrast=zeros(Ns,1);

pre=1;
for i=2:length(UA_time)
    if(UA_time(i)>=IC_time(n))
        
        Angle(n)=min(Angle_AU(pre:i-1));
        Lick(n)=max(Lick_AU(max(pre-12,1):i-1));
        RewardWindow(n)=round(mean(RW_AU(pre:i-1)));
        GoTrials(n)=round(mean(GoTrials_AU(pre:i-1)));
        NogoTrials(n)=round(mean(NoGoTrials_AU(pre:i-1)));
        AirPuff(n)=max(AP_AU(pre:i-1));
        WaterReward(n)=max(WR_AU(pre:i-1));
        XSpeed(n)=mean(XSpeed_AU(pre:i-1));
        YSpeed(n)=mean(YSpeed_AU(pre:i-1));
        TimeOut(n)=round(mean(TO_AU(pre:i-1)));
        ImageFrame(n)=round(mean(ImgFrame_AU(pre:i-1)));
        Contrast(n)=round(mean(Contrast_AU(pre:i-1)));
        
        pre=i;
        n=n+1;
        if(n>Ns)
            break;
        end
    end
    
end

progress=1;

T=1:Ns;

Angle(1)=-1;
S3=sum(diff(Angle)== -1);%%%#go
S4=sum(diff(Angle)== -91);%%%#nogo


GoSE=zeros(S3,2);
NoGoSE=zeros(S4,2);
buf=T(diff(Angle)== 1);
GoSE=[buf(1:S3)',T(diff(Angle)== -1)'];
buf=T(diff(Angle)== 91);
NoGoSE=[buf(1:S4)',T(diff(Angle)== -91)'];


Speed = sqrt(XSpeed.^2 + YSpeed.^2);
Speed=Speed-541;

progress=2;

Hit=zeros(size(GoTrials));
Miss=zeros(size(GoTrials));
CR=zeros(size(GoTrials));
FA=zeros(size(GoTrials));
nH=0;nM=0;nC=0;nF=0;
nHc=0;nMc=0;nCc=0;nFc=0;
HitSE=[];MissSE=[]; FASE=[]; CRSE=[]; HitSE_C=[]; MissSE_C=[]; FASE_C=[]; CRSE_C=[];
for i=1:size(GoSE,1)
    if(max(WaterReward(GoSE(i,1):GoSE(i,2)))==1)
        Hit(GoSE(i,1):GoSE(i,2))=GoTrials(GoSE(i,1):GoSE(i,2));
        nH=nH+1;
        HitSE(nH,:)=GoSE(i,:);
        if(sum(diff(GoTrials(GoSE(i,1):GoSE(i,2)))== -1)==1)
            nHc=nHc+1;
            HitSE_C(nHc,:)=GoSE(i,:);
            
        end
        
        
    else
        Miss(GoSE(i,1):GoSE(i,2))=GoTrials(GoSE(i,1):GoSE(i,2));
        nM=nM+1;
        MissSE(nM,:)=GoSE(i,:);
        if(sum(diff(GoTrials(GoSE(i,1):GoSE(i,2)))== -1)==1)
            nMc=nMc+1;
            MissSE_C(nMc,:)=GoSE(i,:);
            
        end
    end
end
Punish=RewardWindow & AirPuff;
for i=1:size(NoGoSE,1)
    if(max(Punish(NoGoSE(i,1):NoGoSE(i,2)+1))==1)
        FA(NoGoSE(i,1):NoGoSE(i,2))=NogoTrials(NoGoSE(i,1):NoGoSE(i,2));
        nF=nF+1;
        FASE(nF,:)=NoGoSE(i,:);
        if(sum(diff(NogoTrials(NoGoSE(i,1):NoGoSE(i,2)))== -1)==1)
            nFc=nFc+1;
            FASE_C(nFc,:)=NoGoSE(i,:);
            
        end
    else
        CR(NoGoSE(i,1):NoGoSE(i,2))=NogoTrials(NoGoSE(i,1):NoGoSE(i,2));
        nC=nC+1;
        CRSE(nC,:)=NoGoSE(i,:);
        if(sum(diff(NogoTrials(NoGoSE(i,1):NoGoSE(i,2)))== -1)==1)
            nCc=nCc+1;
            CRSE_C(nCc,:)=NoGoSE(i,:);
            
        end
    end
end
progress=3;
Delay_C=zeros(size(GoTrials));
Delay=zeros(size(GoTrials));
AllSE=[GoSE;NoGoSE];
Stimuli_i=GoTrials+NogoTrials;
for i=1:size(AllSE,1)
    TL=AllSE(i,1):AllSE(i,2);
    start=max(TL([1==0;diff(Stimuli_i(TL))==-1]));
    finish=min(TL([diff(RewardWindow(TL))==1;1==0]));
    Delay(start:finish)=1;
    if max(TimeOut(start:finish))==0
        Delay_C(start:finish)=1;
    end
end

progress=4;
AirPuffo{s}=AirPuff;
Angleo{s}=Angle;
GoSEo{s}=GoSE;
GoTrialso{s}=GoTrials;
IC_timeo{s}=IC_time;
Licko{s}=Lick;
NoGoSEo{s}=NoGoSE;
NogoTrialso{s}=NogoTrials;
RewardWindowo{s}=RewardWindow;
TimeOuto{s}=TimeOut;
WaterRewardo{s}=WaterReward;
Delayo{s}=Delay;
Punisho{s}=Punish;
XSpeedo{s}=XSpeed;
YSpeedo{s}=YSpeed;
Speedo{s}=Speed;
Hito{s}=Hit;
Misso{s}=Miss;
CRo{s}=CR;
FAo{s}=FA;
HitSEo{s}=HitSE;
MissSEo{s}=MissSE;
CRSEo{s}=CRSE;
FASEo{s}=FASE;
HitSE_Co{s}=HitSE_C;
MissSE_Co{s}=MissSE_C;
CRSE_Co{s}=CRSE_C;
FASE_Co{s}=FASE_C;
Contrasto{s}=Contrast;
ImageFrameo{s}=ImageFrame;
Delay_Co{s}=Delay_C;




s=s+1;


progress=5;



% 
%     save(strcat(LoadPath{Session+20},'\All_Session'),'ImageFrame','Contrast','X_tile','Y_tile',...
%     'cellCount','IC_time','cellData','Angle','Lick','RewardWindow',...
%     'GoTrials','NogoTrials','AirPuff','WaterReward','XSpeed','YSpeed','Speed','TimeOut'...
%     ,'GoSE','NoGoSE','Hit','Miss','CR','FA','HitSE','MissSE','CRSE','FASE','HitSE_C','MissSE_C','CRSE_C','FASE_C','Delay','Delay_C','Punish','-v7.3');
%     save(strcat(LoadPath{Session+20},'\cellImage'),'cellImage','TileTopLeft','-v7.3')

end




AirPuff=AirPuffo{1};
Angle=Angleo{1};
IC_time=IC_timeo{1};
Lick=Licko{1};
GoSE=GoSEo{1};
GoTrials=GoTrialso{1};
NoGoSE=NoGoSEo{1};
NogoTrials=NogoTrialso{1};
RewardWindow=RewardWindowo{1};
TimeOut=TimeOuto{1};
WaterReward=WaterRewardo{1};
Delay=Delayo{1};
Punish=Punisho{1};
XSpeed=XSpeedo{1};
YSpeed=YSpeedo{1};
Speed=Speedo{1};
Hit=Hito{1};
Miss=Misso{1};
CR=CRo{1};
FA=FAo{1};
HitSE=HitSEo{1};
MissSE=MissSEo{1};
CRSE=CRSEo{1};
FASE=FASEo{1};
HitSE_C=HitSE_Co{1};
MissSE_C=MissSE_Co{1};
CRSE_C=CRSE_Co{1};
FASE_C=FASE_Co{1};
Contrast=Contrasto{1};
ImageFrame=ImageFrameo{1};
Delay_C=Delay_Co{1};



start=SessionLength(1)+1;
for i=2:s-1
    
       
       AirPuff=[AirPuff;AirPuffo{i}];
       Angle=[Angle;Angleo{i}];
       IC_time=[IC_time,IC_timeo{i}+max(IC_time)];
       Lick=[Lick;Licko{i}];
       GoSE=[GoSE;GoSEo{i}+start-1];
       GoTrials=[GoTrials;GoTrialso{i}];
       NoGoSE=[NoGoSE;NoGoSEo{i}+start-1];
       NogoTrials=[NogoTrials;NogoTrialso{i}];
       RewardWindow=[RewardWindow;RewardWindowo{i}];
       TimeOut=[TimeOut;TimeOuto{i}];
       WaterReward=[WaterReward;WaterRewardo{i}];
       Delay=[Delay;Delayo{i}];
       Punish=[Punish;Punisho{i}];
       XSpeed=[XSpeed;XSpeedo{i}];
       YSpeed=[YSpeed;YSpeedo{i}];
       Speed=[Speed;Speedo{i}];
       Hit=[Hit;Hito{i}];
       Miss=[Miss;Misso{i}];
       CR=[CR;CRo{i}];
       FA=[FA;FAo{i}];
       HitSE=[HitSE;HitSEo{i}+start-1];
       MissSE=[MissSE;MissSEo{i}+start-1];
       CRSE=[CRSE;CRSEo{i}+start-1];
       FASE=[FASE;FASEo{i}+start-1];
       HitSE_C=[HitSE_C;HitSE_Co{i}+start-1];
       MissSE_C=[MissSE_C;MissSE_Co{i}+start-1];
       CRSE_C=[CRSE_C;CRSE_Co{i}+start-1];
       FASE_C=[FASE_C;FASE_Co{i}+start-1];
    Contrast=[Contrast;Contrasto{i}];
    ImageFrame=[ImageFrame;ImageFrameo{i}];
    Delay_C=[Delay_C;Delay_Co{i}];
       
       start=start+SessionLength(i);
end

  
    save(strcat(LoadPath{Mouse},'\All_Sessions'),'SessionLength',...
    'cellCount','IC_time','Angle','Lick','RewardWindow',...
    'GoTrials','NogoTrials','AirPuff','WaterReward','XSpeed','YSpeed','Speed','TimeOut'...
    ,'GoSE','NoGoSE','Hit','Miss','CR','FA','HitSE','MissSE','CRSE','FASE','HitSE_C','MissSE_C','CRSE_C','FASE_C','Delay_C','Contrast','ImageFrame','Delay','Punish','-v7.3');
    
    %save(strcat(LoadPath{Mouse},'\cellImage'),'cellImage','TileTopLeft','-v7.3');
     %save(strcat(LoadPath{Mouse},'\cellData'),'cellData_Raw','cellData_Calcium','cellData_Noise','cellData_Baseline','cellIJ','CortexArea','-v7.3');




