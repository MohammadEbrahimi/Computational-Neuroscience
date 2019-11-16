%%%% Manually clear all and open registeration
AddressSetup;

Days=26:32;
ND=length(Days);
TL=0;


for s=1:ND
    
   V=open(strcat(LoadPath{Days(s)},'\All_Session.mat'));

AirPuffo{s}=V.AirPuff;
Angleo{s}=V.Angle;
GoSEo{s}=V.GoSE;
GoTrialso{s}=V.GoTrials;
IC_timeo{s}=V.IC_time;
Licko{s}=V.Lick;
NoGoSEo{s}=V.NoGoSE;
NogoTrialso{s}=V.NogoTrials;
RewardWindowo{s}=V.RewardWindow;
TimeOuto{s}=V.TimeOut;
WaterRewardo{s}=V.WaterReward;
Delayo{s}=V.Delay;
Punisho{s}=V.Punish;
XSpeedo{s}=V.XSpeed;
YSpeedo{s}=V.YSpeed;
Speedo{s}=V.Speed;
Hito{s}=V.Hit;
Misso{s}=V.Miss;
CRo{s}=V.CR;
FAo{s}=V.FA;
HitSEo{s}=V.HitSE;
MissSEo{s}=V.MissSE;
CRSEo{s}=V.CRSE;
FASEo{s}=V.FASE;
HitSE_Co{s}=V.HitSE_C;
MissSE_Co{s}=V.MissSE_C;
CRSE_Co{s}=V.CRSE_C;
FASE_Co{s}=V.FASE_C;
Contrasto{s}=V.Contrast;
ImageFrameo{s}=V.ImageFrame;
Delay_Co{s}=V.Delay_C;


SessionLength{s}=length(V.GoTrials);

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


cellData=zeros(sum(cellMap(:,1)~=0),TL);

%eventBin=zeros(sum(cellMap(:,1)~=0),TL);
TMAP=[cellMap(:,ND),cellMap(:,1:(ND-1))];
I=open(strcat(LoadPath{Days(1)},'\cellImage.mat'));
cellImageo=I.cellImage;
TTL=I.TileTopLeft;
TileTopLeft=TTL((TMAP(:,1)~=0),:);
cellCount=sum(TMAP(:,1)~=0);

start=1;
for i=1:ND
    ncell=1;
    bufcell=cellDatao{i};
    %bufeve=eventBino{i};
    
    
   for j=1:size(TMAP,1)
       if TMAP(j,i)~=0
           if(i==1)
               cellImage{ncell}=cellImageo{TMAP(j,i)};
           end
           
           cellData(ncell,start:(start+SessionLength{i}-1))=bufcell(TMAP(j,i),:);
           %eventBin(ncell,start:(start+SessionLength{i}-1))=bufeve(TMAP(j,i),:);
           %normCellData(ncell,start:(start+SessionLength{i}-1))=(bufcell(TMAP(j,i),:)...
             %  -mean(bufcell(TMAP(j,i),:)))/sqrt(var(bufcell(TMAP(j,i),:)));
           ncell=ncell+1;
           
          
       end
   end
   
   if i>1
       
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
       
       
       
   end
   start=start+SessionLength{i};
end

  
    save('C:\Users\Sadegh\Documents\VLMReborn\L365\Data\allDays\All_Sessions','SessionLength',...
    'cellCount','IC_time','cellData','Angle','Lick','RewardWindow',...
    'GoTrials','NogoTrials','AirPuff','WaterReward','XSpeed','YSpeed','Speed','TimeOut'...
    ,'GoSE','NoGoSE','Hit','Miss','CR','FA','HitSE','MissSE','CRSE','FASE','HitSE_C','MissSE_C','CRSE_C','FASE_C','Delay_C','Contrast','ImageFrame','Delay','Punish','-v7.3');
    
    save('C:\Users\Sadegh\Documents\VLMReborn\L365\Data\allDays\cellImage','cellImage','TileTopLeft','-v7.3')




