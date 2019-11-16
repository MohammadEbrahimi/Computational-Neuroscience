AddressSetup;
Mouse=1;
SpeedMode=1;
ActiveMode=3;
Sessions=4:9;
switch Mouse
        case 1
                Trials={[1:3],[4:6],[7:9],[10:12],[13,14]}; % M1
        case 2
                Trials={[1:3],[4,5],[7:9],[10:12],[13:15],[16:18]}; % M2
        case 3
                Trials={[1:3],[4:6],[7:9],[10:12],[13:15]}; %M3
        case 4
                Trials={[1],[2],[3],[4],[5]}; % M4
        case 5
                Trials={[1],[2],[3],[4,5],[6],[7]}; % M5
        case 6
                Trials={[1],[2],[3],[4,5],[6,7]}; % M6
        case 7
                Trials={[1],[2],[3],[4],[5],[6],[7]}; % M7
end


%     loadpath=LoadPath{Mouse+50};
%     load(strcat(loadpath,'/cellData.mat'));
%     load(strcat(loadpath,'/cellData_ZS.mat'));
%     load(strcat(loadpath,'/SessionLength.mat'));
%     load(strcat(loadpath,'/Datasets/Datasets--5to0.mat'));
    
    TimeVec=(sum(SessionLength(1:(min(Sessions)-1)))+1):sum(SessionLength(1:(max(Sessions))));
    X=cellData_Raw;
    CN=size(X,1);
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   S_H=[];
   for i=unique(HitTrialNumber{SpeedMode,1})
            if sum(HitTrialNumber{SpeedMode,1}==i)==25
                Vec=HitDataset{SpeedMode,1}(HitTrialNumber{SpeedMode,1}==i);
                if max(Vec)<=max(TimeVec) && min(Vec)>=min(TimeVec)
                    S_H=[S_H;X(:,Vec)];
                end
            end
   end
   S_H=reshape(S_H,[CN,size(S_H,1)/CN,25]);
   
   S_C=[];
   for i=unique(CRTrialNumber{SpeedMode,1,ActiveMode})
            if sum(CRTrialNumber{SpeedMode,1,ActiveMode}==i)==25
                Vec=CRDataset{SpeedMode,1,ActiveMode}(CRTrialNumber{SpeedMode,1,ActiveMode}==i);
                if max(Vec)<=max(TimeVec) && min(Vec)>=min(TimeVec)
                    S_C=[S_C;X(:,Vec)];
                end
            end
   end
   S_C=reshape(S_C,[CN,size(S_C,1)/CN,25]);
   
      S_M=[];
   for i=unique(MissTrialNumber{SpeedMode,1,ActiveMode})
            if sum(MissTrialNumber{SpeedMode,1,ActiveMode}==i)==25
                Vec=MissDataset{SpeedMode,1,ActiveMode}(MissTrialNumber{SpeedMode,1,ActiveMode}==i);
                if max(Vec)<=max(TimeVec) && min(Vec)>=min(TimeVec)
                    S_M=[S_M;X(:,Vec)];
                end
            end
   end
   S_M=reshape(S_M,[CN,size(S_M,1)/CN,25]);
   
      S_F=[];
   for i=unique(FATrialNumber{SpeedMode,1})
            if sum(FATrialNumber{SpeedMode,1}==i)==25
                Vec=FADataset{SpeedMode,1}(FATrialNumber{SpeedMode,1}==i);
                if max(Vec)<=max(TimeVec) && min(Vec)>=min(TimeVec)
                    S_F=[S_F;X(:,Vec)];
                end
            end
   end
   S_F=reshape(S_F,[CN,size(S_F,1)/CN,25]);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   D_H=[];
   for i=unique(HitTrialNumber{SpeedMode,2})
            if sum(HitTrialNumber{SpeedMode,2}==i)==5
                Vec=HitDataset{SpeedMode,2}(HitTrialNumber{SpeedMode,2}==i);
                if max(Vec)<=max(TimeVec) && min(Vec)>=min(TimeVec)
                    D_H=[D_H;X(:,Vec)];
                end
            end
   end
   D_H=reshape(D_H,[CN,size(D_H,1)/CN,5]);
   
   D_C=[];
   for i=unique(CRTrialNumber{SpeedMode,2,ActiveMode})
            if sum(CRTrialNumber{SpeedMode,2,ActiveMode}==i)==5
                Vec=CRDataset{SpeedMode,2,ActiveMode}(CRTrialNumber{SpeedMode,2,ActiveMode}==i);
                if max(Vec)<=max(TimeVec) && min(Vec)>=min(TimeVec)
                    D_C=[D_C;X(:,Vec)];
                end
            end
   end
   D_C=reshape(D_C,[CN,size(D_C,1)/CN,5]);
   
      D_M=[];
   for i=unique(MissTrialNumber{SpeedMode,2,ActiveMode})
            if sum(MissTrialNumber{SpeedMode,2,ActiveMode}==i)==5
                Vec=MissDataset{SpeedMode,2,ActiveMode}(MissTrialNumber{SpeedMode,2,ActiveMode}==i);
                if max(Vec)<=max(TimeVec) && min(Vec)>=min(TimeVec)
                    D_M=[D_M;X(:,Vec)];
                end
            end
   end
   D_M=reshape(D_M,[CN,size(D_M,1)/CN,5]);
   
      D_F=[];
   for i=unique(FATrialNumber{SpeedMode,2})
            if sum(FATrialNumber{SpeedMode,2}==i)==5
                Vec=FADataset{SpeedMode,2}(FATrialNumber{SpeedMode,2}==i);
                if max(Vec)<=max(TimeVec) && min(Vec)>=min(TimeVec)
                    D_F=[D_F;X(:,Vec)];
                end
            end
   end
   D_F=reshape(D_F,[CN,size(D_F,1)/CN,5]);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   R_H=[];
   for i=unique(HitTrialNumber{SpeedMode,3})
            if sum(HitTrialNumber{SpeedMode,3}==i)==30
                Vec=HitDataset{SpeedMode,3}(HitTrialNumber{SpeedMode,3}==i);
                if max(Vec)<=max(TimeVec) && min(Vec)>=min(TimeVec)
                    R_H=[R_H;X(:,Vec)];
                end
            end
   end
   R_H=reshape(R_H,[CN,size(R_H,1)/CN,30]);
   
   R_C=[];
   for i=unique(CRTrialNumber{SpeedMode,3,ActiveMode})
            if sum(CRTrialNumber{SpeedMode,3,ActiveMode}==i)==30
                Vec=CRDataset{SpeedMode,3,ActiveMode}(CRTrialNumber{SpeedMode,3,ActiveMode}==i);
                if max(Vec)<=max(TimeVec) && min(Vec)>=min(TimeVec)
                    R_C=[R_C;X(:,Vec)];
                end
            end
   end
   R_C=reshape(R_C,[CN,size(R_C,1)/CN,30]);
   
      R_M=[];
   for i=unique(MissTrialNumber{SpeedMode,3,ActiveMode})
            if sum(MissTrialNumber{SpeedMode,3,ActiveMode}==i)==30
                Vec=MissDataset{SpeedMode,3,ActiveMode}(MissTrialNumber{SpeedMode,3,ActiveMode}==i);
                if max(Vec)<=max(TimeVec) && min(Vec)>=min(TimeVec)
                    R_M=[R_M;X(:,Vec)];
                end
            end
   end
   R_M=reshape(R_M,[CN,size(R_M,1)/CN,30]);
   
      R_F=[];
   for i=unique(FATrialNumber{SpeedMode,3})
            if sum(FATrialNumber{SpeedMode,3}==i)==30
                Vec=FADataset{SpeedMode,3}(FATrialNumber{SpeedMode,3}==i);
                if max(Vec)<=max(TimeVec) && min(Vec)>=min(TimeVec)
                    R_F=[R_F;X(:,Vec)];
                end
            end
   end
   R_F=reshape(R_F,[CN,size(R_F,1)/CN,30]);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   cellType=zeros(1,CN);
   DpVec=zeros(1,CN);
   for cnum=1:CN
        buf=squeeze(abs(mean(D_H(cnum,:,:),2)-mean(D_C(cnum,:,:),2)))'./sqrt(0.5*(var(squeeze(D_H(cnum,:,:)))+var(squeeze(D_C(cnum,:,:)))));
        DpVec(cnum)=max(buf);
   end
   
 cellIndex=find(DpVec'>0.4 & CortexArea<5);
 [~,ind]=sort(DpVec(cellIndex),'descend');
 cellIndex=cellIndex(ind);
 for cnum=cellIndex'
        figHandle = figure(1);
        clf(figHandle);
        returnMap = containers.Map;
        set(figHandle, 'KeyPressFcn', ...
            @(fig_obj , eventDat) readInput(fig_obj, eventDat,returnMap));
        
        td=-0.4:0.1:5.5;
        hold on;%grid on
        title(strcat('  cell Number: ',num2str(cnum)));
        xlabel('Time(s)');
        MF=max([mean(squeeze(S_H(cnum,:,:))),mean(squeeze(D_H(cnum,:,:))),mean(squeeze(R_H(cnum,:,:))),mean(squeeze(S_C(cnum,:,:))),mean(squeeze(D_C(cnum,:,:))),mean(squeeze(R_C(cnum,:,:)))]);
        mF=min([mean(squeeze(S_H(cnum,:,:))),mean(squeeze(D_H(cnum,:,:))),mean(squeeze(R_H(cnum,:,:))),mean(squeeze(S_C(cnum,:,:))),mean(squeeze(D_C(cnum,:,:))),mean(squeeze(R_C(cnum,:,:)))]);
        plot([0,0],[mF,MF],'k--')
        plot([2,2],[mF,MF],'k--')
        plot([2.5,2.5],[mF,MF],'k--')
        shadedErrorBar(td,[mean(squeeze(S_H(cnum,:,:))),mean(squeeze(D_H(cnum,:,:))),mean(squeeze(R_H(cnum,:,:)))],sqrt([var(squeeze(S_H(cnum,:,:))),var(squeeze(D_H(cnum,:,:))),var(squeeze(R_H(cnum,:,:)))]),'b',1);
        shadedErrorBar(td,[mean(squeeze(S_C(cnum,:,:))),mean(squeeze(D_C(cnum,:,:))),mean(squeeze(R_C(cnum,:,:)))],sqrt([var(squeeze(S_C(cnum,:,:))),var(squeeze(D_C(cnum,:,:))),var(squeeze(R_C(cnum,:,:)))]),'k',1);
%          shadedErrorBar(td,[mean(squeeze(S_M(cnum,:,:))),mean(squeeze(D_M(cnum,:,:))),mean(squeeze(R_M(cnum,:,:)))],sqrt([var(squeeze(S_M(cnum,:,:))),var(squeeze(D_M(cnum,:,:))),var(squeeze(R_M(cnum,:,:)))]),'m',1);
%          shadedErrorBar(td,[mean(squeeze(S_F(cnum,:,:))),mean(squeeze(D_F(cnum,:,:))),mean(squeeze(R_F(cnum,:,:)))],sqrt([var(squeeze(S_F(cnum,:,:))),var(squeeze(D_F(cnum,:,:))),var(squeeze(R_F(cnum,:,:)))]),'r',1);
          waitfor(figHandle);
        
        newN= returnMap('newN')
        switch newN
            case 0
              
                cellType(cnum)=0;
            case 1
            
                cellType(cnum)=1;
            case 2
          
                cellType(cnum)=2;
            case 3
     
                cellType(cnum)=3;
                
        end
        
        

 end
 Note='cellData_Raw : cellType=0 : Not Memory 1: constant firing memory 2: decreasing firing memory 3: Perioding response to edge '; 
% save('E:\Reports2\2019_3_1_ShortMemoryInCortex\Mouse5_Day4\cellType','cellType','CortexArea','Note')
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   