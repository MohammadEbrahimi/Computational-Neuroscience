% cellEvents=padding(eventBin,0,2);
% set(0,'defaultlinelinewidth',2);
% set(0,'DefaultAxesFontSize',12);
% 
% Hit=zeros(size(GoTrials));
% Miss=zeros(size(GoTrials));
% CR=zeros(size(GoTrials));
% FA=zeros(size(GoTrials));
% nH=0;nM=0;nC=0;nF=0;
% clear HitSE MissSE FASE CRSE
% for i=1:size(GoSE,1)
%     if(max(WaterReward(GoSE(i,1):GoSE(i,2)))==1)
%         Hit(GoSE(i,1):GoSE(i,2))=GoTrials(GoSE(i,1):GoSE(i,2));
%         nH=nH+1;
%         HitSE(nH,:)=GoSE(i,:);
%     else
%         Miss(GoSE(i,1):GoSE(i,2))=GoTrials(GoSE(i,1):GoSE(i,2));
%         nM=nM+1;
%         MissSE(nM,:)=GoSE(i,:);
%     end
% end
% Punish=RewardWindow & AirPuff;
% for i=1:size(NoGoSE,1)
%     if(max(Punish(NoGoSE(i,1):NoGoSE(i,2)+1))==1)
%         FA(NoGoSE(i,1):NoGoSE(i,2))=NogoTrials(NoGoSE(i,1):NoGoSE(i,2));
%         nF=nF+1;
%         FASE(nF,:)=NoGoSE(i,:);
%     else
%         CR(NoGoSE(i,1):NoGoSE(i,2))=NogoTrials(NoGoSE(i,1):NoGoSE(i,2));
%         nC=nC+1;
%         CRSE(nC,:)=NoGoSE(i,:);
%     end
% end
% % % 
% % % % % 
% % % % % % %%%% significant Mutual Information   Angel=0 is Go
% % %  [pVals MI] = Mutual_Information_Shuffled(cellEvents(:,1:23700)',[GoTrials(1:23700),NogoTrials(1:23700),RewardWindow(1:23700).*(Angle(1:23700)==90),RewardWindow(1:23700).*(Angle(1:23700)==0)],2,1,1000);
% % %  FMI=MI;
% % %  FMI(pVals>0.01)=0;
% % % % [sortedMI cellIndexMI]=sort(FMI(:,2),'descend');
% % % 
% % % %%%Find cell Coordinates and assign area
% cellIJ=zeros(cellCount,2);
% CortexArea=zeros(cellCount,1);
% 
% for i=1:cellCount
%     [row,ind]=max(cellImage{i});
%     [~,indj]=max(row);
%     indi=ind(indj);
%     cellIJ(i,:)=TileTopLeft(i,:)+[indi-1,indj-1]+[5,5];
% end
% 
% for b=1:8
%     a=9-b;
%     vertex=cortexMap{a};
%     [in,on] = inpolygon(cellIJ(:,1),cellIJ(:,2),vertex(:,1),vertex(:,2));
%     CortexArea(in==1)=a;
% end
% 
% 
figure();hold on;
t=cellData(cellIndexMI(2),:)/max(cellData(cellIndexMI(2),:));
p=eventBin(cellIndexMI(2),:)/max(eventBin(cellIndexMI(2),:));

plot(t,'r');
plot(p,'b');
plot(Lick,'m');
plot(NogoTrials,'k--');
plot(GoTrials,'k');




% %%%% Firing Rate
% Lt=36;
% 
% group=RewardWindow;
% SE=HitSE;
% 
% MIVec=ones(cellCount,1);%FMI(:,4);
% 
% dd=zeros(size(SE,1),Lt,8);
% Resp=zeros(8,Lt);
% 
% Cntr=zeros(8,Lt);
% for area=1:8
%     X=mean(cellEvents(( MIVec > 0 & CortexArea==area),:));
%     for k=2:length(SE) 
%         j=0;
%         while group(SE(k,1)+j)==0 && ((SE(k,1)+j) < SE(k,2))
%             j=j+1;
%         end
%         
%         stimL=0;
%         for i=0:(Lt-5)
%             if group(SE(k,1)+j+i)==0 && ((SE(k,1)+j+i) < SE(k,2))
%                 stimL=i;
%                 break;
%             end
%         end
%         stimL
% 
%         Cntr(area,1:5+stimL)=Cntr(area,1:5+stimL)+ones(1,5+stimL);
%         for i=1:5+stimL
%             dd(Cntr(area,i),i,area)=X(SE(k,1)+j-5+i);
%         end
%         
%         
%     end
% end
% 
% CNT_G_R_G=Cntr;
% MFR_G_R_G=dd;



% %%%% have to do during reward for Hit/CR
% dd1=MFR1_R_G;
% dd2=MFR2_R_G;
% dd3=MFR3_R_G;
% Cntr1=CNT1_R_G;
% Cntr2=CNT2_R_G;
% Cntr3=CNT3_R_G;
% Lt=36;
% 
% %%%active inactive
% Cntr=Cntr1+Cntr2;
% %Cntr=Cntr3;
% 
% dd=zeros(max(max(Cntr)),Lt,8);
% for area=1:8
%     for i=1:Lt
%         %dd(1:(Cntr1(area,i)+Cntr2(area,i)+Cntr3(area,i)),i,area)=[squeeze(dd1(1:Cntr1(area,i),i,area));squeeze(dd2(1:Cntr2(area,i),i,area));squeeze(dd3(1:Cntr3(area,i),i,area))];
%         
%         dd(1:(Cntr1(area,i)+Cntr2(area,i)),i,area)=[squeeze(dd1(1:Cntr1(area,i),i,area));squeeze(dd2(1:Cntr2(area,i),i,area))];
%         %dd(1:(Cntr3(area,i)),i,area)=[squeeze(dd3(1:Cntr3(area,i),i,area))];
%         
%     end
% end
% 
Lt=36
dd=MFR_G_R_G;
Cntr=CNT_G_R_G;

avgResp=zeros(8,Lt);
errResp=zeros(8,Lt);
color=['b','r','r','r','g','y','k','m'];
for area=1:8
    avgResp(area,:)=sum(dd(:,:,area))./Cntr(area,:);
    for i=1:Lt
        errResp(area,i)=sqrt(var(dd(1:Cntr(area,i),i,area)))/sqrt(Cntr(area,i));
    end
end

avgResp_g_a_g=avgResp;
errResp_g_a_g=errResp;

% 

save('C:\Users\Sadegh\Documents\VLMReborn\Reports\2015_11_10_Firing rate\L362bis\L362_8_3_2015_duringReward\HitCR\AllCells\AllCells',...
    'avgResp_g_a_g','avgResp_g_a_n','errResp_g_a_g','errResp_g_a_n')
% 
% % % % % % % save('C:\Users\Sadegh\Documents\VLMReborn\Reports\2015_11_10_Firing rate\L347_2015_08_4_duringReward\HitCr\NogoCells',...
% % % % % % %     'avgResp_n_a_g','avgResp_n_a_n','avgResp_n_i_g','avgResp_n_i_n',...
% % % % % % %     'errResp_n_a_g','errResp_n_a_n','errResp_n_i_g','errResp_n_i_n');
% % 
TitleV={'V1','LV','MV','PTLP','A','S','M','RSC'};

Lt=35
td=-0.4:0.1:3;
for area=[1,2,3,4,5,6,7,8];
figure();
 hold on;

shadedErrorBar(td,avgResp_g_a_g(area,1:Lt),errResp_g_a_g(area,1:Lt),'b',1);
shadedErrorBar(td,avgResp_g_a_n(area,1:Lt),errResp_g_a_n(area,1:Lt),'r',1);

% shadedErrorBar(td,avgResp_i_g(area,1:Lt),errResp_i_g(area,1:Lt),'b--',1);
% shadedErrorBar(td,avgResp_i_n(area,1:Lt),errResp_i_n(area,1:Lt),'r--',1);
title(TitleV(area));
end


% for area=[1,2,3,4,5,6,7,8];
% figure();
%  hold on;
% 
% shadedErrorBar(td,avgResp_n_a_g(area,1:Lt),errResp_n_a_g(area,1:Lt),'b',1);
% shadedErrorBar(td,avgResp_n_a_n(area,1:Lt),errResp_n_a_n(area,1:Lt),'r',1);
% % shadedErrorBar(td,avgResp_n_i_g(area,1:Lt),errResp_n_i_g(area,1:Lt),'b--',1);
% % shadedErrorBar(td,avgResp_n_i_n(area,1:Lt),errResp_n_i_n(area,1:Lt),'r--',1);
% title(TitleV(area));
% end



% 
fil=sum(BMAT,2);
cellCount=size(cellData,1);
cellMap=zeros(1017,1017,3);
% areaMap=zeros(1017,1017,3);
CortexArea=zeros(cellCount,1);
RGBarea=[0 3 0 ; 0 0 3 ; 0 1 1.5 ; 0 1 1.5 ; 3 0 0 ; 0 3 0 ; 2 0 2 ;2 1 0; 3 2 0 ];
% for a=1:8
%     areaMap(:,:,1)=areaMap(:,:,1)+Area(:,:,a)*RGBarea(a+1,1);
%     areaMap(:,:,2)=areaMap(:,:,2)+Area(:,:,a)*RGBarea(a+1,2);
%     areaMap(:,:,3)=areaMap(:,:,3)+Area(:,:,a)*RGBarea(a+1,3);
% end


for i=1:cellCount

    if fil(i) ~=0
        i
        clear buf
        ind=TileTopLeft(i,:);
        buf=cellImage{i};
        
        s=size(buf);
        cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,1)=...
            cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,1)+(RGBarea(CortexArea(i)+1,1)*buf);
        
        cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,2)=...
            cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,2)+(RGBarea(CortexArea(i)+1,2)*buf);
        
        cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,3)=...
            cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,3)+(RGBarea(CortexArea(i)+1,3)*buf);
        
    end
    
end

cellMap(:,:,1)=cellMap(:,:,1)-min(min(cellMap(:,:,1)));
cellMap(:,:,1)=floor(250*cellMap(:,:,1)/max(max(cellMap(:,:,1))));

cellMap(:,:,2)=cellMap(:,:,2)-min(min(cellMap(:,:,2)));
cellMap(:,:,2)=floor(250*cellMap(:,:,2)/max(max(cellMap(:,:,2))));

cellMap(:,:,3)=cellMap(:,:,3)-min(min(cellMap(:,:,3)));
cellMap(:,:,3)=floor(250*cellMap(:,:,3)/max(max(cellMap(:,:,3))));


figure()
imshow(uint8(cellMap));


VisualMap=areaMap;
VisualMap=VisualMap-min(min(VisualMap));
VisualMap=VisualMap/max(max(VisualMap));
VisualMap=uint8(255*VisualMap);










