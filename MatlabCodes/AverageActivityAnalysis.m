CodingTable=zeros(8,8,8);
GNCodingTable=zeros(8,3,8,2);
clear DPAll
clear ActivityRatio
clear ActivityRatioArea
ActivityRatio{4}=[];
ActivityRatioArea{4}=[];
DPAll{8,3,8}=[];
DPMMM=zeros(3,7);
MaxDprime=zeros(9,3,7);

AddressSetup;
path='E:\Reports2\2018_03_09_AverageActivity\RawBin\FACR\';

for Mouse=[1,3:7]
    Mouse
    load(strcat(path,'ConcatMouse',num2str(Mouse),'.mat'));
    load(strcat(LoadPath{Mouse+50},'\cellData.mat'));
    %%%%%%%%%%%%%%%%%%%%%
    
    cells_POA=zeros(cellCount,3);
    psh=zeros(cellCount,13);
    timeBin{1}=[3:6];
    timeBin{2}=[7];
    timeBin{3}=[8:13];
    pvalue=0.001;
    num_repeate=20;
    
    Dprime=zeros(cellCount,13);
    
    Dprime_sh=zeros(cellCount,13,100);
    Fstat=zeros(cellCount,13);
    DegFstat=zeros(cellCount,13);
    cellType=zeros(cellCount,3);
    
    for rep=1%1:num_repeate
    Average0_binned=squeeze(Resp0_binned(:,:,rep)./c0_binned(:,:,rep));
    Average1_binned=squeeze(Resp1_binned(:,:,rep)./c1_binned(:,:,rep));
    sem0_binned=((qResp0_binned(:,:,rep)./(c0_binned(:,:,rep)))- Average0_binned.^2)./sqrt(c0_binned(:,:,rep));
    sem1_binned=((qResp1_binned(:,:,rep)./(c1_binned(:,:,rep)))- Average1_binned.^2)./sqrt(c1_binned(:,:,rep));
    num_shuff=1000;
    for i= 1:cellCount
        M0=Average0_binned(i,:);
        C0=mean(c0_binned(i,:,1),3);
        V0=sem0_binned(i,:) .* sqrt(C0);
     
        
        M1=Average1_binned(i,:);
        C1=mean(c1_binned(i,:,1),3);
        V1=sem1_binned(i,:) .* sqrt(C1);
              for bin=1:13
            Dprime(i,bin)=Dprime(i,bin)+(abs(M0(bin)-M1(bin)) / sqrt(0.5*(V0(bin)+V1(bin))));
 
              end
    end
    end
    %Dprime=Dprime/num_repeate;
        
     for i= 1:cellCount   
        for bin=1:13
            for r=1:num_shuff
                M0_sh=Average0_binned_sh(i,bin,r);
                V0_sh=sem0_binned_sh(i,bin,r) .* sqrt(C0(bin));
                
                M1_sh=Average1_binned_sh(i,bin,r);
                V1_sh=sem1_binned_sh(i,bin,r) .* sqrt(C1(bin));
                Dprime_sh(i,bin,r)=abs(M0_sh-M1_sh) / sqrt(0.5*(V0_sh+V1_sh));
            end
            
        end
    end
    
    
    for r=1:num_shuff
        psh(:,:)=psh(:,:)+(Dprime_sh(:,:,r)>Dprime(:,:));
    end
    psh(:,:)=psh(:,:)/num_shuff;
    
    for mode=1:3
        cells_POA(:,mode)=max(psh(:,timeBin{mode})<pvalue,[],2);
        for cnum=1:cellCount
            [~,idxMax]=max(Dprime(cnum,timeBin{mode}));
            if Average0_binned(cnum,min(timeBin{mode})+idxMax-1)> Average1_binned(cnum,min(timeBin{mode})+idxMax-1)
                cellType(cnum,mode)=0;
            else
                cellType(cnum,mode)=1;
            end
           
        end 
    end
    
    if sum(max(psh(:,4:6)<pvalue,[],2)==1)>1
        DPMMM(1,Mouse)=max(max(Dprime(max(psh(:,4:6)<pvalue,[],2)==1,4:6),[],2));
        DPMMM(2,Mouse)=mean(max(Dprime(max(psh(:,4:6)<pvalue,[],2)==1,4:6),[],2));
        DPMMM(3,Mouse)=min(max(Dprime(max(psh(:,4:6)<pvalue,[],2)==1,4:6),[],2));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    S= cells_POA(:,1)==1;
    D= cells_POA(:,2)==1;
    R= cells_POA(:,3)==1;
    
    
    
    for area=1:8
        DPAll{area,1,Mouse}= max(Dprime(S & CortexArea==area,3:6),[],2);
        DPAll{area,2,Mouse}= Dprime(D & CortexArea==area,7);
        DPAll{area,3,Mouse}= max(Dprime(R & CortexArea==area,8:13),[],2);
        MaxDprime(area,1,Mouse)=max([DPAll{area,1,Mouse};0]);
        MaxDprime(area,2,Mouse)=max([DPAll{area,2,Mouse};0]);
        MaxDprime(area,3,Mouse)=max([DPAll{area,3,Mouse};0]);
        
        DPAll{area,1,8}= [DPAll{area,1,8};max(Dprime(S & CortexArea==area,3:6),[],2)];
        DPAll{area,2,8}= [DPAll{area,2,8};Dprime(D & CortexArea==area,7)];
        DPAll{area,3,8}= [DPAll{area,3,8};max(Dprime(R & CortexArea==area,8:13),[],2)];
        %%SDR
        CodingTable(area,1,Mouse)=sum(S & D & R & CortexArea==area);
        %%SD
        CodingTable(area,2,Mouse)=sum(S & D & ~R & CortexArea==area);
        %%SR
        CodingTable(area,3,Mouse)=sum(S & ~D & R & CortexArea==area);
        %%S
        CodingTable(area,4,Mouse)=sum(S & ~D & ~R & CortexArea==area);
        %%DR
        CodingTable(area,5,Mouse)=sum(~S & D & R & CortexArea==area);
        %%D
        CodingTable(area,6,Mouse)=sum(~S & D & ~R & CortexArea==area);
        %%R
        CodingTable(area,7,Mouse)=sum(~S & ~D & R & CortexArea==area);
        %%cells
        CodingTable(area,8,Mouse)=sum( CortexArea==area);
        
        %%%%% Go Nogo Coding 
        %%%%%%%%%Stimuli
        GNCodingTable(area,1,Mouse,1)=sum(S & CortexArea==area & cellType(:,1)==0);
        GNCodingTable(area,1,Mouse,2)=sum(S & CortexArea==area & cellType(:,1)==1);
        
        %%%%%%%%%Delay
        GNCodingTable(area,2,Mouse,1)=sum(D & CortexArea==area & cellType(:,2)==0);
        GNCodingTable(area,2,Mouse,2)=sum(D & CortexArea==area & cellType(:,2)==1);
        
        %%%%%%%%%Reward
        GNCodingTable(area,3,Mouse,1)=sum(R & CortexArea==area & cellType(:,3)==0);
        GNCodingTable(area,3,Mouse,2)=sum(R & CortexArea==area & cellType(:,3)==1);
    end
    
        MaxDprime(9,1,Mouse)=max(MaxDprime(1:8,1,Mouse));
        MaxDprime(9,2,Mouse)=max(MaxDprime(1:8,2,Mouse));
        MaxDprime(9,3,Mouse)=max(MaxDprime(1:8,3,Mouse));
% HC= load(strcat('C:\Users\Sadegh\Documents\VLMReborn\Reports2\2017_01_15_AverageActivity_BalSpeed\Denoised\HitCR\StimResp_Mouse',num2str(Mouse),'.mat'));
% 
% 
% for i=1:4
%     Mask=HC.cellStart(:,i)==1 & HC.cellType==0;
%     PC=max(Average0_binned(Mask,2+i:6),[],2);
%     PE=max(Average1_binned(Mask,2+i:6),[],2);
%     ActivityRatio{i}=[ActivityRatio{i}; (PC-PE)./PC];
%     ActivityRatioArea{i}=[ActivityRatioArea{i}; CortexArea(Mask)];
% end    


end

 CodingTable(:,:,8)=sum(CodingTable(:,:,1:7),3);
 GNCodingTable(:,:,8,:)=sum(GNCodingTable(:,:,1:7,:),3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NBin=5;
mode=3;
percentile=0.90;
meanDP=zeros(8,1);
LHPercentile=zeros(8,2);

xdp=zeros(8,NBin);
pdp=zeros(8,NBin);

for area=[1:6]
    DPX=sort(DPAll{area,mode,8});
    %DPX=DPX(DPX>1);
    BinW=(max(DPX)-min(DPX))/5;
    meanDP(area)=mean(DPX);
%     LHPercentile(area,1)=DPX(ceil((1-percentile)*length(DPX)));
%     LHPercentile(area,2)=DPX(floor((percentile)*length(DPX)));
    
    for nb=1:NBin
       xdp(area,nb)=((nb-0.5)*BinW) + min(DPX) ;
       pdp(area,nb)=sum(DPX<= ((nb*BinW)+ min(DPX)) &  DPX > (((nb-1)*BinW) + min(DPX)))/length(DPX);
 
    end
end

% plot(meanDP,'o')
% hold on
% plot(LHPercentile(:,1),'ro')
% plot(LHPercentile(:,2),'ro')
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% cellStart=zeros(cellCount,4);
% cellType=-1*ones(cellCount,1);
% cellPeak=zeros(cellCount,1);
% 
% cellStart(:,1)= psh(:,3)<0.001;
% cellStart(:,2)= psh(:,3)>0.05 & psh(:,4)<0.001;
% cellStart(:,3)= psh(:,3)>0.05 & psh(:,4)>0.05 & psh(:,5)<0.001;
% cellStart(:,4)= psh(:,3)>0.05 & psh(:,4)>0.05 & psh(:,5)>0.05 & psh(:,6)<0.001;
% 
% for cnum=2566:cellCount
%    
%     [ms,st]=max(cellStart(cnum,:));
%     if ms==1
%          cnum
%         %     cellPeak(cnum)=max([max(Average0_binned(cnum,2+st:6)), max(Average1_binned(cnum,2+st:6))]);
%         %     if max(Average0_binned(cnum,2+st:6)) < max(Average1_binned(cnum,2+st:6))
%         %         cellType(cnum)=1;
%         %     end
%         
%         figHandle = figure(1);
%         clf(figHandle);
%         returnMap = containers.Map;
%         set(figHandle, 'KeyPressFcn', ...
%             @(fig_obj , eventDat) readInput(fig_obj, eventDat,returnMap));
%         
%         td=-0.9:0.1:5.5;
%         hold on;%grid on
%         title(strcat('  cell Number: ',num2str(cnum),' start: ',num2str(st)));
%         shadedErrorBar(td,Average1(cnum,1:length(td)),sem1(cnum,1:length(td)),'b',1);
%         shadedErrorBar(td,Average0(cnum,1:length(td)),sem0(cnum,1:length(td)),'k',1);
%         waitfor(figHandle);
%         
%         newN= returnMap('newN')
%         switch newN
%             case 0
%                 cellStart(cnum,st)=0;
%                 cellType(cnum)=-1;
%             case 1
%                 cellPeak=max(Average1_binned(cnum,2+st:6));
%                 cellType(cnum)=1;
%             case 2
%                 cellPeak=max(Average0_binned(cnum,2+st:6));
%                 cellType(cnum)=0;
%                 
%         end
%         
%         
%         
%     end
% end
% 
% save(strcat(path,'StimResp_Mouse',num2str(Mouse),'.mat'),'cellStart','cellPeak','cellType');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% 
% for i=1:cellCount
%     %if Dprime(i,3)>dpTh && Dprime(i,8)>dpTh && cellType(i,3)~=cellType(i,8)
%     if S(i)==1 
%         figure();
%         title(strcat(num2str(i),'  Type: ',num2str(cellType(i))))
%         %title(strcat('Mouse: ',num2str(Mouse),'  cell Number: ',num2str(i),' Area: ',areaNames(CortexArea(i))));
%         td=-0.9:0.1:5.5;
%         hold on;grid on
%         shadedErrorBar(td,Average1(i,1:length(td)),sem1(i,1:length(td)),'b',1);
%         shadedErrorBar(td,Average0(i,1:length(td)),sem0(i,1:length(td)),'k',1);
%         
%     end
% end
% 


% 
% for i=[21:40]
%     figure();
% plot(squeeze(max(Dprime_sh(i,3:6,:),[],2)),'.')
% hold on
% plot([0,1000],[max(Dprime(i,3:6)),max(Dprime(i,3:6))],'r')
% end


