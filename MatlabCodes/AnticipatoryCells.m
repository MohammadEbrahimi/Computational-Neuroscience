AddressSetup;
HMInfo=zeros(9,7);
FCInfo=zeros(9,7);

 M=1
    load(strcat('E:\Reports2\2018_03_01_FisherInfo\Raw_Bin\Datasets--15to20\HitMissPreStim\Session6',num2str(M),'_mode1.mat'))
    HMInfo(:,M)=max(mean(Info_Val_Fisher,3),[],2);
    W=MaxDecoders{9,1}/norm(MaxDecoders{9,1});
    for r=2:100
        W=W+(MaxDecoders{9,r}/norm(MaxDecoders{9,r}));
    end
    load(strcat('E:\Reports2\2018_03_01_FisherInfo\Raw_Bin\Datasets--15to20\FACRPreStim\Session6',num2str(M),'_mode1.mat'))
    FCInfo(:,M)=max(mean(Info_Val_Fisher,3),[],2);

load(strcat(LoadPath{M+50},'/All_Sessions.mat'));
load(strcat(LoadPath{M+50},'/cellData_ZS.mat'));
load(strcat('E:\Reports2\2018_03_09_AverageActivity\RawBin\HitMiss\ConcatMouse',num2str(M),'.mat'));

cellCount=size(cellData_Raw_bin,1);
cellType=zeros(1,cellCount);
 cellIndex=find(abs(W-mean(W))>3*sqrt(var(W)));
 [~,ind]=sort(W(cellIndex),'descend')
 cellIndex=cellIndex(ind);
 for cnum=cellIndex'
        figHandle = figure(1);
        clf(figHandle);
        returnMap = containers.Map;
        set(figHandle, 'KeyPressFcn', ...
            @(fig_obj , eventDat) readInput(fig_obj, eventDat,returnMap));
        
        td=-0.9:0.1:5.5;
        hold on;%grid on
        title(strcat('  cell Number: ',num2str(cnum)));
        plot([0,0],[0,0.5],'k--')
        plot([2,2],[0,0.5],'k--')
        plot([2.5,2.5],[0,0.5],'k--')
        shadedErrorBar(td,Average1(cnum,1:length(td)),sem1(cnum,1:length(td)).*sqrt(c1(cnum,1:length(td))),'b',1);
        shadedErrorBar(td,Average0(cnum,1:length(td)),sem0(cnum,1:length(td)).*sqrt(c0(cnum,1:length(td))),'k',1);
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
