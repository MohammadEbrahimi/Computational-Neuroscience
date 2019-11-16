%%%%% Score Correlation based on projection on decoder
P1='E:\Reports2\2018_03_02_FisherInstant\ROI_Bin\Datasets--5to0\HitCR';
Areas={'V1','LV','MV','PPC','A','S','M','RSC'};
AddressSetup;
rep=100;


for Mouse=61%61:67
clear ScoreH;
clear ScoreC;
clear W;
W{9}=[];

    
    S=load(strcat(P1,'\Session',num2str(Mouse),'_mode1.mat'));
    D=load(strcat(P1,'\Session',num2str(Mouse),'_mode2.mat'));
    R=load(strcat(P1,'\Session',num2str(Mouse),'_mode3.mat'));
for area=1:9
    for itd=1:25
        SD=S.MaxDecoders{area,itd,1}/rep;
        for r=2:rep
            if length(S.MaxDecoders{area,itd,r})>0
                SD=SD+(S.MaxDecoders{area,itd,r}/rep);
            end
        end
        W{area}=[W{area},(SD/norm(SD))];
    end
    
    for itd=1:5
        SD=D.MaxDecoders{area,itd,1}/rep;
        for r=2:rep
            if length(D.MaxDecoders{area,itd,r})>0
                SD=SD+(D.MaxDecoders{area,itd,r}/rep);
            end
            
        end
        W{area}=[W{area},(SD/norm(SD))];
    end
    
    for itd=1:30
        SD=R.MaxDecoders{area,itd,1}/rep;
        for r=2:rep
            if length(R.MaxDecoders{area,itd,r})>0
                SD=SD+(R.MaxDecoders{area,itd,r}/rep);
            end
        end
        W{area}=[W{area},(SD/norm(SD))];
    end
end





    switch (Mouse-60)
        case 1
            %Trials={[1:3],[4:6],[7:9],[10:12],[13,14]}; % M1
            Trials={[1:3],[4:6],[7:9],[10:12]}; % M1
        case 2
            Trials={[1:3],[4,5],[7:9],[10:12],[13:15],[16:18]}; % M2
        case 3
            Trials={[1:3],[4:6],[7:9],[10:12],[13:15]}; %M3
        case 4
            Trials={[1],[2],[3],[4],[5]}; % M4
        case 5
            Trials={[1],[2],[3],[4],[5],[6],[7]}; % M5
        case 6
            Trials={[1],[2],[3],[4,5],[6,7]}; % M6
        case 7
            Trials={[1],[2],[3],[4],[5],[6],[7]}; % M7
    end
    
    CMtx=zeros(8,8,60,length(Trials)+1);
    CMtx_Shuff=zeros(8,8,60,1000);

    
    
    load(strcat(LoadPath{Mouse-10},'\cellData_ZS.mat'));
    load(strcat(LoadPath{Mouse-10},'\cellData.mat'));
    load(strcat(LoadPath{Mouse-10},'\SessionLength.mat'));
    load(strcat(LoadPath{Mouse-10},'\Datasets\Datasets--5to0'));
    SessionStart=ones(1,length(Trials));
    SessionEnd=sum(SessionLength(1:max(Trials{length(Trials)})))
    for i=2:length(SessionStart)
        SessionStart(i)=sum(SessionLength(1:max(Trials{i-1})))+1;
    end
    for area=1:8
        X=cellData_ROI_bin(CortexArea==area,:);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear HitProj
        clear HitSession
        HitProj{60}=[];
        HitSession{60}=[];
        index=HitDataset{1,1};
        TNum=HitTrialNumber{1,1};
        
        for i=unique(TNum)
            if sum(TNum==i)==25 && SessionEnd>max(index(TNum==i))
                for indt=1:25
                    HitProj{indt}=[HitProj{indt};W{area}(:,indt)'*X(:,min(index(TNum==i))+indt-1)];
                    HitSession{indt}=[HitSession{indt};max(find(SessionStart<=min(index(TNum==i))))];
                end
            end
        end
        
        index=HitDataset{1,2};
        TNum=HitTrialNumber{1,2};
        
        for i=unique(TNum)
            if sum(TNum==i)==5 && SessionEnd>max(index(TNum==i))
                for indt=26:30
                    HitProj{indt}=[HitProj{indt};W{area}(:,indt)'*X(:,min(index(TNum==i))+indt-26)];
                    HitSession{indt}=[HitSession{indt};max(find(SessionStart<=min(index(TNum==i))))];
                end
            end
        end
        
        index=HitDataset{1,3};
        TNum=HitTrialNumber{1,3};
        
        for i=unique(TNum)
            if sum(TNum==i)==30 && SessionEnd>max(index(TNum==i))
                for indt=31:60
                    HitProj{indt}=[HitProj{indt};W{area}(:,indt)'*X(:,min(index(TNum==i))+indt-31)];
                    HitSession{indt}=[HitSession{indt};max(find(SessionStart<=min(index(TNum==i))))];
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear CRProj
        clear CRSession
        CRProj{60}=[];
        CRSession{60}=[];
        index=CRDataset{1,1};
        TNum=CRTrialNumber{1,1};
        
        for i=unique(TNum)
            if sum(TNum==i)==25 && SessionEnd>max(index(TNum==i))
                for indt=1:25
                    CRProj{indt}=[CRProj{indt};W{area}(:,indt)'*X(:,min(index(TNum==i))+indt-1)];
                    CRSession{indt}=[CRSession{indt};max(find(SessionStart<=min(index(TNum==i))))];
                end
            end
        end
        
        index=CRDataset{1,2};
        TNum=CRTrialNumber{1,2};
        
        for i=unique(TNum)
            if sum(TNum==i)==5 && SessionEnd>max(index(TNum==i))
                for indt=26:30
                    CRProj{indt}=[CRProj{indt};W{area}(:,indt)'*X(:,min(index(TNum==i))+indt-26)];
                    CRSession{indt}=[CRSession{indt};max(find(SessionStart<=min(index(TNum==i))))];
                end
            end
        end
        
        index=CRDataset{1,3};
        TNum=CRTrialNumber{1,3};
        
        for i=unique(TNum)
            if sum(TNum==i)==30 && SessionEnd>max(index(TNum==i))
                for indt=31:60
                    CRProj{indt}=[CRProj{indt};W{area}(:,indt)'*X(:,min(index(TNum==i))+indt-31)];
                    CRSession{indt}=[CRSession{indt};max(find(SessionStart<=min(index(TNum==i))))];
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%If you wanne normalize each session separately to avoid inter session correlation
        
        % for i=1:length(SessionStart)
        %    for ntd=1:30
        %     CRProj(HitSession(:,ntd,1)==i,ntd,1)=HitProj(HitSession(:,ntd,1)==i,ntd,1)-mean(HitProj(HitSession(:,ntd,1)==i,ntd,1));
        %     CRProj(CRSession(:,ntd,1)==i,ntd,1)=CRProj(CRSession(:,ntd,1)==i,ntd,1)-mean(CRProj(CRSession(:,ntd,1)==i,ntd,1));
        %    end
        % end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ScoreH{area,1}=HitProj;
        ScoreC{area,1}=CRProj;
        for i=1:length(SessionStart)
            clear b1;
            for indt=1:60
                b1{indt}=HitProj{indt}(HitSession{indt}==i);
            end    
                ScoreH{area,i+1}=b1;
                
            clear b1;
            for indt=1:60
                b1{indt}=CRProj{indt}(CRSession{indt}==i);
            end    
                ScoreC{area,i+1}=b1;
            
        end
    end
    

    
    for sess=1:length(SessionStart)+1
            for td=1:60

                for a1=1:8
                    for a2=1:8
                        if size(ScoreC{a1,sess}{td},1)>0 && size(ScoreC{a2,sess}{td},1)>0
                            x=[(ScoreC{a1,sess}{td}-mean(ScoreC{a1,sess}{td}));(ScoreH{a1,sess}{td}-mean(ScoreH{a1,sess}{td}))];
                            y=[(ScoreC{a2,sess}{td}-mean(ScoreC{a2,sess}{td}));(ScoreH{a2,sess}{td}-mean(ScoreH{a2,sess}{td}))];
                            
                            c=corrcoef(x,y);
                            CMtx(a1,a2,td,sess)=c(1,2);
%                             if sess==1
%                                 for sh=1:1000
%                                     c=corrcoef(x(randperm(length(x))),y(randperm(length(y))));
%                                     CMtx_Shuff(a1,a2,td,sh)=c(1,2);
%                                 end
%                             end
                        end
                        
                    end
                end
            end
    end
    
    
    
end

% Colors=linspecer(10);
% figure();hold on
% MaxInd=zeros(1,length(11:20));
% for i=11:20
%     plot(squeeze(sum(CMtx(1,4,i,:,[1,3:7]),5)),'Color',Colors(i-10,:))
%     [~,MaxInd(i-10)]=max(squeeze(mean(CMtx(1,6,i,:,[1,3:7]),5)));
% end



Areas={'V1','LV','MV','PPC','A','S','M','RSC'};
set(0,'DefaultAxesFontSize',14);
set(0,'defaultAxesFontName','Calibri');
figure();hold on;
xlabel('Time(s)');ylabel('Score Noise Correlations');
Colors=linspecer(8);
for a1=1
for a2=[2:8]
title(strcat(Areas{a1},'  Correlations'));

% 
% for M=1:7
%     plot(-0.4:0.1:5.5,squeeze(CMtx(a1,a2,:,M))','Color',Colors(M,:));
% 
% end
%mask=(squeeze(mean(CMtx(a1,a2,:,:),3)))~=0;
 plot(-0.4:0.1:5.5,squeeze(mean(CMtx(a1,a2,:,2:5),5))','Color',Colors( a2,:),'LineWidth',1.5,'LineStyle','--');
  plot(-0.4:0.1:5.5,squeeze(CMtx(a1,a2,:,1))','Color',Colors( a2,:),'LineWidth',1.5);
  %plot(-0.4:0.1:5.5,squeeze(max(CMtx_Shuff(a1,a2,:,:),[],4))','Color',Colors( a2,:),'LineWidth',1.5,'LineStyle','--');
 plot([0,0],[-0.2,0.5],'Color',[0.6,0.6,0.6],'LineStyle','--');
 plot([2,2],[-0.2,0.5],'Color',[0.6,0.6,0.6],'LineStyle','--');
 plot([2.5,2.5],[-0.2,0.5],'Color',[0.6,0.6,0.6],'LineStyle','--');

end
%legend('V1','LV','MV','PPC','A','S','M','RSC');
end

