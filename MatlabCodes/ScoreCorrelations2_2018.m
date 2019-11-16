%%%%% Score Correlation based on projection on decoder
P1='E:\Reports2\2018_03_01_FisherInfo\ROI_Bin\Datasets--5to0\HitCR';
Areas={'V1','LV','MV','PPC','A','S','M','RSC'};
AddressSetup;
rep=100;


TimeVec=-5:5;

for Mouse=64:67
    clear ScoreH;
    clear ScoreC;
    switch (Mouse-60)
        case 1
            Trials={[1:3],[4:6],[7:9],[10:12],[13,14]}; % M1
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
    
    CMtx=zeros(8,8,60,length(TimeVec),length(Trials)+1);
    CMtx_Shuff=zeros(8,8,60,length(TimeVec),1000);
    CrossMtx=zeros(8,8,length(TimeVec),length(Trials)+1);
    
    S1=load(strcat(P1,'\Session',num2str(Mouse),'_mode1.mat'));
    S1_D=load(strcat(P1,'\Session',num2str(Mouse),'_mode2.mat'));
    S1_R=load(strcat(P1,'\Session',num2str(Mouse),'_mode3.mat'));
    
    load(strcat(LoadPath{Mouse-10},'\cellData_ZS.mat'));
    load(strcat(LoadPath{Mouse-10},'\cellData.mat'));
    load(strcat(LoadPath{Mouse-10},'\SessionLength.mat'));
    load(strcat(LoadPath{Mouse-10},'\Datasets\Datasets--5to0'));
    SessionStart=ones(1,length(Trials));
    for i=2:length(SessionStart)
        SessionStart(i)=sum(SessionLength(1:max(Trials{i-1})))+1;
    end
    for area=1:8
        X=cellData_ROI_bin(CortexArea==area,:);
        D1=S1.MaxDecoders{area,1}/rep;
        D1D=S1_D.MaxDecoders{area,1}/rep;
        D1R=S1_R.MaxDecoders{area,1}/rep;
        for r=2:rep
            D1=D1+(S1.MaxDecoders{area,r}/rep);
            D1D=D1D+(S1_D.MaxDecoders{area,r}/rep);
            D1R=D1R+(S1_R.MaxDecoders{area,r}/rep);
        end
        
        D1=D1/norm(D1);
        D1D=D1D/norm(D1D);
        D1R=D1R/norm(D1R);
        if length(D1)<1
            ScoreH{area}=[];
            ScoreC{area}=[];
            continue;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        HitProj=zeros(1,60,1);
        HitSession=zeros(1,60,1);
        n=1;
        index=HitDataset{1,1};
        TNum=HitTrialNumber{1,1};

        for i=unique(TNum)
            if sum(TNum==i)==25
                HitProj(n,1:25,1)=D1'*X(:,index(TNum==i));
                HitSession(n,1:25,1)=max(find(SessionStart<=min(index(TNum==i))));
                n=n+1;
            end
        end
        nd=1;
        index=HitDataset{1,2};
        TNum=HitTrialNumber{1,2};
        for i=unique(TNum)
            if sum(TNum==i)==5
                HitProj(nd,26:30,1)=D1D'*X(:,index(TNum==i));
                HitSession(nd,26:30,1)=max(find(SessionStart<=min(index(TNum==i))));
                nd=nd+1;
            end
        end
        nr=1;
        index=HitDataset{1,3};
        TNum=HitTrialNumber{1,3};
        for i=unique(TNum)
            if sum(TNum==i)==30
                HitProj(nr,31:60,1)=D1R'*X(:,index(TNum==i));
                HitSession(nr,31:60,1)=max(find(SessionStart<=min(index(TNum==i))));
                nr=nr+1;
            end
        end
        
        HitProj=HitProj(1:min([n,nd,nr]),:,:);
        HitSession=HitSession(1:min([n,nd,nr]),:,:);
       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        CRProj=zeros(1,30,1);
        CRSession=zeros(1,30,1);
        n=1;
        index=CRDataset{1,1};
        TNum=CRTrialNumber{1,1};
        for i=unique(TNum)
            if sum(TNum==i)==25
                CRProj(n,1:25,1)=D1'*X(:,index(TNum==i));
                CRSession(n,1:25,1)=max(find(SessionStart<=min(index(TNum==i))));
                n=n+1;
            end
        end
        nd=1;
        index=CRDataset{1,2};
        TNum=CRTrialNumber{1,2};
        for i=unique(TNum)
            if sum(TNum==i)==5
                CRProj(nd,26:30,1)=D1D'*X(:,index(TNum==i));
                CRSession(nd,26:30,1)=max(find(SessionStart<=min(index(TNum==i))));
                nd=nd+1;
            end
        end
        nr=1;
        index=CRDataset{1,3};
        TNum=CRTrialNumber{1,3};
        for i=unique(TNum)
            if sum(TNum==i)==30
                CRProj(nr,31:60,1)=D1R'*X(:,index(TNum==i));
                CRSession(nr,31:60,1)=max(find(SessionStart<=min(index(TNum==i))));
                nr=nr+1;
            end
        end
        CRProj=CRProj(1:min([n,nd,nr]),:,:);
        CRSession=CRSession(1:min([n,nd,nr]),:,:);
        
        
        %%%%%%%%%%%%%%%%%%%%If you wanne normalize each session separately to avoid inter session correlation
        
        % for i=1:length(SessionStart)
        %    for ntd=1:30
        %     HitProj(HitSession(:,ntd,1)==i,ntd,1)=HitProj(HitSession(:,ntd,1)==i,ntd,1)-mean(HitProj(HitSession(:,ntd,1)==i,ntd,1));
        %     CRProj(CRSession(:,ntd,1)==i,ntd,1)=CRProj(CRSession(:,ntd,1)==i,ntd,1)-mean(CRProj(CRSession(:,ntd,1)==i,ntd,1));
        %    end
        % end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ScoreH{area,1}=HitProj;
        ScoreC{area,1}=CRProj;
        for i=1:length(SessionStart)
            
            b1=HitProj(HitSession(:,1,1)==i,1:25,1);
            b2=HitProj(HitSession(:,26,1)==i,26:30,1);
            b3=HitProj(HitSession(:,31,1)==i,31:60,1);
            sb=min([size(b1,1),size(b2,1),size(b3,1)]);
            ScoreH{area,i+1}=[b1(1:sb,:),b2(1:sb,:),b3(1:sb,:)];
            
            b1=CRProj(CRSession(:,1,1)==i,1:25,1);
            b2=CRProj(CRSession(:,26,1)==i,26:30,1);
            b3=CRProj(CRSession(:,31,1)==i,31:60,1);
            sb=min([size(b1,1),size(b2,1),size(b3,1)]);
            ScoreC{area,i+1}=[b1(1:sb,:),b2(1:sb,:),b3(1:sb,:)];
        end
        
    end
    
    
    
    
    for sess=1:length(SessionStart)+1
        
        for tsind=1:length(TimeVec)
            ts=TimeVec(tsind);
            for td=1:25
                if td+ts <1 || td+ts>25
                    continue;
                end
                for a1=1:8
                    for a2=1:8
                        if size(ScoreC{a1,sess},1)>0 && size(ScoreC{a2,sess},1)>0
                            x=[(ScoreC{a1,sess}(:,td)-mean(ScoreC{a1,sess}(:,td)));(ScoreH{a1,sess}(:,td)-mean(ScoreH{a1,sess}(:,td)))];
                            y=[(ScoreC{a2,sess}(:,td+ts)-mean(ScoreC{a2,sess}(:,td+ts)));(ScoreH{a2,sess}(:,td+ts)-mean(ScoreH{a2,sess}(:,td+ts)))];
                            
                            c=corrcoef(x,y);
                            CMtx(a1,a2,td,tsind,sess)=c(1,2);
                            if sess==1 && ts==0
                                for sh=1:1000
                                    c=corrcoef(x(randperm(length(x))),y(randperm(length(y))));
                                    CMtx_Shuff(a1,a2,td,sh)=c(1,2);
                                end
                            end
                        end
                        
                    end
                end
            end
            for td=26:30
                if td+ts <26 || td+ts>30
                    continue;
                end
                for a1=1:8
                    for a2=1:8
                        if size(ScoreC{a1,sess},1)>0 && size(ScoreC{a2,sess},1)>0
                            x=[(ScoreC{a1,sess}(:,td)-mean(ScoreC{a1,sess}(:,td)));(ScoreH{a1,sess}(:,td)-mean(ScoreH{a1,sess}(:,td)))];
                            y=[(ScoreC{a2,sess}(:,td+ts)-mean(ScoreC{a2,sess}(:,td+ts)));(ScoreH{a2,sess}(:,td+ts)-mean(ScoreH{a2,sess}(:,td+ts)))];
                            c=corrcoef(x,y);
                            CMtx(a1,a2,td,tsind,sess)=c(1,2);
                             if sess==1 && ts==0
                                for sh=1:1000
                                    c=corrcoef(x(randperm(length(x))),y(randperm(length(y))));
                                    CMtx_Shuff(a1,a2,td,sh)=c(1,2);
                                end
                            end
                        end
                        
                    end
                end
            end
            for td=31:60
                if td+ts <31 || td+ts>60
                    continue;
                end
                for a1=1:8
                    for a2=1:8
                        if size(ScoreC{a1,sess},1)>0 && size(ScoreC{a2,sess},1)>0
                            x=[(ScoreC{a1,sess}(:,td)-mean(ScoreC{a1,sess}(:,td)));(ScoreH{a1,sess}(:,td)-mean(ScoreH{a1,sess}(:,td)))];
                            y=[(ScoreC{a2,sess}(:,td+ts)-mean(ScoreC{a2,sess}(:,td+ts)));(ScoreH{a2,sess}(:,td+ts)-mean(ScoreH{a2,sess}(:,td+ts)))];
                            c=corrcoef(x,y);
                            CMtx(a1,a2,td,tsind,sess)=c(1,2);
                             if sess==1 && ts==0
                                for sh=1:1000
                                    c=corrcoef(x(randperm(length(x))),y(randperm(length(y))));
                                    CMtx_Shuff(a1,a2,td,sh)=c(1,2);
                                end
                            end
                        end
                        
                    end
                end
            end
            
            %%%%%% Cross Correlation on
            for a1=1:8
                for a2=1:8
                    a1_s=[];
                    a2_s=[];
                    if size(ScoreC{a1,sess},1)>0 && size(ScoreC{a2,sess},1)>0
                        x=[ScoreC{a1,sess}-repmat(mean(ScoreC{a1,sess}),[size(ScoreC{a1,sess},1),1]);ScoreH{a1,sess}-repmat(mean(ScoreH{a1,sess}),[size(ScoreH{a1,sess},1),1])];
                        y=[ScoreC{a2,sess}-repmat(mean(ScoreC{a2,sess}),[size(ScoreC{a2,sess},1),1]);ScoreH{a2,sess}-repmat(mean(ScoreH{a2,sess}),[size(ScoreH{a2,sess},1),1])];
                        for i=1:size(x,1)
                            a1_s=[a1_s,x(i,max(6,6-ts):min(25,25-ts))];
                            a2_s=[a2_s,y(i,max(6,6+ts):min(25,25+ts))];
                            c=corrcoef(a1_s,a2_s);
                            CrossMtx(a1,a2,tsind,sess)=c(1,2);
                        end
                    end
                end
            end
            
            
        end
        
    end
    Note='Score Correlation based on projection on decoder';
    save(strcat('E:\Reports2\2018_05_03_ScoreCorrelations\ROI_Bin\Mouse',num2str(Mouse)),'CrossMtx','CMtx','CMtx_Shuff');
    
end

Colors=linspecer(15);
figure();hold on
title('PPC to V1 Mi')
xlabel ('Time Shift (s)');ylabel('Cross Correlations');
MaxInd=zeros(1,length(11:20));
for i=1:15
    plot(-0.5:0.1:0.5,squeeze(mean(TSCMtx(4,1,i,:,6),5)),'Color',Colors(i,:))
    
end



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
mask=(squeeze(mean(CMtx(a1,a2,:,6,:),3)))~=0;
 plot(-0.4:0.1:5.5,squeeze(max(CMtx(a1,a2,:,6,2:6),[],5))','Color',Colors( a2,:),'LineWidth',1.5,'LineStyle','--');
  plot(-0.4:0.1:5.5,squeeze(mean(CMtx(a1,a2,:,6,1),5))','Color',Colors( a2,:),'LineWidth',1.5);
  plot(-0.4:0.1:5.5,squeeze(max(CMtx_Shuff(a1,a2,:,:),[],4))','Color',Colors( a2,:),'LineWidth',1.5,'LineStyle','--');
%  plot([0,0],[-0.2,0.5],'Color',[0.6,0.6,0.6],'LineStyle','--');
%  plot([2,2],[-0.2,0.5],'Color',[0.6,0.6,0.6],'LineStyle','--');
%  plot([2.5,2.5],[-0.2,0.5],'Color',[0.6,0.6,0.6],'LineStyle','--');

end
%legend('LV','MV','PPC','A','S','M','RSC');
end

