%Score Correlations based on average projection of the validation set
P1='E:\Reports2\2018_03_01_FisherInfo\Raw_Bin\Datasets--5to0\HitCR';
%{     
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
%}
Sessions={'','_Session1','_Session4','_Session7','_Session10','_Session13'};
AddressSetup;
ActiveMode=1;
rep=100;
TimeVec=-5:5;
CMtx=zeros(8,8,60,length(TimeVec),length(Sessions));
for sess=1:length(Sessions)
for Mouse=63

    
    S=load(strcat(P1,'\Session',num2str(Mouse),'_mode1',Sessions{sess},'.mat'));
    S_D=load(strcat(P1,'\Session',num2str(Mouse),'_mode2',Sessions{sess},'.mat'));
    S_R=load(strcat(P1,'\Session',num2str(Mouse),'_mode3',Sessions{sess},'.mat'));
    clear Scores0
    clear Scores1
    Scores0{9,60}=[];
    Scores1{9,60}=[];
    
    for area=1:8
        
        %%%%%%%%%%%%%%%%%%%%% Stim
        ValProj0=S.ValProj0;
        ValProj1=S.ValProj1;
        Nt0=size(ValProj0,2);
        Nt1=size(ValProj1,2);
        clear mProj0;clear mProj1;
        mProj0{Nt0}=[];
        mProj1{Nt1}=[];
        NP0=zeros(Nt0,1);
        NP1=zeros(Nt1,1);
        for r=51:100
            for i=1:Nt0
                if r==51 mProj0{i}=zeros(1,25);end
                if length(ValProj0{area,i,r})>0
                    mProj0{i}=mProj0{i}+ValProj0{area,i,r};
                    NP0(i)=NP0(i)+1;
                end
                if length(mProj0{i})>0 && r==100
                    mProj0{i}=mProj0{i}/ NP0(i);
                    
                end
                
            end
            for i=1:Nt1
                if r==51 mProj1{i}=zeros(1,25);end
                if length(ValProj1{area,i,r})>0
                    mProj1{i}=mProj1{i}+ValProj1{area,i,r};
                    NP1(i)=NP1(i)+1;
                    
                end
                if length(mProj1{i})>0 && r==100
                    mProj1{i}=mProj1{i}/NP1(i);
                    
                end
            end
        end
        S0=mProj0;
        S1=mProj1;
        %%%%%%%%%%%%%%%%%%%%% Delay
        ValProj0=S_D.ValProj0;
        ValProj1=S_D.ValProj1;
        Nt0=size(ValProj0,2);
        Nt1=size(ValProj1,2);
        clear mProj0;clear mProj1;
        mProj0{Nt0}=[];
        mProj1{Nt1}=[];
        NP0=zeros(Nt0,1);
        NP1=zeros(Nt1,1);
        for r=51:100
            for i=1:Nt0
                if r==51 mProj0{i}=zeros(1,5);end
                if length(ValProj0{area,i,r})>0
                    mProj0{i}=mProj0{i}+ValProj0{area,i,r};
                    NP0(i)=NP0(i)+1;
                end
                if length(mProj0{i})>0 && r==100
                    mProj0{i}=mProj0{i}/ NP0(i);
                    
                end
                
            end
            for i=1:Nt1
                if r==51 mProj1{i}=zeros(1,5);end
                if length(ValProj1{area,i,r})>0
                    mProj1{i}=mProj1{i}+ValProj1{area,i,r};
                    NP1(i)=NP1(i)+1;
                    
                end
                if length(mProj1{i})>0 && r==100
                    mProj1{i}=mProj1{i}/NP1(i);
                    
                end
            end
        end
        D0=mProj0;
        D1=mProj1;
         %%%%%%%%%%%%%%%%%%%%% Reward
        ValProj0=S_R.ValProj0;
        ValProj1=S_R.ValProj1;
        Nt0=size(ValProj0,2);
        Nt1=size(ValProj1,2);
        clear mProj0;clear mProj1;
        mProj0{Nt0}=[];
        mProj1{Nt1}=[];
        NP0=zeros(Nt0,1);
        NP1=zeros(Nt1,1);
        for r=51:100
            for i=1:Nt0
                if r==51 mProj0{i}=zeros(1,30);end
                if length(ValProj0{area,i,r})>0
                    mProj0{i}=mProj0{i}+ValProj0{area,i,r};
                    NP0(i)=NP0(i)+1;
                end
                if length(mProj0{i})>0 && r==100
                    mProj0{i}=mProj0{i}/ NP0(i);
                    
                end
                
            end
            for i=1:Nt1
                if r==51 mProj1{i}=zeros(1,30);end
                if length(ValProj1{area,i,r})>0
                    mProj1{i}=mProj1{i}+ValProj1{area,i,r};
                    NP1(i)=NP1(i)+1;
                    
                end
                if length(mProj1{i})>0 && r==100
                    mProj1{i}=mProj1{i}/NP1(i);
                    
                end
            end
        end
        R0=mProj0;
        R1=mProj1;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        
        
        
        
        for td=1:25
            score=[];
            for k=1:length(S0)
                if ~(isnan(S0{k}(td)))
                  score=[score;S0{k}(td)];
                end
            end
           
            Scores0{area,td}=score;
            score=[];
            for k=1:length(S1)
                if ~(isnan(S1{k}(td)))
                score=[score;S1{k}(td)];
                end
            end
           
            Scores1{area,td}=score;
        end
        
        for td=1:5
            score=[];
            for k=1:length(D0)
                if ~(isnan(D0{k}(td)))
                score=[score;D0{k}(td)];
                end
            end
          
            Scores0{area,25+td}=score;
            score=[];
            for k=1:length(D1)
                if ~(isnan(D1{k}(td)))
                score=[score;D1{k}(td)];
                end
            end
          
            Scores1{area,td+25}=score;
        end
        
        for td=1:30
            score=[];
            for k=1:length(R0)
                if ~(isnan(R0{k}(td)))
                score=[score;R0{k}(td)];
                end
            end
        
            Scores0{area,30+td}=score;
            score=[];
            for k=1:length(R1)
                if ~(isnan(R1{k}(td)))
                score=[score;R1{k}(td)];
                end
            end
     
            Scores1{area,td+30}=score;
        end
 
    end

% for a2=2:8  
%     if length(Scores0{a2,td})>0
%     td=7; 
%     figure();hold on;
%     xlabel('V1 Score');ylabel(strcat(Areas{a2},' Score'));
%     plot((Scores0{1,td}),(Scores0{a2,td}),'r.')
%     plot((Scores1{1,td}),(Scores1{a2,td}),'b.')
%     end
% end  
%     
    
for tsind=1:length(TimeVec)
    ts=TimeVec(tsind);
for td=1:25
    if td+ts <1 || td+ts>25
         continue;
    end
    for a1=1:8
        for a2=1:8
            if length(Scores0{a1,td})>0 && length(Scores0{a2,td})>0
            x=[(Scores0{a1,td}-mean(Scores0{a1,td}));(Scores1{a1,td}-mean(Scores1{a1,td}))];
            y=[(Scores0{a2,td+ts}-mean(Scores0{a2,td+ts}));(Scores1{a2,td+ts}-mean(Scores1{a2,td+ts}))];
            c=corrcoef(x,y);
            CMtx(a1,a2,td,tsind,sess)=c(1,2);
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
            if length(Scores0{a1,td})>0 && length(Scores0{a2,td})>0
            x=[(Scores0{a1,td}-mean(Scores0{a1,td}));(Scores1{a1,td}-mean(Scores1{a1,td}))];
            y=[(Scores0{a2,td+ts}-mean(Scores0{a2,td+ts}));(Scores1{a2,td+ts}-mean(Scores1{a2,td+ts}))];
            c=corrcoef(x,y);
            CMtx(a1,a2,td,tsind,sess)=c(1,2);
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
            if length(Scores0{a1,td})>0 && length(Scores0{a2,td})>0
            x=[(Scores0{a1,td}-mean(Scores0{a1,td}));(Scores1{a1,td}-mean(Scores1{a1,td}))];
            y=[(Scores0{a2,td+ts}-mean(Scores0{a2,td+ts}));(Scores1{a2,td+ts}-mean(Scores1{a2,td+ts}))];
            c=corrcoef(x,y);
            CMtx(a1,a2,td,tsind,sess)=c(1,2);
            end
            
        end
    end
end




end
    
  
    
    
    
    
    
end
end



Colors=linspecer(10);
figure();hold on
MaxInd=zeros(1,length(11:20));
for i=11:20
    plot(squeeze(mean(CMtx(8,1,i,:,1),1)),'Color',Colors(i-10,:),'LineStyle','--')
    plot(squeeze(mean(CMtx(1,8,i,:,1),1)),'Color',Colors(i-10,:))

    [~,MaxInd(i-10)]=max(squeeze(mean(CMtx(1,6,i,:,1),5)));
end
% figure();
% plot(MaxInd);
% 



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
 plot(-0.4:0.1:5.5,squeeze(mean(CMtx(a1,a2,:,6,2:6),5))','Color',Colors( a2,:),'LineWidth',1.5,'LineStyle','--');
 plot(-0.4:0.1:5.5,squeeze(mean(CMtx(a1,a2,:,6,1),5))','Color',Colors( a2,:),'LineWidth',1.5);
 plot([0,0],[-0.2,0.5],'Color',[0.6,0.6,0.6],'LineStyle','--');
 plot([2,2],[-0.2,0.5],'Color',[0.6,0.6,0.6],'LineStyle','--');
 plot([2.5,2.5],[-0.2,0.5],'Color',[0.6,0.6,0.6],'LineStyle','--');

end
%legend('LV','MV','PPC','A','S','M','RSC');
end














