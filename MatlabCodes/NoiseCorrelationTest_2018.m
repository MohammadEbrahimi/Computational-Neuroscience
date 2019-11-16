P1='E:\Reports2\2018_03_01_FisherInfo\ROI_Bin\Datasets--5to0\HitCR';
AddressSetup;
ActiveMode=1;
rep=100;
SMean0=zeros(7,60);
SVar0=zeros(7,60);
SMean1=zeros(7,60);
SVar1=zeros(7,60);
CMtx=zeros(8,8,60,7);
for Mouse=61:67
    S=load(strcat(P1,'\Session',num2str(Mouse),'_mode1.mat'));
    S_D=load(strcat(P1,'\Session',num2str(Mouse),'_mode2.mat'));
    S_R=load(strcat(P1,'\Session',num2str(Mouse),'_mode3.mat'));
    clear Scores0
    clear Scores1
    Scores0{9,60}=[];
    Scores1{9,60}=[];
    
    for area=9
        
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


for td=1:60
    SMean0(Mouse-60,td)=mean(Scores0{area,td});
    SMean1(Mouse-60,td)=mean(Scores1{area,td});
    SVar0(Mouse-60,td)=var(Scores0{area,td});
    SVar1(Mouse-60,td)=var(Scores1{area,td});    
end
    
  

 
    
    
    
    
end

% Areas={'V1','LV','MV','PPC','A','S','M','RSC'};
% set(0,'DefaultAxesFontSize',14);
% set(0,'defaultAxesFontName','Calibri');

figure();hold on;

set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',14);
set(0,'defaultAxesFontName','Calibri');
xlabel('Time(s)');ylabel('Variance');
ScoreVAR=zeros(60,2,6);
nM=0;
for M=[1,3:7]
    nM=nM+1;
   %figure();hold on;
    plot(-0.4:0.1:5.5,0.5*(SVar0(M,1:60)+SVar1(M,1:60)),'Color',[0,0,0.5]);
    plot(-0.35:0.1:5.45,0.5*(diff(SMean0(M,1:60)).^2 +diff(SMean1(M,1:60)).^2)/12,'Color',[1,0.5,0.5]);

  
    ScoreVAR(2:end,1,nM)=0.5*(diff(SMean0(M,1:60)).^2 +diff(SMean1(M,1:60)).^2)/12;
    ScoreVAR(1:end,2,nM)=0.5*(SVar0(M,1:60)+SVar1(M,1:60));
    
%     plot(-0.35:0.1:1.95,0.5*sqrt(diff(SMean0(M,1:25)).^2 ),'b');
%     plot(-0.4:0.1:2,sqrt(SVar0(M,1:25)),'b--');
%     plot(-0.35:0.1:1.95,0.5*sqrt(diff(SMean1(M,1:25)).^2 ),'r');
%     plot(-0.4:0.1:2,sqrt(SVar1(M,1:25)),'r--');


end
legend('Total Variance','State Transition Variance');
plot([0,0],[0,0.3],'Color',[0.6,0.6,0.6],'LineStyle','--');
plot([2,2],[0,0.3],'Color',[0.6,0.6,0.6],'LineStyle','--');
plot([2.5,2.5],[0,0.3],'Color',[0.6,0.6,0.6],'LineStyle','--');




figure();hold on;
shadedErrorBar(-0.35:0.1:5.45,squeeze(mean(ScoreVAR(2:end,1,:),3)),var(squeeze(ScoreVAR(2:end,1,:))'),{'Color',[1,0.5,0.5]},1);
shadedErrorBar(-0.4:0.1:5.5,squeeze(mean(ScoreVAR(:,2,:),3)),var(squeeze(ScoreVAR(:,2,:))'),{'Color',[0,0,0.5]},1);
plot([0,0],[-0.1,0.6],'Color',[0.6,0.6,0.6],'LineStyle','--');
plot([2,2],[-0.1,0.6],'Color',[0.6,0.6,0.6],'LineStyle','--');
plot([2.5,2.5],[-0.1,0.6],'Color',[0.6,0.6,0.6],'LineStyle','--');








