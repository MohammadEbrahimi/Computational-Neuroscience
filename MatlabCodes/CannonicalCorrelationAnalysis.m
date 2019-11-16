SimPath='C:\Users\Sadegh\Documents\VLMReborn\Reports\2016_27_10_Fisher_PLS_TI_TB_ManualMap\DuringStim\';
DelPath='C:\Users\Sadegh\Documents\VLMReborn\Reports\2016_27_10_Fisher_PLS_TI_TB_ManualMap\DuringDelay\';
RewPath='C:\Users\Sadegh\Documents\VLMReborn\Reports\2016_27_10_Fisher_PLS_TI_TB_ManualMap\DuringReward\';

LoadPath{1}= 'C:\Users\Sadegh\Documents\VLMReborn\L347\Data\8_4_2015';
LoadPath{2}= 'C:\Users\Sadegh\Documents\VLMReborn\L347\Data\8_5_2015';
LoadPath{3}= 'C:\Users\Sadegh\Documents\VLMReborn\L347\Data\8_6_2015';
LoadPath{4}= 'C:\Users\Sadegh\Documents\VLMReborn\L347\Data\8_7_2015';
LoadPath{5}= 'C:\Users\Sadegh\Documents\VLMReborn\L347\Data\8_8_2015';

LoadPath{6}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_3_2015';
LoadPath{7}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_4_2015';
LoadPath{8}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_5_2015';
LoadPath{9}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_6_2015';
LoadPath{10}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_7_2015';
LoadPath{11}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_8_2015';

LoadPath{12}= 'C:\Users\Sadegh\Documents\VLMReborn\L362bis\Data\8_3_2015';
LoadPath{13}= 'C:\Users\Sadegh\Documents\VLMReborn\L362bis\Data\8_5_2015';
LoadPath{14}= 'C:\Users\Sadegh\Documents\VLMReborn\L362bis\Data\8_6_2015';
LoadPath{15}= 'C:\Users\Sadegh\Documents\VLMReborn\L362bis\Data\8_7_2015';
LoadPath{16}= 'C:\Users\Sadegh\Documents\VLMReborn\L362bis\Data\8_8_2015';

LoadPath{17}= 'C:\Users\Sadegh\Documents\VLMReborn\L347\Data\allDays';
LoadPath{18}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\allDays';
LoadPath{19}= 'C:\Users\Sadegh\Documents\VLMReborn\L362bis\Data\allDays';

FileNames{1}='L347_8_4_2015_active';
FileNames{2}='L347_8_5_2015_active';
FileNames{3}='L347_8_6_2015_active';
FileNames{4}='L347_8_7_2015_active';
FileNames{5}='L347_8_8_2015_active';

FileNames{6}='L354_8_3_2015_active';
FileNames{7}='L354_8_4_2015_active';
FileNames{8}='L354_8_5_2015_active';
FileNames{9}='L354_8_6_2015_active';
FileNames{10}='L354_8_7_2015_active';
FileNames{11}='L354_8_8_2015_active';%%Short%%

FileNames{12}='L362_8_3_2015_active';
FileNames{13}='L362_8_5_2015_active';
FileNames{14}='L362_8_6_2015_active';
FileNames{15}='L362_8_7_2015_active';
FileNames{16}='L362_8_8_2015_active';

FileNames{17}='L347_all_active';
FileNames{18}='L354_all_active';
FileNames{19}='L362_all_active';


Days=[19];
T=15;
time=-0.4:0.1:5.5;
LowDimVec=1:100;
%Days=[1 2 3 4 5 10 11 12 13 15 16]; %% error days

ND=length(Days);
NT=length(time);
NL=length(LowDimVec);

CrossCorr=zeros(8,8,2*T+1,ND);

clear CanCorrDim
clear CanCorrScores

clear TrialCorrInd
clear TrialCorr





color={'b','b--','k--','r','r--','m','k','g'};
color2={'bo','ro','ko'};
TitleV={'V1','LV','MV','PTLP','A','S','M','RSC'};
set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',16);
set(0,'defaultAxesFontName','Calibri');



for j=1:length(Days)
    i=Days(j);
%     V1=open(strcat(SimPath,'HitCR\',FileNames{i},'_tr.mat'));
    
    [X,CortexArea,HitC,MissC,CRC,FAC,HitSE,MissSE,CRSE,FASE]=RetriveData(LoadPath{Days(j)},1);
   
    SE=[CRSE];
    
    TrialCorrInd{j}=zeros(size(SE,1),8,8);
    TrialCorr{j}=zeros(size(SE,1),8,8);
    group=HitC+CRC;
    L=sum(group);
    ScoreMat=zeros(8,8,2,L);
    for dt=-T:T
        dt
    for a1=1
        for a2=1:8
            X1=X(CortexArea==a1,group==1)';
            X2=X(CortexArea==a2,group==1)';
            X1=X1((max(1+dt,1):min(L,L+dt)),:);
            X2=X2(max(1-dt,1):min(L,L-dt),:);
            
            
            [A,B,r]=canoncorr(X1,X2);
            CrossCorr(a1,a2,dt+T+1,j)=r(1);
            CannonCorrDim{a1,a2,1,dt+T+1,j}=A(:,1);
            CannonCorrDim{a1,a2,2,dt+T+1,j}=B(:,1);


            
            
            
        end
    end
    end
    
%     for a1=1:8
%         a1
%         for a2=1:8
%             [~,maxcorrind]=max(CrossCorr(a1,a2,:,j));
%             CannonCorrScores{a1,a2,1,j}=CannonCorrDim{a1,a2,1,maxcorrind,j}' * X(CortexArea==a1,:);
%             CannonCorrScores{a1,a2,2,j}=CannonCorrDim{a1,a2,2,maxcorrind,j}' * X(CortexArea==a2,:);
%             
%             s1=CannonCorrScores{a1,a2,1,j};
%             s2=CannonCorrScores{a1,a2,2,j};
%             s1(group==0)=0;
%             s2(group==0)=0;
%             
% 
%             
%             
%             for k=1:size(SE,1)
%                 
%                 x1=s1(SE(k,1):SE(k,2));
%                 x2=s2(SE(k,1):SE(k,2));
%                 
%                 [CC,MC,MCind]=CrossCorrelation(x1,x2,T);
%                 TrialCorrInd{j}(k,a1,a2)=MCind;
%                 TrialCorr{j}(k,a1,a2)=MC;
%                 
%                 
%             end
%         end
%     end
    
    
    
    
    
    
end

% hist=zeros(2*T+1,8,8,ND);
% MT=zeros(8,8,ND);
% for j=1:ND
%     for a1=1:8
%         for a2=1:8
%             for i=-T:T
%                 hist(i+T+1,a1,a2,j)=sum(TrialCorrInd{j}(TrialCorr{j}(:,a1,a2)>0,a1,a2)==i);
%             end
%             MT(a1,a2,j)=mean(TrialCorrInd{j}(TrialCorr{j}(:,a1,a2)>0,a1,a2));
%         end
%     end
%     
% end
% 
% 
% 
% 
% 
% a1=4;
% a2=6;
% for i=1:3
% figure();hold on
% xlabel('Time Shift (s)');
% ylabel('Probablity');
% title('<V1,PTLp> Correlation');
% bar((-T:T)/10,hist1(:,a1,a2,i)/sum(hist1(:,a1,a2,i)),1,'b');
% bar((-T:T)/10,hist2(:,a1,a2,i)/sum(hist2(:,a1,a2,i)),0.5,'k');
% legend('Hits','CR');
% % bar(-T:T,hist3(:,a1,a2,i)/sum(hist3(:,a1,a2,i)),0.4,'r');
% % bar(-T:T,hist4(:,a1,a2,i)/sum(hist4(:,a1,a2,i)),0.2,'m');
% end


figure();hold on
for a=1:8
plot((-T:T)/10,squeeze(CrossCorr(1,a,:,1)),color{a});
end







