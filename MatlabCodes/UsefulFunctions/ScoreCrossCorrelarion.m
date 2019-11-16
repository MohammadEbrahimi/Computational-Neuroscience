AddressSetup;
ResultsPath='C:\Users\Sadegh\Documents\VLMReborn\Reports2\2017_05_19_Fisher_PLS\HitCR_alldata\HitCR_alldata';
savePath='C:\Users\Sadegh\Documents\VLMReborn\Reports2\2017_05_20_CrossCorrelationScores\HitCR\'
mode=1;
Sessions=[17,18,19,47,48,49,50];

time=-0.4:0.1:5.5;
LowDimVec=1:50;

ND=length(Sessions);
NT=length(time);
NL=length(LowDimVec);


clear TrialCorrInd
clear TrialCorr




color={'b','b--','k--','r','r--','m','k','g'};
color2={'bo','ro','ko'};
TitleV={'V1','LV','MV','PTLP','A','S','M','RSC'};
set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',16);
set(0,'defaultAxesFontName','Calibri');

Info_Val_Fisher=zeros(8,NL,ND);
MaxInfoDim=zeros(8,ND);
MaxInfo=zeros(8,ND);

for j=1:length(Sessions)
    s=Sessions(j)
    V1=open(strcat(ResultsPath,'Session',num2str(s),'_mode',num2str(mode),'.mat'));

    
    
    
%    figure();
% 
%     maxDim=5;
%     hold on;
%     title(int2str(j))
%     grid on 
%     
    
for area=1:8
    Info_Val_Fisher(area,:,j)=V1.Info_Val_Fisher(area,:);
   %plot(Info_Val_Fisher(area,:,j),color{area})
   [MaxInfo(area,j),MaxInfoDim(area,j)]=max(squeeze(Info_Val_Fisher(area,:,j)));
end


    
    
    
    [X,CortexArea,HitC,MissC,CRC,FAC,HitSE,MissSE,CRSE,FASE]=RetriveData(LoadPath{s},1);
    group=CRC;
    SE=[CRSE];
    
    TrialCorrInd{j}=zeros(size(SE,1),8,8);
    TrialCorr{j}=zeros(size(SE,1),8,8);
    
    L=length(group);
    ScoreMat=zeros(8,sum(group));
    T=15;
    rt=randperm(size(SE,1));
    
    for a1=1:8
                    s1= V1.BMAT(:,a1,MaxInfoDim(a1,j))' * X;
                                s1(group==0)=0;
                                ScoreMat(a1,:)=s1(group==1);
        %  figure();hold on;grid on
        for a2=1:8
            a2

            s2= V1.BMAT(:,a2,MaxInfoDim(a2,j))' * X;
            s2(group==0)=0;

            
            if max(s2)==0
                continue;
            end
            
            
  
            
            for k=1:size(SE,1)
               
                x1=s1(SE(k,1):SE(k,2));
                x2=s2(SE(k,1):SE(k,2));
                
                [CC,MC,MCind]=CrossCorrelation(x1,x2,T);
                %CrossCorr(a1,a2,:,j)=CC;
                TrialCorrInd{j}(k,a1,a2)=MCind;
                TrialCorr{j}(k,a1,a2)=MC;
                
                
                %                 plot(-10:10,squeeze(CrossCorr(a1,a2,:,j)),color{a2})
                %                 plot(MCind,MC,'k*')
                
            end
        end
    end
    data=ScoreMat;
    save(strcat(savePath,'Session',num2str(s),'CRScores_mode',num2str(mode)),'data')
    
    
end


hist=zeros(2*T+1,8,8,ND);
MT=zeros(8,8,ND);
for j=1:ND
    for a1=1:8
        for a2=1:8
            for i=-T:T
                hist(i+T+1,a1,a2,j)=sum(TrialCorrInd{j}(TrialCorr{j}(:,a1,a2)>0.5,a1,a2)==i);
            end
            MT(a1,a2,j)=mean(TrialCorrInd{j}(TrialCorr{j}(:,a1,a2)>0.5,a1,a2));
        end
    end
    
end

save(strcat(savePath,'CRHistogram_mode',num2str(mode)),'hist','MT')


% 
% hist1=histF;
% hist2=histC;
% 
a1=6;
a2=8;
for i=1:7
h=figure();hold on
xlabel('Time Shift (s)');
ylabel('Probablity');
title(strcat('<MV,PPC> Correlation M',num2str(i)));
bar((-T:T)/10,hist2(:,a1,a2,i)/sum(hist2(:,a1,a2,i)),1,'b');
bar((-T:T)/10,hist1(:,a1,a2,i)/sum(hist1(:,a1,a2,i)),0.5,'k');
legend('H','C');
% bar(-T:T,hist3(:,a1,a2,i)/sum(hist3(:,a1,a2,i)),0.4,'r');
% bar(-T:T,hist4(:,a1,a2,i)/sum(hist4(:,a1,a2,i)),0.2,'m');
end







