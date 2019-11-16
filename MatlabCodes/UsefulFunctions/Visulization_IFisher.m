SimPath='C:\Users\Sadegh\Documents\VLMReborn\Reports\2016_10_12_FisherTV\DuringStim\HitCR_alldata_comcells\';
DelPath='C:\Users\Sadegh\Documents\VLMReborn\Reports\2016_28_9_FisherTI_TB_ManualMap\DuringDelay\HitCR\';
RewPath='C:\Users\Sadegh\Documents\VLMReborn\Reports\2016_28_9_FisherTI_TB_ManualMap\DuringReward\HitCR\';

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
LoadPath{11}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_8_2015'; %%Short%%

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
FileNames{11}='L354_8_8_2015_active';

FileNames{12}='L362_8_3_2015_active';
FileNames{13}='L362_8_5_2015_active';
FileNames{14}='L362_8_6_2015_active';
FileNames{15}='L362_8_7_2015_active';
FileNames{16}='L362_8_8_2015_active';

FileNames{17}='L347_all_active';
FileNames{18}='L354_all_active';
FileNames{19}='L362_all_active';


Days=[12:15];
time=-0.4:0.1:2;
%Days=[1 2 3 4 5 10 11 12 13 15 16]; %% error days

ND=length(Days);
NT=length(time);

VC=zeros(ND,NT);
TC=zeros(ND,NT);
FAC=zeros(ND,NT);
MissC=zeros(ND,NT);
EV=zeros(8,NT,ND);
IF=zeros(8,NT,ND);
NIF=zeros(8,NT,ND);
FF=zeros(8,NT,ND);
FC=zeros(8,NT,ND);
MM=zeros(8,NT,ND);
MH=zeros(8,NT,ND);
CellNum=zeros(8,ND);

HitScores=zeros(8,NT,ND);
MissScores=zeros(8,NT,ND);
FAScores=zeros(8,NT,ND);
CRScores=zeros(8,NT,ND);

isTrained=ones(8,NT,ND);

color={'b','b--','k--','r','r--','m','k','g'};

TitleV={'V1','LV','MV','PTLP','A','S','M','RSC'};
set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',16);
set(0,'defaultAxesFontName','Calibri');


for j=1:length(Days)
    i=Days(j);
V1=open(strcat(SimPath,FileNames{i},'_tr.mat'));
V2=open(strcat(DelPath,FileNames{i},'_tr.mat'));
V3=open(strcat(RewPath,FileNames{i},'_tr.mat'));


 V1.IFisher(isnan(V1.IFisher))=0;
 V2.IFisher(isnan(V2.IFisher))=0;
 V3.IFisher(isnan(V3.IFisher))=0;



%  VC(j,:)=2*[V1.ValidationSize*ones(1,25), V2.ValidationSize*ones(1,5) , V3.ValidationSize*ones(1,30)]; 
%  TC(j,:)=2*[V1.TrainingSize*ones(1,25) , V2.TrainingSize*ones(1,5) , V3.TrainingSize*ones(1,30)]; 
%  FAC(j,:)=[V1.FATrialCount(1,:) , V2.FATrialCount(1,6:10) , V3.FATrialCount(1,6:35)]; 
%  MissC(j,:)=[V1.MissTrialCount(1,:), V2.MissTrialCount(1,6:10) , V3.MissTrialCount(1,6:35)]; 

 
EV(:,:,j)=[V1.errorVal];%,V2.errorVal(:,6:10),V3.errorVal(:,6:35)];
%  IF(:,:,j)=[V1.IFisher ,V2.IFisher(:,6:10),V3.IFisher(:,6:35)];
%  FF(:,:,j)=[V1.FAonFA ,V2.FAonFA(:,6:10),V3.FAonFA(:,6:35)];
%  FC(:,:,j)=[V1.FAonCR ,V2.FAonCR(:,6:10),V3.FAonCR(:,6:35)];
%  MM(:,:,j)=[V1.MissonMiss ,V2.MissonMiss(:,6:10),V3.MissonMiss(:,6:35)];
%  MH(:,:,j)=[V1.MissonHit ,V2.MissonHit(:,6:10),V3.MissonHit(:,6:35)];
 

% Data= IF(:,:,j);
% for area=1:8
%  CellNum(area,j)=sum(V1.BMAT(:,area)~=0);
% end
% 


% 
%  Data=EV(:,:,j);
% figure();
% xlabel('Time(s)');
% ylabel('L1-Decoder Error Rate (%)')
% hold on;
% title(int2str(i))
% grid on
% for area=1:8
% plot(time,Data(area,1:length(time)),color{area});
% end


end




% IF(isnan(IF))=0;
% FF(isnan(FF))=0;
% FC(isnan(FC))=0;
% 
L=25;
A=eye(L)/3;
for i=1:L-1
    A(i,i+1)=1/3;
    A(i+1,i)=1/3;
end
A(1:2,1)=0.5;
A(L-1:L,L)=0.5;
A(4:6,5)=[0 ;1 ;0];
% %A(28:30,29)=[0;1;0];
% A(29:31,30)=[0;1;0];
% %A(23:25,24)=[0;1;0];
% A(24:26,25)=[0;1;0];
% 



    figure();
    %title(TitleV{area});
    hold on
    grid on
     plot(zeros(1,80),0.0125:0.0125:1,'k.')
     plot(2*ones(1,80),0.0125:0.0125:1,'k.')
     plot(2.5*ones(1,80),0.0125:0.0125:1,'k.')
    xlabel('Time(s)');
    ylabel('L1 Decoder Error Rate (%)')
%   
maxErrorIndex=zeros(8,4);
accErrorIndex=zeros(8,ND,4);

    
for area=1:8
   
%    figure();
%     title(TitleV{area});
%     hold on
%     grid on
% %      plot(zeros(1,80),0.0125:0.0125:1,'k.')
% %      plot(2*ones(1,80),0.0125:0.0125:1,'k.')
% %      plot(2.5*ones(1,80),0.0125:0.0125:1,'k.')
%     xlabel('Number of Coding Cells');
%     ylabel('L1 Decoder Error Rate (%)')
%     
%  
  
%     filter=squeeze(isTrained(area,:,:))';
     trialCounter=VC;
%     
     Data=squeeze(EV(area,:,:))'; 
%     Data=Data .*filter;
%     trialCounter=trialCounter .*filter;
%     
%     Data =Data .* trialCounter;
%       

      
%     trialCounter2=VC;
%     Data2=squeeze(MH(area,:,:))'; 
%     Data2=Data2 .*filter;
%     trialCounter2=trialCounter2 .*filter;
%     
%     Data2 =Data2 .* trialCounter2;


    x=sum(Data)./ND;%sum(trialCounter);

     plot(time,x,color{area})
%      plot(time,y,'b')

   
%     
%     maxErrorIndex(area,1)=max(x(1:5)./y(1:5));
%     maxErrorIndex(area,2)=max(x(11:25)./y(11:25));
%     maxErrorIndex(area,3)=max(x(26:30)./y(26:30));
%     maxErrorIndex(area,4)=max(x(31:60)./y(31:60));
   
%     accErrorIndex(area,:,1)=sum(x(:,1:5),2)./sum(y(:,1:5),2);
%     accErrorIndex(area,:,2)=sum(x(:,6:25),2)./sum(y(:,6:25),2);
%     accErrorIndex(area,:,3)=sum(x(:,26:30),2)./sum(y(:,26:30),2);
%     accErrorIndex(area,:,4)=sum(x(:,31:60),2)./sum(y(:,31:60),2);
%     

%     shadedErrorBar(time,meanD(1:length(time)),errBar,'m',1);


end

% legend('V1','LV','MV','PTLP','A','S','M','RSC');
% 
% % legend('V1','LV','MV','S'); 
% % 
% % 
% 
% 
% 





