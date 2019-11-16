AddressSetup;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 1 , Prepare the first and second
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% momentum for each session
% T_F=2;
% savingPath='C:\Users\Sadegh\Documents\VLMReborn\Reports\2017_28_03_DprimeMovies\DuringReward\';
% DayVec=[21:46];
% for Day=DayVec
%     hdf5path=hdf5Path{Day};
%     hinfo=hdf5info(hdf5path);
%     %movie=single(hdf5read(hinfo.GroupHierarchy.Datasets(1)));
%     dim=hinfo.GroupHierarchy.Datasets.Dims;
 
%     load(strcat(LoadPath{Day},'\All_Sessions_20hz.mat'));
%     
%     mode=3;
%     
%     SE0=CRSE(2:end,:);
%     SE1=HitSE(2:end,:);
%     SE2=FASE(2:end,:);
%     SE3=MissSE(2:end,:);
%     
%     
%     set(0,'defaultlinelinewidth',2);
%     set(0,'DefaultAxesFontSize',12);
%     SpeedTh=50;
%     ActiveTrialNumber=20;
%     
%     ActiveAnimal=ones(length(Lick),1);
%     % ActiveAnimal(1:204000)=0;
%     % ActiveAnimal(216000:end)=0;
%     for i=1:(length(Lick)-(75*T_F)*ActiveTrialNumber)
%         if max(Lick(i:(i+75*ActiveTrialNumber)))==0
%             ActiveAnimal(i:(i+(75*T_F)*ActiveTrialNumber))=0;
%         end
%     end
%     
%     if mode==1
%         group0=NogoTrials+GoTrials;
%         group1=NogoTrials+GoTrials;
%         LMax=20*T_F;
%         tdVec=-5:(LMax-1);
%     elseif mode==2
%         group0=Delay_C;
%         group1=Delay_C;
%         LMax=5*T_F;
%         tdVec=-5:(LMax-1);
%     elseif mode==3
%         group0=RewardWindow;
%         group1=RewardWindow;
%         LMax=30*T_F;
%         tdVec=-5:(LMax-1);
%     end
%     
%     
%     
%     Status='0 started'
%     
%     dd0=zeros(dim(1),dim(2),length(tdVec));
%     dq0=zeros(dim(1),dim(2),length(tdVec));
%     c0=zeros(1,length(tdVec));
%     for k=1:size(SE0,1)
%         
%         
%         if(max(Speed(SE0(k,1):(SE0(k,1)+55*T_F)))<=SpeedTh && max(ActiveAnimal(SE0(k,1):SE0(k,2)))==1 ...
%                 ) && (max(AirPuff(SE0(k,1):(SE0(k,1)+25*T_F))) ==0 )
%             j=0;
%             while group0(SE0(k,1)+j)==0 && ((SE0(k,1)+j) < SE0(k,2))
%                 j=j+1;
%             end
%             stimL=0;
%             for i=0:LMax
%                 if group0(SE0(k,1)+j+i)==0
%                     stimL=i;
%                     break;
%                 end
%             end
%             
%             if ((SE0(k,1)+j+i) <= SE0(k,2)+1  && stimL>0)
%                 X=readHDF5Subset(hdf5path,[0 0 SE0(k,1)+j-5],[dim(1),dim(2),stimL+5]);
%                 X=single(X);
%                 dd0(:,:,1:5+stimL)=dd0(:,:,1:5+stimL)+X;
%                 dq0(:,:,1:5+stimL)=dq0(:,:,1:5+stimL)+X.^2;
%                 c0(1:5+stimL)=c0(1:5+stimL)+1;
%             end
%         end
%     end
%     Status='0 finished, 1 started'
%     
%     dd1=zeros(dim(1),dim(2),length(tdVec));
%     dq1=zeros(dim(1),dim(2),length(tdVec));
%     c1=zeros(1,length(tdVec));
%     for k=1:size(SE1,1)
%         if(max(Speed(SE1(k,1):(SE1(k,1)+55*T_F)))<=SpeedTh && max(ActiveAnimal(SE1(k,1):SE1(k,2)))==1 ...
%                 )&& (max(AirPuff(SE1(k,1):(SE1(k,1)+25*T_F))) ==0 )
%             j=0;
%             while group1(SE1(k,1)+j)==0 && ((SE1(k,1)+j) < SE1(k,2))
%                 j=j+1;
%             end
%             
%             
%             for i=0:LMax
%                 if group1(SE1(k,1)+j+i)==0
%                     stimL=i;
%                     break;
%                 end
%             end
%             
%             
%             if ((SE1(k,1)+j+i) <= SE1(k,2)+1 && stimL>0)
%                 X=readHDF5Subset(hdf5path,[0 0 SE1(k,1)+j-5],[dim(1),dim(2),stimL+5]);
%                 X=single(X);
%                 dd1(:,:,1:5+stimL)=dd1(:,:,1:5+stimL)+X;
%                 dq1(:,:,1:5+stimL)=dq1(:,:,1:5+stimL)+X.^2;
%                 c1(1:5+stimL)=c1(1:5+stimL)+1;
%             end
%         end
%     end
%     
%     Status='1 finished, 2 started'
%     dd2=zeros(dim(1),dim(2),length(tdVec));
%     dq2=zeros(dim(1),dim(2),length(tdVec));
%     c2=zeros(1,length(tdVec));
%     for k=1:size(SE2,1)
%         
%         
%         if(max(Speed(SE2(k,1):(SE2(k,1)+55*T_F)))<=SpeedTh && max(ActiveAnimal(SE2(k,1):SE2(k,2)))==1 ...
%                 ) && (max(AirPuff(SE2(k,1):(SE2(k,1)+25*T_F))) ==0 )
%             j=0;
%             while group0(SE2(k,1)+j)==0 && ((SE2(k,1)+j) < SE2(k,2))
%                 j=j+1;
%             end
%             stimL=0;
%             for i=0:LMax
%                 if group0(SE2(k,1)+j+i)==0
%                     stimL=i;
%                     break;
%                 end
%             end
%             
%             if ((SE2(k,1)+j+i) <= SE2(k,2)+1  && stimL>0)
%                 X=readHDF5Subset(hdf5path,[0 0 SE2(k,1)+j-5],[dim(1),dim(2),stimL+5]);
%                 X=single(X);
%                 dd2(:,:,1:5+stimL)=dd2(:,:,1:5+stimL)+X;
%                 dq2(:,:,1:5+stimL)=dq2(:,:,1:5+stimL)+X.^2;
%                 c2(1:5+stimL)=c2(1:5+stimL)+1;
%             end
%         end
%     end
%     
%     Status='2finished,3 started'
%     
%     dd3=zeros(dim(1),dim(2),length(tdVec));
%     dq3=zeros(dim(1),dim(2),length(tdVec));
%      c3=zeros(1,length(tdVec));
%      for k=1:size(SE3,1)
%          if(max(Speed(SE3(k,1):(SE3(k,1)+55*T_F)))<=SpeedTh && max(ActiveAnimal(SE3(k,1):SE3(k,2)))==1 ...
%                  )&&( max(AirPuff(SE3(k,1):(SE3(k,1)+25*T_F))) ==0 )
%              j=0;
%              while group1(SE3(k,1)+j)==0 && ((SE3(k,1)+j) < SE3(k,2))
%                  j=j+1;
%              end
%              
%              stimL=0;
%              for i=0:LMax
%                  if group1(SE3(k,1)+j+i)==0
%                      stimL=i;
%                      break;
%                  end
%              end
%              
%              if ((SE3(k,1)+j+i) <= SE3(k,2)+1 && stimL>0)
%                 X=readHDF5Subset(hdf5path,[0 0 SE3(k,1)+j-5],[dim(1),dim(2),stimL+5]);
%                 X=single(X);
%                 dd3(:,:,1:5+stimL)=dd3(:,:,1:5+stimL)+X;
%                 dq3(:,:,1:5+stimL)=dq3(:,:,1:5+stimL)+X.^2;
%                  c3(1:5+stimL)=c3(1:5+stimL)+1;
%              end
%          end
%      end
% %     
% %     %%%%%%Only for the first session
% %     Status='3 finished,Shuffle started'
%     
%     dd0s=zeros(dim(1),dim(2),100);
%     dd1s=zeros(dim(1),dim(2),100);
%     dd2s=zeros(dim(1),dim(2),100);
%     dd3s=zeros(dim(1),dim(2),100);
%     dq0s=zeros(dim(1),dim(2),100);
%     dq1s=zeros(dim(1),dim(2),100);
%     dq2s=zeros(dim(1),dim(2),100);
%     dq3s=zeros(dim(1),dim(2),100);
%     
%     randIndex=floor(rand(100,1)*(dim(3)-500));
%     Sample=zeros(dim(1),dim(2),3000);
%     for i=1:100
%         Sample(:,:,(i-1)*30+1:i*30)=readHDF5Subset(hdf5path,[0 0 randIndex(i)],[dim(1),dim(2),30]);
%     end
%    
%     for iter=1:100
%       for k=1:max([c0(1),c1(1),c2(1),c3(1)])
%  
%              X=Sample(:,:,floor(rand(1)*3000)+1);
%              X=single(X);
%             if k<=c0(1)
%                 dd0s(:,:,iter)=dd0s(:,:,iter)+X;
%                 dq0s(:,:,iter)=dq0s(:,:,iter)+X.^2;
%             end
%             if k<=c1(1)
%                 dd1s(:,:,iter)=dd1s(:,:,iter)+X;
%                 dq1s(:,:,iter)=dq1s(:,:,iter)+X.^2;
%             end
%             if k<=c2(1)
%                 
%                 dd2s(:,:,iter)=dd2s(:,:,iter)+X;
%                 dq2s(:,:,iter)=dq2s(:,:,iter)+X.^2;
%             end
%             if k<=c3(1)
%                 
%                 dd3s(:,:,iter)=dd3s(:,:,iter)+X;
%                 dq3s(:,:,iter)=dq3s(:,:,iter)+X.^2;
%             end
%         end
%     end
%         Status='Done'
%     
%     FirstFrame=readHDF5Subset(hdf5path,[0 0 1],[dim(1),dim(2),1]);
%     
% 
%      save(strcat(savingPath,'DprimeComps_Day',num2str(Day),'_chfm'),'dd0','dd1','dd2','dd3','dq0','dq1','dq2','dq3','dd0s','dd1s','dd2s','dd3s','dq0s','dq1s','dq2s','dq3s','c0','c1','c2','c3','FirstFrame','-v7.3');
%     
%     
%  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PART 2 combine all the first and second
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  momnuntums to find total dprime
% 
% Sessions=40:46
% ResultLoadPath='C:\Users\Sadegh\Documents\VLMReborn\Reports\2017_28_03_DprimeMovies\DuringReward\';
% DataObj=open(strcat(ResultLoadPath,'DprimeComps_Day',int2str(Sessions(1)),'_chfm.mat'));
% meanCR=DataObj.dd0;
% meanHit=DataObj.dd1;
% meanFA=DataObj.dd2;
% meanMiss=DataObj.dd3;
% sqrdCR=DataObj.dq0;
% sqrdHit=DataObj.dq1;
% sqrdFA=DataObj.dq2;
% sqrdMiss=DataObj.dq3;
% meanCR_Shuff=DataObj.dd0s;
% meanHit_Shuff=DataObj.dd1s;
% meanFA_Shuff=DataObj.dd2s;
% meanMiss_Shuff=DataObj.dd3s;
% sqrdCR_Shuff=DataObj.dq0s;
% sqrdHit_Shuff=DataObj.dq1s;
% sqrdFA_Shuff=DataObj.dq2s;
%  sqrdMiss_Shuff=DataObj.dq3s;
% countCR=DataObj.c0;
% countHit=DataObj.c1;
% countFA=DataObj.c2;
% countMiss=DataObj.c3;
% FirstFrame=DataObj.FirstFrame;
% Dim=size(meanCR);
% 
% 
% 
% [optimizer, metric] = imregconfig('multimodal');
% optimizer.InitialRadius = 0.001;
% optimizer.Epsilon = 1.5e-6;
% optimizer.GrowthFactor = 1.01;
% optimizer.MaximumIterations = 1000;
% 
% for s=Sessions(2:end);
%     s
%     DataObj=open(strcat(ResultLoadPath,'DprimeComps_Day',int2str(s),'_chfm.mat'));
%     tform = imregtform(DataObj.FirstFrame, FirstFrame, 'similarity', optimizer, metric);
%     
%     
%     
%     for f=1:Dim(3)
%         meanCR(:,:,f)=meanCR(:,:,f)+ imwarp(DataObj.dd0(:,:,f),tform,'OutputView',imref2d([Dim(1),Dim(2)]));
%         meanHit(:,:,f)=meanHit(:,:,f)+ imwarp(DataObj.dd1(:,:,f),tform,'OutputView',imref2d([Dim(1),Dim(2)]));
%         meanFA(:,:,f)=meanFA(:,:,f)+ imwarp(DataObj.dd2(:,:,f),tform,'OutputView',imref2d([Dim(1),Dim(2)]));
%         meanMiss(:,:,f)=meanMiss(:,:,f)+ imwarp(DataObj.dd3(:,:,f),tform,'OutputView',imref2d([Dim(1),Dim(2)]));
%         
%         sqrdCR(:,:,f)=sqrdCR(:,:,f)+ imwarp(DataObj.dq0(:,:,f),tform,'OutputView',imref2d([Dim(1),Dim(2)]));
%         sqrdHit(:,:,f)=sqrdHit(:,:,f)+ imwarp(DataObj.dq1(:,:,f),tform,'OutputView',imref2d([Dim(1),Dim(2)]));
%         sqrdFA(:,:,f)=sqrdFA(:,:,f)+ imwarp(DataObj.dq2(:,:,f),tform,'OutputView',imref2d([Dim(1),Dim(2)]));
%         sqrdMiss(:,:,f)=sqrdMiss(:,:,f)+ imwarp(DataObj.dq3(:,:,f),tform,'OutputView',imref2d([Dim(1),Dim(2)]));
%         
% 
%         countCR(f)=countCR(f)+DataObj.c0(f);
%         countHit(f)=countHit(f)+DataObj.c1(f);
%         countFA(f)=countFA(f)+DataObj.c2(f);
%         countMiss(f)=countMiss(f)+DataObj.c3(f);
%         
%         
%     end
%     
%     for f=1:100
%         meanCR_Shuff(:,:,f)=meanCR_Shuff(:,:,f)+ imwarp(DataObj.dd0s(:,:,f),tform,'OutputView',imref2d([Dim(1),Dim(2)]));
%         meanHit_Shuff(:,:,f)=meanHit_Shuff(:,:,f)+ imwarp(DataObj.dd1s(:,:,f),tform,'OutputView',imref2d([Dim(1),Dim(2)]));
%         meanFA_Shuff(:,:,f)=meanFA_Shuff(:,:,f)+ imwarp(DataObj.dd2s(:,:,f),tform,'OutputView',imref2d([Dim(1),Dim(2)]));
%         meanMiss_Shuff(:,:,f)=meanMiss_Shuff(:,:,f)+ imwarp(DataObj.dd3s(:,:,f),tform,'OutputView',imref2d([Dim(1),Dim(2)]));
%         
%         sqrdCR_Shuff(:,:,f)=sqrdCR_Shuff(:,:,f)+ imwarp(DataObj.dq0s(:,:,f),tform,'OutputView',imref2d([Dim(1),Dim(2)]));
%         sqrdHit_Shuff(:,:,f)=sqrdHit_Shuff(:,:,f)+ imwarp(DataObj.dq1s(:,:,f),tform,'OutputView',imref2d([Dim(1),Dim(2)]));
%         sqrdFA_Shuff(:,:,f)=sqrdFA_Shuff(:,:,f)+ imwarp(DataObj.dq2s(:,:,f),tform,'OutputView',imref2d([Dim(1),Dim(2)]));
%         sqrdMiss_Shuff(:,:,f)=sqrdMiss_Shuff(:,:,f)+ imwarp(DataObj.dq3s(:,:,f),tform,'OutputView',imref2d([Dim(1),Dim(2)]));
%     end
%     
%     
%     
% end
% 
% pval=0.01;
% dprime=zeros(Dim);
% dprime_shuffle=zeros(Dim(1),Dim(2),100);
% dprime_pval=zeros(Dim(1),Dim(2));
% for i=1:Dim(3)
%     count0=countCR(i);
%     count1=countHit(i);
%     m0=(meanCR(:,:,i)) /(count0);
%     m1=(meanHit(:,:,i))/(count1);
%     
%     q0=(sqrdCR(:,:,i))/(count0);
%     q1=(sqrdHit(:,:,i))/(count1);
%     
%     dprime(:,:,i)=abs(m1-m0) ./ sqrt(0.5*((q1-m1.^2)+...
%         (q0-m0.^2)));
%     
%     if i==1
%         for iter=1:100
%             m0=(meanCR_Shuff(:,:,iter)) /(count0);
%             m1=(meanHit_Shuff(:,:,iter))/(count1);
%             
%             q0=(sqrdCR_Shuff(:,:,iter))/(count0);
%             q1=(sqrdHit_Shuff(:,:,iter))/(count1);
%             
%             dprime_shuffle(:,:,iter)=abs(m1-m0) ./ sqrt(0.5*((q1-m1.^2)+...
%                 (q0-m0.^2)));
%             
%         end
%         for x=1:Dim(1)
%             for y=1:Dim(2)
%                 sortedDP=sort(squeeze(dprime_shuffle(x,y,:)),'descend');
%                 dprime_pval(x,y)=sortedDP(round(100*pval));
%             end
%         end       
%     end
%     
%     Buf=dprime(:,:,i);
%     Buf(Buf<dprime_pval)=0;
%     dprime(:,:,i)=Buf;
% end
% 
% save(strcat(ResultLoadPath,'L368_AllDays'),'meanCR','meanCR_Shuff','meanHit','meanHit_Shuff','meanFA',...
%     'meanFA_Shuff','meanMiss','meanMiss_Shuff','sqrdCR','sqrdHit','sqrdFA','sqrdCR_Shuff','sqrdHit_Shuff','sqrdFA_Shuff',...
%     'sqrdMiss','sqrdMiss_Shuff','countCR','countHit','countFA','countMiss','FirstFrame','-v7.3');

%%%%%%%%%%%%%%%%Part3
Mouse='L362';
load(strcat('C:\Users\Sadegh\Documents\VLMReborn\Reports\2017_28_03_DprimeMovies\DuringStim\',...
    Mouse,'_AllDays.mat'));
Dim=size(meanCR);
Time=[1:5,6:2:45];
dp_stim=zeros(Dim);
for i=1:Dim(3)
    count0=countCR(i);
    count1=countHit(i);
    m0=(meanCR(:,:,i)) /(count0);
    m1=(meanHit(:,:,i))/(count1);
    
    q0=(sqrdCR(:,:,i))/(count0);
    q1=(sqrdHit(:,:,i))/(count1);
    
    dp_stim(:,:,i)=abs(m1-m0) ./ sqrt(0.5*((q1-m1.^2)+...
        (q0-m0.^2)));
end
%dp_stim=dp_stim(:,:,Time);
load(strcat('C:\Users\Sadegh\Documents\VLMReborn\Reports\2017_28_03_DprimeMovies\DuringDelay\',...
    Mouse,'_AllDays.mat'));
Dim=size(meanCR);
Time=[1:5,6:2:15];
dp_delay=zeros(Dim);
for i=1:Dim(3)
    count0=countCR(i);
    count1=countHit(i);
    m0=(meanCR(:,:,i)) /(count0);
    m1=(meanHit(:,:,i))/(count1);
    
    q0=(sqrdCR(:,:,i))/(count0);
    q1=(sqrdHit(:,:,i))/(count1);
    
    dp_delay(:,:,i)=abs(m1-m0) ./ sqrt(0.5*((q1-m1.^2)+...
        (q0-m0.^2)));
end
%dp_delay=dp_delay(:,:,Time);
load(strcat('C:\Users\Sadegh\Documents\VLMReborn\Reports\2017_28_03_DprimeMovies\DuringReward\',...
    Mouse,'_AllDays.mat'));
dp_reward=zeros(Dim);
Dim=size(meanCR);
Time=[1:5,6:2:65];
for i=1:Dim(3)
    count0=countCR(i);
    count1=countHit(i);
    m0=(meanCR(:,:,i)) /(count0);
    m1=(meanHit(:,:,i))/(count1);
    
    q0=(sqrdCR(:,:,i))/(count0);
    q1=(sqrdHit(:,:,i))/(count1);
    
    dp_reward(:,:,i)=abs(m1-m0) ./ sqrt(0.5*((q1-m1.^2)+...
        (q0-m0.^2)));
end
%dp_reward=dp_reward(:,:,Time);
Dim=size(dp_stim);
dprime_hit_cr=zeros(Dim(1),Dim(2),60);
dprime_hit_cr(:,:,1:25)=dp_stim;
dprime_hit_cr(:,:,26:30)=dp_delay(:,:,6:10);
dprime_hit_cr(:,:,31:60)=dp_reward(:,:,6:35);

%dprime_hit_cr_L364=dprime_hit_cr;

dprime_hit_cr=dprime_hit_cr/max(max(max(dprime_hit_cr)));
RGB=uint8(zeros(1017,1017,3,60));
for i=1:60
   % RGB(:,:,:,i)=ind2rgb(dprime_hit_cr(:,:,i),gray(1));
   RGB(:,:,1,i)=uint8(255*dprime_hit_cr(:,:,i));
   RGB(:,:,2,i)=uint8(255*dprime_hit_cr(:,:,i));
   RGB(:,:,3,i)=uint8(255*dprime_hit_cr(:,:,i));
    if (i>5 && i<26)
        RGB(1:50,967:1017,1,i)=255;
    elseif (i>25 && i<31)
        RGB(1:50,967:1017,2,i)=255;
    elseif (i>30 && i<61)
        RGB(1:50,967:1017,3,i)=255;
        
    end
end



writerObj = VideoWriter(strcat('C:\Users\Sadegh\Documents\VLMReborn\Reports\2017_28_03_DprimeMovies\DuringStim\',Mouse,'_hitcr.avi'));
writerObj.FrameRate=2;
open(writerObj)
writeVideo(writerObj,uint8(RGB(1:2:1017,1:2:1017,:,:)));
close(writerObj)

%
%
%
% dprime=zeros(size(dd0));
% for i=1:size(dd0,3)
%     dprime(:,:,i)=abs(dd1(:,:,i)/c1(i)-dd0(:,:,i)/c0(i)) ./ sqrt(0.5*((dq1(:,:,i)/c1(i)-(dd1(:,:,i)/c1(i)).^2)+(dq0(:,:,i)/c0(i)-(dd0(:,:,i)/c0(i)).^2)));
% end

