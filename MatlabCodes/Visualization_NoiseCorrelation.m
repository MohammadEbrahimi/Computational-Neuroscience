Respath='C:\Users\Sadegh\Documents\VLMReborn\Reports2\2018_01_26_NoiseCorrelation\ROI\Datasets-10to5\HitCR';
Days=[51:57];
figure();hold on;
for i=1:length(Days)
    s=Days(i);
   S=load(strcat(Respath,'\Session',num2str(s),'_mode1.mat')); 
   D=load(strcat(Respath,'\Session',num2str(s),'_mode2.mat'));  
   
   IS=max(mean(S.Info_Val_Fisher,3),[],2);
   ID=max(mean(D.Info_Val_Fisher,3),[],2);
   %cellVector=100:100:length(IS)*100;
   cellVector=S.cellVector
   subplot(2,length(Days),i)
   hold on
   plot(cellVector,IS);
   plot(cellVector,ID,'r');
   subplot(2,length(Days),i+length(Days))
   hold on
   plot(cellVector,IS/max(IS));
   plot(cellVector,ID/max(ID),'r'); 
   
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for Mouse=[51:57]
    
    
figure();hold on
xlabel('Num Cells');ylabel('Normalized Fisher Information' );
TP='15to0';
title(strcat(' M',num2str(Mouse-50), '  1.5 to 2'));

P='C:\Users\Sadegh\Documents\VLMReborn\Reports2\2018_01_26_NoiseCorrelation\ROI';
S=load(strcat(P,'\Datasets-',TP,'_allCells\HitCR\Session',num2str(Mouse),'_mode1.mat'));
I1=[0;max(mean(S.Info_Val_Fisher,3),[],2)];
C1=size(S.PLSRot{1,1},1):100:size(S.PLSRot{length(I1)-1,1},1);
[~,index]=max(mean(S.Info_Val_Fisher,3),[],2);
N0=mean(S.CorrNoise0,3);
N1=mean(S.CorrNoise1,3);
N0I=zeros(1,size(N0,1));
N1I=zeros(1,size(N1,1));
for k=1:size(N0,1)
    N0I(k)=N0(k,index(k));
    N1I(k)=N1(k,index(k));
end
% S=load(strcat(P,'\Datasets-0to15_allCells_Sep\HitCR\Session',num2str(Mouse),'_mode2.mat'));
% I1D=[0;max(mean(S.Info_Val_Fisher,3),[],2)];
% C1D=size(S.PLSRot{1,1},1):100:size(S.PLSRot{length(I1D)-1,1},1);

plot([C1],N1I,'b');
plot([C1],N0I,'r');
legend('Go Correlated Noise','Nogo Correlated Noise');

%legend('0 to 0.5','0.5 to 1','1 to 1.5','1.5 to 2')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for Mouse=[61:67]
    
    
figure();hold on

xlabel('Num Cells');ylabel('Normalized Fisher' );
TP='0to15';

P='E:\Reports2\2018_01_26_NoiseCorrelation\ROI_Bin';
 
S=load(strcat(P,'\Datasets-0to15_allCells\MissCR\Session',num2str(Mouse),'_mode1.mat'));
I1=[0;max(mean(S.Info_Val_Fisher,3),[],2)];
Step=size(S.PLSRot{2,1},1)-size(S.PLSRot{1,1},1);
C1=size(S.PLSRot{1,1},1):Step:size(S.PLSRot{length(I1)-1,1},1);

S=load(strcat(P,'\Datasets-0to15_allCells\HitCR\Session',num2str(Mouse),'_mode1.mat'));
I2=[0;max(mean(S.Info_Val_Fisher,3),[],2)];
Step=size(S.PLSRot{2,1},1)-size(S.PLSRot{1,1},1);
C2=size(S.PLSRot{1,1},1):Step:size(S.PLSRot{length(I2)-1,1},1);

S=load(strcat(P,'\Datasets-0to15_allCells\HitCR\Session',num2str(Mouse),'_mode2.mat'));
I3=[0;max(mean(S.Info_Val_Fisher,3),[],2)];
Step=size(S.PLSRot{2,1},1)-size(S.PLSRot{1,1},1);
C3=size(S.PLSRot{1,1},1):Step:size(S.PLSRot{length(I3)-1,1},1);
% 
% S=load(strcat(P,'\Datasets-15to0_allCells_RT\HitCR\Session',num2str(Mouse),'_mode1.mat'));
% I4=[0;max(mean(S.Info_Val_Fisher,3),[],2)];
% C4=size(S.PLSRot{1,1},1):Step:size(S.PLSRot{length(I4)-1,1},1);

S=load(strcat(P,'\Datasets-0to15_allCells\MissCR\Session',num2str(Mouse),'_mode2.mat'));
I1D=[0;max(mean(S.Info_Val_Fisher,3),[],2)];
Step=size(S.PLSRot{2,1},1)-size(S.PLSRot{1,1},1);
C1D=size(S.PLSRot{1,1},1):Step:size(S.PLSRot{length(I1D)-1,1},1);



% S=load(strcat(P,'\Datasets-5to10_allCells\HitCR\Session',num2str(Mouse),'_mode1.mat'));
% I2=[0;max(mean(S.Info_Val_Fisher,3),[],2)];
% C2=size(S.PLSRot{1,1},1):100:size(S.PLSRot{length(I2)-1,1},1);
% 
   plot([0,C2],I2/max(I2),'k');
   plot([0,C3],I3/max(I3),'r');
plot([0,C1],I1/max(I1),'k--');
plot([0,C1D],I1D/max(I1D),'r--');

%  plot([0,C4],I4/max(I4),'r');

%legend('Stimuli','Delay','Stimuli-NoResponse','Delay-NoResponse')
 %legend('0s to 0.5s','0.5s to 1s','1s to 1.5s','1.5s to 2s');
 
title(strcat('M',num2str(Mouse-60)));
end
%%%%%%%%%%%%%%%%%%%%%%%Instant
MouseStartBin=[7,6,6,8,8,9,9];
MouseStartBin1=[6,5,6,8,5,5,7];
CellNum=[5292,4810,2774,2236,4193,3334,3754];
td=-0.4:0.1:2.5;
p='E:\Reports2\2018_01_26_NoiseCorrelation\ROI_Bin\Datasets--5to0_allCells\HitCR';
set(0,'defaultlinelinewidth',1.5);
for Mouse=61%[61,63:67];
Colors=linspecer(60);
S=load(strcat(p,'\Session',num2str(Mouse),'_mode1.mat'));
D=load(strcat(p,'\Session',num2str(Mouse),'_mode2.mat'));
R=load(strcat(p,'\Session',num2str(Mouse),'_mode3.mat'));

CellStep=round((CellNum(Mouse-60)-70)/25);
cellVector=[0,70:CellStep:CellNum(Mouse-60)];
%cellVector=[0,R.cellVector];


HalfPoint=zeros(25,1);
figure();hold on
xlabel('Num Cells');ylabel('Normalized Information');
title(strcat('Mouse ',num2str(Mouse-60)));
for i=MouseStartBin(Mouse-60):60
    if i<=25
    IS=max(mean(S.Info_Val_Fisher(:,i,:,:),4),[],3);
    elseif i<=30
    IS=max(mean(D.Info_Val_Fisher(:,i-25,:,:),4),[],3);
    else
    IS=max(mean(R.Info_Val_Fisher(:,i-30,:,:),4),[],3);
    end
    IS=[0;IS];
     plot(cellVector,IS/max(IS),'Color',Colors(i,:));
    %plot(IS,'Color',Colors(i-4,:));
     

    HalfPoint(i)=min(cellVector(IS/max(IS) >0.5));
    
end

% figure();hold on
% xlabel('Time');ylabel('Num Cells');
% plot(td(MouseStartBin1(Mouse-60):end),HalfPoint(MouseStartBin1(Mouse-60):end),'k');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Compare Instant between areas
MouseStartBin=[7,6,6,8,8,9,9];
MouseStartBin1=[6,5,6,8,5,5,7];
CellNum=[5292,4810,2774,2236,4193,3334,3754];
td=-0.4:0.1:2;
p1='E:\Reports2\2018_01_26_NoiseCorrelation\ROI_Bin\Datasets--5to0_Visual\HitCR';
p2='E:\Reports2\2018_01_26_NoiseCorrelation\ROI_Bin\Datasets--5to0_NonVisual\HitCR';
set(0,'defaultlinelinewidth',1.5);
Colors=linspecer(100);
for Mouse=61:67;

S1=load(strcat(p1,'\Session',num2str(Mouse),'_mode1.mat'));
cellVector1=[0,S1.cellVector];

S2=load(strcat(p2,'\Session',num2str(Mouse),'_mode1.mat'));
cellVector2=[0,S2.cellVector];

NumCell=min(max(cellVector1),max(cellVector2));
n1=sum(cellVector1<=NumCell);
n2=sum(cellVector2<=NumCell);

figure();hold on
xlabel('Num Cells');ylabel('Normalized Information');
title(strcat('Mouse ',num2str(Mouse-60)));
for i=MouseStartBin(Mouse-60):25
    IS=max(mean(S1.Info_Val_Fisher(:,i,:,:),4),[],3);
    IS=[0;IS];
    plot(cellVector1(1:n1),IS(1:n1)/max(IS(1:n1)),'Color',Colors(i-5,:));
     
    IS=max(mean(S2.Info_Val_Fisher(:,i,:,:),4),[],3);
    IS=[0;IS];
    plot(cellVector2(1:n2),IS(1:n2)/max(IS(1:n2)),'Color',Colors(105-i,:));
     
end


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p='C:\Users\Sadegh\Documents\VLMReborn\Reports2\2018_01_26_NoiseCorrelation\ROI\Datasets-0to15_allCells\HitCR\Session57';
for s=[1:7]
 figure();hold on
 xlabel('Num Cells');ylabel('Normalized Fisher Information' );
 S=load(strcat(p,'_mode1_Session',num2str(s),'.mat'));
 IS=max(mean(S.Info_Val_Fisher,3),[],2);
 D=load(strcat(p,'_mode2_Session',num2str(s),'.mat'));
 ID=max(mean(D.Info_Val_Fisher,3),[],2);
 plot(IS/max(IS));
 plot(ID/max(ID),'r');

 legend('Stimuli 0 to 0.5','Delay')
  title(strcat('M',num2str(7),'   Session :',num2str(s),' Training Size S:',num2str(S.TrainingSize),' D:',num2str(D.TrainingSize)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Multi Bin

Colors=linspecer(5);
st=5;
 HalfPoint=zeros(5,25);
figure();hold on;
for i=1:5
    for t=7:(26-i)

    
    IS=squeeze(max(mean(Info_Val_Fisher(:,t,i,:,:),5),[],4));
    HalfPoint(i,t)=min(cellVector(IS/max(IS) >0.5));
    end
 plot(HalfPoint(i,:),'Color',Colors(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%






