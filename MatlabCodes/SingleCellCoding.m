AddressSetup;

for Mouse=[61,63:67]
   
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
    
    
    load(strcat(LoadPath{Mouse-10},'\cellData_ZS.mat'));
    load(strcat(LoadPath{Mouse-10},'\cellData.mat'));
    load(strcat(LoadPath{Mouse-10},'\SessionLength.mat'));
    load(strcat(LoadPath{Mouse-10},'\Datasets\Datasets-5to0'));
    X=cellData_Raw_bin;
    cellCount=size(X,1);
    SessionStart=ones(1,length(Trials)+1);
    for i=2:length(SessionStart)
        SessionStart(i)=sum(SessionLength(1:max(Trials{i-1})))+1;
    end
    SessionStart(i+1)=sum(SessionLength);
    
    Dprime=zeros(cellCount,length(Trials));
    Dprime_Shuff=zeros(cellCount,length(Trials),1000);

    
    for T=1:length(Trials)
        T
        HInd=HitDataset{1,1};
        CInd=CRDataset{1,1,1};
    
        HInd=HInd(HInd>=SessionStart(T) & HInd< SessionStart(T+1));
        CInd=CInd(CInd>=SessionStart(T) & CInd< SessionStart(T+1));
        
        
        Dprime(:,T)=(mean(X(:,HInd),2)-mean(X(:,CInd),2)) ./ sqrt(0.5*(var(X(:,HInd)')+var(X(:,CInd)')))';
    
        for r=1:1000
            Ind=[HInd,CInd];
            Ind=Ind(randperm(length(Ind)));
            HInd=Ind(1:length(HInd));
            CInd=Ind(1+length(HInd):end);
            
            Dprime_Shuff(:,T,r)=(mean(X(:,HInd),2)-mean(X(:,CInd),2)) ./ sqrt(0.5*(var(X(:,HInd)')+var(X(:,CInd)')))';

        end 
    end
    
    save(strcat('H:\Reports2\2019_08_07_SingleCellCoding\HitCR_Raw_bin\Mouse',num2str(Mouse-60)),'Dprime','Dprime_Shuff','-v7.3') 
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Visualiation
pvalue=0.01;
d=3;


nM=1;
DPDiff={};
buf={};
bufAct={};


DPDiff{9}=[];
buf{9}=[];
bufAct{9}=[];

set(0,'defaultlinelinewidth',1.5);
set(0,'DefaultAxesFontSize',14);
set(0,'defaultAxesFontName','Calibri');
for M=[1,7]
  switch M
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
      
load(strcat('H:\Reports2\2019_08_07_SingleCellCoding\HitCR_Raw\Mouse',num2str(M)));
load(strcat(LoadPath{M+50},'\cellData.mat'));
load(strcat(LoadPath{M+50},'\SessionLength.mat'));
    SessionStart=ones(1,length(Trials)+1);
    for i=2:length(Trials)
        SessionStart(i)=sum(SessionLength(1:max(Trials{i-1})))+1;
    end
    SessionStart(i+1)=sum(SessionLength);

% ActiveCell=zeros(length(Trials));
% for t1=1:length(Trials)
%    x= cellData_Raw(:,SessionStart(t1):SessionStart(t1+1));
%    x=remove_baseline(x',1)';
%    for cn=1:size(x,1)
%       no=x(cn,:);
%       no=no(no<0);
%       no=[no,(-1)*no];
%       sd=sqrt(var(no'));
%       if sum(x(cn,:)> 4*sd) > (size(x,2)*0.001)
%           ActiveCell(cn,t1)=1;
%       end
%    
%    end
%     
% end
% save(strcat('H:\Reports2\2019_08_07_SingleCellCoding\HitCR_Raw\ActiveMouse',num2str(M)),'ActiveCell');
load(strcat('H:\Reports2\2019_08_07_SingleCellCoding\HitCR_Raw\ActiveMouse',num2str(M),'.mat'));


 switch M
        case 1
            Trials=1:5; % M1
        case 2
            Trials=1:5; %6 M2
        case 3
            Trials=1:5; %M3
        case 4
            Trials=1:5; % M4
        case 5
            Trials=[1:4,6]; %7 M5
        case 6
            Trials=1:3; % M6, there are too few correctly performed trials for session 4 and 5  
        case 7
            Trials=1:5; % 6,7 M7
    end
    


Dprime=abs(Dprime(:,Trials));
ActiveDays=sum(ActiveCell(:,Trials),2);
Dprime_Shuff=sort(abs(Dprime_Shuff(:,Trials,:)),3,'descend');

PTh=squeeze(Dprime_Shuff(:,:,round(pvalue*1000)));
DP_sig=Dprime>PTh;



NumDayCoding=sum(DP_sig,2);

cellMap=zeros(1017,1017);
for cnum=find(NumDayCoding>0)'
cellMap(max(cellIJ(cnum,1)-d,1):min(cellIJ(cnum,1)+d,1017),max(cellIJ(cnum,2)-d,1):min(cellIJ(cnum,2)+d,1017))=NumDayCoding(cnum);
end


for a=1:8
    Mask=max(Dprime,[],2)>0.5 & CortexArea==a;
    DPDiff{a}=[DPDiff{a};abs(max(Dprime(Mask,:),[],2)-min(Dprime(Mask,:),[],2))];
    
    buf{a}=[buf{a};NumDayCoding(NumDayCoding>0 & CortexArea==a)];
    bufAct{a}=[bufAct{a};NumDayCoding(NumDayCoding>0 & CortexArea==a & ActiveDays==5)];
end
a=9;
    Mask=max(Dprime,[],2)>0.5;
    DPDiff{a}=[DPDiff{a};abs(max(Dprime(Mask,:),[],2)-min(Dprime(Mask,:),[],2))];

    buf{a}=[buf{a};NumDayCoding(NumDayCoding>0)];
    bufAct{a}=[bufAct{a};NumDayCoding(NumDayCoding>0 & ActiveDays==5)];

% figure();
% hist(NumDayCoding,max(NumDayCoding));
% xlim([1 max(NumDayCoding)]);
% title(strcat('Mouse ',num2str(nM),' histogram'));
% xlabel('Number of coding days');
% ylabel('Number of cells');


figure();
imagesc(cellMap);
colormap('jet')
title(strcat('Mouse ',num2str(nM)));
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
nM=nM+1;
end

for a=1:9;
figure();hold on
hist(buf{a},max(buf{a}));
hist(bufAct{a},max(buf{a}));
xlim([1 max(buf{a})]);
title(strcat('Coding cells histogram area',num2str(a)));
xlabel('Number of coding days');
ylabel('Number of cells');

end
%%%%%%%%%%%%%%%%%%Num Coding Across Days
pvalue=0.01;
TotalNumCoding={};
figure();hold on;
for M=[1,3:7]
    
 switch M
        case 1
            Trials=1:5; % M1
        case 2
            Trials=1:5; %6 M2
        case 3
            Trials=1:5; %M3
        case 4
            Trials=1:5; % M4
        case 5
            Trials=[1:4,6]; %7 M5
        case 6
            Trials=1:3; % M6, there are too few correctly performed trials for session 4 and 5  
        case 7
            Trials=1:7; % 6,7 M7
    end
    
    load(strcat('H:\Reports2\2019_08_07_SingleCellCoding\HitCR_Raw\Mouse',num2str(M)));
    
Dprime=abs(Dprime(:,Trials));
Dprime_Shuff=sort(abs(Dprime_Shuff(:,Trials,:)),3,'descend');

PTh=squeeze(Dprime_Shuff(:,:,round(pvalue*1000)));
DP_sig=Dprime>PTh;
plot(sum(DP_sig)/mean(sum(DP_sig)));
    
end






%%%%%%%%%%%%%%%%%%%%%%% Visualization d' and intensity relation
SelCells={};
for M=[7]
     switch M
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
    
    
    load(strcat('H:\Reports2\2019_08_07_SingleCellCoding\HitCR_Raw\ActiveMouse',num2str(M),'.mat'));
    load(strcat('H:\Reports2\2019_08_07_SingleCellCoding\HitCR_Raw_bin\Mouse',num2str(M)));
    load(strcat(LoadPath{M+50},'\cellData.mat'));
    load(strcat(LoadPath{M+50},'\SessionLength.mat'));
    
    SelCells{M}=zeros(size(Dprime,1),1);
    ND = length(Trials);
    
    SessionStart=ones(1,ND+1);
    for i=2:length(Trials)
        SessionStart(i)=sum(SessionLength(1:max(Trials{i-1})))+1;
    end
    SessionStart(i+1)=sum(SessionLength);
    
    
    figure();hold on
    for cnum=1:size(Dprime,1)
       intensity = zeros(1,ND);
       for i=1:ND
           intensity(i)=max(cellData_Raw(cnum,SessionStart(i):SessionStart(i+1)));
       end
       
       mask=ActiveCell(cnum,:)==1;
       if sum(mask)>1
           b1=intensity(mask);
           b2=abs(Dprime(cnum,mask));
           count=1;
           for d1=1:sum(mask)
               for d2=d1:sum(mask)
                   dint=b1(d1)-b1(d2);
                   ddp=b2(d1)-b2(d2);
                   
                   if sign(dint*ddp)==-1
                      if abs(dint)>1 & abs(ddp)>0.2
                          SelCells{M}(cnum)=1;
                          break;
                      end
                   end
               end
               if SelCells{M}(cnum)==1
                   break;
               end
           end
           
       end
      if SelCells{M}(cnum)==1
          [sint,ind]=sort(intensity);
          buf=abs(Dprime(cnum,ind));
        plot(sint(mask(ind)),buf(mask(ind)))
      end  
    end
    
    
    
end



%%%%%%%%%%%%%%%%%%%% dprime for half sessions
AddressSetup;

for Mouse=[61,63:67]
   
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
    
    
    load(strcat(LoadPath{Mouse-10},'\cellData_ZS.mat'));
    load(strcat(LoadPath{Mouse-10},'\cellData.mat'));
    load(strcat(LoadPath{Mouse-10},'\SessionLength.mat'));
    load(strcat(LoadPath{Mouse-10},'\Datasets\Datasets-5to0'));
    X=cellData_Raw;
    cellCount=size(X,1);
      ND = length(Trials);
    
    SessionStart=ones(1,ND+1);
    for i=2:length(Trials)
        SessionStart(i)=sum(SessionLength(1:max(Trials{i-1})))+1;
    end
    SessionStart(i+1)=sum(SessionLength);
    
    
    Dprime=zeros(cellCount,length(Trials),2);
    Dprime_Shuff=zeros(cellCount,length(Trials),2,1000);

    
    for T=1:length(Trials)
        T
        for half=1:2
        HInd=HitDataset{1,1};
        CInd=CRDataset{1,1,1};
        SO=round((half-1)*(SessionStart(T+1)-SessionStart(T)) * 0.5);
        EO=round((half-2)*(SessionStart(T+1)-SessionStart(T)) * 0.5);
        HInd=HInd(HInd>=(SessionStart(T)+SO) & HInd< (SessionStart(T+1)+EO));
        CInd=CInd(CInd>=(SessionStart(T)+SO) & CInd< (SessionStart(T+1)+EO));
        
        
        Dprime(:,T,half)=(mean(X(:,HInd),2)-mean(X(:,CInd),2)) ./ sqrt(0.5*(var(X(:,HInd)')+var(X(:,CInd)')))';
    
        for r=1:1000
            Ind=[HInd,CInd];
            Ind=Ind(randperm(length(Ind)));
            HInd=Ind(1:length(HInd));
            CInd=Ind(1+length(HInd):end);
            
            Dprime_Shuff(:,T,half,r)=(mean(X(:,HInd),2)-mean(X(:,CInd),2)) ./ sqrt(0.5*(var(X(:,HInd)')+var(X(:,CInd)')))';

        end 
        end
    end
    
    save(strcat('H:\Reports2\2019_08_07_SingleCellCoding\HitCR_halfs_Raw\Mouse',num2str(Mouse-60)),'Dprime','Dprime_Shuff','-v7.3') 
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Visualiation for half session
pvalue=0.002;
d=2;
MeanDP=[];
SDDP=[];
HDPDiff=[];
nM=1;
buf=[];
halfTO={};
halfTO{7,5}=[];

DPTH=[0.1,0.3,0.5,0.7,0.9];

set(0,'defaultlinelinewidth',1.5);
set(0,'DefaultAxesFontSize',14);
set(0,'defaultAxesFontName','Calibri');
for M=[1,3:7]
     switch M
        case 1
            Trials=1:5; % M1
        case 2
            Trials=1:6; % M2
        case 3
            Trials=1:5; %M3
        case 4
            Trials=1:5; % M4
        case 5
            Trials=[1:4,6,7]; % M5
        case 6
            Trials=1:3; % M6, there are too few correctly performed trials for session 4 and 5  
        case 7
            Trials=1:7; % M7
    end
    
load(strcat('H:\Reports2\2019_08_07_SingleCellCoding\HitCR_halfs_Raw\Mouse',num2str(M)));
load(strcat(LoadPath{M+50},'\cellData.mat'));

Dprime=abs(Dprime(:,Trials,:));
Dprime_Shuff=sort(abs(Dprime_Shuff(:,Trials,:,:)),4,'descend');

PTh=squeeze(Dprime_Shuff(:,:,:,round(pvalue*1000)));
DP_sig=Dprime>PTh;
for t1=1:length(Trials)
    
    for th=1:5
        maskth=max(Dprime(:,t1,:),[],3)>DPTH(th);
        halfTO{M,th}=[halfTO{M,th},1-(sum(min(DP_sig(maskth,t1,:),[],3))./sum(max(DP_sig(maskth,t1,:),[],3)))];
    end

   Mask=max(Dprime(:,t1,:),[],3)>0.5;
   HDPDiff=[HDPDiff;abs(Dprime(Mask,t1,1)-Dprime(Mask,t1,2))];
end


nM=nM+1;
end


figure();hold on
for th=1:5
    bufhalf=[]
    for M=[1,3:7]
       plot(DPTH(th)*ones(1,length( halfTO{M,th})),halfTO{M,th},'b.')
       bufhalf=[bufhalf,halfTO{M,th}];
    end
    plot(DPTH(th),mean(bufhalf),'r+')
end
%%%%%%%%%%%%%%%%%%%%%%%%%Correlation of half session turn over with between
%%%%%%%%%%%%%%%%%%%%%%%%%session turn over
HTO={};
STO={};
pvalue=0.01;

PS=zeros(6,1);
PSgH=zeros(6,1);
nM=1;
for M=[1,3:7]
     switch M
        case 1
            Trials=1:5; % M1
        case 2
            Trials=1:6; % M2
        case 3
            Trials=1:5; %M3
        case 4
            Trials=1:5; % M4
        case 5
            Trials=[1:4,6,7]; % M5
        case 6
            Trials=1:3; % M6, there are too few correctly performed trials for session 4 and 5  
        case 7
            Trials=1:7; % M7
    end
    
load(strcat('H:\Reports2\2019_08_07_SingleCellCoding\HitCR_halfs_Raw\Mouse',num2str(M)));
Dprime=abs(Dprime(:,Trials,:));
Dprime_Shuff=sort(abs(Dprime_Shuff(:,Trials,:,:)),4,'descend');

PTh=squeeze(Dprime_Shuff(:,:,:,round(pvalue*1000)));
DP_sig=Dprime>PTh;

HTO{M} = zeros(size(Dprime,1),1);

for i=1:length(Trials)
    HTO{M}=HTO{M} | ((min(DP_sig(:,i,:),[],3) ~= max(DP_sig(:,i,:),[],3)) & max(Dprime(:,i,:),[],3)>0.5) ;
end

load(strcat('H:\Reports2\2019_08_07_SingleCellCoding\HitCR_Raw\Mouse',num2str(M)));
Dprime=abs(Dprime(:,Trials));
Dprime_Shuff=sort(abs(Dprime_Shuff(:,Trials,:)),3,'descend');

PTh=squeeze(Dprime_Shuff(:,:,round(pvalue*1000)));
DP_sig=Dprime>PTh;

STO{M} = zeros(size(Dprime,1),1);
STO{M}=((min(DP_sig(:,:),[],2) ~= max(DP_sig(:,:),[],2)) & max(Dprime(:,:),[],2)>0.5) ;

PS(nM)=sum(STO{M})/length(STO{M});
PSgH(nM)= sum(STO{M} & HTO{M}) / sum (HTO{M});


nM=nM+1;
end


figure();hold on
plot(ones(6,1),PSgH./PS,'.')
plot(1,mean(PSgH./PS),'r+')
ylabel('PSgH / PS')

xlim([0,2])
ylim([0,6])






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Single cell visualizations- session to
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% session
%%%%%%%% Mouse 1
3001 -> active, alinged, dprime   /   closeby cells
3181 -> non-active, aligned, *dprime / big vessel 
4451 -> active, dprime / closeby cells, alignment issues 
4466 -> active, dprime, aligned / closeby cells
4563 -> no good 
4592 -> active, dprime / closeby cells, alignment issue
4666 -> active , dprime, alignment / closeby cells

The best active: 1080 (MV) ,2984,3148
The best non-active: 4638,
%%%%%%%%%Mouse 2

The best non-active: 625 (PPC) ,1951,2248,541, 
%%%%%%%%%Mouse 6
686 -> in-active/ dprime, closeby
895 -> in active/drpime, closeby cells
1094 -> in active , dprime / closeby cells
**1124 -> Active, dprime
**2434,1072,1965,1455,2072,2152,2753,3377 -> Active, dprime
**793,726 -> Active , dpirme, /closeby cells
*719,2469, 1258,53 -> in-Active, dprime
*1404,996 -> in-active, dprime, closeby cells
1462 ->  in-active , dprime
**1861 -> in active ,dprime,closeby cells, 
The best active: 2434,1462
The best non-active: 3225,3378
 



Dprime=abs(Dprime);
candidates=find(max(Dprime,[],2) > 0.5 & min(Dprime,[],2) < 0.1 & ActiveDays<5);



cnum=625;


w=20;
for i=[1,4,7,10,13]
    X=squeeze(sum(RegisteredFirstFrame(cellIJ(cnum,1)-w:cellIJ(cnum,1)+w,cellIJ(cnum,2)-w:cellIJ(cnum,2)+w,i:i+2),3));

    
    figure();
    colormap('gray');
  imagesc(X);

end
figure();hold on;
plot((1:size(cellData_Raw,2))/600,cellData_Raw(cnum,:))
plot(ones(2,1)*SessionStart(1:7)/600,[zeros(1,7);10*ones(1,7)],'k--')
ylim([0 max(cellData_Raw(cnum,:))+1])
title(num2str(Dprime(cnum,:)))
abs(Dprime(cnum,:))
ActiveCell(cnum,:)


buf=cellImage{cnum};
figure();
imagesc(buf(max(cellIJ(cnum,1)-TileTopLeft(cnum,1)-w,1):min(cellIJ(cnum,1)-TileTopLeft(cnum,1)+w,250),max(cellIJ(cnum,2)-TileTopLeft(cnum,2)-w,1):cellIJ(cnum,2)-TileTopLeft(cnum,2)+w));

%%%%%%%%%%%%%%%%%%%%%%half session turn over
%%%%%%%%%%%%%%%%Mouse 6 
Session1: 
3276 -> 1 inactive (LV) 
session2:
1509 -> 1, active first half coding during stimulus and then during delay
2115 -> 1, active (V1)
Session5:
1310 -> 2, inactive (not completely)
1321 -> 2, inactive


set(0,'defaultlinelinewidth',1);
set(0,'DefaultAxesFontSize',14);
set(0,'defaultAxesFontName','Calibri');

Candidates=max(abs(Dprime),[],3)>0.8 & min(abs(Dprime),[],3)<0.1;
sum(Candidates)
find(Candidates(:,1)>0)

cnum=2115;Sess=2;
ind=SessionStart(Sess):SessionStart(Sess+1);

figure();hold on;
plot((1:length(ind))/10,cellData_Raw(cnum,ind),'b')
title(num2str(squeeze(abs(Dprime(cnum,Sess,:)))'))
xlabel('Time (s)');
ylabel('Ca2+ Zscore');

w=32;
figure();hold on;

buf=cellImage{cnum};
imagesc(buf(max(cellIJ(cnum,1)-TileTopLeft(cnum,1)-w,1):min(cellIJ(cnum,1)-TileTopLeft(cnum,1)+w,250),max(cellIJ(cnum,2)-TileTopLeft(cnum,2)-w,1):cellIJ(cnum,2)-TileTopLeft(cnum,2)+w));

wt=40
[~,tes]=sort(cellData_Raw(cnum,median(ind):max(ind)),'descend')
test=tes(2)+round(length(ind)/2)

tind=8877+SessionStart(Sess)-1;
tind=tind-wt:tind+wt;
hdf5In='T:\L368_Concat\dfof_single.h5';
Mov=single(readHDF5Subset(hdf5In,[cellIJ(cnum,1)-w cellIJ(cnum,2)-w min(tind)-1 ],[2*w,2*w,2*wt]));

Mov=Mov-min(min(min(Mov)));
Mov=Mov/max(max(max(Mov)));
implay(reshape(Mov,[2*w 2*w 1 2*wt]))

% outputfilename=strcat('C:\Users\Sadegh\Documents\VLMReborn\NCPaper\SienceDraft1_Supplementary\SingleCells\EM_M7_C',num2str(cnum),'_S',num2str(Sess),'H2');
% newVid = VideoWriter(outputfilename, 'MPEG-4'); % New
% newVid.FrameRate = 20;
% newVid.Quality = 100;
% open(newVid);
% 
% writeVideo(newVid,reshape(Mov,[2*w 2*w 1 2*wt]))

% close(newVid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


hind=find(diff(Hit(tind)==1));
cind=find(diff(CR(tind)==1));
RH=[];
RC=[];
for i=1:length(hind)
    RH=[RH;cellData_Raw(1509,tind(hind(i)):tind(hind(i))+35)];
    RC=[RC;cellData_Raw(1509,tind(cind(i)):tind(cind(i))+35)];
end
    
figure();
imagesc([RH;-1*ones(20,36);RC]);
figure();
imagesc(RC);