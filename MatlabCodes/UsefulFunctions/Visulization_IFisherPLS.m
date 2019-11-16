
ResultsPath='C:\Users\Sadegh\Documents\VLMReborn\Reports2\2018_01_08_FisherInfo_ConcatDays_SpeedBal_EqualSets\Raw\Datasets-10to0\HitMiss\';
ResultsPath2='C:\Users\Sadegh\Documents\VLMReborn\Reports2\2018_01_08_FisherInfo_ConcatDays_SpeedBal_EqualSets\Raw\Datasets-10to0\FACR\';
mode=2;
TSize=zeros(7,1);
Sessions=[51:54,57];
Sessions2=[51:54,57];
areaVec=9;
LowDimVec=1:50;
if mode==1
    time=-0.4:0.1:2;
elseif mode==2
    time=-0.4:0.1:0.5;
elseif mode==3
    time=-0.4:0.1:3;
end



NT=length(time);
NL=length(LowDimVec);

c1={'b','--b','-b*','-.b',':b','--bo','--bs','--m'};
r1={'r','--r','-r*','-.r',':r','--ro','--rs','--m'};
color={'b','b--','k--','r','r--','m','k','g','m--'};
TitleV={'V1','LV','MV','PTLP','A','S','M','RSC','All'};
set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',16);
set(0,'defaultAxesFontName','Calibri');


ND=length(Sessions);
Info_Val_Fisher=zeros(areaVec,NL,ND);
MaxInfoDim=zeros(areaVec,ND);
MaxInfo=zeros(areaVec,ND);
MaxInfoMask=zeros(areaVec,ND);


for s=Sessions
    j=s-50;
    V1=open(strcat(ResultsPath,'Session',num2str(s),'_mode',num2str(mode),'.mat'));

    figure();
    
    maxDim=20;
    hold on;

    title(strcat('M',int2str(j)))
    grid on
    
    for area=1:areaVec
        Info_Val_Fisher(area,:,j)=mean(V1.Info_Val_Fisher(area,:,:),3);
        [MaxInfo(area,j),MaxInfoDim(area,j)]=max(squeeze(Info_Val_Fisher(area,:,j)));
        if MaxInfo(area,j)>0%max(V1.Info_Val_Fisher_Sh(area,MaxInfoDim(area,j),:))
            MaxInfoMask(area,j)=1;
        end

        plot(Info_Val_Fisher(area,:,j),color{area});

    end
    
    
end


Res=MaxInfo';
for a=1:areaVec
    Res(8,a)=mean(MaxInfo(a,(MaxInfoMask(a,:)~=0)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ND=length(Sessions2);
Info_Val_Fisher2=zeros(areaVec,NL,ND);
MaxInfoDim=zeros(areaVec,ND);
MaxInfo=zeros(areaVec,ND);
MaxInfoMask2=zeros(areaVec,ND);

for s=Sessions2
    j=s-50;
    V1=open(strcat(ResultsPath2,'Session',num2str(s),'_mode',num2str(mode),'.mat'));

    figure();
    
    maxDim=20;
    hold on;

    title(strcat('M',int2str(j)))
    grid on
    
    for area=1:areaVec
        Info_Val_Fisher2(area,:,j)=mean(V1.Info_Val_Fisher(area,:,:),3);
        [MaxInfo(area,j),MaxInfoDim(area,j)]=max(squeeze(Info_Val_Fisher2(area,:,j)));
        if MaxInfo(area,j)>0%max(V1.Info_Val_Fisher_Sh(area,MaxInfoDim(area,j),:))
            MaxInfoMask2(area,j)=1;
        end

        plot(Info_Val_Fisher2(area,:,j),color{area});

    end
    
    
end


Res2=MaxInfo';
for a=1:areaVec
    Res2(8,a)=mean(MaxInfo(a,(MaxInfoMask2(a,:)~=0)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F=zeros(areaVec,2);
for a=1:areaVec
    mask= MaxInfoMask(a,:) | MaxInfoMask(a,:) ;
    F(a,2)=sum(mask);
    F(a,1)= signrank(Res([mask,false],a),Res2([mask,false],a));
end
    





