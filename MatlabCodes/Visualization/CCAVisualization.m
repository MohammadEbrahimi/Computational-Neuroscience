%%%%%%for mixed stimuli every thing is stored in stimuli 0 -> CR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Areas={'V1','LV','MV','PPC','A','S','M','RSC'};
buf={};
TimeP='-15to20';

for Mouse=[61,63:64]


    loadpath=LoadPath{Mouse-10};
    load(strcat('E:\Reports2\2018_03_14_CannonicalCorrelations\ROI_Bin\Datasets-',TimeP,'\HitCR_L2_MixedStim\Session',num2str(Mouse),'_mode1'));
    load(strcat(loadpath,'/cellData.mat'));
    load(strcat(loadpath,'/cellData_ZS.mat'));
    load(strcat(loadpath,'/SessionLength.mat'));
    load(strcat(loadpath,'/Datasets/Datasets--20to15.mat'));
    cellEvents=cellData_ROI_bin;
    NRepeate=size(Wx0,3);
    
TimeVec=1:25;
    CMtx_avg=zeros(8,8,length(TimeVec),2,20); 
CMtx=zeros(8,8,length(TimeVec),2,20,NRepeate);
dprime=zeros(8,8,length(TimeVec),2,20,NRepeate);
dprime_avg=zeros(8,8,length(TimeVec),2,20);
dprime_sh=zeros(8,8,length(TimeVec),2,2,20,NRepeate);
meanModeProj=zeros(8,8,length(TimeVec),2,2,20);
sdModeProj=zeros(8,8,length(TimeVec),2,2,20);
    
    
    MW1={};MW2={};
    for a1=1:8%x
        for a2=a1+1:8%y
            if sum(CortexArea==a1)>5 && sum(CortexArea==a2)>5

                W1=zeros(size(Wx0{a1,a2,1},1),20,NRepeate);
                W2=zeros(size(Wy0{a1,a2,1},1),20,NRepeate);
                for nr=1:NRepeate
                    W1(:,:,nr)= (Wx0{a1,a2,nr}.* repmat(sign(mean(Wx0{a1,a2,nr})),[size(Wx0{a1,a2,nr},1),1]));
                    W2(:,:,nr)= (Wy0{a1,a2,nr}.* repmat(sign(mean(Wy0{a1,a2,nr})),[size(Wy0{a1,a2,nr},1),1]));
                end
                
                MW1{a1,a2}=squeeze(mean(W1,3));
                MW2{a1,a2}=squeeze(mean(W2,3));
                
                
                HitTrialMtx=zeros(1,length(CortexArea),length(TimeVec));
                n=1;
                index=HitDataset{1,1};
                TNum=HitTrialNumber{1,1};
                
                for i=unique(TNum)
                    if sum(TNum==i)==length(TimeVec)
                        HitTrialMtx(n,:,TimeVec)= cellEvents(:,index(TNum==i));
                        n=n+1;
                    end
                end
                
                CRTrialMtx=zeros(1,length(CortexArea),length(TimeVec));
                n=1;
                index=CRDataset{1,1};
                TNum=CRTrialNumber{1,1};
                
                for i=unique(TNum)
                    if sum(TNum==i)==length(TimeVec)
                        CRTrialMtx(n,:,TimeVec)= cellEvents(:,index(TNum==i));
                        n=n+1;
                    end
                end
                
                
                
                for ccnum=1:20
                    ccnum

                    D1=squeeze(MW1{a1,a2}(:,ccnum)/norm(MW1{a1,a2}(:,ccnum)));
                    D2=squeeze(MW2{a1,a2}(:,ccnum)/norm(MW2{a1,a2}(:,ccnum)));
                    
                    
                    HitProj=zeros(size(HitTrialMtx,1),length(TimeVec),2);
                    CRProj=zeros(size(CRTrialMtx,1),length(TimeVec),2);
                    for td=1:length(TimeVec)
                        

                        HitProj(:,td,1)=D1'*squeeze(HitTrialMtx(:,CortexArea==a1,td))';
                        HitProj(:,td,2)=D2'*squeeze(HitTrialMtx(:,CortexArea==a2,td))';
                        
                        
                        CRProj(:,td,1)=D1'*squeeze(CRTrialMtx(:,CortexArea==a1,td))';
                        CRProj(:,td,2)=D2'*squeeze(CRTrialMtx(:,CortexArea==a2,td))';

                        
                        meanModeProj(a1,a2,td,1,1,ccnum)=mean(CRProj(:,td,1));
                        sdModeProj(a1,a2,td,1,1,ccnum)=sqrt(var(squeeze(CRProj(:,td,1))));
                        
                        meanModeProj(a1,a2,td,2,1,ccnum)=mean(CRProj(:,td,2));
                        sdModeProj(a1,a2,td,2,1,ccnum)=sqrt(var(squeeze(CRProj(:,td,2))));
                        
                        meanModeProj(a1,a2,td,1,2,ccnum)=mean(HitProj(:,td,1));
                        sdModeProj(a1,a2,td,1,2,ccnum)=sqrt(var(squeeze(HitProj(:,td,1))));
                        
                        meanModeProj(a1,a2,td,2,2,ccnum)=mean(HitProj(:,td,2));
                        sdModeProj(a1,a2,td,2,2,ccnum)=sqrt(var(squeeze(HitProj(:,td,2))));
                        
                        
                        
                        
                        cval=corrcoef(CRProj(:,td,1),CRProj(:,td,2));
                        CMtx_avg(a1,a2,td,1,ccnum)=cval(1,2);
                        
                         cval=corrcoef(HitProj(:,td,1),HitProj(:,td,2));
                         CMtx_avg(a1,a2,td,2,ccnum)=cval(1,2);
                    end
                    CMtx_avg(a1,a2,:,1,ccnum)=CMtx_avg(a1,a2,:,1,ccnum)*sign(mean(CMtx_avg(a1,a2,:,1,ccnum)));
                    CMtx_avg(a1,a2,:,2,ccnum)=CMtx_avg(a1,a2,:,2,ccnum)*sign(mean(CMtx_avg(a1,a2,:,2,ccnum)));
                    %%%%%%%%%%dprime_avg(a1,a2,td,a1/a2,Hit/CR,ccnum)
                    %%%%%%%%%% a1
                    dprime_avg(a1,a2,:,1,ccnum)=((mean(HitProj(:,:,1))-mean(CRProj(:,:,1))).^2) ./  (0.5*(var(squeeze(HitProj(:,:,1)))+var(squeeze(CRProj(:,:,1)))));%%%%CR vector on Hit/CR
                    %%%%%%%%%% a2
                    dprime_avg(a1,a2,:,2,ccnum)=((mean(HitProj(:,:,2))-mean(CRProj(:,:,2))).^2) ./  (0.5*(var(squeeze(HitProj(:,:,2)))+var(squeeze(CRProj(:,:,2)))));%%%%CR vector on Hit/CR
                    
                    
                    for sh=1:NRepeate
                        
                        
                        D1=squeeze(W1(:,ccnum,sh)/norm(W1(:,ccnum,sh)));
                        D2=squeeze(W2(:,ccnum,sh)/norm(W2(:,ccnum,sh)));
                        
                        
                        HitProj=zeros(size(HitTrialMtx,1),length(TimeVec),2);
                        CRProj=zeros(size(CRTrialMtx,1),length(TimeVec),2);
                        for td=1:length(TimeVec)
                            
                            HitProj(:,td,1)=D1'*squeeze(HitTrialMtx(:,CortexArea==a1,td))';
                            HitProj(:,td,2)=D2'*squeeze(HitTrialMtx(:,CortexArea==a2,td))';
                            
                            
                            CRProj(:,td,1)=D1'*squeeze(CRTrialMtx(:,CortexArea==a1,td))';
                            CRProj(:,td,2)=D2'*squeeze(CRTrialMtx(:,CortexArea==a2,td))';

                            
                            
                            
                            cval=corrcoef(CRProj(:,td,1),CRProj(:,td,2));
                            CMtx(a1,a2,td,1,ccnum,sh)=cval(1,2);
                            
                             cval=corrcoef(HitProj(:,td,1),HitProj(:,td,2));
                             CMtx(a1,a2,td,2,ccnum,sh)=cval(1,2);
                        end
                        CMtx(a1,a2,:,1,ccnum,sh)=CMtx(a1,a2,:,1,ccnum,sh)*sign(max(CMtx(a1,a2,:,1,ccnum,sh)));
                        CMtx(a1,a2,:,2,ccnum,sh)=CMtx(a1,a2,:,2,ccnum,sh)*sign(max(CMtx(a1,a2,:,2,ccnum,sh)));
                        %%%%%%%%%%dprime(a1,a2,td,a1/a2,Hit/CR,ccnum,Repeate)
                        %%%%%%%%%% a1
                        dprime(a1,a2,:,1,ccnum,sh)=((mean(HitProj(:,:,1))-mean(CRProj(:,:,1))).^2) ./  (0.5*(var(squeeze(HitProj(:,:,1)))+var(squeeze(CRProj(:,:,1)))));%%%%CR vector on Hit/CR
                        %%%%%%%%%% a2
                        dprime(a1,a2,:,2,ccnum,sh)=((mean(HitProj(:,:,2))-mean(CRProj(:,:,2))).^2) ./  (0.5*(var(squeeze(HitProj(:,:,2)))+var(squeeze(CRProj(:,:,2)))));%%%%CR vector on Hit/CR
                        
                        
                        
                      
                        D1=D1(randperm(size(D1,1)));
                        D2=D2(randperm(size(D2,1)));
                        
                        
                        
                        HitProj=zeros(size(HitTrialMtx,1),length(TimeVec),2);
                        CRProj=zeros(size(CRTrialMtx,1),length(TimeVec),2);
                        for td=1:length(TimeVec)
                            
                            HitProj(:,td,1)=D1'*squeeze(HitTrialMtx(:,CortexArea==a1,td))';
                            HitProj(:,td,2)=D2'*squeeze(HitTrialMtx(:,CortexArea==a2,td))';

                            
                            
                            CRProj(:,td,1)=D1'*squeeze(CRTrialMtx(:,CortexArea==a1,td))';
                            CRProj(:,td,2)=D2'*squeeze(CRTrialMtx(:,CortexArea==a2,td))';
         
                        end
                        
                        %%%%%%%%%%dprime(a1,a2,td,a1/a2,Hit/CR,ccnum,Repeate)
                        %%%%%%%%%% a1
                        dprime_sh(a1,a2,:,1,ccnum,sh)=((mean(HitProj(:,:,1))-mean(CRProj(:,:,1))).^2) ./  (0.5*(var(squeeze(HitProj(:,:,1)))+var(squeeze(CRProj(:,:,1))))); %%%%Hit vector on Hit/CR
                        %%%%%%%%%% a2
                        dprime_sh(a1,a2,:,2,ccnum,sh)=((mean(HitProj(:,:,2))-mean(CRProj(:,:,2))).^2) ./  (0.5*(var(squeeze(HitProj(:,:,2)))+var(squeeze(CRProj(:,:,2)))));%%%%CR vector on Hit/CR
                        
                    end
                end
               
            end
        end
    end
        save(strcat('E:\Reports2\2018_03_14_CannonicalCorrelations\ROI_Bin\Datasets-',TimeP,'\HitCR_L2_MixedStim\Session',num2str(Mouse),'_mode1_Results'),...
            'dprime','dprime_avg','dprime_sh','CMtx','CMtx_avg','MW1','MW2','meanModeProj','sdModeProj','-v7.3')
end




%%%%%%%%%%%%%%%%%%%%%%%Draw Map for CCA vector
ccnum=1;Mouse=67;
load(strcat(LoadPath{Mouse-10},'\cellData.mat'));
load(strcat('E:\Reports2\2018_03_14_CannonicalCorrelations\ROI_Bin\Datasets--5to0\HitCR_L2_MixedStim\Session',num2str(Mouse),'_mode1_Results.mat'));

a1=1;a2=6;
Vect=zeros(size(CortexArea));
Vect(CortexArea==a1)=MW1C{a1,a2}(:,ccnum);
%th=1*sqrt(var(MW1C{a1,a2}(:,ccnum)));
buf=sort(abs(Vect(CortexArea==a1)),'descend');th=buf(20);
se = strel('disk',4,4);
cellMap1=zeros(1017,1017);

for cnum=1:length(Vect)
    x=zeros(1017,1017);
    x(cellIJ(cnum,1),cellIJ(cnum,2))=abs(Vect(cnum))>th;
    x=imdilate(x,se);
    x=x*Vect(cnum);
    cellMap1(x>0)=x(x>0);
end
Vect=zeros(size(CortexArea));
Vect(CortexArea==a2)=MW2C{a1,a2}(:,ccnum);
%th=1*sqrt(var(MW2C{a1,a2}(:,ccnum)));
buf=sort(abs(Vect(CortexArea==a2)),'descend');th=buf(20);
se = strel('disk',4,4);
cellMap2=zeros(1017,1017);

for cnum=1:length(Vect)
    x=zeros(1017,1017);
    x(cellIJ(cnum,1),cellIJ(cnum,2))=abs(Vect(cnum))>th;
    x=imdilate(x,se);
    x=x*Vect(cnum);
    cellMap2(x>0)=x(x>0);
end

%imagesc([cellMap1,max(max(cellMap1))*ones(1017,5),cellMap2])

imagesc([cellMap1+cellMap2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Average CCA
TimeVec=1:25;
MeanCCA=zeros(8,8,20,6);
TimeCCA=zeros(8,8,max(TimeVec),2,20,6);
ActivityCCA=zeros(8,8,max(TimeVec),2,2,20,6);

nM=0;
for Mouse=[61,63:67]
    nM=nM+1;
    load(strcat('E:\Reports2\2018_03_14_CannonicalCorrelations\ROI_Bin\Datasets--15to20\HitCR_L2_MixedStim\Session',num2str(Mouse),'_mode1.mat'));
    load(strcat('E:\Reports2\2018_03_14_CannonicalCorrelations\ROI_Bin\Datasets--15to20\HitCR_L2_MixedStim\Session',num2str(Mouse),'_mode1_Results.mat'));
    
    MeanCCA(:,:,:,nM)=squeeze(mean(CanCorr0_val,3));
    TimeCCA(:,:,:,1,:,nM)=squeeze(mean(CMtx(:,:,:,1,:,:),6));
    TimeCCA(:,:,:,2,:,nM)=squeeze(mean(CMtx(:,:,:,2,:,:),6));
    ActivityCCA(:,:,:,:,:,:,nM)=meanModeProj;
end

AverageAnimalCCA=zeros(8,8,20);
for cc=1:20
for a1=1:8
    for a2=a1+1:8
        
        Mask=squeeze(MeanCCA(a1,a2,cc,:)>0);
        AverageAnimalCCA(a1,a2,cc)=mean(MeanCCA(a1,a2,cc,Mask));
        AverageAnimalCCA(a2,a1,cc)=mean(MeanCCA(a1,a2,cc,Mask));
    end
end
end


AreaColors=double([0,0,128;0,191,255;100,149,237;151,13,255;0,250,154;230,20,60;255,192,185;255,150,50;0,0,0])/256;
set(0,'defaultlinelinewidth',1.5);
set(0,'DefaultAxesFontSize',14);
set(0,'defaultAxesFontName','Calibri');
%%%%%%%CCA dimensions
for ba1=1:8
    figure();hold on;
    xlabel('CCA Mode Number');
    ylabel('Correlation Coefficient');
    
    for ba2=1:8
        if ba1~=ba2
            a1=min(ba1,ba2);
            a2=max(ba1,ba2);
            
            mask=squeeze(MeanCCA(a1,a2,1,:))>0;
            
            
            hold on;
            shadedErrorBar(1:20,squeeze(mean(MeanCCA(a1,a2,:,mask),4))',sqrt(var(squeeze(MeanCCA(a1,a2,:,mask))')/6),{'Color',AreaColors(ba2,:),'Marker','.','MarkerSize',15},1);
        end
    end
end
legend('V1','LV','MV','PPC','A','S','M','RSC');
%%%%%%%%
%%%%%%%%CCA over Time
TimeCCA(isnan(TimeCCA))=0;
TimeVec1=-1.4:0.1:0.5;
TimeVec=-0.4:0.1:2;
ccn=1;
TimeCCAMax=max(TimeCCA05,TimeCCA_50);
TimeCCAMax(2,5,:,:,:,:)=TimeCCA_50(2,5,:,:,:,:);
for ba1=1:8
    figure();hold on;
    xlabel('Time(s)');
    ylabel('Correlation Coefficient');
    for ba2=1:8
            %title(strcat('A',num2str(ba1),',A',num2str(ba2),':CC',num2str(ccn)));
        if ba1~=ba2
            a1=min(ba1,ba2);
            a2=max(ba1,ba2);
            mask=squeeze(TimeCCA_1520(a1,a2,end,1,ccn,:))~=0;
            
            
            plot([0,0],[0,1],'Color',[0.6,0.6,0.6],'LineStyle','--');
            CH=1:2;
            shadedErrorBar(TimeVec1,squeeze(mean(mean(TimeCCA_1520(a1,a2,6:25,CH,ccn,mask),4),6))',sqrt(var(squeeze(mean(TimeCCA_1520(a1,a2,6:25,CH,ccn,mask),4))')/sum(mask*100)),{'Color',AreaColors(ba2,:),'LineStyle','--'},1);
           
           % mask=squeeze(TimeCCAMax(a1,a2,end,1,ccn,:))~=0;
%             CH=1;
%             shadedErrorBar(TimeVec,squeeze(mean(mean(TimeCCAMax(a1,a2,:,CH,ccn,mask),4),6))',sqrt(var(squeeze(mean(TimeCCAMax(a1,a2,:,CH,ccn,mask),4))')/sum(mask)),{'Color',AreaColors(ba2,:),'LineStyle','--'},1);
            CH=1:2;
            shadedErrorBar(TimeVec,squeeze(mean(mean(TimeCCAMax(a1,a2,:,CH,ccn,mask),4),6))',sqrt(var(squeeze(mean(TimeCCAMax(a1,a2,:,CH,ccn,mask),4))')/sum(mask*100)),{'Color',AreaColors(ba2,:)},1);
            ylim([0,0.8]);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pre/Post change of CCA
PreCCA=TimeCCA_1520;
PostCCA=max(TimeCCA05,TimeCCA_50);
PostCCA(2,5,:,:,:,:)=TimeCCA_50(2,5,:,:,:,:);
PreTime=5:20;
PostTime=8:13;

CCAChange=zeros(9,9,20,2);
TimeCCAMax=max(TimeCCA05,TimeCCA_50);
Free=0.05;

for CH=1:2;
for ccn=1:20;

for ba1=1:8
    
    for ba2=1:8
        if ba1~=ba2
            a1=min(ba1,ba2);
            a2=max(ba1,ba2);
            mask=squeeze(PreCCA(a1,a2,end,1,ccn,:))~=0;
            
            PreVec=squeeze(mean(mean(PreCCA(a1,a2,PreTime,:,ccn,mask),4),3));
            PreVecSD=sqrt(var(squeeze(mean(PreCCA(a1,a2,PreTime,:,ccn,mask),4))));
            PostVec=squeeze(mean(mean(PostCCA(a1,a2,PreTime,CH,ccn,mask),4),3));
            PostVecSD=sqrt(var(squeeze(mean(PostCCA(a1,a2,PreTime,CH,ccn,mask),4))));
          %pval=ttest2(PreVec,PostVec);
            if sum(PreVec>=PostVec-Free)>=sum(mask) || sum(PreVec-Free<=PostVec)>=sum(mask)
            
                if (mean(PostVec)>0.1 || mean(PreVec)>0.1) && abs(mean(PostVec)-mean(PreVec))>(Free+0.05)
                %if pval==1
                    CCAChange(ba1,ba2,ccn,CH)=mean(PostVec)-mean(PreVec);
                end
                
            end
        end
    end
end
end
end
CCAChange(9,:,:,:)=min(min(min(min(CCAChange))));
CCAChange(:,9,:,:)=min(min(min(min(CCAChange))));
figure();colormap('jet')
imagesc(reshape(permute(CCAChange,[1,3,2,4]),[9*20,9*2]))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Info In CC Modes
AreaColors=double([0,0,128;0,191,255;100,149,237;151,13,255;0,250,154;230,20,60;255,192,185;255,150,50;0,0,0])/256;
TimeVec=-0.4:0.1:2;
for Mouse=[61,63:67]
    load(strcat('E:\Reports2\2018_03_14_CannonicalCorrelations\ROI_Bin\Datasets-0to5\HitCR_L2_MixedStim\Session',num2str(Mouse),'_mode1_Results.mat'));
    NRepeate=size(dprime,6);
    ccnumVec=1:6;
    figure();hold on;
    adisp=1;
    a1p=1;
    for a2p=1:8
        a1=min(a1p,a2p);
        a2=max(a1p,a2p);
        if a1p<a2p
            adisp=1;
        else
            adisp=2;
        end
        for ccnum=ccnumVec
            subplot(length(ccnumVec),1,ccnum);hold on;
            %shadedErrorBar(-0.4:0.1:2,dprime_avg(a1,a2,:,1,1,ccnum),sqrt(var(squeeze(dprime(a1,a2,:,1,1,ccnum,:))')/100),{'Color',AreaColors(a2,:)},1);
            shadedErrorBar(TimeVec,mean(dprime(a1,a2,:,adisp,ccnum,:),6),sqrt(var(squeeze(dprime(a1,a2,:,adisp,ccnum,:))')/100),{'Color',AreaColors(a2p,:)},1);
            shadedErrorBar(TimeVec,mean(dprime_sh(a1,a2,:,adisp,ccnum,:),6),sqrt(var(squeeze(dprime_sh(a1,a2,:,adisp,ccnum,:))')/100),{'Color',[0.5,0.5,0.5],'LineStyle','--'},1);
            
        end
    end
end


NormalizedDprime=zeros(8,8,20,6);
Dprime=zeros(8,8,20,25,6);
nM=1;
for Mouse=[61,63:67]
    load(strcat('E:\Reports2\2018_03_14_CannonicalCorrelations\ROI_Bin\Datasets-0to5\HitCR_L2_MixedStim\Session',num2str(Mouse),'_mode1_Results.mat'));
    adisp=1;
    for a1p=1:8
        for a2p=1:8
            if a1p~=a2p
                a1=min(a1p,a2p);
                a2=max(a1p,a2p);
                if a1p<a2p
                    adisp=1;
                else
                    adisp=2;
                end
                for ccnum=1:20
                    NormalizedDprime(a1p,a2p,ccnum,nM)=max(mean(dprime(a1,a2,:,adisp,ccnum,:),6),[],3);
                    Dprime(a1p,a2p,ccnum,:,nM)=squeeze(mean(dprime(a1,a2,:,adisp,ccnum,:),6));
                end
                
            end
        end
    end
    nM=nM+1;
end
figure();
imagesc(reshape(permute(NormalizedDprime,[1,3,2,4]),[8*20,8*6]))


MeanNormalizedDprime=zeros(9,8,20);

for a1=1:8
    for a2=1:8
        if a1~=a2
            Mask=squeeze(NormalizedDprime(a1,a2,1,:)>0);
            for cc=1:20
                MeanNormalizedDprime(a1,a2,cc)=mean(NormalizedDprime(a1,a2,cc,Mask));
               % MeanDprime(a1,a2,ccnum,:)=mean(Dprime(a1,a2,cc,:,Mask));
            end
            
        end
    end
end
 MeanNormalizedDprime(9,:,:)=0;

figure();
imagesc(reshape(permute(MeanNormalizedDprime,[1,3,2]),[9*20,8]))

figure();
a1=1;
for cc=1:8
    for a2=2:8
        subplot(1,8,cc);hold on;
        xlim([-1 2])
        ylim([0 5])
        Mask=squeeze(NormalizedDprime(a1,a2,1,:)>0);
        shadedErrorBar(-0.4:0.1:2,mean(Dprime(a1,a2,cc,:,Mask),5),sqrt(var(squeeze(Dprime(a1,a2,cc,:,Mask))')/100),{'Color',AreaColors(a2,:)},1);
        
        
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pre/Post Stim Simillarity
ModeSimCC=zeros(8,8,20,5);
nM=1;
for Mouse=[61,63:67]
    Mouse
ModeSim=zeros(8,8,20,20);

    load(strcat('E:\Reports2\2018_03_14_CannonicalCorrelations\ROI_Bin\Datasets-0to5\HitCR_L2_MixedStim\Session',num2str(Mouse),'_mode1'));
    MW1_post={};MW2_post={};
 NRepeate=size(Wx0,3);
 for a1=1:8%x
        for a2=a1+1:8%y
            if length(Wx0{a1,a2,1})>1 
                W1=zeros(size(Wx0{a1,a2,1},1),20,NRepeate);
                W2=zeros(size(Wy0{a1,a2,1},1),20,NRepeate);
                for nr=1:NRepeate
                    W1(:,:,nr)= (Wx0{a1,a2,nr}.* repmat(sign(mean(Wx0{a1,a2,nr})),[size(Wx0{a1,a2,nr},1),1]));
                    W2(:,:,nr)= (Wy0{a1,a2,nr}.* repmat(sign(mean(Wy0{a1,a2,nr})),[size(Wy0{a1,a2,nr},1),1]));
                end
                MW1_post{a1,a2}=squeeze(mean(W1,3));
                MW2_post{a1,a2}=squeeze(mean(W2,3));
            end
        end
    end
    
    load(strcat('E:\Reports2\2018_03_14_CannonicalCorrelations\ROI_Bin\Datasets--15to20\HitCR_L2_MixedStim\Session',num2str(Mouse),'_mode1'));
    MW1_pre={};MW2_pre={};
    NRepeate=size(Wx0,3);
    for a1=1:8%x
 
        for a2=a1+1:8%y
  
            if length(Wx0{a1,a2,1})>1 

                W1=zeros(size(Wx0{a1,a2,1},1),20,NRepeate);
                W2=zeros(size(Wy0{a1,a2,1},1),20,NRepeate);
                for nr=1:NRepeate
                    W1(:,:,nr)= (Wx0{a1,a2,nr}.* repmat(sign(mean(Wx0{a1,a2,nr})),[size(Wx0{a1,a2,nr},1),1]));
                    W2(:,:,nr)= (Wy0{a1,a2,nr}.* repmat(sign(mean(Wy0{a1,a2,nr})),[size(Wy0{a1,a2,nr},1),1]));
                end
                
                MW1_pre{a1,a2}=squeeze(mean(W1,3));
                MW2_pre{a1,a2}=squeeze(mean(W2,3));
            end
        end
    end

    for cc1=1:20
        for cc2=1:20
            for a1=1:8
                for a2=1:8
                    if a1<a2 && length(Wx0{a1,a2,1})>1 
                        ModeSim(a1,a2,cc1,cc2)=similarity([MW1_pre{a1,a2}(:,cc1),MW1_post{a1,a2}(:,cc2)]);
                    elseif a2<a1 && length(Wx0{a2,a1,1})>1 
                        ModeSim(a1,a2,cc1,cc2)=similarity([MW2_pre{a2,a1}(:,cc1),MW2_post{a2,a1}(:,cc2)]);
                    end
                    
                end
            end
        end
    end

    for cc=1:20
        ModeSimCC(:,:,cc,nM)=abs(ModeSim(:,:,cc,cc));
    end
    
    ModeSim=permute(ModeSim,[1,3,2,4]);
     figure();
     imagesc(reshape(abs(ModeSim),[8*20,8*20]));
    nM=nM+1;
    
end



Res=zeros(8,9,20);

for a1=1:8
    for a2=1:8
        for cc=1:20
            Res(a1,a2,cc)=mean(ModeSimCC(a1,a2,cc,(ModeSimCC(a1,a2,cc,:)~=0)));
            Res(a1,9,cc)=0;
        end
    end
end

        
 figure()
    imagesc(reshape(abs(Res),[8,9*20])); 
                




%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Mode Simillarities
AllRes=[]
SimMtx=zeros(8,8,8,20,6);
nM=0;
for Mouse=[61,63:67]
    nM=nM+1;
    load(strcat('E:\Reports2\2018_03_14_CannonicalCorrelations\ROI_Bin\Datasets--15to20\HitCR_L2_MixedStim\Session',num2str(Mouse),'_mode1_Results.mat'));
    load(strcat('E:\Reports2\2018_03_01_FisherInfo\ROI_Bin\Datasets--5to0\HitCR\Session',num2str(Mouse),'_mode1.mat'));

    
    Decoder={};
    for a=1:8
        Decoder{a}=MaxDecodersFix{a,51};
        for r=52:100
            Decoder{a}= Decoder{a}+MaxDecodersFix{a,51};
        end
        Decoder{a}=Decoder{a}/50;
    end
    
    DecoderSimMtx=zeros(8,8,2,20);
    DecoderSimMtx_Sh=zeros(8,8,2,20,100);
    
    
    GoNogoSimMtx=zeros(8,8,20);
    for ccnum=1:20
        
        for a1=1:8;
            for a2=1:8
                for a3=1:8
                    check=0;
                   
                        if a1<a2 && min(size(MW1)>=[a1,a2]) && length(MW1{a1,a2}) >0
                            Vec1=MW1{a1,a2}(:,ccnum);
                            check=check+1;
                        elseif a2<a1 && min(size(MW2)>=[a2,a1]) && length(MW2{a2,a1}) >0
                            Vec1=MW2{a2,a1}(:,ccnum);
                            check=check+1;
                        end
                        if a1<a3 && min(size(MW1)>=[a1,a3]) && length(MW1{a1,a3}) >0
                            Vec2=MW1{a1,a3}(:,ccnum);
                            check=check+1;
                        elseif a3<a1 && min(size(MW2)>=[a3,a1]) && length(MW2{a3,a1}) >0
                            Vec2=MW2{a3,a1}(:,ccnum);
                            check=check+1;
                        end
                        if check==2
                            SimMtx(a1,a2,a3,ccnum,nM)=similarity([Vec1,Vec2]);
                        end
                    
                end
            end
        end
        for a1=1:8
            for a2=1:8
                if a1<a2 && min(size(MW1)>=[a1,a2]) && length(MW1{a1,a2}) > 0
                    DecoderSimMtx(a1,a2,1,ccnum)=similarity([MW1{a1,a2}(:,ccnum),Decoder{a1}]);
          
                    for sh=1:100
                        bufd=MW1{a1,a2}(:,ccnum);
                        bufd=bufd(randperm(length(bufd)));
                        
                        DecoderSimMtx_Sh(a1,a2,1,ccnum,sh)=similarity([bufd,Decoder{a1}]);
                  

                    end
                elseif a2<a1 && min(size(MW1)>=[a2,a1]) && length(MW1{a2,a1}) > 0
                    DecoderSimMtx(a1,a2,1,ccnum)=similarity([MW2{a2,a1}(:,ccnum),Decoder{a1}]);

                    for sh=1:100
                        bufd=MW2{a2,a1}(:,ccnum);
                        bufd=bufd(randperm(length(bufd)));
                        
                        DecoderSimMtx_Sh(a1,a2,1,ccnum,sh)=similarity([bufd,Decoder{a1}]);

                    end
                end
                
            end
        end
        
        
    end
    DecoderSimMtx=abs(DecoderSimMtx);
    DecoderSimMtx(DecoderSimMtx==0)=nan;
    
    DecoderSimMtx_Sh=abs(DecoderSimMtx_Sh);
    DecoderSimMtx_Sh(DecoderSimMtx_Sh==0)=nan;

    
    ccnumVec=1:10;
    
    Result=[]
    for i=ccnumVec
        buf=DecoderSimMtx(:,:,1,i);
        buf(max(DecoderSimMtx_Sh(:,:,1,i,:),[],5) >DecoderSimMtx(:,:,1,i))=nan;
        Result=[Result;buf;nan(1,8)];
    end
    figure();imagesc(Result)
    AllRes=[AllRes,nan(90,1),Result];
end

meanSim=zeros(8,8,8);
for ccnum=1:20
    for a1=1:8
        for a2=1:8
            for a3=1:8
                SVec=abs(SimMtx(a1,a2,a3,ccnum,:));
                meanSim(a1,a2,a3,ccnum)=mean(SVec(SVec>0));
                
            end
        end
    end
end
figure();
for a=1:8
    subplot(8,1,a);
    imagesc(squeeze(meanSim(a,:,:)));
    set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])


end





%%%%%%%%%%%%%%%%%%%%%%%%%% Example Video
Mouse=67;
%TimeInt=53500:56500;

TimeInt=51500:54500;
AddressSetup;
load(strcat('E:\Reports2\2018_03_14_CannonicalCorrelations\ROI_Bin\Datasets--5to0\HitCR_L2\Session',num2str(Mouse),'_mode1_Results.mat'))
load(strcat(LoadPath{Mouse-10},'\cellData_ZS.mat'))
load(strcat(LoadPath{Mouse-10},'\cellData.mat'))

cbuf=zeros(size(CortexArea));
cbuf(CortexArea==1)=MW1C{1,6}(:,1);
[~,ind1]=sort(cbuf,'descend');
cbuf=zeros(size(CortexArea));
cbuf(CortexArea==6)=MW2C{1,6}(:,1);
[~,ind2]=sort(cbuf,'descend');


Buffer=readHDF5Subset('T:\L368_Concat\dfof_single.h5',[0 0 TimeInt(1)-1],[1017,1017,length(TimeInt)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncell=[4,5,6,9,12,13,17,19];
Vid1=zeros(21*length(ncell),21,3001);
for j=1:length(ncell)
    i=ncell(j);
    X=Buffer(cellIJ(ind1(i),1)-10:cellIJ(ind1(i),1)+10,cellIJ(ind1(i),2)-10:cellIJ(ind1(i),2)+10,:);
    X=X-min(min(min(X)));
    X=X/max(max(max(X)));
    M=mean(X,3);
    for f=1:3001
        X(:,:,f)=X(:,:,f)-M;
        X(:,:,f)=X(:,:,f)./M;
    end
    Vid1(((j-1)*21) +1: (j*21),:,:)=X;
end
time [552,744,2400,2493,2537];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncell=[5,7,8,38,18,19,20,39];
%ncell=21:40;
Vid2=zeros(21*length(ncell),21,3001);
for j=1:length(ncell)
    i=ncell(j);
    X=Buffer(cellIJ(ind2(i),1)-10:cellIJ(ind2(i),1)+10,cellIJ(ind2(i),2)-10:cellIJ(ind2(i),2)+10,:);
    X=X-min(min(min(X)));
    X=X/max(max(max(X)));
    M=mean(X,3);
    for f=1:3001
        X(:,:,f)=X(:,:,f)-M;
        X(:,:,f)=X(:,:,f)./M;
    end
    Vid2(((j-1)*21) +1: (j*21),:,:)=X;
end

X=zeros(size(Vid2,1),size(Vid2,2),1,101);
X(:,:,1,:)=Vid2(:,:,2350:2450);
X(X<0)=0;
X(X>1)=1;

writerObj = VideoWriter(strcat('C:\Users\Sadegh\Documents\VLMReborn\Presentations\P13_08_15_2018\CCA\A2_2350_2450.avi'));
writerObj.FrameRate=10;
open(writerObj)
writeVideo(writerObj,X);
close(writerObj)


P1=MW1C{1,6}(:,1)' *cellData_ROI_bin(CortexArea==1,:);
P2=MW2C{1,6}(:,1)' *cellData_ROI_bin(CortexArea==6,:);


Colors=linspecer(40);
figure();hold on;
for i=1:10
    plot((1:length(TimeInt))/10,cellData_ROI(ind1(i),TimeInt)+(8*i),'Color',Colors(i,:));
    plot((1:length(TimeInt))/10,cellData_ROI(ind2(i),TimeInt)+(8*i)+90,'Color',Colors(41-i,:));
end
plot((1:length(TimeInt))/10,4*P1(TimeInt)+190,'b')
plot((1:length(TimeInt))/10,4*P2(TimeInt)+210,'r')

