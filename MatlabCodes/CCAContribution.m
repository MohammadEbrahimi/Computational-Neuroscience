Areas={'V1','LV','MV','PPC','A','S','M','RSC'};
buf={};
latentOrig=zeros(8,8,25,2,20,6);
latentExcl=zeros(8,8,25,2,20,6);
latentExcl_sh=zeros(8,8,25,2,20,6,100);
nM=1;
for Mouse=[61,63:67]



    loadpath=LoadPath{Mouse-10};
    load(strcat('E:\Reports2\2018_03_14_CannonicalCorrelations\ROI_Bin\Datasets-0to5\HitCR_L2_MixedStim\Session',num2str(Mouse),'_mode1'));
    load(strcat(loadpath,'/cellData.mat'));
    load(strcat(loadpath,'/cellData_ZS.mat'));
    load(strcat(loadpath,'/SessionLength.mat'));
    load(strcat(loadpath,'/Datasets/Datasets--5to0.mat'));
    cellEvents=cellData_ROI_bin;
    MW1={};MW2={};
    NRepeate=size(Wx0,3);
    
    CMtx_avg=zeros(8,8,25,2,20); 
    CMtx=zeros(8,8,25,2,20,NRepeate);
    dprime=zeros(8,8,25,2,20,NRepeate);
    dprime_avg=zeros(8,8,25,2,20);
    dprime_sh=zeros(8,8,25,2,2,20,NRepeate);
    
    for a1=1:8%x
        for a2=a1+1:8%y
            if sum(CortexArea==a1)>20 && sum(CortexArea==a2)>20

                W1=zeros(size(Wx0{a1,a2,1},1),20, NRepeate);
                W2=zeros(size(Wy0{a1,a2,1},1),20, NRepeate);
                for nr=1:NRepeate
                    W1(:,:,nr)= (Wx0{a1,a2,nr}.* repmat(sign(mean(Wx0{a1,a2,nr})),[size(Wx0{a1,a2,nr},1),1]));
                    W2(:,:,nr)= (Wy0{a1,a2,nr}.* repmat(sign(mean(Wy0{a1,a2,nr})),[size(Wy0{a1,a2,nr},1),1]));
                end
                
                MW1{a1,a2}=squeeze(mean(W1,3));
                MW2{a1,a2}=squeeze(mean(W2,3));
                
                
                HitTrialMtx=zeros(1,length(CortexArea),25);
                n=1;
                index=HitDataset{1,1};
                TNum=HitTrialNumber{1,1};
                
                for i=unique(TNum)
                    if sum(TNum==i)==25
                        HitTrialMtx(n,:,1:25)= cellEvents(:,index(TNum==i));
                        n=n+1;
                    end
                end
                nH=n-1;
                
                CRTrialMtx=zeros(1,length(CortexArea),25);
                n=1;
                index=CRDataset{1,1};
                TNum=CRTrialNumber{1,1};
                
                for i=unique(TNum)
                    if sum(TNum==i)==25
                        CRTrialMtx(n,:,1:25)= cellEvents(:,index(TNum==i));
                        n=n+1;
                    end
                end
                nC=n-1;
                
                if nH>nC
                    mask=zeros(1,nH);
                    mask(randperm(nH,nC))=1;
                    HitTrialMtx=HitTrialMtx(mask==1,:,:);
                elseif nH<nC
                    mask=zeros(1,nC);
                    mask(randperm(nC,nH))=1;
                    CRTrialMtx=CRTrialMtx(mask==1,:,:);
                end
                
                nTrial=min(nH,nC);
                
                
                for ccnum=1:20
                    ccnum

                    D1=squeeze(MW1{a1,a2}(:,ccnum)/norm(MW1{a1,a2}(:,ccnum)));
                    D2=squeeze(MW2{a1,a2}(:,ccnum)/norm(MW2{a1,a2}(:,ccnum)));
                    
                    
                    HitProj=zeros(size(HitTrialMtx,1),25,2);
                    CRProj=zeros(size(CRTrialMtx,1),25,2);
                    for td=1:25
                        

%                         HitProj(:,td,1)=D1'*squeeze(HitTrialMtx(:,CortexArea==a1,td))';
%                         HitProj(:,td,2)=D2'*squeeze(HitTrialMtx(:,CortexArea==a2,td))';
%                         
%                         
%                         CRProj(:,td,1)=D1'*squeeze(CRTrialMtx(:,CortexArea==a1,td))';
%                         CRProj(:,td,2)=D2'*squeeze(CRTrialMtx(:,CortexArea==a2,td))';

                            dd0=CRTrialMtx(:,CortexArea==a1,td);
                            dd1=HitTrialMtx(:,CortexArea==a1,td);
                            X=[dd0 - repmat(mean(dd0),[nTrial,1]);dd1 - repmat(mean(dd1),[nTrial,1])];
                            [~,~,lat]=pca(X);
                            latentOrig(a1,a2,td,1,ccnum,nM)=lat(ccnum);
                            
                            latentExcl(a1,a2,td,1,ccnum,nM)=var(X*D1);
                            for sh=1:100 
                                latentExcl_sh(a1,a2,td,1,ccnum,nM,sh)=var(X*D1(randperm(length(D1))));
                            end
                           
                            dd0=CRTrialMtx(:,CortexArea==a2,td);
                            dd1=HitTrialMtx(:,CortexArea==a2,td);
                            X=[dd0 - repmat(mean(dd0),[nTrial,1]);dd1 - repmat(mean(dd1),[nTrial,1])];
                            [~,~,lat]=pca(X);
                            latentOrig(a1,a2,td,2,ccnum,nM)=lat(ccnum);
                            
                            latentExcl(a1,a2,td,2,ccnum,nM)=var(X*D2);
                            for sh=1:100 
                                latentExcl_sh(a1,a2,td,2,ccnum,nM,sh)=var(X*D2(randperm(length(D2))));
                            end
                      
                        
           
                    end
                    
                   
                end
               
            end
        end
    end
         save(strcat('E:\Reports2\2018_03_14_CannonicalCorrelations\ROI_Bin\Datasets-0to5\HitCR_L2_MixedStim\Session',num2str(Mouse),'_mode1_Results'),...
             'dprime','dprime_avg','dprime_sh','CMtx','CMtx_avg','MW1','MW2','meanModeProj','sdModeProj','-v7.3')
nM=nM+1;
end


%%%%%%%%%%%%%%%%%CCA Contribution
figure();hold on
set(0,'defaultlinelinewidth',1.5);
set(0,'DefaultAxesFontSize',14);
set(0,'defaultAxesFontName','Calibri');
a1=1;
td=15;
plot(1:20,squeeze(mean(mean(latentOrig(a1,a1+1:8,td,1,:,:),2),6)))
plot(1:20,squeeze(mean(mean(latentExcl(a1,a1+1:8,td,1,:,:),2),6)),'r')
plot(1:20,squeeze(mean(mean(max(latentExcl_sh(a1,a1+1:8,td,1,:,:,:),[],7),2),6)),'k')

%
figure();
set(0,'DefaultAxesFontSize',6);
for cc=1:10
    subplot(1,10,cc);
    plot(-0.4:0.1:2,squeeze(mean(mean(latentExcl(a1,a1+1:8,:,1,cc,:),2),6)),'k')
end


