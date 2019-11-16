AddressSetup;
OutputPath='J:\GoNogo\';
nM=1;
for Mouse=[61,63:67]
    nM
    loadpath=LoadPath{Mouse-10};

    
    load(strcat(loadpath,'/cellData.mat'));
    load(strcat(loadpath,'/cellData_ZS.mat'));
    load(strcat(loadpath,'/Datasets/Datasets--5to0.mat'));
    
    index=HitDataset{1,1};
    TN=HitTrialNumber{1,1};
    X=zeros(length(unique(TN)),size(cellData_ROI_bin,1));
    
  for tind=1:25
    for n=unique(TN)   
        buf=index(TN==n);
        X(n,:)=cellData_ROI_bin(:,buf(tind));
    end
    csvwrite(strcat(OutputPath,'M',num2str(nM),'\Hit_',num2str(tind),'.cvs'),X)
  end
  
  
    index=CRDataset{1,1,1};
    TN=CRTrialNumber{1,1,1};
    X=zeros(length(unique(TN)),size(cellData_ROI_bin,1));
    
  for tind=1:25
    for n=unique(TN)   
        buf=index(TN==n);
        X(n,:)=cellData_ROI_bin(:,buf(tind));
    end
    csvwrite(strcat(OutputPath,'M',num2str(nM),'\CR_',num2str(tind),'.cvs'),X)
  end
  
    
    
    
    
    
    
    
    nM=nM+1;
    
end


