AddressSetup;
MeanFRtable=zeros(8,6);
MedFRtable=zeros(8,6);
nM=1;
for Mouse=[61,63:67]
    nM
    loadpath=LoadPath{Mouse-10};

    
    load(strcat(loadpath,'/cellData.mat'));
    load(strcat(loadpath,'/cellData_ZS.mat'));
    load(strcat(loadpath,'/Datasets/Datasets-0to0.mat'));
    TimeVec=[HitDataset{1,1},CRDataset{1,1,3}];
    
    clear V;
    clear P;
     V.T       = size(cellData_ROI,2); % # of time steps
    V.dt    = 1/10;  % time step size
    V.smc_iter_max = 1;
    
    % initialize params
    P.a     = 1;    % observation scale
    P.b     = 0;    % observation bias
    tau     = 2;    % decay time constant
    P.gam   = 1-V.dt/tau; % C(t) = gam*C(t-1)
    P.lam   =0.05;
    P.sig   = 0.5;  % standard deviation of observation noise
    
    cellData_oopsi=zeros(size(cellData_ROI));
    
    
    parfor i=1:size(cellData_ROI,1)
        i
        cellData_oopsi(i,:)=(fast_oopsi(double(cellData_ROI(i,:)),V,P));
        
    end
    save(strcat('J:\VLMReborn\FNNDTraces\M',num2str(nM)),'cellData_oopsi','-v7.3')
    
    for a=1:8
        if sum(CortexArea==a)>5
            MeanFRtable(a,nM)=mean( mean(cellData_oopsi(CortexArea==a,TimeVec)));
            MedFRtable(a,nM)=median( mean(cellData_oopsi(CortexArea==a,TimeVec),2));

        end
                      
            
    end
    nM=nM+1;
end
   
for a=1:8
    MeanFR(a)=10*mean(MeanFRtable(a,MeanFRtable(a,:)>0));
    MedFR(a)=10*mean(MedFRtable(a,MedFRtable(a,:)>0));
end