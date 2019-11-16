MovieAddress='E:\dfof_single.h5';
load('C:\Users\Sadegh\Documents\VLMReborn\L368\Data\ConcatDays\Datasets\Datasets--5to0.mat');
Hit_Avg_Movie=zeros(1017,1017,60);
CR_Avg_Movie=zeros(1017,1017,60);
FA_Avg_Movie=zeros(1017,1017,60);
Miss_Avg_Movie=zeros(1017,1017,60);

Hit_Sqr_Movie=zeros(1017,1017,60);
CR_Sqr_Movie=zeros(1017,1017,60);
FA_Sqr_Movie=zeros(1017,1017,60);
Miss_Sqr_Movie=zeros(1017,1017,60);


Hit_Var_Movie=zeros(1017,1017,60);
CR_Var_Movie=zeros(1017,1017,60);
FA_Var_Movie=zeros(1017,1017,60);
Miss_Var_Movie=zeros(1017,1017,60);

CHUNK_SIZE=10000;

hinfo=hdf5info(MovieAddress);
Dim=hinfo.GroupHierarchy.Datasets.Dims;

nH=[0,0,0];
nC=[0,0,0];
nM=[0,0,0];
nF=[0,0,0];

for f=1:CHUNK_SIZE:Dim(3)
    f
    if Dim(3)-f+1<CHUNK_SIZE
        CHUNK_SIZE=Dim(3)-f+1;
    end
    Buffer=readHDF5Subset(MovieAddress,[0 0 f-1],[Dim(1),Dim(2),CHUNK_SIZE]);
    
    
    for i=1:max(HitTrialNumber{1,1})
        TimeVec=HitDataset{1,1}(HitTrialNumber{1,1}==i);
        if min(TimeVec)>=f && (min(TimeVec)+25)<f+CHUNK_SIZE
            Hit_Avg_Movie(:,:,1:25)=Hit_Avg_Movie(:,:,1:25) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+25);
            Hit_Sqr_Movie(:,:,1:25)=Hit_Sqr_Movie(:,:,1:25) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+25).^2;
            nH(1)=nH(1)+1;
        end
    end
    
    
    for i=1:max(HitTrialNumber{1,2})
        TimeVec=HitDataset{1,2}(HitTrialNumber{1,2}==i);
        if min(TimeVec)>=f && (min(TimeVec)+5)<f+CHUNK_SIZE
            Hit_Avg_Movie(:,:,26:30)=Hit_Avg_Movie(:,:,26:30) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+5);
            Hit_Sqr_Movie(:,:,26:30)=Hit_Sqr_Movie(:,:,26:30) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+5).^2;
            nH(2)=nH(2)+1;
        end
    end
    
    
    for i=1:max(HitTrialNumber{1,3})
        TimeVec=HitDataset{1,3}(HitTrialNumber{1,3}==i);
        if min(TimeVec)>=f && (min(TimeVec)+30)<f+CHUNK_SIZE
            Hit_Avg_Movie(:,:,31:60)=Hit_Avg_Movie(:,:,31:60) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+30);
            Hit_Sqr_Movie(:,:,31:60)=Hit_Sqr_Movie(:,:,31:60) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+30).^2;
            nH(3)=nH(3)+1;
        end
    end
    
  %%%%%%%%%%%%%%%%%%  
    
    for i=1:max(CRTrialNumber{1,1})
        TimeVec=CRDataset{1,1}(CRTrialNumber{1,1}==i);
        if min(TimeVec)>=f && (min(TimeVec)+25)<f+CHUNK_SIZE
            CR_Avg_Movie(:,:,1:25)=CR_Avg_Movie(:,:,1:25) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+25);
            CR_Sqr_Movie(:,:,1:25)=CR_Sqr_Movie(:,:,1:25) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+25).^2;
            nC(1)=nC(1)+1;
        end
    end
    
    
    for i=1:max(CRTrialNumber{1,2})
        TimeVec=CRDataset{1,2}(CRTrialNumber{1,2}==i);
        if min(TimeVec)>=f && (min(TimeVec)+5)<f+CHUNK_SIZE
            CR_Avg_Movie(:,:,26:30)=CR_Avg_Movie(:,:,26:30) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+5);
            CR_Sqr_Movie(:,:,26:30)=CR_Sqr_Movie(:,:,26:30) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+5).^2;
            nC(2)=nC(2)+1;
        end
    end
    
    
    for i=1:max(CRTrialNumber{1,3})
        TimeVec=CRDataset{1,3}(CRTrialNumber{1,3}==i);
        if min(TimeVec)>=f && (min(TimeVec)+30)<f+CHUNK_SIZE
            CR_Avg_Movie(:,:,31:60)=CR_Avg_Movie(:,:,31:60) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+30);
            CR_Sqr_Movie(:,:,31:60)=CR_Sqr_Movie(:,:,31:60) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+30).^2;
            nC(3)=nC(3)+1;
        end
    end
    
   %%%%%%%%%%%%%%%%%%% 
    
    for i=1:max(FATrialNumber{1,1})
        TimeVec=FADataset{1,1}(FATrialNumber{1,1}==i);
        if min(TimeVec)>=f && (min(TimeVec)+25)<f+CHUNK_SIZE
            FA_Avg_Movie(:,:,1:25)=FA_Avg_Movie(:,:,1:25) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+25);
            FA_Sqr_Movie(:,:,1:25)=FA_Sqr_Movie(:,:,1:25) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+25).^2;
            nF(1)=nF(1)+1;
        end
    end
    
    
    for i=1:max(FATrialNumber{1,2})
        TimeVec=FADataset{1,2}(FATrialNumber{1,2}==i);
        if min(TimeVec)>=f && (min(TimeVec)+5)<f+CHUNK_SIZE
            FA_Avg_Movie(:,:,26:30)=FA_Avg_Movie(:,:,26:30) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+5);
            FA_Sqr_Movie(:,:,26:30)=FA_Sqr_Movie(:,:,26:30) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+5).^2;
            nF(2)=nF(2)+1;
        end
    end
    
    
%     for i=1:max(FATrialNumber{1,3})
%         TimeVec=FADataset{1,3}(FATrialNumber{1,3}==i);
%         if min(TimeVec)>=f && (min(TimeVec)+30)<f+CHUNK_SIZE
%             FA_Avg_Movie(:,:,31:60)=FA_Avg_Movie(:,:,31:60) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+30);
%             FA_Sqr_Movie(:,:,31:60)=FA_Sqr_Movie(:,:,31:60) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+30).^2;
%             nF(3)=nF(3)+1;
%         end
%     end
%     
    %%%%%%%%%%%%%%%%
    
    for i=1:max(MissTrialNumber{1,1})
        TimeVec=MissDataset{1,1}(MissTrialNumber{1,1}==i);
        if min(TimeVec)>=f && (min(TimeVec)+25)<f+CHUNK_SIZE
            Miss_Avg_Movie(:,:,1:25)=Miss_Avg_Movie(:,:,1:25) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+25);
            Miss_Sqr_Movie(:,:,1:25)=Miss_Sqr_Movie(:,:,1:25) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+25).^2;
            nM(1)=nM(1)+1;
        end
    end
    
    
    for i=1:max(MissTrialNumber{1,2})
        TimeVec=MissDataset{1,2}(MissTrialNumber{1,2}==i);
        if min(TimeVec)>=f && (min(TimeVec)+5)<f+CHUNK_SIZE
            Miss_Avg_Movie(:,:,26:30)=Miss_Avg_Movie(:,:,26:30) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+5);
            Miss_Sqr_Movie(:,:,26:30)=Miss_Sqr_Movie(:,:,26:30) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+5).^2;
            nM(2)=nM(2)+1;
        end
    end
    
    
    for i=1:max(MissTrialNumber{1,3})
        TimeVec=MissDataset{1,3}(MissTrialNumber{1,3}==i);
        if min(TimeVec)>=f && (min(TimeVec)+30)<f+CHUNK_SIZE
            Miss_Avg_Movie(:,:,31:60)=Miss_Avg_Movie(:,:,31:60) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+30);
            Miss_Sqr_Movie(:,:,31:60)=Miss_Sqr_Movie(:,:,31:60) + Buffer(:,:,min(TimeVec)-f+1:min(TimeVec)-f+30).^2;
            nM(3)=nM(3)+1;
        end
    end
    
    
    
end

    Hit_Avg_Movie(:,:,1:25)=Hit_Avg_Movie(:,:,1:25)/nH(1);
    Hit_Avg_Movie(:,:,26:30)=Hit_Avg_Movie(:,:,26:30)/nH(2);
    Hit_Avg_Movie(:,:,31:60)=Hit_Avg_Movie(:,:,31:60)/nH(3);
    
    Hit_Var_Movie(:,:,1:25)=(Hit_Sqr_Movie(:,:,1:25)/nH(1)) - Hit_Avg_Movie(:,:,1:25).^2 ;
    Hit_Var_Movie(:,:,26:30)=(Hit_Sqr_Movie(:,:,26:30)/nH(2))- Hit_Avg_Movie(:,:,26:30).^2;
    Hit_Var_Movie(:,:,31:60)=(Hit_Sqr_Movie(:,:,31:60)/nH(3)) - Hit_Avg_Movie(:,:,31:60).^2;
    %
    CR_Avg_Movie(:,:,1:25)=CR_Avg_Movie(:,:,1:25)/nC(1);
    CR_Avg_Movie(:,:,26:30)=CR_Avg_Movie(:,:,26:30)/nC(2);
    CR_Avg_Movie(:,:,31:60)=CR_Avg_Movie(:,:,31:60)/nC(3);
    
    CR_Var_Movie(:,:,1:25)=(CR_Sqr_Movie(:,:,1:25)/nC(1)) - CR_Avg_Movie(:,:,1:25).^2 ;
    CR_Var_Movie(:,:,26:30)=(CR_Sqr_Movie(:,:,26:30)/nC(2))- CR_Avg_Movie(:,:,26:30).^2;
    CR_Var_Movie(:,:,31:60)=(CR_Sqr_Movie(:,:,31:60)/nC(3)) - CR_Avg_Movie(:,:,31:60).^2;
    %
    FA_Avg_Movie(:,:,1:25)=FA_Avg_Movie(:,:,1:25)/nF(1);
    FA_Avg_Movie(:,:,26:30)=FA_Avg_Movie(:,:,26:30)/nF(2);
    %FA_Avg_Movie(:,:,31:60)=FA_Avg_Movie(:,:,31:60)/nF(3);
    
    FA_Var_Movie(:,:,1:25)=(FA_Sqr_Movie(:,:,1:25)/nF(1)) - FA_Avg_Movie(:,:,1:25).^2 ;
    FA_Var_Movie(:,:,26:30)=(FA_Sqr_Movie(:,:,26:30)/nF(2))- FA_Avg_Movie(:,:,26:30).^2;
    %FA_Var_Movie(:,:,31:60)=(FA_Sqr_Movie(:,:,31:60)/nF(3)) - FA_Avg_Movie(:,:,31:60).^2;
    %
    Miss_Avg_Movie(:,:,1:25)=Miss_Avg_Movie(:,:,1:25)/nM(1);
    Miss_Avg_Movie(:,:,26:30)=Miss_Avg_Movie(:,:,26:30)/nM(2);
    Miss_Avg_Movie(:,:,31:60)=Miss_Avg_Movie(:,:,31:60)/nM(3);
    
    Miss_Var_Movie(:,:,1:25)=(Miss_Sqr_Movie(:,:,1:25)/nM(1)) - Miss_Avg_Movie(:,:,1:25).^2 ;
    Miss_Var_Movie(:,:,26:30)=(Miss_Sqr_Movie(:,:,26:30)/nM(2))- Miss_Avg_Movie(:,:,26:30).^2;
    Miss_Var_Movie(:,:,31:60)=(Miss_Sqr_Movie(:,:,31:60)/nM(3)) - Miss_Avg_Movie(:,:,31:60).^2;
    
    
save('E:\Reports2\2018_03_02_AverageMovies\L368','Hit_Avg_Movie','CR_Avg_Movie','FA_Avg_Movie','Miss_Avg_Movie','Hit_Var_Movie','CR_Var_Movie','FA_Var_Movie','Miss_Var_Movie','-v7.3')



% X=Hit_Avg_Movie;
% X=X-min(min(min(X)));
% X=X/max(max(max(X)));
% X(1:30,990:1017,5:25)=1;
% X(1:15,1002:1017,26:30)=1;
% implay(X,1)
% 



