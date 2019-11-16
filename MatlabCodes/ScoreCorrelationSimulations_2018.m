p='E:\Reports2\2018_05_03_ScoreCorrelations\ROI_Bin\';
ConCMtx=zeros(8,8,60,6);
TSCMtx=zeros(8,8,20,11,6);
CRCMtx=zeros(8,8,11,6);
Mask=zeros(8,8,6);
SesCMtx=zeros(8,8,60,6);
nM=0;
for M=[61,63:67]
    nM=nM+1;
   load(strcat(p,'Mouse',num2str(M))); 
   for a1=1:8
       for a2=1:8
           ConCMtx(a1,a2,:,nM)=squeeze(CMtx(a1,a2,:,6,1));
           TSCMtx(a1,a2,:,:,nM)=squeeze(CMtx(a1,a2,6:25,:,1));
           CRCMtx(a1,a2,:,nM)=squeeze(CrossMtx(a1,a2,:,1));
           Mask(a1,a2,nM)=(min(ConCMtx(a1,a2,:,nM))~=max(ConCMtx(a1,a2,:,nM)));
           SesCMtx(a1,a2,:,nM)=squeeze(mean(CMtx(a1,a2,:,6,2:end),5));
           if nM==6
                SesCMtx(a1,a2,:,nM)=squeeze(mean(CMtx(a1,a2,:,6,[2,3,4,6]),5));
           elseif nM==7
                SesCMtx(a1,a2,:,nM)=squeeze(mean(CMtx(a1,a2,:,6,[2,3,4,8]),5));
           end
       end
   end
end

MeanConCMtx=zeros(8,8,60);
VarConCMtx=zeros(8,8,60);
MeanSesCMtx=zeros(8,8,60);
MaxSesCMtx=zeros(8,8,60);
for a1=1:8
    for a2=1:8
        
        MeanConCMtx(a1,a2,:)=mean(ConCMtx(a1,a2,:,Mask(a1,a2,:)==1),4);
        VarConCMtx(a1,a2,:)=var(squeeze(ConCMtx(a1,a2,:,Mask(a1,a2,:)==1))');
        MeanSesCMtx(a1,a2,:)=mean(SesCMtx(a1,a2,:,Mask(a1,a2,:)==1),4);
        MaxSesCMtx(a1,a2,:)=max(SesCMtx(a1,a2,:,Mask(a1,a2,:)==1),[],4);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
FisherInfo=zeros(9,60);
SimRat=zeros(56,6);
for M=1:6
load(strcat(p,'\ScoreMouse',num2str(61))); 
     
M0=zeros(8,60);
M1=zeros(8,60);
SD0=zeros(8,60);
SD1=zeros(8,60);
for a=1:8
        M0(a,:)=mean(ScoreC{a,1});
        M1(a,:)=mean(ScoreH{a,1});
        SD0(a,:)=sqrt(var(ScoreC{a,1}));
        SD1(a,:)=sqrt(var(ScoreH{a,1}));
end
 


for td=1:60
SD=0.5*(SD0(:,td)+SD1(:,td));

    
S=squeeze(ConCMtx(:,:,td,M)) .* (SD*SD');
mask=diag(S)>0;
    for a=find(mask)'
            FisherInfo(a,td)=((M1(a,td)-M0(a,td))/(SD(a)))^2;  
    end
    
    W=inv(squeeze(S(mask,mask)))*(M1(mask,td)-M0(mask,td));
    
    FisherInfo(9,td)=(W'*(M1(mask,td)-M0(mask,td)))^2/(W'*(squeeze(S(mask,mask)))*W);     
end


SimRat(:,M)=squeeze(sum(FisherInfo(1:8,5:end)))./FisherInfo(9,5:end);
end

plot(Time,mean(SimRat,2),'r--')
%plot(Time,squeeze(sum(FisherInfo(1:8,5:end)))./FisherInfo(9,5:end));
























