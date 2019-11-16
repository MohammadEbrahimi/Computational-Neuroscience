global StorageArray
hinfo=hdf5info('A:\Mapping\L354_Preproccessed\Mapping0_orig-Objects\Obj_1 - L354Mapping0deg(1).h5');
%%% movC-> mov for motion removal 
movC=single(hdf5read(hinfo.GroupHierarchy.Datasets(1)));
movC(isnan(movC))=0;
obNum=1;
timeBin=2;
dim=size(movC);
%%%%%%%%%%%
l=dim(3);
clear hinfo

StimFr=StorageArray.Children{1, obNum}.Object.Children{1, 5}.Object.Data;
ImFr=StorageArray.Children{1, obNum}.Object.Children{1, 6}.Object.Data;
movX=StorageArray.Children{1, obNum}.Object.Children{1, 7}.Object.Data;
movY=StorageArray.Children{1, obNum}.Object.Children{1, 8}.Object.Data;
mot=movX .^2 +movY.^2;

StimPos=zeros(1,l);
Motion=zeros(1,l);
G = fspecial('gaussian',[10 10],5);


Step=1
for i=1:l
    %mov(:,:,i)=imfilter(mov(:,:,i),G,'same');
    StimPos(i)=mean(StimFr(ImFr== (i-1)*timeBin));
    Motion(i)=mean(mot(ImFr== (i-1)*timeBin));
end
n=1;

%%%%%%%should be trialEnd1 for Motiona removal
clear trialEnd
for i=2:l-1
    if(StimPos(i) < StimPos(i-1) && StimPos(i) < StimPos(i+1))
        trialEnd(n)=i;
        n=n+1;
    end
end
trialCount=length(trialEnd);
p=round(mean(diff(trialEnd)));


%%%%% un-comment for motion removal
n=0;
for i=1:trialCount
    if(max(Motion(trialEnd1(i)-p+1:trialEnd1(i))) < 500)
        n=n+1;
    end
end
trialEnd=trialEnd1(1:n);
movC=single(zeros(dim(1),dim(2),trialEnd(n)));
n=1;
for i=1:trialCount
    if(max(Motion(trialEnd1(i)-p+1:trialEnd1(i))) < 500)
        movC(:,:,trialEnd(n)-p+1:trialEnd(n))=mov(:,:,trialEnd1(i)-p+1:trialEnd1(i));
        n=n+1;
    end
end
trialCount=length(trialEnd);
clear mov
l=size(movC,3);
%%%%%%%%%%%%%%%%%Denoise the Data
Step=2
sigma = 3;
filtSize = 10;
x = linspace(-filtSize / 2, filtSize / 2, filtSize);
gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
gaussFilter = gaussFilter / sum (gaussFilter); % normalize


sigma2 = 2*p;
filtSize2 = 6*p;
x2 = linspace(-filtSize2 / 2, filtSize2 / 2, filtSize2);
gaussFilter2 = exp(-x2 .^ 2 / (2 * sigma2 ^ 2));
gaussFilter2 = gaussFilter2 / sum (gaussFilter2); % normalize



Nh=32;%number of harmonies
Fs=10;
%NFFT=2^nextpow2(l);
NFFT=l;
f=Fs/2*linspace(0,1,NFFT/2+1);
fj=NFFT/p;
Fp=zeros(1,2*Nh+2);
Fp(1)=2;
Fp(2*Nh+2)=NFFT;
for i=1:Nh
    Fp(i+1)=Fp(i)+fj;
    Fp(2*Nh+2-i)=Fp(2*Nh+3-i)-fj;
end
clear H
H(1)=0;
for k=2:NFFT
    H(k)=1 / (max(min(abs(Fp(1,:)-k)),1));
end



clear y
clear Y
clear Y2
clear y2
FFTweight=zeros(dim(1),dim(2));
avemov=single(zeros(dim(1),dim(2),p));
FFTphase=zeros(dim(1),dim(2));

for i=1:dim(1)
    i
    for j=1:dim(2)
        
        for k=1:(trialCount-1)
            avemov(i,j,:)=squeeze(avemov(i,j,:))+squeeze(movC(i,j,trialEnd(k)-p+1:trialEnd(k)));
        end
        
        y=squeeze(movC(i,j,:));
        Y=fft(y,NFFT)/l;
        
        Y2=Y .* H';
        
        FFTphase(i,j)=angle(Y2(104));
        y2=squeeze(real(ifft(Y2,NFFT)));
        movC(i,j,:)= y2 - conv (squeeze(y2), gaussFilter2, 'same');
        
        FFTweight(i,j)=sum(abs(Y2).^2)/sum(abs(Y).^2);
        
    end
end



[~,AVGphase]=min(avemov(:,:,30:end),[],3);
AVGphase=AVGphase/Fs;
AVGweight=max(avemov,[],3)-min(avemov,[],3);

% %%%%%%%%%%%%%% correlation phase
% Step=3
% xr=squeeze(mov(40,300,:));
% phase=zeros(dim(1),dim(2));
% pVec=(0:p-1)-floor(p/2);
% 
% for i=1:dim(1)
%     i
%     for j=1:dim(2)
%         x2=squeeze(mov(i,j,:));
% %         [~,phase(i,j)]=min(x2(30:end));
% %         phase(i,j)=phase(i,j)+30;
%         corr=zeros(160,1);
%         
%         for k=1:p 
%             corr(k)=x2(max(-pVec(k),1):min(l,l-pVec(k)+1))' * xr(max(1,pVec(k)):min(l,l+pVec(k)+1));
%         end
%         
%         [~,ind]=max(corr);
%   
%         phase(i,j)=pVec(ind);
%         
%         
%     end
% end

%figure();imagesc(phase);colormap('Jet');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Averaging
Step=3
Nt=10;
Ns=trialCount-Nt;
trialPhase=zeros(dim(1),dim(2),Ns-1);

for i=1:Ns-1
    i
    ATrials=zeros(dim(1),dim(2),p);
    for j=1:Nt
        ATrials=ATrials+movC(:,:,trialEnd((i-1)+j)-p+1:trialEnd((i-1)+j));
    end
    ATrials=ATrials/Nt;
    [~,trialPhase(:,:,i)]=min(ATrials(:,:,30:end),[],3);
    
    trialPhase(:,:,i)=trialPhase(:,:,i)+30;
    
end
% mTrial=mTrial/Ns;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 varMap=zeros(dim(1),dim(2));
 mPhase=zeros(dim(1),dim(2));
%  mAmp=zeros(dim(1),dim(2));
% 
%  mAmp=max(mTrial(:,:,25:end),[],3)-min(mTrial(:,:,25:end),[],3);
% % 
for k=1:Ns-1
    mPhase(:,:)=mPhase(:,:)+trialPhase(:,:,k);
end
 mPhase(:,:)=mPhase(:,:)/((Ns-1)*Fs);

    for i=1:dim(1)
        for j=1:dim(2)
            varMap(i,j)=var(squeeze(trialPhase(i,j,1:Ns-1)));
        end
    end
 
%  varMap(varMap>1500)=1500;   
%     %
  figure();imagesc(mPhase);colormap('Jet');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%assuming mPhase0 and mPhase90
% for i=1:Ns
%     trialPhase0(:,:,i)=trialPhase0(:,:,i)-min(min(trialPhase0(:,:,i)));
%     trialPhase90(:,:,i)=trialPhase90(:,:,i)-min(min(trialPhase90(:,:,i))); 
% end
% 
% 
% 
% 
% 
%mPhase0=mPhase0/10;
% mPhase90=mPhase90/10;
% G1 = fspecial('gaussian',[36 36],12);
%  smp0=imfilter(mPhase0,G1,'same','symmetric');
%  smp0=imfilter(smp0,G1,'same','symmetric');
%  smp0=imfilter(smp0,G1,'same','symmetric');
% smp0=imfilter(smp0,G1,'same','symmetric');
%  smp0=imfilter(smp0,G1,'same','symmetric');
% smp0=imfilter(smp0,G1,'same','symmetric');
% % smp0=imfilter(smp0,G1,'same','symmetric');
% 
%  smp90=imfilter(mPhase90,G1,'same','symmetric');
% smp90=imfilter(smp90,G1,'same','symmetric');
% smp90=imfilter(smp90,G1,'same','symmetric');
% smp90=imfilter(smp90,G1,'same','symmetric');
% smp90=imfilter(smp90,G1,'same','symmetric');
% smp90=imfilter(smp90,G1,'same','symmetric');
% % smp90=imfilter(smp90,G1,'same','symmetric');
% % figure();imagesc(smp0);colormap('Jet');
% % figure();imagesc(smp90);colormap('Jet');
% 
% 
% [Gmag90, Gdir90] = imgradient(smp90,'prewitt');
% [Gmag0, Gdir0] = imgradient(smp0,'prewitt');
% 
% map=sind(Gdir0-Gdir90);
% figure();imagesc(map);colormap('Jet');
% % % % % % 
% % % % 
% G3 = fspecial('gaussian',[10 10],5);
% 
% Mask1=((FFTweight90+FFTweight0));
% Mask1=Mask1/max(max(Mask1));
% Mask1(Mask1>0.03)=0.03;
% 
% Mask2=(1./(varMap90+varMap0));
% Mask2(Mask2>0.005)=0.005;
% Mask2=imfilter(Mask2,G3,'same','symmetric');
% 
% Mask3=((AVGweight90+AVGweight0));
% Mask3=Mask3/max(max(Mask3));
% Mask3(Mask3>0.15)=0.15;
%     
%  
%   figure();imagesc(Mask1 .* Mask2 .* Mask3.*map);colormap('Jet'); 
% % % % 
%  i=100;
%  j=100;
% % % 
%  figure();
%  plot(StimPos);
%  hold on
% % %poi=squeeze(ATrials(i,j,:));
%   s=squeeze(mov(i,j,:));
%   plot(s*60/max(s),'k');
%   plot(Motion/1000,'r');
% %poi=poi*60/max(poi);
%plot(poi,'r');
% 



