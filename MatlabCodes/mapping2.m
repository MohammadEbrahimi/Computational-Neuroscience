% hinfo=hdf5info('Q:\L364\Mapping\L364_90deg_original-Objects\Obj_1 - L364901_00001.h5');
% movie=single(hdf5read(hinfo.GroupHierarchy.Datasets(1)));
% movie(isnan(movie))=0;
warning('off');
path='I:\20170412\L368\180deg\';
prefix='L368deg1801_'
dim=[1024 1024];
Stim=StimPos0;
ImgFrame=ImgFrame0;

St=[1;diff(Stim)< -70];%& Stim>=0;
En=[diff(Stim)< -70;0];%& [Stim(2:end);0]>=0;

StFrame=ImgFrame(St==1);
StFrame=StFrame(1:(length(StFrame)-1));
EnFrame=ImgFrame(En==1);
L=length(EnFrame);

P=min(EnFrame-StFrame);
G2 = fspecial('gaussian',[6 6],2);

AvgMovie=zeros(dim(1),dim(2),P);
%Movie=zeros(dim(1),dim(2),P*L);
for i=1:P
    for j=1:L
        t=Tiff(strcat(path,prefix,num2str(StFrame(j)+i,'%05d'),'.tif'),'r');
        imageframe=t.read();
        t.close();
        AvgMovie(:,:,i)=AvgMovie(:,:,i)+ double((1/L)*imageframe);
        %Movie(:,:,(j-1)*P+i)=imageframe;%-imfilter(imageframe,G2,'same','symmetric');
        %AvgMovie(:,:,i)=AvgMovie(:,:,i)+ ((1/L)*Movie(:,:,(j-1)*P+i));
    end
end
DFF=zeros([dim(1),dim(2),P],'double');
for i=1:dim(1)
    for j=1:dim(2)
        x=squeeze(AvgMovie(i,j,:));
        M=mean(x);
        x=x-M;
        x=x/M;
        DFF(i,j,:)=x;
    end
end
% M_F=DFF;
% for i=1:P
%     M_F(:,:,i)=imfilter(DFF(:,:,i),G2,'same','symmetric');
% end
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% [coeff,score,latent]=pca(reshape(DFF,dim(1)*dim(2),dim(3))');
% PC=reshape(coeff(:,1),1024,1024);
% dim=size(DFF90);
% [U,S,V]=svd(reshape(DFF90,dim(1)*dim(2),dim(3)),0);
% k=20;
% LR=U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
% LR=reshape(LR,dim(1),dim(2),dim(3));
% figure();view_movie(LR,40);
% figure();imagesc(PC)
% 

% phasemap=zeros(dim(1),dim(2));
% for i=1:dim(1)
%     for j=1:dim(2)
%         phi=angle(fft(squeeze(DFFMovie(i,j,:))));
%         phasemap(i,j)=phi(57);
%     end
% end
% figure();imagesc(phasemap);
%   
% G1 = fspecial('gaussian',[12 12],4);
% pm=imfilter(phasemap,G1,'same','symmetric');
% pm=imfilter(pm,G1,'same','symmetric');
% pm=imfilter(pm,G1,'same','symmetric');
% figure();imagesc(pm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AvgMovieTrial=uint16(zeros(dim(1)/2,dim(2)/2,floor(L/10)*P));
%     for j=1:floor(L/10)
% for k=1:10
% for i=1:P
% 
%         t=Tiff(strcat(path,prefix,num2str(StFrame(j)+i,'%05d'),'.tif'),'r');
%         imageframe=t.read();
%         t.close();
%         imageframe_s=downsamp2d(imageframe,[2 2]);
%         AvgMovieTrial(:,:,(j-1)*P+i)=AvgMovieTrial(:,:,(j-1)*P+i)+ ((1/10)*imageframe_s);
%     end
% end
%     end
% 
% F0=mean(AvgMovieTrial,3);
% Dim=size(AvgMovieTrial);
% DFFMovie=zeros(Dim,'single');
% for i=1:Dim(3)
%     
% DFFMovie(:,:,i)=(single(AvgMovieTrial(:,:,i))-F0)./F0;
% end












