%w = 1536; h = 1024;
FilterSize=16;
stride=8;
image_size=FilterSize+(2*stride);
segmented_dataset=zeros(FilterSize,FilterSize,3,3,1);


%img_dataset=zeros(h-4,w-4,100);
img_dataset={};
for i=1:151
    
    %f1 = fopen(strcat('E:\Reports2\2018_07_17_PredictiveCoding_RaoModel\Vanhateren_Dataset\imk',num2str(i,'%05d'),'.imc'), 'rb', 'ieee-be');
    %buf = fread(f1, [w, h], 'uint16');
    %img_dataset(:,:,i)=buf(3:w-2,3:h-2)';
    
    buf = imread(strcat('E:\Reports2\2018_07_17_PredictiveCoding_RaoModel\NaturalScenes_files\Img',num2str(i),'.jpg'));
    [h,w,~]=size(buf);
    
    img_dataset{i}=mean(buf(5:h-5,5:w-5,:),3);
    
end

% We use 6 times sigma as coefficient beyond that are close to zero

AppliedFilter=fspecial('gaussian',24,8);
RFFilter=fspecial('gaussian',FilterSize,round(FilterSize/4));
RFFilter=RFFilter/norm(RFFilter);


ni=1;
crop_dataset=zeros(image_size,image_size,1);
imageIndex=0;
for i=1:151
    [h,w]=size(img_dataset{i});
    x_seg=1:image_size:(h-image_size);
    y_seg=1:image_size:(w-image_size);
    for ix=1:length(x_seg)
        for iy=1:length(y_seg)
            crop_dataset(:,:,ni)=img_dataset{i}(x_seg(ix):x_seg(ix)+(2*FilterSize)-1,y_seg(iy):y_seg(iy)+(2*FilterSize)-1);
            imageIndex(ni)=i;
            ni=ni+1;
        end
    end
end
ni=ni-1;
meanDataset=mean(crop_dataset,3);
sdDataset=reshape(sqrt(var(reshape(crop_dataset,[image_size^2,ni])')),[image_size,image_size]);

crop_dataset=crop_dataset-repmat(mean(crop_dataset,3),[1,1,ni]);
crop_dataset=crop_dataset./repmat(reshape(sqrt(var(reshape(crop_dataset,[image_size^2,ni])')),[image_size,image_size]),[1,1,ni]);

for i=1:ni
    crop_dataset(:,:,i)=crop_dataset(:,:,i)-imfilter(crop_dataset(:,:,i),AppliedFilter,'symmetric');
    crop_dataset(:,:,i)= crop_dataset(:,:,i)/norm( crop_dataset(:,:,i));
    for ixf=1:3
        for iyf=1:3
            segmented_dataset(:,:,ixf,iyf,i)=crop_dataset(1+((ixf-1)*stride):((ixf-1)*stride)+FilterSize,1+((iyf-1)*stride):((iyf-1)*stride)+FilterSize,i);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_iter_r=100;
batchRun=1;
trainSize=3000;
num_r=64;
num_rh=256;
K1=0.1;
K2=1;
var1=5;
vartd=10;
alpha=0.7;
alphah=0.05;
lambda=0.02;
Thresh1=0.001;
Thresh2=0.0006;



ErrorT1=zeros(3,3,max_iter_r);
ErrorT2=zeros(1,max_iter_r);
Error1=zeros(3,3,max_iter_r);
Error2=zeros(1,max_iter_r);

U=randn(FilterSize^2,num_r);
Uh=randn(num_r*9,num_rh);

U=normCol(U);
Uh=normCol(Uh);

for br=1:batchRun
    for i=1:trainSize
        r=randn(num_r*9,1);
        rh=randn(num_rh,1);
        r=0.1*r/norm(r);
        rh=0.1*rh/norm(rh);
        
        
        for iter=1:max_iter_r
            rtd=Uh*rh;
            rtd(abs(rtd)<Thresh2)=0;
            
            for ixf=1:3
                for iyf=1:3
                    FilterInd=(((((ixf-1)*3)+(iyf-1))*num_r)+1):((((ixf-1)*3)+iyf)*num_r);
                    I=reshape( segmented_dataset(:,:,ixf,iyf,i),[FilterSize^2,1]);
                    
                    Fx=tanh(U*r(FilterInd));
                    dFx=1-Fx.^2;
%                      dFx(abs(Fx)<Thresh1)=0;
%                      Fx(abs(Fx)<Thresh1)=0;
                    dFx=diag(dFx);
                    
                    
                      dr=(((K1/var1)*U'*dFx *(I-(Fx)))+((K1/vartd)*(rtd(FilterInd)-r(FilterInd)))-((K1*alpha)*r(FilterInd))) ;
                 %dr(abs(dr)<Thresh1)=0;
                     r(FilterInd)=r(FilterInd)+dr;
                
                                        
                    Error1(ixf,iyf,iter)=(norm(I-Fx)^2 / norm(I)^2);
                    
                end
                
                
            end
            
            Fx=tanh(Uh*rh);
            dFx=1-Fx.^2;
%              dFx(abs(Fx)<Thresh2)=0;
%              Fx(abs(Fx)<Thresh2)=0;
            dFx=diag(dFx);
            
                    drh=(((K1/var1)*Uh'*dFx*(r-(Fx)))-(K1*alphah*rh));
        %drh(abs(drh)<Thresh2)=0;
        rh=rh+drh;
        
            
            Error2(iter)=(norm(r-Fx)^2 / norm(r)^2) ;
            
        end
        
        
        for ixf=1:3
            for iyf=1:3
                FilterInd=(((((ixf-1)*3)+(iyf-1))*num_r)+1):((((ixf-1)*3)+iyf)*num_r);
                I=reshape( segmented_dataset(:,:,ixf,iyf,i),[FilterSize^2,1]);
                
                Fx=tanh(U*r(FilterInd));
                dFx=1-Fx.^2;
%                  dFx(abs(Fx)<Thresh1)=0;
%                  Fx(abs(Fx)<Thresh1)=0;
                dFx=diag(dFx);
                
                U=U+((K2/var1)*((dFx*(I-Fx)*r(FilterInd)')-(K2*lambda*U)));
                
                ErrorT1(ixf,iyf,((br-1)*trainSize)+i)=(norm(I-Fx)^2 / norm(I)^2);
                
            end
        end
        
        Fx=tanh(Uh*rh);
        dFx=1-Fx.^2;
%          dFx(abs(Fx)<Thresh2)=0;
%          Fx(abs(Fx)<Thresh2)=0;
        dFx=diag(dFx);
        
        Uh=Uh+((K2/var1)*((dFx*(r-Fx)*rh')-(K2*lambda*Uh)));
        
        
        Uh=normCol(Uh);
        U=normCol(U);
        
        
        ErrorT2(((br-1)*trainSize)+i)=(norm(r-Fx)^2 / norm(r)^2);
        
        
        if rem(i,100)==0
            i
            K2=K2/1.02;
        end
        
        
    end
    K2=K2/2;
end

Filters1=reshape(U,[FilterSize,FilterSize,num_r]);

save('E:\Reports2\2018_07_17_PredictiveCoding_RaoModel\Models\Model2_tanh_orig','max_iter_r','batchRun','trainSize',...
    'num_r','num_rh','K1','K2','var1','vartd','alpha','alphah','lambda','ErrorT1','ErrorT2','Error1','Error2',...
    'U','Uh','Filters1','FilterSize','stride','segmented_dataset','AppliedFilter');



%colormap(gray);
%imagesc(img_dataset(:,:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%testing
K1=0.1;
ti=1;
max_iter_r=100;
Thresh1=0.001;
Thresh2=0.0006;
testImg=img_dataset{ti};
figure();
subplot(1,3,1)
colormap('gray');
imagesc(testImg)
subplot(1,3,2)
colormap('gray');
imagesc(testImg-imfilter(testImg,AppliedFilter,'symmetric'))

[h,w]=size(testImg);
x_seg=1:image_size:(h-image_size);
y_seg=1:image_size:(w-image_size);

Reconstructed=zeros(max(x_seg)+image_size-1,max(y_seg)+image_size-1);
testIndex=find(imageIndex==ti);
TestError1=zeros(3,3,max_iter_r+1,length(testIndex));

for si=1:length(testIndex)
    
    
    %TestError1=zeros(3,3,max_iter_r+1);
    TestError2=zeros(1,max_iter_r+1);
    
    r=randn(num_r*9,1);
    rh=randn(num_rh,1);
    r=r/norm(r);
    rh=rh/norm(rh);
    
    
    
    for ixf=1:3
        for iyf=1:3
            FilterInd=(((((ixf-1)*3)+(iyf-1))*num_r)+1):((((ixf-1)*3)+iyf)*num_r);
            I=reshape( segmented_dataset(:,:,ixf,iyf,testIndex(si)),[FilterSize^2,1]);
            Fx=U*r(FilterInd);
%             Fx(Fx<0)=0;
            TestError1(ixf,iyf,1)=(norm(I-Fx)^2) / (norm(I)^2);
        end
    end
    
    for iter=1:max_iter_r
        rtd=Uh*rh;
        
        for ixf=1:3
            for iyf=1:3
                FilterInd=(((((ixf-1)*3)+(iyf-1))*num_r)+1):((((ixf-1)*3)+iyf)*num_r);
                I=reshape( segmented_dataset(:,:,ixf,iyf,testIndex(si)),[FilterSize^2,1]);
                Fx=U*r(FilterInd);
                dFx=ones(size(Fx));
                 dFx(abs(Fx)<Thresh1)=0;
                 Fx(abs(Fx)<Thresh1)=0;
                dFx=diag(dFx);
                
                
                r(FilterInd)=r(FilterInd)+(((K1/var1)*U'*dFx *(I-(Fx)))+((K1/vartd)*(rtd(FilterInd)-r(FilterInd)))-((K1*alpha)*r(FilterInd))) ;
                
                
                
                TestError1(ixf,iyf,iter+1,si)=(norm(I-Fx)^2) / (norm(I)^2);
                
            end
            
            
        end
        Fx=Uh*rh;
        dFx=ones(size(Fx));
         dFx(abs(Fx)<Thresh2)=0;
         Fx(abs(Fx)<Thresh2)=0;
        dFx=diag(dFx);
        
        rh=rh+(((K1/var1)*Uh'*dFx*(r-(Fx)))-(K1*alphah*rh));
        
        
        TestError2(iter)=(norm(r-Fx)^2)/(norm(r)^2);
        
    end
    
    segImg=zeros(FilterSize,FilterSize,3,3);
    recImg=zeros(image_size,image_size);
    recImgOL=zeros(image_size,image_size);
    for ixf=1:3
        for iyf=1:3
            FilterInd=(((((ixf-1)*3)+(iyf-1))*num_r)+1):((((ixf-1)*3)+iyf)*num_r);
            recImg(1+((ixf-1)*stride):1+((ixf-1)*stride)+FilterSize-1,1+((iyf-1)*stride):1+((iyf-1)*stride)+FilterSize-1)=recImg(1+((ixf-1)*stride):1+((ixf-1)*stride)+FilterSize-1,1+((iyf-1)*stride):1+((iyf-1)*stride)+FilterSize-1)+reshape(U*r(FilterInd),[FilterSize,FilterSize]);
            recImgOL(1+((ixf-1)*stride):1+((ixf-1)*stride)+FilterSize-1,1+((iyf-1)*stride):1+((iyf-1)*stride)+FilterSize-1)=recImgOL(1+((ixf-1)*stride):1+((ixf-1)*stride)+FilterSize-1,1+((iyf-1)*stride):1+((iyf-1)*stride)+FilterSize-1)+1;
            
        end
    end
    recImg=recImg./recImgOL;
    
    y_seg_Index=rem(si,length(y_seg));
    x_seg_Index=floor(si/length(y_seg))+1;
    if y_seg_Index==0
        y_seg_Index=length(y_seg);
    end
    if floor(si/length(y_seg))==ceil(si/length(y_seg))
        x_seg_Index=x_seg_Index-1;
    end
    
    Reconstructed(x_seg(x_seg_Index):x_seg(x_seg_Index)+image_size-1,y_seg(y_seg_Index):y_seg(y_seg_Index)+image_size-1)=recImg;
end

subplot(1,3,3);colormap('gray');imagesc(Reconstructed)








