
load('E:\Reports2\2018_07_17_PredictiveCoding_RaoModel\Models\Model2_Relu.mat');
GoStim = double(imread(strcat('E:\Reports2\2018_07_17_PredictiveCoding_RaoModel\GoNogo\GoStim2_32_32.bmp')))/255;
Screen = double(imread(strcat('E:\Reports2\2018_07_17_PredictiveCoding_RaoModel\GoNogo\Blank_32_32.bmp')))/255;

Screen=mean(Screen,3);
GoStim=mean(GoStim,3);
ScrImg=Screen-imfilter(Screen,AppliedFilter,'symmetric');
ScrImg=ScrImg/norm(ScrImg);
testImg=GoStim-imfilter(GoStim,AppliedFilter,'symmetric');
testImg=testImg/norm(testImg);

K1=0.1;
NRepeat=100;
Thresh1=0.001;
Thresh2=0.0006;
NoiseP=0.03;
PreStimTime=100;
StimTime=100;
max_iter_r=PreStimTime+StimTime;
image_size=FilterSize+(2*stride);

Video=zeros(image_size,image_size,max_iter_r);
recVideo=zeros(image_size,image_size,max_iter_r);
rBuf=zeros(num_r*9,max_iter_r,NRepeat);
rhBuf=zeros(num_rh,max_iter_r,NRepeat);


for Nr=1:NRepeat
    Nr
    Video(:,:,1:PreStimTime)=repmat(ScrImg,[1,1,PreStimTime])+ (NoiseP*randn([image_size,image_size,PreStimTime]));
    Video(:,:,PreStimTime+1:end)=repmat(testImg,[1,1,StimTime])+ (NoiseP*randn([image_size,image_size,StimTime]));
    
    
    
    
    segTestImg=zeros(FilterSize,FilterSize,3,3,max_iter_r);
    for i=1:max_iter_r
        for ixf=1:3
            for iyf=1:3
                segTestImg(:,:,ixf,iyf,i)=Video(1+((ixf-1)*stride):((ixf-1)*stride)+FilterSize,1+((iyf-1)*stride):((iyf-1)*stride)+FilterSize,i);
            end
        end
    end
    
    
    
    
    TestError1=zeros(3,3,max_iter_r+1);
    TestError2=zeros(1,max_iter_r);
    
    r=randn(num_r*9,1);
    rh=randn(num_rh,1);
    r=0.1*r/norm(r);
    rh=0.1*rh/norm(rh);
    
    
    
    for ixf=1:3
        for iyf=1:3
            FilterInd=(((((ixf-1)*3)+(iyf-1))*num_r)+1):((((ixf-1)*3)+iyf)*num_r);
            I=reshape( segTestImg(:,:,ixf,iyf,1),[FilterSize^2,1]);
            Fx=U*r(FilterInd);
            
            Fx(Fx<Thresh1)=0;
            
            TestError1(ixf,iyf,1)=(norm(I-Fx)^2) / (norm(I)^2);
        end
    end
    
    for iter=1:max_iter_r
        rtd=Uh*rh;
        rtd(abs(rtd)<Thresh2)=0;
        
        for ixf=1:3
            for iyf=1:3
                FilterInd=(((((ixf-1)*3)+(iyf-1))*num_r)+1):((((ixf-1)*3)+iyf)*num_r);
                I=reshape( segTestImg(:,:,ixf,iyf,iter),[FilterSize^2,1]);
                
                Fx=U*r(FilterInd);
                dFx=ones(size(Fx));
                dFx(abs(Fx)<Thresh1)=0;
                Fx(abs(Fx)<Thresh1)=0;
                dFx=diag(dFx);
                
                
                    dr=(((K1/var1)*U'*dFx *(I-(Fx)))+((K1/vartd)*(rtd(FilterInd)-r(FilterInd)))-((K1*alpha)*r(FilterInd))) ;
                dr(abs(dr)<Thresh1)=0;
                r(FilterInd)=r(FilterInd)+dr;
                
               % TestError1(ixf,iyf,iter+1)=norm(dr);
                
                 TestError1(ixf,iyf,iter+1)=(norm(I-U*r(FilterInd))^2) / (norm(I)^2);
                
                
                
            end
            
            
        end
        Fx=Uh*rh;
        dFx=ones(size(Fx));
         dFx(abs(Fx)<Thresh2)=0;
         Fx(abs(Fx)<Thresh2)=0;
        dFx=diag(dFx);
        drh=(((K1/var1)*Uh'*dFx*(r-(Fx)))-(K1*alphah*rh));
        drh(abs(drh)<Thresh2)=0;
        rh=rh+drh;
        %TestError2(iter)=norm(drh);
        TestError2(iter)=(norm(r-Fx)^2)/(norm(r)^2);
        
        
        rBuf(:,iter,Nr)=r;
        rhBuf(:,iter,Nr)=rh;
        
        
        
        
        recImg=zeros(image_size,image_size);
        recImgOL=zeros(image_size,image_size);
        for ixf=1:3
            for iyf=1:3
                FilterInd=(((((ixf-1)*3)+(iyf-1))*num_r)+1):((((ixf-1)*3)+iyf)*num_r);
                recImg(1+((ixf-1)*stride):1+((ixf-1)*stride)+FilterSize-1,1+((iyf-1)*stride):1+((iyf-1)*stride)+FilterSize-1)=recImg(1+((ixf-1)*stride):1+((ixf-1)*stride)+FilterSize-1,1+((iyf-1)*stride):1+((iyf-1)*stride)+FilterSize-1)+reshape(U*r(FilterInd),[FilterSize,FilterSize]);
                recImgOL(1+((ixf-1)*stride):1+((ixf-1)*stride)+FilterSize-1,1+((iyf-1)*stride):1+((iyf-1)*stride)+FilterSize-1)=recImgOL(1+((ixf-1)*stride):1+((ixf-1)*stride)+FilterSize-1,1+((iyf-1)*stride):1+((iyf-1)*stride)+FilterSize-1)+1;
                
            end
        end
        recVideo(:,:,iter)=recImg./recImgOL;
    end
end

meanCorrR=zeros(1,max_iter_r);
meanCorrRh=zeros(1,max_iter_r);
CMtx=zeros(num_r*9,num_r*9,max_iter_r);
CMtxh=zeros(num_rh,num_rh,max_iter_r);

for iter=1:max_iter_r
    buf=corrcoef(squeeze(rBuf(:,iter,:))');
    buf(eye(num_r*9)==1)=0;
    CMtx(:,:,iter)=buf;
    meanCorrR(iter)=mean(abs(buf(eye(num_r*9)==0)));
    
    buf=corrcoef(squeeze(rhBuf(:,iter,:))');
    buf(eye(num_rh)==1)=0;
    CMtxh(:,:,iter)=buf;
     meanCorrRh(iter)=mean(abs(CMtxh(eye(num_rh)==0)));
    
end
%recVideo=Video;
NormalVideo=recVideo-min(min(min(recVideo)));
NormalVideo=NormalVideo/max(max(max(NormalVideo)));
implay(NormalVideo,10)

% 
% X=reshape(NormalVideo(:,:,85:130),[32,32,1,46]);
% writerObj = VideoWriter(strcat('C:\Users\Sadegh\Documents\VLMReborn\Presentations\P13_08_15_2018\PC\rec_gonogo.avi'));
% writerObj.FrameRate=3;
% open(writerObj)
% writeVideo(writerObj,X);
% close(writerObj)
%     

