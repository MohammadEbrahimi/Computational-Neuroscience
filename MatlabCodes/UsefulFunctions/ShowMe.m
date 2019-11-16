function [out]=ShowMe(X,normalize,fps)
%%%normalize: 1=whole data , 2= frame by frame
dim=size(X);
X= double(X);
if normalize==1
    
    X=X-min(min(min(X)));
    X=X/max(max(max(X)));
end
if normalize==2
    for frame=1:dim(3)
        X(:,:,frame)=X(:,:,frame)-min(min(X(:,:,frame)));
        X(:,:,frame)=X(:,:,frame)/max(max(X(:,:,frame)));
    end
end
    X=uint8(X*255);

    if length(size(X))==3
        out=zeros(dim(1),dim(2),3,dim(3));
        for i=1:dim(3)
            out(:,:,:,i)=ind2rgb(X(:,:,i),jet(255));
        end
        implay(out,fps);
        
    end
%     if length(size(X))==2
%         
%         out=ind2rgb(X,jet(255));
%         imshow(out)
%     end
%     
    
end



% 
% inp=L354_linnue;
% out=zeros(500,500,25);
% for i=1:5
%     for j-1:5
%         for f=1:25
%             out((i-1)*100+1:i*100,(j-1)*100+1:j*100,f)=inp(i,j,f);
%         end
%     end
% end
% RGB=ShowMe(out);
% RGB(1:20,1:20,1,6:25)=1;
% writerObj = VideoWriter('C:\Users\Sadegh\Documents\VLMReborn\Reports\2016_18_10_AvgMovie\DuringStim\HitCR\L362_b3.avi');
% writerObj.FrameRate=2;
% open(writerObj)
% writeVideo(writerObj,RGB);
% close(writerObj)
