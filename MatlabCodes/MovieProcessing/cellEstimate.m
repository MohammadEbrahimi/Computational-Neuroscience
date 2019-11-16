%function [Ncells]=cellEstimate(TileProj)
TileProj=MaxProjOut{10};
Dim=size(TileProj);
Centroids=TileProj;
Ncells=0;
w=2;
D1=reshape(TileProj,[1,Dim(1)*Dim(2)]);
D1=sort(D1(D1>0));
M=D1(floor(length(D1)/10));
maxColor=max(max(TileProj));

for i=1+w:Dim(1)-w
    for j=1+w:Dim(2)-w
        filt=TileProj(i-w:i+w,j-w:j+w);
        filt(1+w,1+w)=min(min(filt));
        if TileProj(i,j)>= max(max(filt))&&  min(min(TileProj(i-1:i+1,j-1:j+1)))> M  && max(max(Centroids(i-w:i+w,j-w:j+w)))<maxColor
            Centroids(i,j)=maxColor;
            Ncells=Ncells+1;
        end
    end
end
Ncells=round(Ncells*1.1);
%end 


imagesc(Centroids)