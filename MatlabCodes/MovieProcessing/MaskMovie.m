function MaskMovie(MaskPath,hdf5inPath,hdf5outPath)
% hdf5inPath='Q:\L364\Calcium\2017_03_11\CaMovie.h5';
% hdf5outPath='Q:\L364\Calcium\2017_03_11\tiles\CaMovie_masktiled_';
% MaskPath='Q:\L364\Calcium\2017_03_11\S1_Mask.png';


CACHE_SIZE=100;
clear CachedMovieIn CachedMovieOut
CachedMovieIn=MovieCache(CACHE_SIZE,hdf5inPath);
dim=CachedMovieIn.Dims;
[~,~,Mask]=imread(MaskPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Tiles
TileX=250;
TileY=250;
OverlapTile=5;

MovieSize=dim;

BeginningX=1:(TileX-OverlapTile):MovieSize(1);
EndX=min(MovieSize(1),BeginningX+TileX-1);

if (EndX(end)-BeginningX(end))<TileX/2
    EndX(end)=[];
    BeginningX(end)=[];
    EndX(end)=MovieSize(1);
end

BeginningY=1:TileY-OverlapTile:MovieSize(2);
EndY=min(MovieSize(2),BeginningY+TileY-1);

if EndY(end)-BeginningY(end)<TileY/2
    EndY(end)=[];
    BeginningY(end)=[];
    EndY(end)=MovieSize(2);
end

[ListStartX,ListStartY]=meshgrid(BeginningY,BeginningY);
[ListEndX,ListEndY]=meshgrid(EndX,EndY);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:numel(ListStartX)
    CachedMovieOut{j}=MovieQueue(CACHE_SIZE,[ListEndX(j)-ListStartX(j)+1 ListEndY(j)-ListStartY(j)+1 MovieSize(3)],...
        strcat(hdf5outPath,int2str(j),'.h5'));
    
end


for i=1:dim(3)
    [Buf,CachedMovieIn]=getFrame(CachedMovieIn,i);
    Buf(Mask~=0)=0;

    for j=1:numel(ListStartX)
        CachedMovieOut{j}=PushFrame(CachedMovieOut{j},Buf(ListStartX(j):ListEndX(j),ListStartY(j):ListEndY(j)));
    end
end


for j=1:numel(ListStartX)
    FlushQueue(CachedMovieOut{j});
end
 end