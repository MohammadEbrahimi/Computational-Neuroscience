function SpatialHighPass(hdf5inPath,hdf5outPath,sigmaValue,BaselineCheck)
%%%%%%% read a hdf5 dataset and pass it through high pass filter
% hdf5inPath='A:\L362bis\Calcium\5_8_2015\concatMovie-Objects\Obj_1 - Concatenated Movie.h5';
% hdf5outPath='A:\L362bis\Calcium\5_8_2015\concatMovie-Objects\highpass.h5';
% hinfo=hdf5info(hdf5inPath);
CACHE_SIZE=100;
clear CachedMovieIn CachedMovieOut
CachedMovieIn=MovieCache(CACHE_SIZE,hdf5inPath);
dim=CachedMovieIn.Dims;
CacheMovieOut=MovieQueue(CACHE_SIZE,dim(1:3),hdf5outPath);

DataSize=dim;


% BaselineCheck=1;
% sigmaValue=20;

% We use 6 times sigma as coefficient beyond that are close to zero
FilterSize=round(6*sigmaValue);
AppliedFilter=fspecial('gaussian',FilterSize,sigmaValue);

if FilterSize/2==round(FilterSize/2)
    Identity=padarray([0.25 0.25;0.25 0.25],[(FilterSize)/2-1 (FilterSize)/2-1]);
else
    Identity=padarray(1,[(FilterSize-1)/2 (FilterSize-1)/2]);
end

% We create high pass from low pass. 
% 2 times identity to keep constant brightness
AppliedFilter=Identity-AppliedFilter;


if BaselineCheck==1
    MeanImage=zeros(DataSize(1:2),'single');
    for i=1:DataSize(3)
        [Buf,CachedMovieIn]=getFrame(CachedMovieIn,i);
        %Buf=readHDF5Subset(hdf5inPath,[0 0 i-1],[dim(1),dim(2),1]);
        MeanImage=MeanImage+single(Buf)/DataSize(3);
        multiWaitbar('Calculating baseline',i/DataSize(3));
    end
    multiWaitbar('Calculating baseline','Close');
end





for i=1:DataSize(3)
    [Buf,CachedMovieIn]=getFrame(CachedMovieIn,i);
    %Buf=readHDF5Subset(hdf5inPath,[0 0 i-1],[dim(1),dim(2),1]);
    data=single(Buf);
    if BaselineCheck==0
        AppOutput=imfilter(data,AppliedFilter);
    else
        AppOutput=imfilter(data-MeanImage,AppliedFilter)+MeanImage;
    end
    CacheMovieOut=PushFrame(CacheMovieOut,AppOutput);
%    frame=i
end


FlushQueue(CacheMovieOut);
end

