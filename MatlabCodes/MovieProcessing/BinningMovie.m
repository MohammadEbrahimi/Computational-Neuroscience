function BinningMovie(hdf5inPath,hdf5outPath,BinningFactorSpace,BinningFactorTime)

CACHE_SIZE=2048;
clear CachedMovieIn CachedMovieOut
CachedMovieIn=MovieCache(CACHE_SIZE,hdf5inPath);
dim=CachedMovieIn.Dims;

% BinningFactorSpace=2;
% BinningFactorTime=2;

DataSize=dim;
if isempty(BinningFactorSpace) || BinningFactorSpace<1 || round(BinningFactorSpace)~=BinningFactorSpace
    error('MOSAIC:Binning_Movie:NotABug:BinFact','Please select a valid binning factor');
end

if isempty(BinningFactorTime) || BinningFactorTime<1 || round(BinningFactorTime)~=BinningFactorTime
    error('MOSAIC:Binning_Movie:NotABug:BinFact','Please select a valid binning factor');
end

[FirstFrame,CachedMovieIn]=getFrame(CachedMovieIn,1);
LocalImage=downsamp2d(FirstFrame,[BinningFactorSpace BinningFactorSpace]);
FinalSize=size(LocalImage);

% We preallocate the output.
TotalFrame=DataSize(3);
FinalNumberFrame=floor(TotalFrame/BinningFactorTime);
CacheMovieOut=MovieQueue(CACHE_SIZE,[FinalSize(1) FinalSize(2) FinalNumberFrame],hdf5outPath);

multiWaitbar('Changing resolution...',0);

TimeBinningCounter=0;
OldClass=class(LocalImage);
LocalImage=zeros([FinalSize(1) FinalSize(2) 1],'single');
kOutput=1;
LocalTime=0;
kInput=1;

NewImage=zeros([FinalSize(1) FinalSize(2) 1],'single');



for i=1:TotalFrame
        [Buf,CachedMovieIn]=getFrame(CachedMovieIn,kInput);
        NewImage=downsamp2d(single(Buf),[BinningFactorSpace BinningFactorSpace]);

        kInput=kInput+1;
        
    if TimeBinningCounter==0
        LocalImage=NewImage;
    else
        LocalImage=LocalImage+NewImage;
    end
    TimeBinningCounter=TimeBinningCounter+1;

    if TimeBinningCounter==BinningFactorTime % || i==DataSize(3), we trash the end to make sure of accurate timing of all frames
        AppOutput=cast(LocalImage/TimeBinningCounter,OldClass);
        CacheMovieOut=PushFrame(CacheMovieOut,AppOutput);
        kOutput=kOutput+1;
        TimeBinningCounter=0;
    end
    multiWaitbar('Changing resolution...',i/DataSize(3));
end

FlushQueue(CacheMovieOut);
multiWaitbar('Changing resolution...','Close');
end
