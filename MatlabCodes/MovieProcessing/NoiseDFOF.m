function NoiseDFOF(hdf5inPath,hdf5outPath,FrameBlock)

CACHE_SIZE=200;
clear CachedMovieIn CachedMovieOut
CachedMovieIn=MovieCache(CACHE_SIZE,hdf5inPath);
dim=CachedMovieIn.Dims;
CacheMovieOut=MovieQueue(CACHE_SIZE,dim(1:3),hdf5outPath);

DataSize=dim;
NbFrame=DataSize(3);
%%%10 s
% FrameBlock=100
% We block process this for speed, especially in the HD case
BlockProc=20;
IndexMatrix=1:BlockProc:NbFrame;
NumberBlock=numel(IndexMatrix);
F0Image=0;
k=0;

multiWaitbar('Calculating Edges',0);
[FirstFrame,CachedMovieIn]=getFrame(CachedMovieIn,1);
F0Min=FirstFrame;
F0Max=FirstFrame;

% We first calculate FO image
for i=IndexMatrix
    k=k+1;
    
    currentIndexBlock=i:min(NbFrame,(i+BlockProc-1));
    
    % Mean is weighted as last block can have different number of
    % frames
    [Buf,CachedMovieIn]=getFrame(CachedMovieIn,currentIndexBlock);
    F0Min=min(min(Buf,[],3),F0Min);
    F0Max=max(max(Buf,[],3),F0Max);

    multiWaitbar('Calculating Edges',k/numel(IndexMatrix));
end
multiWaitbar('Calculating Edges','Close');

F0Min=single(F0Min);
F0Max=single(F0Max);

% This is to get rid of saturation
NumberBinsOut=1;

% Edge matrix is normalized
nBin=round(1+log2(DataSize(3)));
Edge=0:1/nBin:(1-NumberBinsOut/nBin);

multiWaitbar('Calculating frame distribution',0);

% Initiate
CountedBin=zeros(size(histc(FirstFrame,Edge,3)));
k=0;

% We first calculate FO image
for i=IndexMatrix
    k=k+1;
    
    currentIndexBlock=i:min(NbFrame,(i+BlockProc-1));
    [Buf,CachedMovieIn]=getFrame(CachedMovieIn,currentIndexBlock);
    currentBlock=single(Buf);
    currentBlock=(currentBlock-repmat(F0Min,[1 1 numel(currentIndexBlock)]))./...
        repmat(F0Max-F0Min,[1 1 numel(currentIndexBlock)]);
    
    CountedBin=CountedBin+histc(currentBlock,Edge,3);
    
    multiWaitbar('Calculating frame distribution',k/numel(IndexMatrix));
end
multiWaitbar('Calculating frame distribution','Close');

multiWaitbar('Calculating F0',0);

[~,FoundIndex]=max(CountedBin,[],3);

% initiate
F0Image=double(FirstFrame);

% We calculate the mean of both edges.
% Last edge is a particular case
CorrespondingValue=Edge/2;
CorrespondingValue(1:(end-1))=CorrespondingValue(1:(end-1))+CorrespondingValue(2:end);
CorrespondingValue(end)=Edge(end);

F0Image(:)=F0Min(:)+(F0Max(:)-F0Min(:)).*CorrespondingValue(FoundIndex(:))';
F0Block=single(F0Image);

multiWaitbar('Calculating F0','Close');


multiWaitbar('Calculating variance',0);

% Initiate
CountedBin=zeros(size(histc(FirstFrame,Edge,3)));
k=0;
VarianceImage=Inf*ones(size(F0Image),'single');

IndexMatrixNoise=1:FrameBlock:NbFrame;

% We estimate the noise for each pixel
for i=IndexMatrixNoise
    k=k+1;
        k=k+1;
    currentBlock=i:min(NbFrame,(i+FrameBlock-1));
    
    if size(F0Block,3)~=numel(currentBlock)
        F0Block=repmat(F0Image,[1,1,length(currentBlock)]);
    end
    
    [Buf,CachedMovieIn]=getFrame(CachedMovieIn,currentBlock);
    LocalBlock=single(Buf);
    
    if numel(currentBlock)==FrameBlock
        VarianceImage=min(VarianceImage,sum((LocalBlock-F0Block).^2/(size(LocalBlock,3)-1),3));
    end
    
    multiWaitbar('Calculating variance',k/numel(IndexMatrix));
end
VarianceImage=sqrt(VarianceImage);

multiWaitbar('Calculating variance','Close');

multiWaitbar('Normalizing noise',0);

% We preallocate
F0Block=single(F0Image);
VarianceBlock=single(VarianceImage);

k=0;

for i=IndexMatrix
    k=k+1;
    currentBlock=i:min(NbFrame,(i+BlockProc-1));
    
    if size(F0Block,3)~=numel(currentBlock)
        F0Block=repmat(F0Image,[1,1,length(currentBlock)]);
        VarianceBlock=repmat(VarianceImage,[1,1,length(currentBlock)]);
    end
    
    [LocalBlock,CachedMovieIn]=getFrame(CachedMovieIn,currentBlock);
    LocalBlock=~(F0Block==0).*(single(LocalBlock)-F0Block)./VarianceBlock;
    %AppOutput.Data(1:DataSize(1),1:DataSize(2),currentBlock,1)=LocalBlock;
    CacheMovieOut=PushFrame(CacheMovieOut,LocalBlock);
    multiWaitbar('Normalizing noise',k/NumberBlock);
end

multiWaitbar('Normalizing noise','Close');

FlushQueue(CacheMovieOut);
end
