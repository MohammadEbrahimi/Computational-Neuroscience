 hdf5In='E:\L364\Calcium\2017_03_11\CaMovie_highpass.h5'
 hdf5Out='E:\L364\Calcium\2017_03_11\CaMovie_highpass_motcorr.h5'
 TraceOut='E:\L364\Calcium\2017_03_11\MotCorrOutput'
 
 


hinfo=hdf5info(hdf5In);
Dim=hinfo.GroupHierarchy.Datasets.Dims;
Dim(3)=500;
MovieInput=single(readHDF5Subset(hdf5In,[0 0 0],[Dim(1),Dim(2),Dim(3)]));;
ImageInput=MovieInput(:,:,1);

% Motion Correction setup--------------------------------------------------
% Update handles
%handles=guidata(handles.output);
% We estimate the pyramid size based on the picture size
MIN_SIZE=12;
pyramidDepth = 1;
% This code is taken from ImageJ TurboReg plugin
sw=Dim(1);
sh=Dim(2);
while (((2 * MIN_SIZE) <= sw) && ((2 * MIN_SIZE) <= sh))
    sw=sw/2;
    sh=sh/2;
    pyramidDepth=pyramidDepth+1;
end

DataSize=Dim;
%MaxFrame=DataSize(3);
%MatrixMotCorrDispl=zeros(3,MaxFrame);
InvertImage=1;

% if isempty(ImageInput)
%     RefPic=MovieInput.Data(:,:,1,1);
% else
% Force user to choose a reference image
    RefPic=ImageInput*double(1.0);
% end

% Default values if no selection
BW_ROI=ones(size(RefPic),'single');

% We gather region selection from interface.
PolyPos=[357,157;...
    757,157;...
    840,512;...
    738,875;...
    333,900;...
    128,512];

% We convert into pixels coordinates
% PolyPos(:,1) = axes2pix(DataSize(2), MovieInput.Xposition, PolyPos(:,1));
% PolyPos(:,2) = axes2pix(DataSize(1), MovieInput.Yposition, PolyPos(:,2));

if ~iscell(PolyPos)
    BW_ROI = poly2mask(PolyPos(:,1), PolyPos(:,2), DataSize(1), DataSize(2));
end

% single then double is to force copy on write to lose reference to feed
% into mex files (that works by reference).
Mask1=single(double(BW_ROI));
Mask2=single(double(BW_ROI));

TurboRegOptions.RegisType=1;
TurboRegOptions.SmoothX=2;
TurboRegOptions.SmoothY=2;

% This is just to provide some more flexibility to adjust computation speed
SpeedTradeOff=double(0.1);
if SpeedTradeOff>0.5
    Lastlevels=2;
    minGain=2*SpeedTradeOff-0.5;
else
    Lastlevels=1;
    minGain=2*SpeedTradeOff;
end

TurboRegOptions.minGain=minGain;
TurboRegOptions.Levels=pyramidDepth;
TurboRegOptions.Lastlevels=Lastlevels;
TurboRegOptions.Epsilon=1.192092896E-07;
TurboRegOptions.zapMean=0;
TurboRegOptions.Interp=1;

OriginalClass=class(MovieInput);

if (TurboRegOptions.RegisType==1 || TurboRegOptions.RegisType==2)
    TransformationType='affine';
else
    TransformationType='projective';
end

ClassMovie=class(MovieInput);

SubsSpatialMeanCheck=1;

SpatialMeanCheck=1;
ApplySpatMeanPixels=str2num('5');

RefPic=single(RefPic);

SpatMeanPixels=str2num('20');

% PARFOR loops to see all variables no matter if they are used or not.
hDisk=fspecial('disk', SpatMeanPixels);

if SubsSpatialMeanCheck
    RefPic=RefPic-imfilter(RefPic,hDisk);
end

if InvertImage
    RefPic=-RefPic;
end

hDiskApply=fspecial('disk', ApplySpatMeanPixels);

if SpatialMeanCheck
    RefPic=imfilter(RefPic,hDiskApply);
end

CLim(1)=str2double('-inf');
CLim(2)=str2double('inf');

MinPixelThreshold=min(RefPic(:));
MaxPixelThreshold=max(RefPic(:));

MinPixelThreshold=max(MinPixelThreshold,CLim(1));
MaxPixelThreshold=min(MaxPixelThreshold,CLim(2));

RefPic(MinPixelThreshold>RefPic)=MinPixelThreshold;
RefPic(MaxPixelThreshold<RefPic)=MaxPixelThreshold;

%------------------Crop setup---------------------------------------------
%Crop at least 4 pixels on each side to get rid of motion correction
%artifacts
posImg(1)=4;
posImg(2)=4;
posImg(3)=(DataSize(2)-8);
posImg(4)=(DataSize(1)-8);
Xposition=1:DataSize(2);
Yposition=1:DataSize(1);
%Converting to frame internal coordinates
posImgAbs(1) = pix2axes(DataSize(2),Xposition,posImg(1));
posImgAbs(2) = pix2axes(DataSize(1),Yposition,posImg(2));
posImgAbs(3) = pix2axes(DataSize(2),Xposition,posImg(1)+posImg(3))-posImgAbs(1);
posImgAbs(4) = pix2axes(DataSize(1),Yposition,posImg(2)+posImg(4))-posImgAbs(2);

rect=posImgAbs;

TestImage=imcrop(Xposition,Yposition,MovieInput(:,:,1),rect);
SizeFinalCropped=size(TestImage);
% 
% Xposition=MovieInput.Xposition;
% Yposition=MovieInput.Yposition;

% We preallocate the output.
%AppOutput=MovieObject([SizeFinal(1) SizeFinal(2) MovieInput.DataSize(3) MovieInput.DataSize(4)],AppInput.DataClass);
%-------------------------------------------------------------------------

%---------------Binning (subsampling) setup-------------------------------
BinningFactorSpace=1;
BinningFactorTime=1;

if isempty(BinningFactorSpace) || BinningFactorSpace<1 || round(BinningFactorSpace)~=BinningFactorSpace
    error('MOSAIC:Binning_Movie:NotABug:BinFact','Please select a valid binning factor');
end

if isempty(BinningFactorTime) || BinningFactorTime<1 || round(BinningFactorTime)~=BinningFactorTime
    error('MOSAIC:Binning_Movie:NotABug:BinFact','Please select a valid binning factor');
end

LocalImage=downsamp2d(TestImage,[BinningFactorSpace BinningFactorSpace]);
FinalSizeSubsSpace=size(LocalImage);

FinalNumberFrame=floor(DataSize(3)/BinningFactorTime);

TimeBinningCounter=0;
OldClass=class(LocalImage);
LocalImage=zeros([FinalSizeSubsSpace(1) FinalSizeSubsSpace(2) 1],'single');
kOutput=1;
LocalTime=0;
kInput=1;
Period=0.05;
NewImage=zeros([FinalSizeSubsSpace(1) FinalSizeSubsSpace(2) 1],'single');
NewTime=Period/2;


%-------------------------------------------------------------------------

% We preallocate the output.
% FullSize movie
MovieOutFull=single(zeros([SizeFinalCropped(1) SizeFinalCropped(2) DataSize(3)]));
ParallelLoop=1;

if matlabpool('size')==0
    matlabpool('open');
end
    
AllowedSpacePerBlock=4e9; % Gigabytes
SingleFrame=MovieInput(:,:,1);
FrameProp=whos('SingleFrame');
FrameSize=FrameProp.bytes;
Block=floor(AllowedSpacePerBlock/FrameSize);
Block=BinningFactorTime*floor(Block/BinningFactorTime);
% We may trash several last frames but timimng of full movie and susampoled in time would be the same.
ActualMaxFrame=BinningFactorTime*floor(DataSize(3)/BinningFactorTime);

IterBlock=1:Block:ActualMaxFrame; 

MatrixMotCorrDispl=zeros(3,ActualMaxFrame);

k=0;

for i_block=IterBlock
    i_block
    k=k+1;
    EndBlock=min(i_block-1+Block,ActualMaxFrame);
    DataBlock=single(MovieInput(:,:,i_block:EndBlock));
    
    [~,~,L]=size(DataBlock);
    ResultBlock=zeros([SizeFinalCropped(1), SizeFinalCropped(2), L],'single');
    MatrixMotCorrDisplTmp=zeros(3,L);
    
    parfor i=1:L
    
        %------------ Motion Correction --------------------
        ToAlign=DataBlock(:,:,i);
        ToAlignMotion=ToAlign;
        if SubsSpatialMeanCheck
            ToAlignMotion=ToAlignMotion-imfilter(ToAlignMotion,hDisk);
        end
        if InvertImage
            ToAlignMotion=-ToAlignMotion;
        end
        if SpatialMeanCheck
            ToAlignMotion=imfilter(ToAlignMotion,hDiskApply);
        end
        ToAlignMotion(MinPixelThreshold>ToAlignMotion)=MinPixelThreshold;
        ToAlignMotion(MaxPixelThreshold<ToAlignMotion)=MaxPixelThreshold;
        ResultsOut=turboreg(RefPic,single(ToAlignMotion),Mask1,Mask2,TurboRegOptions);
        MatrixMotCorrDisplTmp(:,i)=[ResultsOut.Translation(1) ResultsOut.Translation(2) ResultsOut.Rotation];
        Frame=single(transfturboreg(ToAlign,ones(size(ToAlign),'single'),ResultsOut));
        %------------ Motion Correction --------------------
     
        %------------ Crop --------------------
        FrameCropped=imcrop(Xposition,Yposition,Frame,rect);
        
        %------------ Crop --------------------
        ResultBlock(:,:,i)=FrameCropped;
    end
    
    
    MovieOutFull(:,:,i_block:EndBlock)=ResultBlock;
    MatrixMotCorrDispl(:,i_block:EndBlock)=MatrixMotCorrDisplTmp;
    
end





try % use try / catch here, since delete(struct) will raise an error.
        delete(ppm);
    catch me
end

