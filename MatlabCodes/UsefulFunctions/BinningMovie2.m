function [AppOutput]=BinningMovie2(In,BinningFactorSpace,BinningFactorTime)
DataSize=size(In);
FirstFrame=In(:,:,1);
LocalImage=downsamp2d(FirstFrame,[BinningFactorSpace BinningFactorSpace]);
FinalSize=size(LocalImage);

% We preallocate the output.
TotalFrame=DataSize(3);
FinalNumberFrame=floor(TotalFrame/BinningFactorTime);

TimeBinningCounter=0;
OldClass=class(LocalImage);
LocalImage=zeros([FinalSize(1) FinalSize(2) 1]);
kOutput=1;
kInput=1;



AppOutput=zeros([FinalSize(1) FinalSize(2) FinalNumberFrame]);

for i=1:TotalFrame
        Buf=In(:,:,i);
        NewImage=downsamp2d(single(Buf),[BinningFactorSpace BinningFactorSpace]);

        kInput=kInput+1;
        
    if TimeBinningCounter==0
        LocalImage=NewImage;
    else
        LocalImage=LocalImage+NewImage;
    end
    TimeBinningCounter=TimeBinningCounter+1;

    if TimeBinningCounter==BinningFactorTime % || i==DataSize(3), we trash the end to make sure of accurate timing of all frames
        AppOutput(:,:,kOutput)=cast(LocalImage/TimeBinningCounter,OldClass);
        
        kOutput=kOutput+1;
        TimeBinningCounter=0;
    end
   
end

end