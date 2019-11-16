function MotionCorrection(hdf5In,hdf5Out,TraceOut)
% hdf5In='E:\L364\Calcium\2017_03_11\CaMovie_highpass.h5'
% hdf5Out='E:\L364\Calcium\2017_03_11\CaMovie_highpass_motcorr.h5'
% TraceOut='E:\L364\Calcium\2017_03_11\MotCorrOutput'
CHUNK_SIZE=1000;
ioptions.inputDatasetName = '/Object';
ioptions.turboregRotation = 0;
ioptions.RegisType = 1;
ioptions.parallel = 1;
ioptions.meanSubtract = 1;
ioptions.normalizeType = 'bandpass'; % matlabDisk is alternative input. Done on input to turboreg but NOT on final movie.
ioptions.registrationFxn = 'transfturboreg';
ioptions.normalizeBeforeRegister = ''; % set to blank if don't want any filtering on output movie
ioptions.imagejFFTLarge = 10000;
ioptions.imagejFFTSmall = 80;
ioptions.saveNormalizeBeforeRegister = [];
ioptions.cropCoords = [256,256,768,768];
ioptions.closeMatlabPool = 0;
ioptions.refFrame = 1;
ioptions.refFrameMatrix = [];

hinfo=hdf5info(hdf5In);
Dim=hinfo.GroupHierarchy.Datasets.Dims;
Buffer=single(zeros(Dim(1),Dim(1),CHUNK_SIZE+1));
Buffer(:,:,1)=single(readHDF5Subset(hdf5In,[0 0 0],[Dim(1),Dim(2),1]));
clear Transformation


for f=1:CHUNK_SIZE:Dim(3)
    if Dim(3)-f+1<CHUNK_SIZE
        CHUNK_SIZE=Dim(3)-f+1;
    end
        
 Buffer(:,:,2:CHUNK_SIZE+1)=single(readHDF5Subset(hdf5In,[0 0 f-1],[Dim(1),Dim(2),CHUNK_SIZE])); 
[Buffer ResultsOutOriginal] = turboregMovie(Buffer(:,:,1:(CHUNK_SIZE+1)),'options',ioptions);
for k=f:(f+CHUNK_SIZE-1)
    Transformation{k}=ResultsOutOriginal{k-f+2};
end

if f==1
    createHdf5File(hdf5Out,'/Object', Buffer(:,:,2:end));
else
    appendDataToHdf5(hdf5Out,'/Object', Buffer(:,:,2:end));
    
end


end
save(TraceOut,'Transformation','-v7.3');
end

