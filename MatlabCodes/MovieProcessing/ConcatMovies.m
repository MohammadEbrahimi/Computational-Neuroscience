%function ConcatMovies(inPath,hdf5outPath,Sessions,NumFrames)
% hdf5inPath='Q:\L364\Calcium\2017_03_11\CaMovie.h5';
AddressSetup
 hdf5outPath='E:\L362\ConcatDays.h5';
% MaskPath='Q:\L364\Calcium\2017_03_11\S1_Mask.png';

Sessions=12:16;
NumFrames=[35992,35994,35995,35995,23997,...
    35994,34982,35996,35994,35994,23998,...
    35995,35994,35996,35995,23996];


CACHE_SIZE=100;
clear CachedMovieIn CachedMovieOut
load(strcat(LoadPath{Sessions(1)},'\FirstFrame.mat'));
Mask_Ref=FirstFrame;
CachedMovieOut=MovieQueue(CACHE_SIZE,[1017,1017,sum(NumFrames(Sessions))],hdf5outPath);


for s=Sessions
    s
	clear CachedMovieIn
	CachedMovieIn=MovieCache(CACHE_SIZE,hdf5Path{s});
	dim=CachedMovieIn.Dims;
	
	load(strcat(LoadPath{s},'\FirstFrame.mat'));
	Mask=FirstFrame;
	
        [optimizer, metric] = imregconfig('multimodal');
    
        optimizer.InitialRadius = 0.001;
    	optimizer.Epsilon = 1.5e-6;
    	optimizer.GrowthFactor = 1.01;
    	optimizer.MaximumIterations = 1000;
    
    	tform = imregtform(Mask, Mask_Ref, 'similarity', optimizer, metric);


	for i=1:dim(3)
        i
    		[Buf,CachedMovieIn]=getFrame(CachedMovieIn,i);
		
		if s~=Sessions(1)
			Buf = imwarp(Buf,tform,'OutputView',imref2d([1017,1017]));
		end
		CachedMovieOut=PushFrame(CachedMovieOut,Buf);
	end
end
	
FlushQueue(CachedMovieOut);
%end


