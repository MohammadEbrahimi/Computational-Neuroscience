function Tiff2hdf5(tiffPath,hdf5outPath,prefix,Dim)
CACHE_SIZE=100;


CacheMovieOut=MovieQueue(CACHE_SIZE,Dim,hdf5outPath);

for i=1:Dim(3)
    t=Tiff(strcat(tiffPath,prefix,num2str(i,'%05d'),'.tif'),'r');
    imageframe=t.read();
    t.close();
    CacheMovieOut=PushFrame(CacheMovieOut,imageframe);
end

FlushQueue(CacheMovieOut);
end