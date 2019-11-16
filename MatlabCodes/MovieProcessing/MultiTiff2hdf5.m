AddressSetup;
Sessions=45;

for s=Sessions
hdf5Out= strcat('E:\Old_Sessions\Session-',num2str(s),'\CaMovie.h5')
nofile=1;
    for chunk=NumTiff{s}
        chunk
        Nf=size(imfinfo(strcat(tiffInOld{s},'\trial(',num2str(chunk),').tif')),1);
        Imds=single(zeros(1024,1024,Nf));
        Im=single(zeros(2048,2048,Nf));
        for f=1:Nf
            f
            Im(:,:,f)=imread(strcat(tiffInOld{s},'\trial(',num2str(chunk),').tif'),f);
            %Imds=downsamp2d(single(Im),[2 2]);
            Imds(:,:,f)=single(0.5*(Im(1:2:2048,1:2:2048,f)+Im(2:2:2048,2:2:2048,f)));
            
        end
                if nofile==1
                    createHdf5File(hdf5Out,'/Object', Imds);
                    nofile=0;
                else
                    appendDataToHdf5(hdf5Out,'/Object', Imds);

                end
        
        
        
    end

    
    
end
