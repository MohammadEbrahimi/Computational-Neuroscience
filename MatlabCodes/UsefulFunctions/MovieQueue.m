classdef MovieQueue
    properties
        Size
        Nframe
        Data
        Dims
        Hdf5path
        FileCreated
    end
    methods
        function obj=MovieQueue(size,dims,hdf5path)
            obj.Size=size;
            obj.Nframe=0;
            obj.Hdf5path=hdf5path;
            obj.Dims=dims;
            obj.FileCreated=0;
        end
        function [obj]=PushFrame(obj,Frame)
            InDim=size(Frame);
            if length(InDim)==2
                Nf=1;
            else
                Nf=InDim(3);
            end
                
                
            if obj.Nframe+Nf>obj.Size
                FlushQueue(obj);
                obj.Nframe=0;
                obj.FileCreated=1;
            end
            
            obj.Nframe=obj.Nframe+Nf;
            obj.Data(:,:,(obj.Nframe-Nf+1):obj.Nframe)=single(Frame);
            if obj.Nframe==obj.Size
                FlushQueue(obj)
                obj.Nframe=0;
                obj.FileCreated=1;
            end
        end
        function FlushQueue(obj)
            if obj.Nframe>0
                if obj.FileCreated==0
                    createHdf5File(obj.Hdf5path,'/Object', obj.Data(:,:,1:obj.Nframe));
                else
                    appendDataToHdf5(obj.Hdf5path,'/Object' , obj.Data(:,:,1:obj.Nframe));
                end                             
            end
            
        end
    end
end
