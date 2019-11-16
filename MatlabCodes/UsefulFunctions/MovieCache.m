classdef MovieCache
    properties 
        StartFrame
        Size
        Data
        Dims
        Hdf5path
        HitCount
        MissCount
    end
    methods
        function obj=MovieCache(size,hdf5path)
            obj.Size=size;
            obj.StartFrame=0;
            obj.Hdf5path=hdf5path;
            obj.HitCount=0;
            obj.MissCount=0;
            
            hinfo=hdf5info(obj.Hdf5path);
            obj.Dims=hinfo.GroupHierarchy.Datasets.Dims;
            obj.Data=zeros(obj.Dims(1),obj.Dims(2),obj.Size);
        end
        function [Frame,obj]=getFrame(obj,i)
            if obj.StartFrame==0 || min(i)< obj.StartFrame || max(i)>(obj.StartFrame+obj.Size-1)
                obj.Data=readHDF5Subset(obj.Hdf5path,[0 0 min(i)-1 ],[obj.Dims(1),obj.Dims(2),min([obj.Size,obj.Dims(3)-min(i)+1])]);
                obj.StartFrame=min(i);
                obj.MissCount=obj.MissCount+1;
            else
                obj.HitCount=obj.HitCount+1;
            end
            
                index=i-obj.StartFrame+1;
                Frame=obj.Data(:,:,index);
        end
       
    end
end
                