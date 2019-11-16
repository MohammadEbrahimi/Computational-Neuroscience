tile=1;
load(strcat('C:\Users\Sadegh\Documents\VLMReborn\L367\Data\ConcatDays\CellSorting\ExtractTraces\AllCells_extract_',num2str(tile),'.mat'));
%load(strcat('E:\L347_tiles\MaximumProjections.mat'));
movieAddress=strcat('E:\L367_tiles\CaMovie_masktiled_',num2str(tile),'.h5');
hinfo=hdf5info(movieAddress);
dims=hinfo.GroupHierarchy.Datasets.Dims;
movietile=readHDF5Subset(movieAddress,[0,0,0],dims);

% cellImages=zeros(dims(1),dims(2),Nc(tile));
% Start=1;
% if tile>1
%     Start=sum(Nc(1:tile-1))+1;
% end
% End=sum(Nc(1:tile));
% 
% for i=1:Nc(tile)
%    cellImages(:,:,i)=ImageOut{Start+i-1}; 
% end
timeSample=[1];
for j=1:1000:dims(3)
    timeSample=[timeSample,j:j+200];
end

%[~,~,choice]= signalSorter(cellImages,TraceOut(Start:End,timeSample),'inputMovie',movietile(:,:,timeSample))%,'valid',choice)
cellImages=output.spatial_weights;
TraceOut=double(output.temporal_weights)';

[~,~,choice]= signalSorter(cellImages,TraceOut(:,timeSample),'inputMovie',movietile(:,:,timeSample));%,'valid',choice)

% clear movietile
% sortedCells{tile}=choice;
% save('C:\Users\Sadegh\Documents\VLMReborn\L367\Data\ConcatDays\CellSorting\ExtractTraces\Manual_extract_sorted','sortedCells');