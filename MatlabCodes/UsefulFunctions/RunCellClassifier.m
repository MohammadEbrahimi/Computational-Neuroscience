Mouse='L367';
ExPath=strcat('C:\Users\Sadegh\Documents\VLMReborn\',Mouse,'\Data\ConcatDays\CellSorting\ExtractTraces\');

%%%%Training
 ManualTiles=[9,10,13];
% load(strcat('C:\Users\Sadegh\Documents\VLMReborn\',Mouse,'\Data\ConcatDays\AllCells_pcaica.mat'));
load(strcat(ExPath,'Manual_extract_sorted.mat'));
InImages={};
InSignals={};
InTargets={};
InMovie={};
% 
 for i=1:length(ManualTiles)
     tile=ManualTiles(i)
     load(strcat(ExPath,'AllCells_extract_',num2str(tile)));
%     
%     Start=1;
%     if tile>1
%         Start=sum(Nc(1:tile-1))+1;
%     end
%     End=sum(Nc(1:tile));
%     dim=size(ImageOut{Start});
%     cellImages=zeros(dim(1),dim(2),Nc(tile));
%     for k=1:Nc(tile)
%         cellImages(:,:,k)=ImageOut{Start+k-1}; 
%     end
%  
%     InImages{1,i}=cellImages;
%     InSignals{1,i}=TraceOut(Start:End,:);
    InImages{1,i}=output.spatial_weights;
    InSignals{1,i}=double(output.temporal_weights)';
    InTargets{1,i}=sortedCells{tile};
    InMovie{1,i}=strcat('E:\',Mouse,'_tiles\CaMovie_masktiled_',num2str(tile),'.h5');
end
[classifyStruct] = classifySignals(InImages,InSignals,'classifierType','all','trainingOrClassify','training'...
    ,'inputMovieList',InMovie,'inputTargets',InTargets);

save(strcat(ExPath,'CellClassifier'),'classifyStruct','ManualTiles');
%%%%%%%
AutoTiles=[1:16];
load(strcat(ExPath,'\CellClassifier.mat'));
Nc=zeros(16,1);
AutoSortedCells={};


% for i=1:length(AutoTiles)
%     tile=AutoTiles(i)
%     
%     Start=1;
%     if tile>1
%         Start=sum(Nc(1:tile-1))+1;
%     end
%     End=sum(Nc(1:tile));
%     dim=size(ImageOut{Start});
%     cellImages=zeros(dim(1),dim(2),Nc(tile));
%     for k=1:Nc(tile)
%         cellImages(:,:,k)=ImageOut{Start+k-1}; 
%     end
%  
%     InImages{1,i}=cellImages;
%     InSignals{1,i}=TraceOut(Start:End,:);
%     InMovie{1,i}=strcat('E:\',Mouse,'_tiles\CaMovie_masktiled_',num2str(tile),'.h5');
% 
% end

for i=1:length(AutoTiles)
    InImages={};
    InSignals={};
    InMovie={};
    tile=AutoTiles(i)
    load(strcat(ExPath,'AllCells_extract_',num2str(tile)));
    Nc(tile)=size(output.temporal_weights,2);
    InImages{1,1}=output.spatial_weights;
    InSignals{1,1}=double(output.temporal_weights)';
    InMovie{1,1}=strcat('E:\',Mouse,'_tiles\CaMovie_masktiled_',num2str(tile),'.h5');

    
[OutClassifyStruct] = classifySignals(InImages,InSignals,...
    'classifierType','all',...
    'trainingOrClassify','classify',...
    'inputMovieList',InMovie,...
    'inputStruct',classifyStruct);


save(strcat(ExPath,'Auto_extract_sorted_',num2str(tile)),'OutClassifyStruct');

end
















