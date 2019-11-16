function [CellTraces,CellImages]=PCA_ICA_Tiled(hdf5In,OutPath,Nc)
  NCells{1}=round([25,85,25,5,...
    40,85,90,20,...
    70,200,150,35,...
    15,75,30,16]*4/3);%%L364
NCells{2}=round([65,125,140,60,...
    450,500,500,100,...
    300,650,400,400,...
    35,150,150,15]*11/10);%%L365
NCells{3}=round([40,130,65,10,...
    150,175,200,35,...
    150,350,200,35,...
    35,150,75,20]*11/10);%%L367
NCells{4}=round([175,250,250,150,...
    450,500,550,400,...
    470,550,500,350,...
    110,475,500,100]*11/10);%%L368
      

CellTraces=[];
CurCell=1;

for T=1:length(Nc)
    
hinfo=hdf5info(strcat(hdf5In,T,'.h5'));
movie=single(hdf5read(hinfo.GroupHierarchy.Datasets(1)));

[IcaTraces,IcaFilters]=PCA_ICA(movie,2*Nc(T),Nc(T));

if T==1
    CellTraces=IcaTraces;
else
    CellTraces=[CellTraces;IcaTraces];
end
for n=1:Nc(T)
    CellImages{CurCell}=squeeze(IcaFilters(n,:,:));
    CurCell=CurCell+1;
end
end
save(strcat(OutPath,'AllCells'),'CellTraces','CellImages','-v7.3')

end


end