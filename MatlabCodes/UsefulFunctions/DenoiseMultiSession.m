clear all
AddressSetup;
AllDays=17;
Th=0.5;
A=open(strcat(LoadPath{AllDays},'\All_Sessions.mat'));
B=open(strcat(LoadPath{AllDays},'\registerationUnion.mat'));
C=open(strcat(LoadPath{AllDays},'\AllCellsUnion.mat'));
%load(strcat(LoadPath{AllDays},'\Calcium.mat'));
load(strcat(LoadPath{AllDays},'\areas.mat'));


AR=C.CellsActiveRatio;
cellNumbers=find(AR>Th);
cellData=C.AllCells;
SessionLength=A.SessionLength;
cellIndex=B.cellIndex;
cellIm=B.cellImage;


Dim=[length(cellNumbers),size(cellData,2)];
cellData_Calcium=zeros(Dim);
cellData_Raw=zeros(Dim);
cellData_Noise=zeros(Dim);
cellData_Baseline=zeros(Dim);
cellData_Mask=zeros(Dim);
cellImage={};




start=1;
for session=1:length(SessionLength)
    session
    timePeriod=start:start+SessionLength{session}-1;
    for cellind=1:length(cellNumbers)
        cellind
        cnum=cellNumbers(cellind);
        if cellIndex(cnum,session)>0
            [S,N,B]=CaDecomposition(cellData(cnum,timePeriod));
            cellData_Calcium(cellind,timePeriod)=S;
            cellData_Raw(cellind,timePeriod)=cellData(cnum,timePeriod);
            cellData_Noise(cellind,timePeriod)=N;
            cellData_Baseline(cellind,timePeriod)=B;
            cellData_Mask(cellind,timePeriod)=1;
        end
        
    end
    start=start+SessionLength{session};
end



for cellind=1:length(cellNumbers)
    cellind
    cnum=cellNumbers(cellind);
    cellImage{cellind}=cellIm{cnum};
    start=1;
    Noise=cellData_Noise(cellind,cellData_Mask(cellind,:)==1);
    for session=1:length(SessionLength)
        timePeriod=start:start+SessionLength{session}-1; 
        if cellIndex(cnum,session)==0
            Noise1=Noise(randperm(length(Noise),length(timePeriod)));
            cellData_Raw(cellind,timePeriod)=Noise1;
        end
       start=start+SessionLength{session}; 
    end
    
end



TileTopLeft=B.TileTopLeft(cellNumbers,:);
se=strel('disk',10,8);
cellCount=size(cellData_Raw,1);
cellIJ=zeros(cellCount,2);
CortexArea=zeros(cellCount,1);
for i=1:cellCount
    [row,ind]=max(cellImage{i});
    [~,indj]=max(row);
    indi=ind(indj);
    cellIJ(i,:)=TileTopLeft(i,:)+[indi-1,indj-1];
    for a=[8,7,5,1,6,2,3,4]
        mask=imdilate(Area(:,:,a),se);
        if mask(cellIJ(i,1),cellIJ(i,2))
            CortexArea(i)=a;
        end
    end
end
CortexArea(CortexArea==0)=9;


ActiveRatio=AR(cellNumbers);
save(strcat(LoadPath{AllDays},'\Calcium.mat'),'cellData_Calcium','cellData_Raw','cellData_Noise','cellData_Baseline','CortexArea','cellNumbers','ActiveRatio','TileTopLeft','-v7.3')