%%%% after reapplying pca ica filters to the concat movie 
%%%% this code normalizes the traces and also cleans the baseline
AddressSetup;
Win=1000;
Mouse=51;
load(strcat(LoadPath{Mouse},'\CellSorting\cellData_ROI.mat'));
load(strcat(LoadPath{Mouse},'\cellImage.mat'));
load(strcat(LoadPath{Mouse},'\SessionLength.mat'));
cellData_ROINorm=zeros(size(cellData_ROI,1),sum(SessionLength));
cellV=1:size(cellData_ROI,1);

% for cnum=cellV
%     cellData_ROI(cnum,:)=cellData_ROI(cnum,:)/...
%         (sum(sum(cellImage{cnum})) *size(cellImage{cnum},1)*size(cellImage{cnum},2));
% end
baseline={};
startS=1;
for s=1:length(SessionLength)
    s
    endS=sum(SessionLength(1:s));
    parfor f=startS:endS
        baseline{f}=median(cellData_ROI(:,max(startS,f-Win):min(endS,f+Win)),2);
        cellData_ROINorm(:,f)=cellData_ROI(:,f)-baseline{f};
    end
    startS=startS+SessionLength(s);
end

save(strcat(LoadPath{Mouse},'\cellData_ROI'),'cellData_ROINorm','-v7.3');

