Mouse={'L347','L354','L362bis','L364','L365','L367','L368'};
for i=1:7
p=strcat('C:\Users\Sadegh\Documents\VLMReborn\',Mouse{i},'\Data\ConcatDays\');
load(strcat(p,'cellData_ROI.mat'));
load(strcat(p,'cellData.mat'));
load(strcat(p,'cellData_RT.mat'));
load(strcat(p,'cellData_Raw_RT.mat'));
load(strcat(p,'SessionLength.mat'));

cellData_Raw=remove_baseline(cellData_Raw,1);
cellData_ROINorm=remove_baseline(cellData_ROINorm,1);

cellData_Raw_PerSession=zeros(size(cellData_Raw),'single');
cellData_ROI_PerSession=zeros(size(cellData_Raw),'single');
cellData_Raw_bin=zeros(size(cellData_Raw),'single');
cellData_ROI_bin=zeros(size(cellData_Raw),'single');

for s=1:length(SessionLength)
    s
    STime=sum(SessionLength(1:s-1))+1:sum(SessionLength(1:s));
    cellData_Raw_PerSession(:,STime)=zscore(cellData_Raw(:,STime)')';
    cellData_ROI_PerSession(:,STime)=zscore(cellData_ROINorm(:,STime)')';
  
end
cellData_Raw=zscore(cellData_Raw')';
cellData_ROI=zscore(cellData_ROINorm')';

cellData_Raw_bin=cellData_Raw_Events>0;
cellData_ROI_bin=cellData_Events>0;

Note='in per sessions Individual Sessions are zscored separately';
save(strcat(p,'cellData_ZS'),'cellData_Raw','cellData_ROI','cellData_Raw_PerSession','cellData_Raw_bin','cellData_ROI_PerSession','cellData_ROI_bin','-v7.3');
end