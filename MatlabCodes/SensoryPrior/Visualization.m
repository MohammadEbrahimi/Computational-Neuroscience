p='E:\Reports2\2018_03_01_FisherInfo\Raw_Bin\Datasets--15to20\';
FACR=zeros(9,7);
HitMiss=zeros(9,7);

for Mouse=61:67
    
load(strcat(p,'FACRPreStim\Session',num2str(Mouse),'_mode1.mat'));
[res,ind]=max(mean(Info_Val_Fisher,3),[],2);
res_sh=zeros(9,1);
for i=1:9
    res_sh(i)=max(Info_Val_Fisher_Sh(i,ind(i),:));
end
FACR(res>res_sh,Mouse-60)=res(res>res_sh);

load(strcat(p,'HitMissPreStim\Session',num2str(Mouse),'_mode1.mat'));
[res,ind]=max(mean(Info_Val_Fisher,3),[],2);
res_sh=zeros(9,1);
for i=1:9
    res_sh(i)=max(Info_Val_Fisher_Sh(i,ind(i),:));
end
HitMiss(res>res_sh,Mouse-60)=res(res>res_sh);


end