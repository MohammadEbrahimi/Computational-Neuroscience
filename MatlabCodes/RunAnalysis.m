AddressSetup;
savePath='C:\Users\Sadegh\Documents\VLMReborn\Reports2\2017_05_19_Fisher_PLS_HoldOut\GoNogo_alldata\';

Sessions=[17,18,19,47,48,49,50];

for s=Sessions
    Fisher_PLS_HoldOut(LoadPath{s},strcat(savePath,'Session',num2str(s)),1,1,0);
    Fisher_PLS_HoldOut(LoadPath{s},strcat(savePath,'Session',num2str(s)),2,1,0);
    %Fisher_PLS_HoldOut(LoadPath{s},strcat(savePath,'Session',num2str(s)),3,1,0);
    
end
