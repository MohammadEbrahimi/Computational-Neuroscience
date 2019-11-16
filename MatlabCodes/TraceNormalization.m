
AddressSetup;
for Mouse=2:7
load(strcat(LoadPath{50+Mouse},'\cellData.mat'));
load(strcat(LoadPath{50+Mouse},'\SessionLength.mat'));
start=1;
cellCount=size(cellData_Raw,1);
ActiveRatio=zeros(cellCount,1);
cellData_Normalized=zeros(size(cellData_Calcium));

for s=1:length(SessionLength)
    s
    TP=start:sum(SessionLength(1:s));
    for cnum=1:cellCount
        noise=cellData_Noise(cnum,:);
        cstd=sqrt(var(noise(cellData_Noise(cnum,:)~=0)));
        Amp=max(cellData_Calcium(cnum,:));
        if sum(cellData_Calcium(cnum,TP)> 3*cstd)/length(TP) > 0.003
            ActiveRatio(cnum)=ActiveRatio(cnum)+length(TP);
            signal=cellData_Calcium(cnum,TP);
            M=mean(signal(cellData_Calcium(cnum,TP)> 3*cstd));
            cellData_Normalized(cnum,TP)=cellData_Calcium(cnum,TP) * (Amp / M);
        else
            %cellData_Normalized(cnum,TP)=cellData_Calcium(cnum,TP);
        end
        
        
    end
    start=start+SessionLength(s);
end

ActiveRatio=ActiveRatio/sum(SessionLength);
save(strcat(LoadPath{50+Mouse},'\cellData_Normalized'),'cellData_Normalized','ActiveRatio','-v7.3');
end