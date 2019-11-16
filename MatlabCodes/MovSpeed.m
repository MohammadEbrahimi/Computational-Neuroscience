cellEvents=padding(eventBin,0,2);
set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',12);

%%%%Find cell Coordinates and assign area
cellIJ=zeros(cellCount,2);
CortexArea=zeros(cellCount,1);

for i=1:cellCount
    [row,ind]=max(cellImage{i});
    [~,indj]=max(row);
    indi=ind(indj);
    cellIJ(i,:)=TileTopLeft(i,:)+[indi-1,indj-1];
end

for a=1:8
    vertex=cortexMap{a};
    [in,on] = inpolygon(cellIJ(:,1),cellIJ(:,2),vertex(:,1),vertex(:,2));
    CortexArea(in==1)=a;
end


%%%%%%%%%%%%

Speed3=sqrt((XSpeed.^2)+(YSpeed.^2));
Resp3=zeros(8,length(IC_time));

for area=1:8
    Resp3(area,:)=mean(cellEvents( (CortexArea==area),:));
 
end

Resp=[Resp1,Resp2,Resp3];
Speed=[Speed1;Speed2;Speed3];

RHO=zeros(8,1);
PVAL=zeros(8,1);
r=zeros(8,1);
m=zeros(8,1);
b=zeros(8,1);
for area=1:8
    [RHO(area),PVAL(area)] = corr(Resp(area,:)',Speed,'Type','Spearman');
end

plot(Resp(8,:),Speed,'.')
title('RSC Spearsman Correlation Rho=-0.057 P<1e-26') 
xlabel('Average Firing Rate')
ylabel('Movement Speed')