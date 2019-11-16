cellEvents=padding(eventBin,0,2);
Tone=zeros(length(RewardWindow),1);
Tone(:,1)=(diff([RewardWindow;0])==1);
Tone=(padding(Tone',0,6))';
WaitPeriod=padding((diff([RewardWindow;0]')==1),5,0)';
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

n=1;
group=zeros(sum(GoTrials+NogoTrials),1);
for i=1:length(GoTrials)
    if(GoTrials(i)==1)
        group(n)=1;
        n=n+1;
    end
    
    if(NogoTrials(i)==1)
        group(n)=0;
        n=n+1;
    end
    
    
end


[B,intercept,MissCurve,FalseAlarmCurve,Lambda]=lassoglmcv(cellEvents((CortexArea==1),(GoTrials+NogoTrials==1)),group,5)

% [pVals MI] = Mutual_Information_Shuffled(cellEvents',Tone,2,1,1000);
% FMI=MI;
% FMI(pVals>0.01)=0;
% [sortedMI cellIndexMI]=sort(FMI,'descend');
% 
% cellMap=zeros(1017,1017);
% for i=1:100%sum(sortedMI>0)
%     s=size(cellImage{cellIndexMI(i)});
%     ind=TileTopLeft(cellIndexMI(i),:);
%     cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1)=...
%         cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1)+(sortedMI(i)*cellImage{cellIndexMI(i)});
% end
% figure();imagesc(cellMap);