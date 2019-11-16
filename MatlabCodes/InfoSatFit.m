function [Nhalf,errorVec]=InfoSatFit(Cells,Info)
Im=max(Info);
errorVec=zeros(max(Cells),1);
for c=1:max(Cells)
    eps=1/double(c);
   for i=1:length(Cells)
        errorVec(c)=errorVec(c)+((((eps*Im*Cells(i))/(1+eps*Cells(i)))-Info(i))^2);    
   end

end
[~,Nhalf]=min(errorVec);
end