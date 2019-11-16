% p='C:\Users\Sadegh\Documents\VLMReborn\L365\Data\ConcatDays\';
% load(strcat(p,'cellData.mat'));
% load(strcat(p,'SessionLength.mat'));
% cellData=cellData_Raw;
% cellData_Raw_Events=zeros(size(cellData),'single');
RWB=3;
RWA=7;
for s=1:length(SessionLength)
    s
    cellData(:,sum(SessionLength(1:s-1))+1:sum(SessionLength(1:s)))=...
        remove_baseline(cellData(:,sum(SessionLength(1:s-1))+1:sum(SessionLength(1:s))),1);
   for cnum=1:size(cellData,1)
       x=cellData(cnum,sum(SessionLength(1:s-1))+1:sum(SessionLength(1:s)));
       N=[x(x<=0),-x(x<0)];
       sd=sqrt(var(N));
       y=x;
       y(x<2*sd)=0;
       z=zeros(size(y),'single');
       
       
       LastMin=1;
       LastMax=1;
       Detected=0;
       for i=1:length(y)
           if y(i)<=min(y(max(i-RWB,1):min(i+RWB,length(y))) )
               if (Detected==1 && y(LastMin)>=y(i)) || Detected==0
                 LastMin=i;
                 Detected=1;
               end
           end
           if y(i)>=max(y(max(i-RWA,1):min(i+RWA,length(y))) )
               LastMax=i;
           end
           if LastMax>LastMin && Detected==1 &&  x(LastMax)-x(LastMin)>sd
               Detected=0;
               z(LastMin:LastMax)=y(LastMin:LastMax);
           end
       end
       cellData_Raw_Events(cnum,sum(SessionLength(1:s-1))+1:sum(SessionLength(1:s)))=z;

               

   end
end
save(strcat(p,'cellData_Raw_RT'),'cellData_Raw_Events','-v7.3');