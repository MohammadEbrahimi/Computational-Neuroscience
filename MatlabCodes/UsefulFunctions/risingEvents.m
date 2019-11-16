function [A]=risingEvents(cellTrace,cellEvents,alpha)



dim =size(cellTrace);
A=zeros(dim);
for i=1:dim(1)
    for j=1:dim(2)
        if (cellEvents(i,j)==1)
            tr=1;
            while(tr < j-1 && ((cellTrace(i,j-tr)-cellTrace(i,j-tr-1))/...
                    (cellTrace(i,j-tr+1)-cellTrace(i,j-tr))) > alpha)
                tr=tr+1;
            end
            
            A(i,j-tr:j)=1;
        end
    end
end
        
end
