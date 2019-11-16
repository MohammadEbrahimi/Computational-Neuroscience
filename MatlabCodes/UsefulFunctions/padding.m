function [A]=padding(B,s,e)
dim=size(B);
A=zeros(dim(1),dim(2));
for i=1:dim(1)
    for j=1:dim(2)
        if (B(i,j)==1)
            A(i,max(1,j-s):min(j+e,dim(2)))= A(i,max(1,j-s):min(j+e,dim(2)))+1;
        end
    end
end
