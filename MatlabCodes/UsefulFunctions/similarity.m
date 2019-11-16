%% Correlation coefficient matix S for columns in the matrix X
function [S]=similarity(X)
S=zeros(size(X,2));
for i=1:size(X,2)
    for j=1:size(X,2)
        S(i,j)=(X(:,i)' / norm(X(:,i))) * (X(:,j) / norm(X(:,j)));
    end
end
if size(X,2)==2
    S=S(1,2);
end

end
