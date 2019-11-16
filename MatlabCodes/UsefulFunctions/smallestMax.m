function [ind]=smallestMax(X,a)
X=squeeze(X);
ind=1;
absMax=max(X);
for i=1:length(X)
if (abs(X(i)-absMax)/absMax) <a
    ind=i;
    break;
end

end
end