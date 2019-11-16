function [X]=clipp(I,a)
if length(size(I))==2
    L=reshape(I,1,size(I,1)*size(I,2));
end
if length(size(I))==3
    L=reshape(I,1,size(I,1)*size(I,2)*size(I,3));
end
L=sort(L);
NT=round(length(L)*a);
LowTh=L(NT);
HighTh=L(length(L)-NT);
X=I;
X(I>HighTh)=HighTh;
X(I<LowTh)=LowTh;
end

