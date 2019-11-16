
function M=downsamp2d(M,bindims)
%DOWNSAMP2D - simple tool for 2D downsampling
%
%  M=downsamp2d(M,bindims)
%
%in:
%
% M: a matrix
% bindims: a vector [p,q] specifying pxq downsampling
%
%out:
%
% M: the downsized matrix

OriginalClass=class(M);
p=bindims(1); q=bindims(2);

[m,n]=size(M); %M is the original matrix
newm=p*floor(m/p);
newn=q*floor(n/q);
M=M(1:newm,1:newn);
m=newm;
n=newn;

M=sum(reshape(M,p,[]) ,1 );
M=reshape(M,m/p,[]).'; %Note transpose

M=sum(reshape(M,q,[]) ,1);
M=reshape(M,n/q,[]).'; %Note transpose

M=cast(M/(p*q),OriginalClass);
end