%%%% Quantize each row of X=p by L in a way that each of the NQ bins has equal number of samples in them
%%%% it returns the quantized version of the data NQX=b by L and the quantization thresholds binTh=p by NQ-1
function [NQX,binTh]=QuantEqual(X,NQ)
dim=size(X);
bs=floor(dim(2)/NQ);
binTh=zeros(dim(1),NQ-1);
NQX=zeros(size(X));

for i=1:dim(1)
    S=sort(X(i,:));
    
    for b=1:NQ
        binTh(i,b)=S((b-1)*bs + 1);
        NQX(i,:)=NQX(i,:)+(X(i,:)>=binTh(i,b));
    end
    

end
end
