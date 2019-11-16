function module=findmodule(ss_stats,lociw,lociv,thresh)
n1=length(lociw);
n2=length(lociv);
msize=0;
for i=1:n1
    for j=1:n2
        fdr=sum(sum(ss_stats(1:i,1:j)>0.05))/i/j;
        if fdr<thresh && i*j>msize
            msize=i*j;
            module.w=lociw(1:i);
            module.v=lociv(1:j);
        end
    end
end