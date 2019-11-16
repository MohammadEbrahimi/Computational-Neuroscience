function [P,Xp ] = Distribution(x,n)
P=zeros(n,1);
Xp=zeros(n,1);
bin=(max(x)-min(x))/n
for i=1:n
    P(i)=sum(x( x> (min(x)+((i-1)*bin)) & x <(min(x)+((i)*bin))));
    Xp(i)= (min(x)+((i-1)*bin)) + (bin/2);
end



end

