function [CC,MC,MCind]=CrossCorrelation(x1,x2,T)
CC=zeros(1,(2*T)+1);
L=length(x1);

for dt= (-1*T):T
                bcor=corrcoef(x1(max(1+dt,1):min(L,L+dt)),x2(max(1-dt,1):min(L,L-dt)));
                CC(T+1+dt)= bcor(1,2);
end

[MC,ind]=max(CC);
MCind=ind-T-1;


end
