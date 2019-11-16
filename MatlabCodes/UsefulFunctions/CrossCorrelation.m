%%%%%Compute the Cross correlation between input x1 of 1 by L and x2 of 1 by L over the window of -T to T
%%%%Inputs x1=1byL data, x2=1byL data, T scalar time window of -T to T for cross correlations
%%%%output CC= 1by2T+1 cross correlation coefficient, MC= scalar maximum cross correlation coeeficient, MCind= scalar time shift for maximum cross correlations
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
