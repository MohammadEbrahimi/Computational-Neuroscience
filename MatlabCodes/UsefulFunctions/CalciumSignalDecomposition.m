%%%This code gets a single cell trace on dimension 1,n and outputs
%%%outputs S=calcium events N=Noise and B=Baseline
%%%you can change the mWindow parameter for baseline computation 
%%%Sadegh Ebrahimi
function [S,N,B]=CalciumSignalDecomposition(x)
tol=0.001;
mWindow=100;

L=length(x);
B=zeros(1,L);%Baseline
N=zeros(1,L);%Noise
S=zeros(1,L);%Signal
pth=0;
th=3*sqrt(var(x));
iter=1;
sd=0;
rate=0;
while abs(th-pth)>tol
    sd(iter)=th/3;
    
    S = x-B;
    S(S<th) = 0;
    rate(iter)=sum(S>0)/length(S);
    
    T=x-S;
    for i=1:L
        %%%%mean or median
        B(i)=mean(T(max(1,i-(mWindow/2)):min(L,i+(mWindow/2))));
    end
    N=(x-B)-S;
    pth=th;
    th=3*sqrt(var(N));
    
    iter=iter+1;
end
end   
