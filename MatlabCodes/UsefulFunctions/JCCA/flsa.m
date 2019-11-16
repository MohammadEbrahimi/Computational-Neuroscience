function [y,yf]=flsa(A,lam1)
%% Initialization
A=permute(A,[3,1,2]);
[K,p,nl]=size(A);
A=reshape(A,K,p*nl,1);

[A_sort,A_order]=sort(A);
pos=0:1:(p*nl-1);
rppos=repmat(pos,K,1)*K;
A_order=A_order+rppos;
trueA=A_sort;
newc=1:K;
newc=repmat(newc',1,p*nl)*2-(K+1);
fusions=false(K,p*nl);
y=zeros(size(A));
%% Begin loop
for iter=1:K-1
    [~,ordermats]=sort(A_sort);
    betasg=A_sort-lam1*newc;
    [~,ordermats_new]=sort(betasg);
    ordermats_new(ordermats_new+rppos)=ordermats;
    fusions(1:K-1,:)=((diff(ordermats_new)-diff(ordermats))<0)|(diff(A_sort)==0)|fusions(1:K-1,:);
%     fusions
    A_sort=fmean(trueA,double(fusions));
   newc=fmean(ordermats,double(fusions))*2-(K+1);
end
y(A_order)=A_sort-lam1*newc;
y=reshape(y,K,p,nl);
y=permute(y,[2,3,1]);
yf=y;
end