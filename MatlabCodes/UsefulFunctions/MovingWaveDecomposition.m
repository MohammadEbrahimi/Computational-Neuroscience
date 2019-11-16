function [A,f,p,Xh]=MovingWaveDecomposition(X)
%%% Decomposing X into the form of X~ A .* f(t-p) 
Dim=size(X);
A=rand(Dim(1),Dim(2))*0.01;
f=rand(1,2*Dim(3))-0.5;
p=uint16(rand(Dim(1),Dim(2))*20-10);
NCore=16;
tol=0.02;
maxIter=20;
alpha=0.01;

G = fspecial('gaussian',[30 30],10);


Loss=zeros(maxIter,1);
Xh=zeros(size(X));
error=2*tol;
iter=0;
while error>tol && iter<maxIter
iter=iter+1;
df=zeros(1,2*Dim(3));
dA=zeros(Dim(1),Dim(2));
fprime=diff(f);
fprime=[fprime,fprime(end)];
for i=1:Dim(1)
    parfor (j=1:Dim(2),NCore)
        Xh(i,j,:)=(A(i,j)*f(Dim(3)-p(i,j):(2*Dim(3))-p(i,j)-1));
        dA(i,j)=sum( squeeze((X(i,j,:)-Xh(i,j,:)).*((Xh(i,j,:))./ A(i,j))));
    end 
        for j=1:Dim(2)
            df(Dim(3)-p(i,j):(2*Dim(3))-p(i,j)-1) = df(Dim(3)-p(i,j):(2*Dim(3))-p(i,j)-1)-(squeeze(Xh(i,j,:)-X(i,j,:))' *A(i,j));
        end
    
end

f=f+ (alpha * df);
A=A+ (alpha * dA);
MCInd=zeros(Dim(2),1);
for i=1:Dim(1)
    parfor (j=1:Dim(2),NCore)
        [~,~,MCInd(j)]=CrossCorrelation((A(i,j)*f(Dim(3)-p(i,j):(2*Dim(3))-p(i,j)-1)),squeeze(X(i,j,:)),10);
        p(i,j)=p(i,j)-MCInd(j);
    end
end

p=imfilter(p,G,'same','symmetric');

p(p>=Dim(3))=Dim(3)-1;
error=sqrt(sum(sum(sum((X-Xh).^2)))) / sqrt(sum(sum(sum((X).^2))))
Loss(iter)=error;
end

%%%Reconstruction Xh = A .* f(t-p) 
 Dim=[size(A),length(f)/2];
 Xh=zeros(Dim);
 for i=1:Dim(1)
     for j=1:Dim(2)
         Xh(i,j,:)=(A(i,j)*f(Dim(3)-p(i,j):(2*Dim(3))-p(i,j)-1));
     end 
 end
 
end
