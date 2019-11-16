MapTimeL364={[100:340],[100:300],[80:290],[1:427]};
MapTimeL365={[100:300],[100:300],[90:300],[1:420]};
MapTimeL367={[100:340],[1:280],[90:350],[120,420]};
MapTimeL368={[100:330],[1:330],[90:370],[120:420]};
% 
% X=BinningMovie2(DFF(:,:,100:340),2,1);
% Dim=size(X);
% A=rand(Dim(1),Dim(2))*0.01;
% f=rand(1,2*Dim(3))-0.5;
% p=uint16(rand(Dim(1),Dim(2))*20-10);
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

%%%Reconstruction
% Dim=[size(A),length(f)/2];
% Xh=zeros(Dim);
% for i=1:Dim(1)
%     for j=1:Dim(2)
%         Xh(i,j,:)=(A(i,j)*f(Dim(3)-p(i,j):(2*Dim(3))-p(i,j)-1));
%     end 
% end
% 



%   smp0=imfilter(p0,G,'same','symmetric');
% smp90=imfilter(p270,G,'same','symmetric');
%   smp0=imfilter(smp0,G,'same','symmetric');
% smp90=imfilter(smp90,G,'same','symmetric');
% [Gmag90, Gdir90] = imgradient(smp90,'prewitt');
% [Gmag0, Gdir0] = imgradient(smp0,'prewitt');
% 
% map=sind(Gdir0-Gdir90);
% figure();imagesc(map);colormap('Jet');
%     