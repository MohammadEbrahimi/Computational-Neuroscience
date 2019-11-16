function [wxMat,wyMat,rVec]=SparseCCA(X,Y,cx,cy,reg,numCC)
% X=dd0(CortexArea==af & NearCortexArea(:,am)==0,:)';
% Y=dd0(CortexArea==am & NearCortexArea(:,af)==0 ,:)';
%  cx=par.cxVec{af,am}(indx);
%  cy=par.cyVec{af,am}(indx);
%%% X= n x p1   Y=n x p2
dx=size(X);
dy=size(Y);
wxMat=zeros(dx(2),numCC);
wyMat=zeros(dy(2),numCC);
rVec=zeros(numCC,1);



xmask=max(X)>min(X);
ymask=max(Y)>min(Y);

X=X(:,xmask);
Y=Y(:,ymask);

X=X-repmat(mean(X),[dx(1),1]);
X=X./repmat(sqrt(var(X)),[dx(1),1]);
Y=Y-repmat(mean(Y),[dy(1),1]);
Y=Y./repmat(sqrt(var(Y)),[dy(1),1]);
tol=1e-3;
X(isnan(X))=0;
Y(isnan(Y))=0;


dx=size(X);
dy=size(Y);






for icc=1:numCC

wx=randn(dx(2),1);
wy=randn(dy(2),1);
wx=wx/norm(wx);
wy=wy/norm(wy);

wxp=zeros(dx(2),1);
wyp=zeros(dy(2),1);
iter=0;

while (norm(wx-wxp)/norm(wx)) > tol && (norm(wy-wyp)/norm(wy)) > tol
    wxp=wx;
    wyp=wy;

    wx = X' * Y *wy;
    if reg==1
        wx=softTreshold(wx,cx);
    elseif reg==2
        wx=wx/norm(wx);
    end
    
    wy=(wx' * X' * Y)';
    wy(isnan(wy))=0;
    if reg==1
        wy=softTreshold(wy,cy);
    elseif reg==2
        wy=wy/norm(wy);
    end

    iter=iter+1;
    if iter>500
        break;
    end

end
c=corrcoef([X * wx , Y * wy]);
r=c(1,2);

wxMat(xmask,icc)=wx;
wyMat(ymask,icc)=wy;
rVec(icc)=r;

X=X-(X*wx*wx');
Y=Y-(Y*wy*wy');

end  


end