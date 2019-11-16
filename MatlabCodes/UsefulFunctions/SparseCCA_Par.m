function [maxwx,maxwy,cxVec,cyVec,rVal,rTrain]=SparseCCA_Par(X,Y,XVal,YVal,reg)
%  X=dd0(CF,:)';
%  Y=dd0(CM,:)';
%  XVal=dd0v(CF,:)';
%  YVal=dd0v(CM,:)';
%%%%%%%% X= n x p1   Y=n x p2
nSearch=20;
dx=size(X);
dy=size(Y);
cxVec=logspace(log10(sqrt(5)),log10(sqrt(dx(2))),nSearch);
cyVec=logspace(log10(sqrt(5)),log10(sqrt(dy(2))),nSearch);
%cxVec=linspace((sqrt(5)),(sqrt(dx(2))),nSearch);
%cyVec=linspace((sqrt(5)),(sqrt(dy(2))),nSearch);
rVal=zeros(nSearch,nSearch);
rTrain=zeros(nSearch,nSearch);

X=X-repmat(mean(X),[dx(1),1]);
X=X./repmat(sqrt(var(X)),[dx(1),1]);
Y=Y-repmat(mean(Y),[dy(1),1]);
Y=Y./repmat(sqrt(var(Y)),[dy(1),1]);
X(isnan(X))=0;
Y(isnan(Y))=0;


wx={};
wy={};
wx{20,20}=[];
wy{20,20}=[];
c={};
for ix=1:nSearch
    parfor iy=1:nSearch

        [wx{ix,iy},wy{ix,iy},rTrain(ix,iy)]=SparseCCA(X,Y,cxVec(ix),cyVec(iy),reg);
        c{iy}=corrcoef(XVal*wx{ix,iy},YVal*wy{ix,iy});
        rVal(ix,iy)=c{iy}(1,2);

    end
end

[~,indx]=max(max(rVal,[],2));
[~,indy]=max(max(rVal,[],1));
maxwx=wx{indx,indy};
maxwy=wy{indx,indy};


end