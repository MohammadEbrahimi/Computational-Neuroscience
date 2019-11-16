function [pv,max_w,max_v]=cal_diff_set(Data,w,v,vhat,i1,i2)
label=Data.label;
ntrain=size(Data.m1,1);
nperumu=10000;
ss=calss(w,v,vhat,i1,i2);

%% prefiltering
[pro1,l]=sort(ss.w,'descend');
indt=pro1>0.1;
max_w=l(indt);
[pro2,m]=sort(ss.v,'descend');
indt=pro2>0.1;
max_v=m(indt);
m1=Data.m1(:,max_w);m2=Data.m2(:,max_v);
%% permutation test
temp=find(label==i1);
corr1=cal_corr(Data.m1(temp,max_w),Data.m2(temp,max_v));
temp=find(label==i2);
corr2=cal_corr(Data.m1(temp,max_w),Data.m2(temp,max_v));
pv=0;
for i=1:nperumu
    indb=randperm(ntrain);
    m1t=m1(indb,:);
    m2t=m2(indb,:);
    temp=find(label==i1);
    m1t(temp,:)=zscore(m1t(temp,:));
    m2t(temp,:)=zscore(m2t(temp,:));
    corr1t=cal_corr(m1t(temp,:),m2t(temp,:));
    temp=find(label==i2);
    m1t(temp,:)=zscore(m1t(temp,:));
    m2t(temp,:)=zscore(m2t(temp,:));
    corr2t=cal_corr(m1t(temp,:),m2t(temp,:));
    pv=pv+double((corr1t-corr2t)>(corr1-corr2));
end
pv=pv/nperumu;
end

function corr=cal_corr(m1,m2)
corr=abs(m1'*m2./(sqrt(sum(m1.^2))'*sqrt(sum(m2.^2))));
end

function ss=calss(w,v,vhat,i1,i2)
npermu=size(w,2);
ss.w=sum(abs(w)>0,2)/npermu;
ss.v=sum((abs(v(:,:,i1))>abs(v(:,:,i2)))&(vhat(:,:,i1)~=vhat(:,:,i2)),2)/npermu;
end