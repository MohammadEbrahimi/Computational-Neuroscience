%% Initialization
%Generate Data
p1=50;p2=50;pt=5000;nclass=2;pd=20;
n=100;
sig1=1;sig2=1;%noise level
a=zeros(pt,1);
a(1:p1)=randu(p1);
b=zeros(pt,1);
b(1:p2)=randu(p2);
Data.m1=[];Data.m2=[];
for ic=1:nclass
bt=b;
bt(p2+(ic-1)*pd+1:p2+ic*pd)=randu(pd);
ind=false(pt,1);
ind(1:p2)=true;
ind(p2+(ic-1)*pd+1:p2+ic*pd)=true;
mu=randn(n,1);
mu=sign(mu).*(abs(mu)+0.1);
mu=zscore(mu);
mu=mu*sig1;% Common feature
m1=mu*a'+randn(n,pt)*sig2;
m2=mu*bt'+randn(n,pt)*sig2;
h_ind((ic-1)*n+1:ic*n)=ic;
Data.m1=[Data.m1;m1];Data.m2=[Data.m2;m2];
end
seq=1:nclass*n;
seq=reshape(seq,n,nclass);
seq=reshape(seq',n*nclass,1);
Data.m1=Data.m1(seq,:);
Data.m2=Data.m2(seq,:);
Data.label=h_ind(seq)';
%% Run JSCCA
opts.nss=100;opts.lamtv=0.1;
opts.ss_th(1)=300;opts.ss_th(2)=300;
[~,opts.v0]=svds_initial(Data.m1,Data.m2,opts.ss_th);
[w,v,vhat]=jscca_ss(Data,opts);
[ss_stats,lociw,lociv]=cal_diff_set(Data,w,v,vhat,1,2);
module=findmodule(ss_stats,lociw,lociv,0.05);



%% Differential Module Detection
%%


