function [w,v,vhat]=solvejscca(Data,ss_th,lamtv, ss_mtx, v0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the joint sparse canonical correlation analysis with stability
% selection
% Data      - The data contains two types of data Data.m1 and Data.m2: m X p1 for Data.m1 and m X p2 for Data.m2  where m is the number of samples and p1,p2 are the
%          dimension of features. The data should be normalized for each
%          class
% Data.label - A m X 1 vector with positive integers indicates the class label
% of each sample
% ss_th- parameter for sparse level, default is the number of samples.
% lamtv- parameter for class fusion
% ss_mtx- monto carlo sampling matrix, which can avoid 'for' loop to accelerate runing speed.
% v0- initialization [Fixed in stability selection to avoid mismatching].
% The algorithms are summarized in the following papers: 
%
%
% params.max_iter[reccomended: 400]
% Maximal number of iterations allowed.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (2014): Jian Fang, Yu-ping Wang.
% 
% This code is distributed under the terms of the GNU General Public
% License 2.0.
% 
% Permission to use, copy, modify, and distribute this software
% for any purpose without fee is hereby granted, provided that
% this entire notice is included in all copies of any software
% which is or includes a copy or modification of this software
% and in all copies of the supporting documentation for such
% software. This software is being provided "as is", without any
% express or implied warranty. In particular, the authors do not
% make any representation or warranty of any kind concerning the
% merchantability of this software or its fitness for any
% particular purpose."
%
% Contact email: jfang3@tulane.edu, wyp@tulane.edu
% Biomedical Engineering Department, Tulane University, New Orleans, LA 70118, USA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
p1=size(Data.m1,2); % Number of variables in the first class
p2=size(Data.m2,2); % Number of variables in the second class
[~,nss]=size(ss_mtx); % Number of subsampling in stability selection
label=Data.label;
nclass=length(unique(label)); % Number of classes
cv_train_1=cell(nclass,1); 
cv_train_2=cell(nclass,1);
ptw=cell(nclass,1); % Indicater matrix for the selected samples during each subsampling
nc=zeros(nclass,1); % Number of total samples in each class
ntr=zeros(nclass,nss); % Number of samples in each subsampling
for i=1:nclass
    ind_temp=(label==i);
    cv_train_1{i}=Data.m1(ind_temp,:);
    cv_train_2{i}=Data.m2(ind_temp,:);
    ptw{i}=double(ss_mtx(ind_temp,:));
    nc(i)=sum(ind_temp);
    ntr(i,:)=sum(repmat(ind_temp,1,nss)&ss_mtx);
end
% Initialize the loading vector
v=repmat(v0,1,nss,nclass); % Use the same initialization to avoid mismatching.
%% Block coordinate decent for optimization
for ii=1:400
    lam1=max(p1*0.5/1.1^ii,ss_th(1)); 
    lam2=max(p2*0.5/1.1^ii,ss_th(2));
    w=zeros(p1,nss);
    for ic=1:nclass
        w=w+cv_train_1{ic}'*((cv_train_2{ic}*v(:,:,ic)).*ptw{ic}./repmat(ntr(ic,:),nc(ic),1));
    end
    w=soft(w,fix(lam1));
    for ic=1:nclass
        v(:,:,ic)=cv_train_2{ic}'*((cv_train_1{ic}*w).*ptw{ic}./repmat(ntr(ic,:),nc(ic),1));
    end
	if nclass>2
% The solver for multi-class fused lasso is slow in 2-class cases.
    [v,vhat]=flsa(v,lamtv);
	else
	[v,vhat]=flsa_2c(v,lamtv);
	end
    v=vsoft(v,fix(lam2));
end
end


% function to solve fused lasso in two-class cases.
function [vf,vhat]=flsa_2c(v,lam)
dx=(v(:,:,1)-v(:,:,2))/2;
dx=sign(dx).*min(abs(dx),lam);
vf(:,:,1)=v(:,:,1)-dx;vf(:,:,2)=v(:,:,2)+dx;
vhat=vf;
end

% Soft-thresholding with normalization
function y=soft(x,lambda)
n=size(x,1);
temp=sort(abs(x),'descend');
th=temp(lambda+1,:);
y=sign(x).*max(abs(x)-repmat(th,n,1),0);
ny=sqrt(sum(y.^2));
y=y./repmat(ny,n,1);
end

% vSoft-thresholding with normalization
function v=vsoft(x,lam)
[n,k,nc]=size(x);
th_pos=fix(lam);
th_pos=repmat(th_pos,1,k*nc)+(1:n:k*n*nc);
x=reshape(x,n,k*nc,1);
temp=sort(abs(x),'descend');
th=temp(th_pos);
y=sign(x).*max(abs(x)-repmat(th,n,1),0);
v=reshape(y,n,k,nc);
fnorm=sqrt(sum(sum(v.^2,3),1));
v=v./repmat(fnorm,n,1,nc);
end

