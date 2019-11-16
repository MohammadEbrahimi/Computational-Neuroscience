function [w,v,vhat]=jscca_ss(Data,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the joint sparse canonical correlation analysis with stability
% selection
% Data      - The data contains two types of data Data.m1 and Data.m2: m X p1 for Data.m1 and m X p2 for Data.m2  where m is the number of samples and p1,p2 are the
%          dimension of features. The data should be normalized for each
%          class
% Data.label - A m X 1 vector with positive integers indicates the class label
% of each sample
% opts.ss_th- parameter for sparse level, default is the number of samples.
% opts.lamtv- parameter for class fusion
% opts.v0- initialization [Fixed in stability selection to avoid mismatching].
% The algorithms are summarized in the following papers: 

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
ntrain=size(Data.m1,1);
label=Data.label;
%% JSCCA with stability selection
%% Normalization
for i=1:length(unique(label))
    tind=(label==i);
    Data.m1(tind,:)=zscore(Data.m1(tind,:));
    Data.m2(tind,:)=zscore(Data.m2(tind,:));
end
% get permutation
ss_mtx=rand(ntrain,opts.nss);
temp=sort(ss_mtx,'descend');
th=temp(fix(ntrain/2),:);
ss_mtx=ss_mtx>repmat(th,ntrain,1);
[w,v,vhat]=solvejscca(Data,opts.ss_th,opts.lamtv, ss_mtx,opts.v0);
end




