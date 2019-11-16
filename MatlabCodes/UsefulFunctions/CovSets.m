%%%%Computing average covariance matrix for Set0 of p by n0 and and Set1 of p by n1
function [S]=CovSets(Set0,Set1)
 
            m0=mean(Set0,2);
            m1=mean(Set1,2);
            d0=size(Set0,2);
            d1=size(Set1,2);
            
            S0=((Set0-m0*ones(1,d0))*(Set0-m0*ones(1,d0))')/(d0-1);
            S1=((Set1-m1*ones(1,d1))*(Set1-m1*ones(1,d1))')/(d1-1);
            
            S=(S0+S1)/2;
end
