### Linear estimation of Fisher information between to distributions samples by Set0=p by N0 and Set1=p by N1 
function [IF]=FisherInformation(Set0,Set1)
 
            m0=mean(Set0,2);
            m1=mean(Set1,2);
            d0=size(Set0,2);
            d1=size(Set1,2);
            
            S0=((Set0-m0*ones(1,d0))*(Set0-m0*ones(1,d0))')/(d0-1);
            S1=((Set1-m1*ones(1,d1))*(Set1-m1*ones(1,d1))')/(d1-1);
            
            S=(S0+S1)/2;
            dm=m1-m0;
            
            
            IF= dm' * inv(S) * dm;
end
