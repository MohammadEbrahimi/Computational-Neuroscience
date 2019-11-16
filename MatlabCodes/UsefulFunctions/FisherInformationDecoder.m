%%% Linear estimation of Binary Fisher information between two validation set samples specified by Set0=p by N0 and Set1=p by N1 
%%% Through a Linear decoder specified by B = p by 1 
%%% It returns Fisher Information in IF , And the Signal and Noise part separately 
function [IF,Signal,Noise]=FisherInformationDecoder(Set0,Set1,B)

            m0=mean(Set0,2);
            m1=mean(Set1,2);
            d0=size(Set0,2);
            d1=size(Set1,2);

            S0=cov(Set0');
            S1=cov(Set1');

            S=(S0+S1)/2;
            dm=m1-m0;

           Signal= (B'* dm)^2;
           Noise= B' * S * B;

           IF= Signal / Noise;
end
