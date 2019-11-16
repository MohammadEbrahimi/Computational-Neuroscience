function [IF,W,Signal,Noise]=FisherInformationVal(TSet0,TSet1,VSet0,VSet1)

            m0=mean(TSet0,2);
            m1=mean(TSet1,2);
            d0=size(TSet0,2);
            d1=size(TSet1,2);

            S0=cov(TSet0');
            S1=cov(TSet1');

            S=(S0+S1)/2;

            dm=m1-m0;

            W= inv(S) * dm;
            [IF,Signal,Noise]=FisherInformationDecoder(VSet0,VSet1,W);

end