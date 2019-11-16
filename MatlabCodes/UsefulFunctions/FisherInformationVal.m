%%%%% Computing Fisher Information (IF) on the Validation set (VSet0=pbyN0,VSet1=pbyN1) through a Fisher decoder (W=pby1) optimized on the training Set (TSet0=pbyN00,TSet1=pbyN11)
%%%%% The Signal and Noise part is returned separately 
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
