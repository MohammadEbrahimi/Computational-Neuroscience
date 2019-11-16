 function [movC,FFTweight]= filterVisStim(movie,frameRate,stimLength)
dim = size(movie)
l=dim(3)
p=stimLength*frameRate;
Nh=32;%number of harmonies
Fs=frameRate;
NFFT=l;


f=Fs/2*linspace(0,1,NFFT/2+1);
fj=NFFT/p;
Fp=zeros(1,2*Nh+2);
Fp(1)=2;
Fp(2*Nh+2)=NFFT;
for i=1:Nh
    Fp(i+1)=Fp(i)+fj;
    Fp(2*Nh+2-i)=Fp(2*Nh+3-i)-fj;
end

H(1)=0;
for k=2:NFFT
    H(k)=1 / (max(min(abs(Fp(1,:)-k)),1));
end

sigma2 = 2*p;
filtSize2 = 6*p;
x2 = linspace(-filtSize2 / 2, filtSize2 / 2, filtSize2);
gaussFilter2 = exp(-x2 .^ 2 / (2 * sigma2 ^ 2));
gaussFilter2 = gaussFilter2 / sum (gaussFilter2); % normalize


FFTweight=zeros(dim(1),dim(2));

for i=1:dim(1)
    percentProgress = round( 100* i/dim(1))
    for j=1:dim(2)
  
        y=squeeze(movie(i,j,:));
        Y=fft(y,NFFT)/l;
        
        Y2=Y .* H';
        
        y2=squeeze(real(ifft(Y2,NFFT)));
        movC(i,j,:)= y2 - conv (squeeze(y2), gaussFilter2, 'same');
        
        FFTweight(i,j)=sum(abs(Y2).^2)/sum(abs(Y).^2);
        
    end
end

end

