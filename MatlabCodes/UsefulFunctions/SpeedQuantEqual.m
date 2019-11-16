function [SQ,binTh]=SpeedQuantEqual(Speed,ValidTime,Nbin)
bs=floor(sum(ValidTime)/Nbin);
binTh=zeros(1,Nbin-1);
SQ=zeros(size(Speed));
sortedSpeed=sort(Speed(ValidTime));
for b=2:Nbin-1
    binTh(b)=sortedSpeed((b-1)*bs+1);
end

binTh=binTh(diff(binTh)>0);
for b=1:length(binTh)
    SQ=SQ+(Speed>=binTh(b));
end
end