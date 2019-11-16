function [TimeStamp,Er]=correctImageFrame(IF,T,Nf)
step=max(T)/Nf;
TimeStamp=0:step:(Nf-1)*step;

RiseTime=T([diff(IF),0]>0);

Er=zeros(size(TimeStamp));
for i=1:length(TimeStamp)
    [Er(i),ind]=min(abs(TimeStamp(i)-RiseTime));
    TimeStamp(i)=RiseTime(ind);
    RiseTime(ind)=max(TimeStamp)*2;
end


end
