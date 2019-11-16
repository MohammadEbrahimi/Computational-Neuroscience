clear all
AddressSetup;
Sessions=50;
savePath='C:\Users\Sadegh\Documents\VLMReborn\Reports2\2017_05_20_BehavioralAnalysis\';



load(strcat(LoadPath{Sessions},'\All_Sessions.mat'));
SpeedTh=50;
ActiveTrialNumber=20;

ActiveAnimal=ones(length(Lick),1);
for i=1:(length(Lick)-75*ActiveTrialNumber)
    if max(Lick(i:(i+75*ActiveTrialNumber)))==0
        ActiveAnimal(i:(i+75*ActiveTrialNumber))=0;
    end
end


Result=zeros(4,3);
clear ActiveMask
%        N-Clean     N-Running    N-Disengage 
%Hit    
%CR
%Miss
%FA    
er=0;
for c=1:4
    switch c
        case 1
            SE=HitSE;
        case 2
            SE=CRSE;
        case 3
            SE=MissSE;
        case 4
            SE=FASE;
    end
mask=ones(size(SE,1),1);
for k=1:(size(SE,1)-1)
    if max(ActiveAnimal(SE(k,1):SE(k,2)))==0
        mask(k)=0;
        Result(c,3)=Result(c,3)+1;
    elseif max (Speed(SE(k,1):SE(k,2)+55)) > SpeedTh
        Result(c,2)=Result(c,2)+1;
    else
        Result(c,1)=Result(c,1)+1;
    end
    
end

ActiveMask{c}=mask;
end
W=30*(60*10);  
HFTraces=zeros(2,length(Lick)-W);
for i=1:length(Lick)-W
    nh=sum(HitSE(ActiveMask{1}==1,1)> i & HitSE(ActiveMask{1}==1,1)<i+W);
    nc=sum(CRSE(ActiveMask{2}==1,1)> i & CRSE(ActiveMask{2}==1,1)<i+W);
    nm=sum(MissSE(ActiveMask{3}==1,1)> i & MissSE(ActiveMask{3}==1,1)<i+W);
    nf=sum(FASE(ActiveMask{4}==1,1)> i & FASE(ActiveMask{4}==1,1)<i+W);
    
    HFTraces(:,i)=[nh/(nh+nm); nf/(nf+nc)];
    
    
end

HF=zeros(2,1);
HF(1)=sum(Result(1,1:2))/(sum(Result(1,1:2))+sum(Result(3,1:2)));
HF(2)=sum(Result(4,1:2))/(sum(Result(4,1:2))+sum(Result(2,1:2)));

Disengage=sum(Result(:,3))/sum(sum(Result));

save(strcat(savePath,'Session',num2str(Sessions)),'HF','HFTraces','Disengage','Result');

%%%%%%%%%%%%%%%%%%%%%%%
Sessions=[17:19,47:50];
Res=zeros(7,3);
figure();hold on;grid on
c1={'b','--b','-b*','-.b',':b','--bo','--bs'};
r1={'r','--r','-r*','-.r',':r','--ro','--rs'};
for i= 1:length(Sessions)
    s=Sessions(i);
    load(strcat(savePath,'Session',num2str(s),'.mat'))
    Res(i,1)=HF(1);
    Res(i,2)=HF(2);
    Res(i,3)=Disengage;
    L=length(HFTraces);
    HFTraces=HFTraces(:,1:3000:L);
    plot(1:5:(5*length(HFTraces)),HFTraces(1,:),c1{i})
    plot(1:5:(5*length(HFTraces)),HFTraces(2,:),r1{i})
end








