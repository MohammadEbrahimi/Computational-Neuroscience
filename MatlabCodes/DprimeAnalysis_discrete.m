

LoadPath{1}= 'C:\Users\Sadegh\Documents\VLMReborn\L347\Data\8_4_2015';
LoadPath{2}= 'C:\Users\Sadegh\Documents\VLMReborn\L347\Data\8_5_2015';
LoadPath{3}= 'C:\Users\Sadegh\Documents\VLMReborn\L347\Data\8_6_2015';
LoadPath{4}= 'C:\Users\Sadegh\Documents\VLMReborn\L347\Data\8_7_2015';
LoadPath{5}= 'C:\Users\Sadegh\Documents\VLMReborn\L347\Data\8_8_2015';

LoadPath{6}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_3_2015';
LoadPath{7}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_4_2015';
LoadPath{8}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_5_2015';
LoadPath{9}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_6_2015';
LoadPath{10}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_7_2015';
LoadPath{11}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_8_2015';

LoadPath{12}= 'C:\Users\Sadegh\Documents\VLMReborn\L362bis\Data\8_3_2015';
LoadPath{13}= 'C:\Users\Sadegh\Documents\VLMReborn\L362bis\Data\8_5_2015';
LoadPath{14}= 'C:\Users\Sadegh\Documents\VLMReborn\L362bis\Data\8_6_2015';
LoadPath{15}= 'C:\Users\Sadegh\Documents\VLMReborn\L362bis\Data\8_7_2015';
LoadPath{16}= 'C:\Users\Sadegh\Documents\VLMReborn\L362bis\Data\8_8_2015';

LoadPath{17}= 'C:\Users\Sadegh\Documents\VLMReborn\L347\Data\allDays';
LoadPath{18}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\allDays';
LoadPath{19}= 'C:\Users\Sadegh\Documents\VLMReborn\L362bis\Data\allDays';


Day=19;
sll=4;

% %SessionLength{1}=[5999 5999 5999 5998 5999 5998 5999 5999 5999 5999 5999 5999 5999 5999 5999 5999 5999 6000 5999 5999 5999 5999 5999 6000 5999 5999 5999 6000];
% SessionLength{1}=[11998 11997 11997 11998 11998 11998 11998 11998 11999 11998 11998 11999 11998 11999];
% %SessionLength{3}=[23995 23995 23996 23996 23997 23997 23997];
% SessionLength{2}=[35992 35994 35995 35995 23997];
% SessionLength{3}=[71986 71990];
% SessionLength{4}=[167973];

%%%%% SessionLength{1}=[5999 5999 5999 5999 5999 5999 5999 5999 5999 5999 5999 4987 5999 5999 5999 6000 5999 6000 5999 5999 5999 5999 5999 5999 5999 5999 5999 5999 5999 5999 5999 6000 5999 6000];
% SessionLength{1}=[11998 11998 11998 11998 11998 10986 11998 11999 11999 11998 11998 11998 11998 11998 11998 11999 11999]; 
% SessionLength{2}=[35994 34982 35996 35994 35994 23998];
% SessionLength{3}=[70976 71990 59992];
% SessionLength{4}=[202958];
% 
%%%%% SessionLength{1}=[5999 5999 5999 5999 5999 6000 5999 5999 5999 5999 5999 5999 5999 5999 5999 6000 5999 6000 5999 5999 5999 5999 5999 6000 5999 5999 5999 5999]; 
SessionLength{1}=[11998 11998 11999 11998 11998 11998 11998 11999 11999 11998 11998 11999 11998 11998];  
SessionLength{2}=[35995 35994 35996 35995 23996];
SessionLength{3}=[71989 71991];
SessionLength{4}=[167976];


path=LoadPath{Day};

mode=1;
[cellEvents,CortexArea,HitC,MissC,CRC,FAC,Delay,RewardWindow,Speed]=RetriveData(path,mode);


group0=CRC;
group1=HitC;


time=1:length(group0);

cellEvents_sh=cellEvents;

tvec=time(group0==1);
cellEvents_sh(:,tvec)=cellEvents(:,tvec(randperm(length(tvec))));
tvec=time(group1==1);
cellEvents_sh(:,tvec)=cellEvents(:,tvec(randperm(length(tvec))));


clear DPMAT

clear DPMAT_TRUESHUFF
clear DPMAT_SHUFF


for area=1:8
    area
    X=cellEvents(CortexArea==area,:);
    dim=size(X);
    
    for s=1:sll
        L=SessionLength{s};
        DPbuf=zeros(dim(1),length(L));
        start=1;
        for it=1:length(L)
            trace=X(:,start:sum(L(1:it)))';
            m0=mean(trace(group0(start:sum(L(1:it)))==1,:));
            v0=var(trace(group0(start:sum(L(1:it)))==1,:));
            m1=mean(trace(group1(start:sum(L(1:it)))==1,:));
            v1=var(trace(group1(start:sum(L(1:it)))==1,:));
            DPbuf(:,it)=(abs(m1-m0)./sqrt(0.5*(v1+v0)))';
            start=sum(L(1:it));
        end
        
        DPMAT{area,s}=DPbuf;
        
    end

end


for area=1:8
    area
    X=cellEvents_sh(CortexArea==area,:);
    dim=size(X);
    
    for s=1:sll
        L=SessionLength{s};
        DPbuf=zeros(dim(1),length(L));
        start=1;
        for it=1:length(L)
            trace=X(:,start:sum(L(1:it)))';
            m0=mean(trace(group0(start:sum(L(1:it)))==1,:));
            v0=var(trace(group0(start:sum(L(1:it)))==1,:));
            m1=mean(trace(group1(start:sum(L(1:it)))==1,:));
            v1=var(trace(group1(start:sum(L(1:it)))==1,:));
            DPbuf(:,it)=(abs(m1-m0)./sqrt(0.5*(v1+v0)))';
            start=sum(L(1:it));
        end
        
        DPMAT_TRUESHUFF{area,s}=DPbuf;
        
    end

end





 X=cellEvents;
    dim=size(X);
    
    X_shuff=X(randperm(dim(1),1000),:);
    for shi=1:1000
        X_shuff(shi,:)=X_shuff(shi,randperm(dim(2)));
    end
    
    for s=1:sll
        L=SessionLength{s};
        DPbuf=zeros(1000,length(L));
        start=1;
        for it=1:length(L)
            trace=X_shuff(:,start:sum(L(1:it)))';
            m0=mean(trace(group0(start:sum(L(1:it)))==1,:));
            v0=var(trace(group0(start:sum(L(1:it)))==1,:));
            m1=mean(trace(group1(start:sum(L(1:it)))==1,:));
            v1=var(trace(group1(start:sum(L(1:it)))==1,:));
            DPbuf(:,it)=(abs(m1-m0)./sqrt(0.5*(v1+v0)))';
            start=sum(L(1:it));
        end
        
        DPMAT_SHUFF{s}=DPbuf;
        
    end


    

area=1;
meandp=zeros(3,sll);

DPO=DPMAT; 
figure();hold on;    
for s=1:sll
    V=mean(DPO{area,s});
    V=V(~isnan(V));
    x=s*ones(1,length(V));
    plot(x,V,'b.');
    plot(s,mean(V),'r*');
    meandp(1,s)=mean(V);
end

DPO=DPMAT_TRUESHUFF;  
for s=1:sll
    V=mean(DPO{area,s});
    V=V(~isnan(V));
    meandp(2,s)=mean(V);
end


DPO=DPMAT_SHUFF;  
for s=1:sll
    V=mean(DPO{s});
    V=V(~isnan(V));
    meandp(3,s)=mean(V);
end


figure();hold on
plot(meandp(1,:),'b');
plot(meandp(2,:),'k');
plot(meandp(3,:),'r');
legend('normal','true shuff','shuff');






