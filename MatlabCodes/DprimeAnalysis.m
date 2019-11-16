

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


Day=17;
ML=1200;
MXL=110000;



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
clear DPMAT_time

clear DPMAT_TRUESHUFF
clear DPMAT_SHUFF
clear LMAT
clear LMAT_SHUFF


for area=1:8
    area
    X=cellEvents(CortexArea==area,:);
    dim=size(X);
   
    
    L=ML:2000:MXL;
    DPbuf=zeros(dim(1),length(L));
    DPbuf_t=zeros(dim(1),length(L));
    
    for s=1:length(L)
        k=L(s);
        %progress=round(100*k/MXL)
        dbuf=zeros(dim(1),dim(2)-k+1);
        for it=1:100:dim(2)-k+1
            trace=X(:,it:it+k-1)';
            m0=mean(trace(group0(it:it+k-1)==1,:));
            v0=var(trace(group0(it:it+k-1)==1,:));
            m1=mean(trace(group1(it:it+k-1)==1,:));
            v1=var(trace(group1(it:it+k-1)==1,:));
            dbuf(:,it)=(abs(m1-m0)./sqrt(0.5*(v1+v0)))';

        end
        
        [DPbuf(:,s),DPbuf_t(:,s)]=max(dbuf,[],2);
        
    end
    DPMAT_time{area}=DPbuf_t;
    DPMAT{area}=DPbuf;

    LMAT{area}=L;
end


for area=1:8
    area
    X=cellEvents_sh(CortexArea==area,:);
    dim=size(X);
   
    
    L=ML:2000:MXL;
    DPbuf=zeros(dim(1),length(L));

    for s=1:length(L)
        k=L(s);
        %progress=round(100*k/MXL)
        dbuf=zeros(dim(1),dim(2)-k+1);
        for it=1:100:dim(2)-k+1
            trace=X(:,it:it+k-1)';
            m0=mean(trace(group0(it:it+k-1)==1,:));
            v0=var(trace(group0(it:it+k-1)==1,:));
            m1=mean(trace(group1(it:it+k-1)==1,:));
            v1=var(trace(group1(it:it+k-1)==1,:));
            dbuf(:,it)=(abs(m1-m0)./sqrt(0.5*(v1+v0)))';

        end
        
        DPbuf(:,s)=max(dbuf,[],2);
        
    end
   
    DPMAT_TRUESHUFF{area}=DPbuf;
end


for shuffle=1:1
    shuffle
    X=cellEvents;
    dim=size(X);
    
    X_shuff=X(randperm(dim(1),1000),:);
    for shi=1:1000
        X_shuff(shi,:)=X_shuff(shi,randperm(dim(2)));
    end
    
    L=ML:2000:MXL;
    DPbufsh=zeros(1000,length(L));
    
    for s=1:length(L)
        k=L(s);
        progress=round(100*k/MXL)
        dbufsh=zeros(1000,dim(2)-k+1);
        for it=1:100:dim(2)-k+1

            trace=X_shuff(:,it:it+k-1)';
            m0=mean(trace(group0(it:it+k-1)==1,:));
            v0=var(trace(group0(it:it+k-1)==1,:));
            m1=mean(trace(group1(it:it+k-1)==1,:));
            v1=var(trace(group1(it:it+k-1)==1,:));
            dbufsh(:,it)=(abs(m1-m0)./sqrt(0.5*(v1+v0)))';
            
        
        end
        
        DPbufsh(:,s)=max(dbufsh,[],2);
    end
    DPMAT_SHUFF{shuffle}=DPbufsh;
    LMAT_SHUFF{shuffle}=L;
end


%test=DPMAT{1};

DPMEAN_S=zeros(8,length(LMAT{1}));
y=max(DPMAT_SHUFF{1});
for area =1:8
    for td=1:length(LMAT{1})
        x=DPMAT{area};
        x(isinf(x))=0;
        x(isnan(x))=0;
        DPMEAN_S(area,td)=mean(x((x(:,td)>y(td)),td));
    end
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cellCount=size(BMAT,1);
% 
% cellIJ=zeros(cellCount,2);
% CortexArea=zeros(cellCount,1);
% 
% for i=1:cellCount
%     [row,ind]=max(cellImage{i});
%     [~,indj]=max(row);
%     indi=ind(indj);
%     cellIJ(i,:)=TileTopLeft(i,:)+[indi-1,indj-1];
%     for a=1:8
%         
%         if Area(cellIJ(i,1),cellIJ(i,2),a)
%             CortexArea(i)=a;
%             
%         end
%         
%         
%     end
%     
% end
% 
% 
% color={'b','b--','k--','r','r--','m','k','g'};
% TitleV={'V1','LV','MV','PTLP','A','S','M','RSC'};
% set(0,'defaultlinelinewidth',2);
% set(0,'DefaultAxesFontSize',16);
% set(0,'defaultAxesFontName','Calibri');
% 
% figure();hold on
% for a=1:8
%     fx=BMAT(CortexArea==a,a);
%     x=DPMAT{a};
%             x(isinf(x))=0;
%         x(isnan(x))=0;
%     plot(LMAT{a}/600,mean(x(fx~=0,:)),color{a})
% end
% plot(LMAT_SHUFF{1}/600,max(DPMAT_SHUFF{1}),'y')
%                 
    
    
    
    
