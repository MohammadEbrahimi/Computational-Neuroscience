SimPath='C:\Users\Sadegh\Documents\VLMReborn\Reports\2016_27_10_Fisher_PLS_TI_TB_ManualMap\DuringStim\';
DelPath='C:\Users\Sadegh\Documents\VLMReborn\Reports\2016_27_10_Fisher_PLS_TI_TB_ManualMap\DuringDelay\';
RewPath='C:\Users\Sadegh\Documents\VLMReborn\Reports\2016_27_10_Fisher_PLS_TI_TB_ManualMap\DuringReward\';

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

FileNames{1}='L347_8_4_2015_active';
FileNames{2}='L347_8_5_2015_active';
FileNames{3}='L347_8_6_2015_active';
FileNames{4}='L347_8_7_2015_active';
FileNames{5}='L347_8_8_2015_active';

FileNames{6}='L354_8_3_2015_active';
FileNames{7}='L354_8_4_2015_active';
FileNames{8}='L354_8_5_2015_active';
FileNames{9}='L354_8_6_2015_active';
FileNames{10}='L354_8_7_2015_active';
FileNames{11}='L354_8_8_2015_active';%%Short%%

FileNames{12}='L362_8_3_2015_active';
FileNames{13}='L362_8_5_2015_active';
FileNames{14}='L362_8_6_2015_active';
FileNames{15}='L362_8_7_2015_active';
FileNames{16}='L362_8_8_2015_active';

FileNames{17}='L347_all_active';
FileNames{18}='L354_all_active';
FileNames{19}='L362_all_active';


Days=[17];
time=-0.4:0.1:5.5;
LowDimVec=1:100;
%Days=[1 2 3 4 5 10 11 12 13 15 16]; %% error days

ND=length(Days);
NT=length(time);
NL=length(LowDimVec);


IFPLSC=zeros(8,NL,ND);
MaxIFPLSC=zeros(ND,8,3);
PLSDim=zeros(ND,8,3);
TimeData=zeros(ND,8,NT);

pEgE=zeros(8,8,ND);
pEgC=zeros(8,8,ND);
pCgC=zeros(8,8,ND);
pCgE=zeros(8,8,ND);
pC=zeros(8,ND);
pE=zeros(8,ND);
Pt=zeros(8,8,4);

pMgA=zeros(8,4,ND);
pMgAC=zeros(8,2,ND);


color={'b','b--','k--','r','r--','m','k','g'};
color2={'bo','ro','ko'};
TitleV={'V1','LV','MV','PTLP','A','S','M','RSC'};
set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',16);
set(0,'defaultAxesFontName','Calibri');




for j=1:length(Days)
    i=Days(j);
    %V1=open(strcat(SimPath,'GoNogo_alldata\',FileNames{i},'_tr.mat'));

    
    maxDim=10;

    
    for area=1:8
        
        [MaxIFPLSC(j,area,1),PLSDim(j,area,1)]=max(V1.IFisherPLSC_ValNT(area,1:maxDim));

    end
    

    %[X,CortexArea,HitC,MissC,CRC,FAC,HitSE,MissSE,CRSE,FASE]=RetriveData(LoadPath{Days(j)},1);
    
    SE=[CRSE];
    
    group=(CRC)';
    class=0;

    L=length(group);
    ScoreMat=zeros(8,L);

    
    for a=1:8
        ScoreMat(a,:)=LRClassify(X' , V1.BMAT(:,a,PLSDim(j,a,1)) ,V1.InterceptMAT(a,1,PLSDim(j,a,1)));
        ScoreMat(a,group==0)=0.5;
    end
    
    
    for a1=1:8
        for a2=1:8
            if a1~=a2
               
            
            C1=(ScoreMat(a1,:)==class & group==1);% | (ScoreMat(a1,:)==1 & group1==1);
            E1=(ScoreMat(a1,:)==(1-class) & group==1);% | (ScoreMat(a1,:)==0 & group1==1);
            C2=(ScoreMat(a2,:)==class & group==1);% | (ScoreMat(a2,:)==1 & group1==1);
            E2=(ScoreMat(a2,:)==(1-class) & group==1);% | (ScoreMat(a2,:)==0 & group1==1);
            
            pEgE(a1,a2,j)=sum(E1 & E2) / sum(E2);
            pEgC(a1,a2,j)=sum(E1 & C2) / sum(C2);
            pCgE(a1,a2,j)=sum(C1 & E2) / sum(E2);
            pCgC(a1,a2,j)=sum(C1 & C2) / sum(C2);
            
            pE(a1,j)=sum(E1==1)/sum(group==1);
            pC(a1,j)=sum(C1==1)/sum(group==1);
            
            Pt(a1,a2,1)= (pEgE(a1,a2,j)*pC(a1,j))/((pEgE(a1,a2,j)*pC(a1,j))+...
                pE(a1,j)-(pEgE(a1,a2,j)*pE(a1,j)));
            
            Pt(a1,a2,4)= (pCgC(a1,a2,j)*pE(a1,j))/((pCgC(a1,a2,j)*pE(a1,j))+...
                pC(a1,j)-(pCgC(a1,a2,j)*pC(a1,j)));
            
            end
        end
        
    %         pMgA(a1,1,j)=sum(FAC'==1 & ScoreMat(a1,:)==1)/...
    %            ( sum(ScoreMat(a1,:)==1 & FAC'==1)+sum(ScoreMat(a1,:)==1 & CRC'==1));
    %         pMgA(a1,2,j)=sum(MissC'==1 & ScoreMat(a1,:)==0)/...
    %             (sum(ScoreMat(a1,:)==0 & MissC'==1)+sum(ScoreMat(a1,:)==0 & HitC'==1));
    %         pMgA(a1,3,j)=sum(FAC'==1)/sum(FAC+CRC==1);
    %         pMgA(a1,4,j)=sum(MissC'==1)/sum(MissC+HitC==1);
    %         
    %         pMgAC(a1,1,j)=sum(CRC'==1 & ScoreMat(a1,:)==1)/...
    %            ( sum(ScoreMat(a1,:)==1 & FAC'==1)+sum(ScoreMat(a1,:)==1 & CRC'==1));
    %         pMgAC(a1,2,j)=sum(HitC'==1 & ScoreMat(a1,:)==0)/...
    %             (sum(ScoreMat(a1,:)==0 & MissC'==1)+sum(ScoreMat(a1,:)==0 & HitC'==1));
        
    end
    
    
    
end





Pt(1,1,:)=0.5;Pt(3,3,:)=0.5;Pt(4,4,:)=0.5;Pt(6,6,:)=0.5;Pt(8,8,:)=0.5;


f=[1 3 4 6 8];
figure();
imagesc([Pt(f,f,1),Pt(f,f,4)]);








