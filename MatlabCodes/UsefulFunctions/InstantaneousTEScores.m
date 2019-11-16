
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

SavePath='C:\Users\Sadegh\Documents\VLMReborn\Reports\2016_14_11_TransferEntropy\DuringStim\HitCR_inst\L362\L362_hitcr_';
Days=[19];
maxDim=10;
maxHist=5;
time=-4:20;
LowDimVec=1:100;

ND=length(Days);
NT=length(time);
NL=length(LowDimVec);



IFPLSC=zeros(8,NL,ND);
MaxIFPLSC=zeros(ND,8,3);
PLSDim=zeros(ND,8,3);
TimeData=zeros(ND,8,NT);


color={'b','b--','k--','r','r--','m','k','g'};
TitleV={'V1','LV','MV','PTLP','A','S','M','RSC'};
set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',16);
set(0,'defaultAxesFontName','Calibri');


for j=1:length(Days)
    i=Days(j);
    V1=open(strcat(SimPath,'HitCR_alldata\',FileNames{i},'_tr.mat'));
    
    
    
    
    
    for area=1:8
        
        [MaxIFPLSC(j,area,1),PLSDim(j,area,1)]=max(V1.IFisherPLSC_ValNT(area,1:maxDim));
        
    end
    
    
    
    
    [X,CortexArea,HitC,MissC,CRC,FAC,HitSE,MissSE,CRSE,FASE]=RetriveData(LoadPath{Days(j)},1);
    
    SE=[HitSE;CRSE];
    group=HitC+CRC;
    
    L=length(group);
    ScoreMat=zeros(8,L);
    
    shuff=randperm(size(SE,1));
    
    for a1=1:8
        ScoreMat(a1,:)= V1.BMAT(:,a1,PLSDim(j,a1,1))' * X;
    end
    
    %             if Days(j)==17 && a2==5
    %                 continue;
    %             end
    
    Label=1;
    for td=time
        c0=1;
        
        dd=zeros(8,L);
        for si=1:length(shuff)
            
            k=shuff(si);
            
            jj=0;
            while group(SE(k,1)+jj)==0 && ((SE(k,1)+jj) < SE(k,2))
                jj=jj+1;
            end
            
     
            stimL=0;
            for ii=0:max(time)
                if group(SE(k,1)+jj+ii)==0
                    stimL=ii;
           
                    break;
                end
            end
            if ( stimL>0 && td<=stimL)
                dd(:,c0:c0+maxHist)=ScoreMat(:,(SE(k,1)+jj+td-1-maxHist):(SE(k,1)+jj+td-1) );
                c0=c0+maxHist+1;
            end
        end
   
    c0=c0-1;
    data=dd(:,1:c0);
    save(strcat(SavePath,int2str(Label)),'data');
    Label=Label+1;
     end
    
end





















