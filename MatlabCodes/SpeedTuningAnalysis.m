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


Days=[17:19];

ND=length(Days);
Bins=20;

SpeedTuneHist=zeros(8,Bins,ND,3);
SpeedTuneBins=zeros(8,Bins,ND,3);
MST=zeros(8,ND,3);
TunedRatio=zeros(8,ND,3);
DataLength=zeros(ND,3);
for mode=1
for j=1:ND
        [X,CortexArea,HitC,MissC,CRC,FAC,Delay,Reward,Speed]=RetriveData(LoadPath{Days(j)},1);
        if mode==1
        group=HitC+CRC+MissC+FAC;
        elseif mode ==2
            group=Delay;
        elseif mode==3
            group=Reward;
        end
       
        length(group)
        Speed(group==0)=0;
        Speed=Speed-min(Speed);
        Speed=Speed/max(Speed);
        q=quantile(Speed(Speed>0.001),20);
        [binCounts,S_index]=histc(Speed,q);
        group(Speed<0.001)=0;
        
        cellCorr=zeros(size(X,1),1);
        cellCorrP=zeros(size(X,1),1);
        for i=1:size(X,1)
             [cellCorr(i),cellCorrP(i)]=corr(X(i,group==1)',S_index(group==1),'Type','Spearman');
             DataLength(j,mode)=sum(group==1);
%             mdl=fitlm(X(i,group==1),Speed(group==1),'linear');
%             cellCorr(i)=mdl.Coefficients.Estimate(2);
        end
        
        
        
        
        for a=1:8
            SCVEC=cellCorr(CortexArea==a & cellCorrP<0.05 );
            if length(SCVEC)>0
            TunedRatio(a,j,mode)=sum(CortexArea==a & cellCorrP<0.05 )/sum(CortexArea==a );
            MST(a,j,mode)=mean(SCVEC);
            SpeedTuneBins(a,:,j,mode)=linspace(min(SCVEC),max(SCVEC),Bins);
            SpeedTuneHist(a,:,j,mode)=histc(SCVEC, SpeedTuneBins(a,:,j,mode));
            end
        end
end
end         
    



color={'b','y','r','r','r--','m','k','g'};        
mode=1;

for j=1:ND
    figure();hold on
    width=1;
    for area=[1,3]
        if max(SpeedTuneBins(area,:,j,mode)>0)
      bar(SpeedTuneBins(area,:,j,mode),SpeedTuneHist(area,:,j,mode)/sum(SpeedTuneHist(area,:,j,mode)),width,color{area}) ;
      width=width-0.4;
        end
    end
    legend('V1','PM,AM')
end
    


