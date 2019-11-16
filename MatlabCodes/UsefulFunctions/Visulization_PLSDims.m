ResultsPath='C:\Users\Sadegh\Documents\VLMReborn\Reports2\2017_05_19_Fisher_PLS\HitCR_alldata\HitCR_alldata';
mode=1;
Sessions=[17,18,19,47,48,49,50];

Days=[17];
time=-0.4:0.1:5.5;
LowDimVec=1:50;

ND=length(Sessions);
NT=length(time);
NL=length(LowDimVec);


clear TrialCorrInd
clear TrialCorr




color={'b','b--','k--','r','r--','m','k','g'};
color2={'bo','ro','ko'};
TitleV={'V1','LV','MV','PTLP','A','S','M','RSC'};
set(0,'defaultlinelinewidth',2);
set(0,'DefaultAxesFontSize',16);
set(0,'defaultAxesFontName','Calibri');

Info_Val_Fisher=zeros(8,NL,ND);
MaxInfoDim=zeros(8,ND);
MaxInfo=zeros(8,ND);

for j=1:length(Sessions)
    s=Sessions(j)
    V1=open(strcat(ResultsPath,'Session',num2str(s),'_mode',num2str(mode),'.mat'));

    
    
    
%    figure();
% 
%     maxDim=5;
%     hold on;
%     title(int2str(j))
%     grid on 
%     
    
for area=1:8
    Info_Val_Fisher(area,:,j)=V1.Info_Val_Fisher(area,:);
   plot(Info_Val_Fisher(area,:,j),color{area})
   [MaxInfo(area,j),MaxInfoDim(area,j)]=max(squeeze(Info_Val_Fisher(area,:,j)));
end


    
    
    [X,CortexArea,HitC,MissC,CRC,FAC,HitSE,MissSE,CRSE,FASE]=RetriveData(LoadPath{Sessions(j)},1);
    
    SE=[CRSE];
    
    group=FAC;
    L=sum(group);
    ScoreMat=zeros(8,L);
    T=15;
    rt=randperm(size(SE,1));
    
    for a1=1
        
        %  figure();hold on;grid on
        for a2=3
            a2
            s1= V1.BMAT(:,a1,PLSDim(j,a1,1))' * X;
            s2= V1.BMAT(:,a2,PLSDim(j,a2,1))' * X;
%             s1(group==0)=0;
%             s2(group==0)=0;
            
%             if Days(j)==17 && a2==5
%                 continue;
%             end
            
            
            ScoreMat(a2,:)=s2(group==1);
            
%             n=1;
%             for k=1:size(SE,1)
%               i=rt(k);
%               x2=s2(SE(i,1):SE(i,2));
%               x2=x2(group(SE(i,1):SE(i,2))==1);
%              ScoreMat(a2,n:n+length(x2)-1)=x2;
%              n=n+length(x2);
%              if n>4279
%                  break;
%              end
%             end 
            
            
     

        end
    end
    
    
    
    
end

% 
% %[~,Qdata]=Jdist(ScoreMat',[8 8 8 8 8 8 8 8]);
% Qdata=QuantEqual(ScoreMat,8);
% dim=size(Qdata);
% mdt=0;
% V=zeros(mdt+1,dim(2)-mdt);
% MV=zeros(mdt+1,dim(2)-mdt);
% Pt=zeros(mdt+1,dim(2)-mdt);
% S=zeros(mdt+1,dim(2)-mdt);
% RSC=zeros(mdt+1,dim(2)-mdt);
% for dt=1:mdt+1
%     V(dt,:)=Qdata(1,dt:(dim(2)+dt-1-mdt));
%     MV(dt,:)=Qdata(3,dt:(dim(2)+dt-1-mdt));
%     Pt(dt,:)=Qdata(4,dt:(dim(2)+dt-1-mdt));
%     S(dt,:)=Qdata(6,dt:(dim(2)+dt-1-mdt));
%     RSC(dt,:)=Qdata(8,dt:(dim(2)+dt-1-mdt));
% end
% s=(HitC*2)+CRC;
% s=s(s>0);
% s=s(mdt+1:dim(2));
% score = score_dags([V;MV;Pt;S;RSC], [8 8 8 8 8], dags,'discrete',[1 2 3 4 5], 'scoring_fn', 'bic');
% %ns=8*ones(1,5);%ns(1)=2;
% 
% %[sampled_graphs, accept_ratio] = learn_struct_mcmc([V;MV;Pt;S;RSC], ns, 'nsamples', 1000, 'burnin', 200);
% %mcmc_post = mcmc_sample_to_hist(sampled_graphs, dags);
% 
% %post = normalise(exp(score));
% labels={'V1','MV','PTLp','S','RSC'};
% %labels={'Stim','V5','V4','V3','V2','V1','V0','MV5','MV4','MV3','MV2','MV1','MV0',...
%  %   'P5','P4','P3','P2','P1','P0','S5','S4','S3','S2','S1','S0'};
% S=exp(score-max(score));
% S=S/sum(S);
% [b]=max(S);
% f=(S>(0.95*b));
% 
% figure();
% %  subplot(1,sum(f)+1,1)
%   plot(S)
% 
% %bar(mcmc_post)
% for i=1:sum(f)
% %subplot(sum(f),1,i)
% figure();
% title('m3 hitcr bic');
% [~,b]=max(S);
% draw_graph(dags{b},labels);
% S(b)=min(S);
% end
% 
% 
% 
% % figure();plot(score)
% % [a b]=max(score);
% % f=(score==a);
% % sum(f)
% % hold on;plot(f*a,'k*')
% % 







