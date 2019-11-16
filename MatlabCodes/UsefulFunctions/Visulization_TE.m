AddressSetup;
ResultsPath='C:\Users\Sadegh\Documents\VLMReborn\Reports2\2017_05_20_TE\HitCR_alldata\';
mode=1;

Sessions=[17,18,19,47,48,49,50];
ND=length(Sessions);
A=[1,3,4,6,8];
% TEHit=zeros(5,5,ND);
% TECR=zeros(5,5,ND);


for j=[1:7]
    s=Sessions(j);
   V1=open(strcat(ResultsPath,'Session',num2str(s),'CRScores_mode',num2str(mode),'_\results\linnue_meanReshapeMtx.mat'));
A=[1,3,4,6,8];
if (j==2 || j==5 || j==7)  A=1:5; end
TECR(:,:,j)=V1.meanRes(A,A);
% figure();
% 
% imagesc(V1.meanRes(A,A));
% title(strcat('M',num2str(j)));
end

% fval=zeros(5,5);
% for a1=1:5
%     for a2=1:5
%         fval(a1,a2)=repeatedMeasures([squeeze(TECR(a1,a2,:))';squeeze(TEHit(a1,a2,:))']);
%     end
% end
% 
