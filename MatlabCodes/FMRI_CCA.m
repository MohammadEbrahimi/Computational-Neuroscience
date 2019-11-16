Human={'restingState_kgs_071512','restingState_knk_071812','restingState_kw_071512'};
LV1=load(strcat('E:\FMRI\',Human{3},'\LV1\Scan3.mat'));
RV1=load(strcat('E:\FMRI\',Human{3},'\RV1\Scan3.mat'));

trialTime=size(LV1.tSeries,1);
ML=repmat(mean(LV1.tSeries),[size(LV1.tSeries,1),1]);
MR=repmat(mean(RV1.tSeries),[size(RV1.tSeries,1),1]);

LV1.tSeries=(LV1.tSeries - ML);
RV1.tSeries=(RV1.tSeries - MR);
LV1.tSeries=(LV1.tSeries ./ ML);
RV1.tSeries=(RV1.tSeries ./ MR);

RVCells=(RV1.coords(2,:)<130);
LVCells=(LV1.coords(2,:)>145);
rVal=zeros(20,100);
rTrain=zeros(20,100);
rVal_sh=zeros(20,100);
rTrain_sh=zeros(20,100);

wxMat=zeros(sum(LVCells),20,100);
wyMat=zeros(sum(RVCells),20,100);

for sh=1:100

tind=sort(randperm(trialTime,trialTime/2),'ascend');
vind=[];
for i=1:trialTime
    if max(tind==i)==0
        vind=[vind,i];
    end
end
X=LV1.tSeries(tind,LVCells);
Y=RV1.tSeries(tind,RVCells);
[wxMat(:,:,sh),wyMat(:,:,sh),rTrain(:,sh)]=SparseCCA(X,Y,1,1,2,20);

X=LV1.tSeries(tind(randperm(length(tind))),LVCells);
Y=RV1.tSeries(tind(randperm(length(tind))),RVCells);
[wxMat_sh,wyMat_sh,rTrain_sh(:,sh)]=SparseCCA(X,Y,1,1,2,20);





X=LV1.tSeries(vind,LVCells);
Y=RV1.tSeries(vind,RVCells);
for cc=1:20
    if sign(mean(wxMat(:,cc,sh)))==sign(mean(wyMat(:,cc,sh)))
wxMat(:,cc,sh)=wxMat(:,cc,sh)*sign(mean(wxMat(:,cc,sh)));    
wyMat(:,cc,sh)=wyMat(:,cc,sh)*sign(mean(wyMat(:,cc,sh)));   
    else
        wxMat(:,cc,sh)=wxMat(:,cc,sh)*sign(mean(wxMat(:,cc,sh)));
        wyMat(:,cc,sh)=wyMat(:,cc,sh)*sign(mean(wyMat(:,cc,sh)))*-1;
    end
    
corrbuf=corrcoef([X*wxMat(:,cc,sh),Y*wyMat(:,cc,sh)]);
rVal(cc,sh)=corrbuf(1,2);
corrbuf=corrcoef([X(randperm(length(vind)),:)*wxMat_sh(:,cc),Y(randperm(length(vind)),:)*wyMat_sh(:,cc)]);
rVal_sh(cc,sh)=corrbuf(1,2);
end

end

figure();hold on
xlabel('Mode Number')
ylabel('Correlation Coefficient')
plot(mean(rTrain,2))
plot(max(rTrain_sh,[],2),'b--')
plot(mean(rVal,2),'r*-')
plot(max(rVal_sh,[],2),'r--')

legend('Training','Training Shuffled','Validation','Validation Shuffled (max of 100)'); 

figure();
for cc=1:1;
subplot(1,1,cc);hold on
title(strcat('Mode ',num2str(cc)));
ylabel('Projection Value');
xlabel('Time(s)');
LProj=LV1.tSeries(:,LVCells)*mean(wxMat(:,cc,:),3);
LProj=LProj/norm(LProj);
plot(1:2:288,LProj);
RProj=RV1.tSeries(:,RVCells)*mean(wyMat(:,cc,:),3)
RProj=RProj/norm(RProj);
plot(1:2:288,RProj,'r');

end
legend('Left V1 Projection','Right V1 Projection');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cc=1
% D=zeros(length([LVCells,RVCells]),1);
% D([LVCells,RVCells])=[wxMat(:,cc);wyMat(:,cc)];
% cord=[LV1.coords,RV1.coords];
% cord=cord-repmat(min(cord,[],2)-1,[1,size(cord,2)]);
% cord=round(cord);
% z=unique(cord(3,:));
% fmap=zeros(max(cord(1,:)),(max(cord(2,:))+1)*max(cord(3,:)));
% 
% for i=z
%     index=find(cord(3,:)==i);
%    for j=index
%        if D(j)~=0
%         fmap(cord(1,j),cord(2,j)+((i-1)*(max(cord(2,:))+1)))=D(j);
%        end
%    end
%    fmap([1:max(cord(1,:))],(max(cord(2,:))+1)*i)=min(D);
%    
% end
% imagesc(fmap);
% 
