function [choice]=sadeghSorter(inputImages,inputSignals,MaxProj)
    Dim=size(inputSignals);
    choice=(zeros(1,Dim(1))==ones(1,Dim(1)));
    topEvents=zeros(20,41);
    ImDim=size(MaxProj);
    BackGround=zeros([ImDim,3],'uint8');
    MaxProj(MaxProj>20)=20;
    
    
    BG=MaxProj-min(min(MaxProj));
    BG=BG/max(max(BG));
    BG=round(255*BG);
    BackGround(:,:,1)=uint8(BG);
    BackGround(:,:,2)=uint8(BG);
    BackGround(:,:,3)=uint8(BG);
    
    for cnum=1:Dim(1)
x=inputSignals(cnum,:);
[~,sortedIndex]=sort(x,'descend');
n=1;
eventList=zeros(20,1);
for k=1:Dim(2)
    if min(abs(eventList-sortedIndex(k)))>100
        eventList(n)=sortedIndex(k);
        topEvents(n,:)=inputSignals(cnum,eventList(n)-20:eventList(n)+20);
        n=n+1;
        if n==21
            break;
        end
    end
    
    
end

z=inputImages(:,:,cnum)-min(min(inputImages(:,:,cnum)));
z=z/max(max(z));
w=round(255*edge(z));
z=round(255*z);

maxJ=max(z);
[~,indJ]=max(maxJ);
[~,indI]=max(z(:,indJ));


cellImage=BackGround;
cellImage(:,:,1)=cellImage(:,:,1)+uint8(z);
buf=cellImage(:,:,1); 
buf(w>0)=0;
cellImage(:,:,1)=buf;
buf=cellImage(:,:,3); 
buf(w>0)=0;
cellImage(:,:,3)=buf;
buf=cellImage(:,:,2); 
buf(w>0)=0;
cellImage(:,:,2)=buf+uint8(w);
   


figHandle = figure('units','normalized','outerposition',[0.01 0.4 0.7 0.5]);



clf(figHandle);
returnMap = containers.Map;
set(figHandle, 'KeyPressFcn', ...
    @(fig_obj , eventDat) readInput(fig_obj, eventDat,returnMap));


subplot(1,3,1);
hold on 
set(0,'defaultlinelinewidth',0.5);
title(num2str(cnum));
plot(topEvents','b');
set(0,'defaultlinelinewidth',3);
plot(mean(topEvents),'r');
subplot(1,3,2);
imshow(cellImage(max(indI-30,1):min(indI+30,ImDim(1)),max(indJ-30,1):min(indJ+30,ImDim(2)),:)); 
subplot(1,3,3);
BufBG=BackGround;
BufBG(indI,indJ,3)=255;
BufBG(indI,indJ,1)=0;
BufBG(indI,indJ,2)=0;
imshow(BufBG(max(indI-30,1):min(indI+30,ImDim(1)),max(indJ-30,1):min(indJ+30,ImDim(2)),:)); 


waitfor(figHandle);

 newN= returnMap('newN')
 if newN==3 %Up
    choice(cnum)=(1==1);
 elseif newN==0 %Down
    choice(cnum)=(1==0);
 end


    end
end
