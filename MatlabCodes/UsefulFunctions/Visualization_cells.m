AddressSetup;

Days=[1,6,12,21,26,33,40];
AllDays=[17,18,19,47,48,49,50]



cellMap=zeros(1017,1017,3);
RGBarea=[2 0 0 ; -2 0 2 ; 0 1.5 1.5 ; 1.5 0 1.5; 0 3 0 ; 2 0 2 ;2 2 2; 3 2 0 ];
cellCount=zeros(7,1);
area=1;
for j=7:7%1:length(Days)
 


%A=open(strcat(LoadPath{i},'\All_Sessions.mat'));
% M=open(strcat(LoadPath{i},'\areas.mat'));
 I=open(strcat(LoadPath{Days(j)},'\cellImage.mat'));
 I2=open(strcat(LoadPath{AllDays(j)},'\cellImage.mat'));
% cellData=A.cellData;
 cellImage=I.cellImage;
 TileTopLeft=I.TileTopLeft;
cellCount=size(TileTopLeft,1);



for i=1:cellCount


        clear buf
        ind=TileTopLeft(i,:);
        buf=cellImage{i};
        Dim=size(buf);
        x=sort(reshape(buf,Dim(1)*Dim(2),1),'descend');
        th=x(20);
        buf(buf<th)=0;
        buf(buf>=th)=1;
        
        
        
        
        buf=buf(1:size(buf,1)-5,1:size(buf,2)-5);
        s=size(buf);
        cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,1)=...
            cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,1)+(RGBarea(area,1)*buf);
        
        cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,2)=...
            cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,2)+(RGBarea(area,2)*buf);
        
        cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,3)=...
            cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,3)+(RGBarea(area,3)*buf);
        

    
end



 cellImage=I2.cellImage;
 TileTopLeft=I2.TileTopLeft;
cellCount=size(TileTopLeft,1);


area=2
for i=1:cellCount


        clear buf
        ind=TileTopLeft(i,:);
        buf=cellImage{i};
        Dim=size(buf);
        x=sort(reshape(buf,Dim(1)*Dim(2),1),'descend');
        th=x(20);
        buf(buf<th)=0;
        buf(buf>=th)=1;
        
        
        
        
        buf=buf(1:size(buf,1)-5,1:size(buf,2)-5);
        s=size(buf);
        cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,1)=...
            cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,1)+(RGBarea(area,1)*buf);
        
        cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,2)=...
            cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,2)+(RGBarea(area,2)*buf);
        
        cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,3)=...
            cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,3)+(RGBarea(area,3)*buf);
        

    
end

end

SAT=900;
cellMap(:,:,1)=cellMap(:,:,1)-min(min(cellMap(:,:,1)));
cellMap(:,:,1)=floor(SAT*cellMap(:,:,1)/max(max(cellMap(:,:,1))));

cellMap(:,:,2)=cellMap(:,:,2)-min(min(cellMap(:,:,2)));
cellMap(:,:,2)=floor(SAT*cellMap(:,:,2)/max(max(cellMap(:,:,2))));

cellMap(:,:,3)=cellMap(:,:,3)-min(min(cellMap(:,:,3)));
cellMap(:,:,3)=floor(SAT*cellMap(:,:,3)/max(max(cellMap(:,:,3))));


 
 figure()
imshow(uint8(cellMap));



% 
% 
% dim=size(cellMapVidPLS);
% rgb=zeros(dim(1),dim(1),3,dim(3));
% for i=1:dim(3)
%     rgb(:,:,:,i)=ind2rgb(uint8(cellMapVidPLS(:,:,i)),jet(255));
% end
% rgb(1:30,987:end,1,5:25)=1;
% rgb(1:30,987:end,2,5:25)=0;
% rgb(1:30,987:end,3,5:25)=0;
% implay(rgb,2)
% 
% Video=VideoWriter('C:\Users\Sadegh\Documents\VLMReborn\Reports\2016_17_10_FisherTV_TB_ManualMap\DuringStim\HitCR_alldata\L362_PLS.avi');
% Video.FrameRate=2;
% open(Video);
% writeVideo(Video,rgb)
% close(Video);







