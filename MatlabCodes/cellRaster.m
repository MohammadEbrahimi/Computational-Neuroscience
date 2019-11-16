Day=22;
area=1;
ActiveTrialNumber=20;
SpeedTh=2;
CleanPeriod=19;

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
LoadPath{20}= 'C:\Users\Sadegh\Documents\VLMReborn\L364\Data\allDays';
LoadPath{21}= 'C:\Users\Sadegh\Documents\VLMReborn\L365\Data\allDays';
LoadPath{22}= 'C:\Users\Sadegh\Documents\VLMReborn\L367\Data\allDays';
LoadPath{23}= 'C:\Users\Sadegh\Documents\VLMReborn\L368\Data\allDays';



path=LoadPath{Day};
load(strcat(path,'\All_Sessions.mat'));
load(strcat(path,'\Calcium.mat'));
load(strcat(path,'\areas.mat'));
load(strcat(path,'\cellImage.mat'));




% ts=1;
% te=0;
% for i=1:0
%    ts=ts+ SessionLength{i};
% end
% for i=1:5
%     te=te+SessionLength{i};
% end
% td=ts:te;


X=cellData;%padding(eventBin,3,0);;
%X=[Speed';Speed'];


SE0=CRSE(2:end,:);
SE1=HitSE(2:end,:);
SE2=FASE(2:end,:);
SE3=MissSE(2:end,:);
group0=CR+FA;
group1=Hit+Miss;

Resp0=zeros(cellCount,65);
Resp1=zeros(cellCount,65);
Resp2=zeros(cellCount,65);
Resp3=zeros(cellCount,65);
sem0=zeros(cellCount,65);
sem1=zeros(cellCount,65);
sem2=zeros(cellCount,65);
sem3=zeros(cellCount,65);
qResp0=zeros(cellCount,65);
qResp1=zeros(cellCount,65);
qResp2=zeros(cellCount,65);
qResp3=zeros(cellCount,65);
c0=zeros(cellCount,65);
c1=zeros(cellCount,65);
c2=zeros(cellCount,65);
c3=zeros(cellCount,65);




ActiveAnimal=ones(length(Lick),1);
for i=1:(length(Lick)-75*ActiveTrialNumber)
    if max(Lick(i:(i+75*ActiveTrialNumber)))==0
        ActiveAnimal(i:(i+75*ActiveTrialNumber))=0;
    end
end
% Speed(1:ts)=100;
% Speed(te:end)=100;



%%%%%%%%%%%%%%%%%%%%%
%%%%Find cell Coordinates and assign area
cellIJ=zeros(cellCount,2);
CortexArea=zeros(cellCount,1);

for i=1:cellCount
    [row,ind]=max(cellImage{i});
    [~,indj]=max(row);
    indi=ind(indj);
    cellIJ(i,:)=TileTopLeft(i,:)+[indi-1,indj-1];
    for a=1:8
        
        if Area(cellIJ(i,1),cellIJ(i,2),a)
            CortexArea(i)=a;
            
        end
        
        
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lick(Lick>0)=0;

for k=1:size(SE0,1)
     
        if(max(Speed(SE0(k,1):SE0(k,1)+CleanPeriod))<=SpeedTh  && max(ActiveAnimal(SE0(k,1):SE0(k,2)))==1 && max(Lick(SE0(k,1):SE0(k,1)+CleanPeriod))==0)
            j=0;
            while group0(SE0(k,1)+j)==0 && ((SE0(k,1)+j) < SE0(k,2))
                j=j+1;
            end
            stimL=0;
            for i=0:20
                if group0(SE0(k,1)+j+i)==0 && ((SE0(k,1)+j+i) <= SE0(k,2)+2)
                    stimL=i;
                    break;
                end
            end

            if (stimL>0)
                Resp0(:,1:10+stimL)=Resp0(:,1:10+stimL)+X(:,SE0(k,1)+j-10:SE0(k,1)+j+stimL-1);
                qResp0(:,1:10+stimL)=qResp0(:,1:10+stimL)+(X(:,SE0(k,1)+j-10:SE0(k,1)+j+stimL-1).^2);
                c0(:,1:10+stimL)=c0(:,1:10+stimL)+1;
            end
            
             j=0;
            while RewardWindow(SE0(k,1)+j)==0 && ((SE0(k,1)+j) < SE0(k,2))
                j=j+1;
            end
            rewL=0;
            for i=0:30
                if RewardWindow(SE0(k,1)+j+i)==0 && ((SE0(k,1)+j+i) <= SE0(k,2)+2)
                    rewL=i;
                    break;
                end
            end
            
            if rewL>0
                Resp0(:,31:35+rewL)=Resp0(:,31:35+rewL)+X(:,SE0(k,1)+j-5:SE0(k,1)+j+rewL-1);
                qResp0(:,31:35+rewL)=qResp0(:,31:35+rewL)+(X(:,SE0(k,1)+j-5:SE0(k,1)+j+rewL-1).^2);
                c0(:,31:35+rewL)=c0(:,31:35+rewL)+1;
            end
            
            
            
       end
end
Raster0 = Resp0 ./ c0;
sem0=((qResp0./(c0-1))- Raster0.^2)./sqrt(c0);



for k=1:size(SE1,1)
     
        if(max(Speed(SE1(k,1):SE1(k,1)+CleanPeriod))<=SpeedTh && max(ActiveAnimal(SE1(k,1):SE1(k,2)))==1 && max(Lick(SE1(k,1):SE1(k,1)+CleanPeriod))==0)
            j=0;
            while group1(SE1(k,1)+j)==0 && ((SE1(k,1)+j) < SE1(k,2))
                j=j+1;
            end
            stimL=0;
            for i=0:20
                if group1(SE1(k,1)+j+i)==0 && ((SE1(k,1)+j+i) <= SE1(k,2)+2)
                    stimL=i;
                    break;
                end
            end

            if (stimL>0)
                Resp1(:,1:10+stimL)=Resp1(:,1:10+stimL)+X(:,SE1(k,1)+j-10:SE1(k,1)+j+stimL-1);
                qResp1(:,1:10+stimL)=qResp1(:,1:10+stimL)+(X(:,SE1(k,1)+j-10:SE1(k,1)+j+stimL-1).^2);
                c1(:,1:10+stimL)=c1(:,1:10+stimL)+1;
            end
            
             j=0;
            while RewardWindow(SE1(k,1)+j)==0 && ((SE1(k,1)+j) < SE1(k,2))
                j=j+1;
            end
            rewL=0;
            for i=0:30
                if RewardWindow(SE1(k,1)+j+i)==0 && ((SE1(k,1)+j+i) <= SE1(k,2)+2)
                    rewL=i;
                    break;
                end
            end
            
            if rewL>0
                Resp1(:,31:35+rewL)=Resp1(:,31:35+rewL)+X(:,SE1(k,1)+j-5:SE1(k,1)+j+rewL-1);
                qResp1(:,31:35+rewL)=qResp1(:,31:35+rewL)+(X(:,SE1(k,1)+j-5:SE1(k,1)+j+rewL-1).^2);
                c1(:,31:35+rewL)=c1(:,31:35+rewL)+1;
            end
            
            
            
       end
end
Raster1 = Resp1 ./ c1;
sem1=((qResp1./(c1-1))- Raster1.^2)./sqrt(c1);



for k=1:size(SE2,1)
     
        if(max(Speed(SE2(k,1):SE2(k,1)+CleanPeriod))<=SpeedTh && max(ActiveAnimal(SE2(k,1):SE2(k,2)))==1 && max(Lick(SE2(k,1):SE2(k,1)+CleanPeriod))==0)
            j=0;
            while group0(SE2(k,1)+j)==0 && ((SE2(k,1)+j) < SE2(k,2))
                j=j+1;
            end
            stimL=0;
            for i=0:20
                if group0(SE2(k,1)+j+i)==0 && ((SE2(k,1)+j+i) <= SE2(k,2)+2)
                    stimL=i;
                    break;
                end
            end

            if (stimL>0)
                Resp2(:,1:10+stimL)=Resp2(:,1:10+stimL)+X(:,SE2(k,1)+j-10:SE2(k,1)+j+stimL-1);
                qResp2(:,1:10+stimL)=qResp2(:,1:10+stimL)+(X(:,SE2(k,1)+j-10:SE2(k,1)+j+stimL-1).^2);
                c2(:,1:10+stimL)=c2(:,1:10+stimL)+1;
            end
            
             j=0;
            while RewardWindow(SE2(k,1)+j)==0 && ((SE2(k,1)+j) < SE2(k,2))
                j=j+1;
            end
            rewL=0;
            for i=0:30
                if RewardWindow(SE2(k,1)+j+i)==0 && ((SE2(k,1)+j+i) <= SE2(k,2)+2)
                    rewL=i;
                    break;
                end
            end
            
            if rewL>0
                Resp2(:,31:35+rewL)=Resp2(:,31:35+rewL)+X(:,SE2(k,1)+j-5:SE2(k,1)+j+rewL-1);
                qResp2(:,31:35+rewL)=qResp2(:,31:35+rewL)+(X(:,SE2(k,1)+j-5:SE2(k,1)+j+rewL-1).^2);
                c2(:,31:35+rewL)=c2(:,31:35+rewL)+1;
            end
            
            
            
       end
end
Raster2 = Resp2 ./ c2;
sem2=((qResp2./(c2-1))- Raster2.^2)./sqrt(c2);



for k=1:size(SE3,1)
     
        if(max(Speed(SE3(k,1):SE3(k,1)+CleanPeriod))<=SpeedTh  && max(ActiveAnimal(SE3(k,1):SE3(k,2)))==1 && max(Lick(SE3(k,1):SE3(k,1)+CleanPeriod))==0)
            j=0;
            while group1(SE3(k,1)+j)==0 && ((SE3(k,1)+j) < SE3(k,2))
                j=j+1;
            end
            stimL=0;
            for i=0:20
                if group1(SE3(k,1)+j+i)==0 && ((SE3(k,1)+j+i) <= SE3(k,2)+2)
                    stimL=i;
                    break;
                end
            end

            if (stimL>0)
                Resp3(:,1:10+stimL)=Resp3(:,1:10+stimL)+X(:,SE3(k,1)+j-10:SE3(k,1)+j+stimL-1);
                qResp3(:,1:10+stimL)=qResp3(:,1:10+stimL)+(X(:,SE3(k,1)+j-10:SE3(k,1)+j+stimL-1).^2);
                c3(:,1:10+stimL)=c3(:,1:10+stimL)+1;
            end
            
             j=0;
            while RewardWindow(SE3(k,1)+j)==0 && ((SE3(k,1)+j) < SE3(k,2))
                j=j+1;
            end
            rewL=0;
            for i=0:30
                if RewardWindow(SE3(k,1)+j+i)==0 && ((SE3(k,1)+j+i) <= SE3(k,2)+2)
                    rewL=i;
                    break;
                end
            end
            
            if rewL>0
                Resp3(:,31:35+rewL)=Resp3(:,31:35+rewL)+X(:,SE3(k,1)+j-5:SE3(k,1)+j+rewL-1);
                qResp3(:,31:35+rewL)=qResp3(:,31:35+rewL)+(X(:,SE3(k,1)+j-5:SE3(k,1)+j+rewL-1).^2);
                c3(:,31:35+rewL)=c3(:,31:35+rewL)+1;
            end
            
            
            
       end
end
Raster3 = Resp3 ./ c3;
sem3=((qResp3./(c3-1))- Raster3.^2)./sqrt(c3);


cellNum=2;
index=find(CortexArea==4)

for cellNum=index'

figure();
td=-0.9:0.1:2;
hold on;grid on
shadedErrorBar(td,Raster1(cellNum,1:length(td)),sem1(cellNum,1:length(td)),'b',1);
shadedErrorBar(td,Raster3(cellNum,1:length(td)),sem3(cellNum,1:length(td)),'m',1);

% figure();
% hold on;grid on
shadedErrorBar(td,Raster0(cellNum,1:length(td)),sem0(cellNum,1:length(td)),'k',1); 
shadedErrorBar(td,Raster2(cellNum,1:length(td)),sem2(cellNum,1:length(td)),'r',1);
hold off
end

% 
% 
% C0=CountH;
% SS0=SEMH .* C0;
% M0=MeanH;
% C1=CountM;
% SS1=SEMM .* C1;
% M1=MeanM;
% F=zeros(1,length(td));
% DegF=zeros(1,length(td));
% for i=1:length(td)
%     SSW=(SS0(i)+SS1(i)) / (C1(i)+C0(i)-2);
%     MT=((C0(i)*M0(i))+(C1(i)*M1(i)))/(C0(i)+C1(i));
%     SSB=(C0(i)*((M0(i)-MT)^2))  +  (C1(i)*((M1(i)-MT)^2));
%     
%     F(i)=SSB/SSW;
%     DegF(i)=(C1(i)+C0(i)-2);    
% end
% 
% figure();plot(td,F)
% figure();plot(td,DegF)


%%%%%%%%%%%%%%%%%%%%Show cell Map

% cellMap=zeros(1017,1017,3);
% RGBarea=[0 0 3 ; 0 1.5 1.5 ; 0 1.5 1.5 ; 3 0 0 ; 0 3 0 ; 2 0 2 ;2 2 2; 3 2 0 ];
% for area =5:5
%     Mask=zeros(1,cellCount);
%     Mask(oIndex([sc]))=1;
%     %Mask(BMAT(:,area)~=0)=1;
%     
% for i=1:cellCount
% 
%     if Mask(i) ~=0
%         clear buf
%         ind=TileTopLeft(i,:);
%         buf=cellImage{i};
%         
%         s=size(buf);
%         cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,1)=...
%             cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,1)+(RGBarea(area,1)*buf);
%         
%         cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,2)=...
%             cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,2)+(RGBarea(area,2)*buf);
%         
%         cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,3)=...
%             cellMap(ind(1):ind(1)+s(1)-1,ind(2):ind(2)+s(2)-1,3)+(RGBarea(area,3)*buf);
%         
%     end
%     
% end
% end
% 
% cellMap(:,:,1)=cellMap(:,:,1)+ 0.001*Area(:,:,farea);
% cellMap(:,:,1)=cellMap(:,:,1)-min(min(cellMap(:,:,1)));
% cellMap(:,:,1)=floor(250*cellMap(:,:,1)/max(max(cellMap(:,:,1))));
% 
% cellMap(:,:,2)=cellMap(:,:,2)-min(min(cellMap(:,:,2)));
% cellMap(:,:,2)=floor(250*cellMap(:,:,2)/max(max(cellMap(:,:,2))));
% 
% cellMap(:,:,3)=cellMap(:,:,3)-min(min(cellMap(:,:,3)));
% cellMap(:,:,3)=floor(250*cellMap(:,:,3)/max(max(cellMap(:,:,3))));
% 
% 
% 
% figure()
% imshow(uint8(cellMap));
