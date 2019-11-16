hinfo=hdf5info('A:\Mapping\L347_Preproccessed\L347Mapping90_orig-Objects\Obj_1 - L347Mapping90deg(1).h5');
movie=single(hdf5read(hinfo.GroupHierarchy.Datasets(1)));
movie(isnan(movie))=0;

% global StorageArray
% Stim=StorageArray.Children{1, 1}.Object.Children{1, 5}.Object.Data;
% Time=StorageArray.Children{1, 1}.Object.Children{1, 5}.Object.XVector;
% Frames=StorageArray.Children{1, 1}.Object.Children{1, 6}.Object.Data;
% 
% L=length(Stim);
% t1=[diff(Stim)<-30 & Stim(1:L-1) <=0;0];
% t2=[diff(Stim)<-30 & Stim(1:L-1) >0;0];
% StartT=Frames(t1==1);
% EndT=Frames(t2==1);

% AvgMovie=single(zeros(size(movie,1),size(movie,2),216));
% for i=1:min([length(StartT),length(EndT)]);
%     AvgMovie=AvgMovie+movie(:,:,StartT(i)-25:StartT(i)+190);
% end
% AvgMovie=AvgMovie/min([length(StartT),length(EndT)]);
% 
% ShowMe(AvgMovie,1,20);



 dim=size(movie);
p=157;
nt=0;

avgMovie=zeros(dim(1),dim(2),p);
for i=1:p:(dim(3)-p+1)
    avgMovie=avgMovie+movie(:,:,i:i+p-1);
    nt=nt+1;
end


avgdff=uint8(zeros(dim(1),dim(2),p));
for i=1:dim(1)
    for j=1:dim(2)
        x=squeeze(avgMovie(i,j,:));
        x=(x-mean(x))/mean(x);
        %x=x/max(x);
        
        avgdff(i,j,:)=uint8(255*x);
    end
end


% p=size(AvgMovie,3);
% dim=size(AvgMovie);
% rgb=zeros(dim(1),dim(1),3,p);
% for i=1:p
%     rgb(:,:,:,i)=ind2rgb(uint8(AvgMovie(:,:,i)),jet(255));
% end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LoadPath{1}= 'C:\Users\Sadegh\Documents\VLMReborn\L347\Data\8_4_2015';
% LoadPath{2}= 'C:\Users\Sadegh\Documents\VLMReborn\L347\Data\8_5_2015';
% LoadPath{3}= 'C:\Users\Sadegh\Documents\VLMReborn\L347\Data\8_6_2015';
% LoadPath{4}= 'C:\Users\Sadegh\Documents\VLMReborn\L347\Data\8_7_2015';
% LoadPath{5}= 'C:\Users\Sadegh\Documents\VLMReborn\L347\Data\8_8_2015';
% 
% LoadPath{6}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_3_2015';
% LoadPath{7}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_4_2015';
% LoadPath{8}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_5_2015';
% LoadPath{9}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_6_2015';
% LoadPath{10}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_7_2015';
% LoadPath{11}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\8_8_2015';
% 
% LoadPath{12}= 'C:\Users\Sadegh\Documents\VLMReborn\L362bis\Data\8_3_2015';
% LoadPath{13}= 'C:\Users\Sadegh\Documents\VLMReborn\L362bis\Data\8_5_2015';
% LoadPath{14}= 'C:\Users\Sadegh\Documents\VLMReborn\L362bis\Data\8_6_2015';
% LoadPath{15}= 'C:\Users\Sadegh\Documents\VLMReborn\L362bis\Data\8_7_2015';
% LoadPath{16}= 'C:\Users\Sadegh\Documents\VLMReborn\L362bis\Data\8_8_2015';
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MapPath='C:\Users\Sadegh\Documents\VLMReborn\L347\Mapping\';
% Days=1:5;
% ND=length(Days);

% MapDay=open(strcat(MapPath,'\FirstFrame.mat'));
% moving=imrotate(MapDay.Mapping0,90);
% 
% 
% for k=1:ND
%     ExpDay=open(strcat(LoadPath{Days(k)},'\FirstFrame.mat'));
%     fixed=ExpDay.FirstFrame;
%     
%     [optimizer, metric] = imregconfig('multimodal');
%     
%     optimizer.InitialRadius = 0.001;
%     optimizer.Epsilon = 1.5e-6;
%     optimizer.GrowthFactor = 1.01;
%     optimizer.MaximumIterations = 1000;
%     
%     tform = imregtform(moving, fixed, 'similarity', optimizer, metric);
%     
% %      movingRegistered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)));
% % 
% %      figure
% %      imshowpair(fixed, movingRegistered,'Scaling','joint');
% %      
%     
%     
%     for i=1:8
%         [~,~,mask]=imread(strcat(MapPath,'ManualMapMasks\a',int2str(i),'.png'));
%         mask=imrotate((mask ~=0),90);
%         maskRegistered(:,:,i) = imwarp(mask,tform,'OutputView',imref2d([1017,1017]));
%         imwrite(maskRegistered(:,:,i),strcat(LoadPath{Days(k)},'\areaMasks\','Manual_a',int2str(i),'.png'));
%     end
%     
% end
%%%%% reload after correction and save as matlab array
% Area=zeros(1017,1017,8);
% for k=1:ND
% for i=1:8
%     x=imread(strcat(LoadPath{Days(k)},'\areaMasks\Manual_a',int2str(i),'.png'));
%     Area(:,:,i)=x;%(4:1020,4:1020);
% end
% save(strcat(LoadPath{Days(k)},'\areas'),'Area');
% 
% 
% end
% 



