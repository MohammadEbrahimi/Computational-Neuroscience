AddressSetup;

Days=1:5;
AllDays=17;
ND=length(Days);
DCent=5;
corrTh=0.5;

clear tform
clear cellImages
clear tcorr
%%%%%%%Load the first day
preDay=imread(strcat(LoadPath{Days(1)},'\areaMasks\a4.png'));
preIm=open(strcat(LoadPath{Days(1)},'\cellImage.mat'));
fixed=uint8(preDay);
cellImPre=preIm.cellImage;
TTLPre=preIm.TileTopLeft;
CentroidPre=zeros(size(TTLPre));
cellIndex=zeros(length(cellImPre),length(Days));
for i=1:length(cellImPre)
    x=cellImPre{i};
    cellImages{i,1}=x;
    [max1,ind1]= max(x,[],1);
    [~,loc2]= max(max1);
    [~,loc1]= max(x(:,loc2));
    CentroidPre(i,:)=[loc1,loc2];
    cellIndex(i,1)=i;
end
cellNumber=length(cellImPre)+1;




for pdind=2:ND
    pdind
    %postDay=imread(strcat(CaMovies{Days(pdind)-20},'S',num2str(Days(pdind)-20),'_Mask.png'));
    postDay=imread(strcat(LoadPath{Days(pdind)},'\areaMasks\a4.png'));
    postIm=open(strcat(LoadPath{Days(pdind)},'\cellImage.mat'));
    
    
    
    moving=uint8(postDay);
    cellImPost=postIm.cellImage;
    TTLPost=postIm.TileTopLeft;
    
    [optimizer, metric] = imregconfig('multimodal');
    
    optimizer.InitialRadius = 0.001;
    optimizer.Epsilon = 1.5e-6;
    optimizer.GrowthFactor = 1.01;
    optimizer.MaximumIterations = 1000;
    
    tform{pdind} = imregtform(moving, fixed, 'similarity', optimizer, metric);
    
    %     movingRegistered = imwarp(moving,tform{i},'OutputView',imref2d(size(fixed)));
    %
    %     figure
    %     imshowpair(fixed, movingRegistered,'Scaling','joint');
    
    clear cellImPostReg
    TTLPostReg=zeros(size(TTLPost));
    CentroidPostReg=zeros(size(TTLPost));
    TTLVal=[1,246,491,736];
    for ipost=1:length(cellImPost) 
     
        y=zeros(1017,1017);
        indy=TTLPost(ipost,:);
        bufy=cellImPost{ipost};
        sy=size(bufy);
        
        y(indy(1):indy(1)+sy(1)-1,indy(2):indy(2)+sy(2)-1)= ...
            y(indy(1):indy(1)+sy(1)-1,indy(2):indy(2)+sy(2)-1)+bufy;
        
        yReg = imwarp(y,tform{pdind},'OutputView',imref2d([1017,1017]));
        
        [max1,ind1]=max(yReg,[],1);
        [~,loc2]=max(max1);
        [~,loc1]=max(yReg(:,loc2));
        
        TTLPostReg(ipost,:)=[(floor((loc1-1)/245)*245)+1 , (floor((loc2-1)/245)*245)+1] ;
        s1=249;
        s2=249;
        if TTLPostReg(ipost,1)>=736
            TTLPostReg(ipost,1)=736;
            s1=281;
        end
        if TTLPostReg(ipost,2)>=736
            TTLPostReg(ipost,2)=736;
            s2=281;
        end
        cellImPostReg{ipost}=yReg(TTLPostReg(ipost,1):(TTLPostReg(ipost,1)+s1),...
            TTLPostReg(ipost,2):(TTLPostReg(ipost,2)+s2));
        
        [max2,ind2]= max(cellImPostReg{ipost},[],1);
        [~,loc22]= max(max2);
        [~,loc12]= max(cellImPostReg{ipost}(:,loc22));
        CentroidPostReg(ipost,:)=[loc12,loc22];

    end
    


    
    corrMat=zeros(length(cellImPre),length(cellImPost));
    
    for ipre=1:length(cellImPre)
        x=cellImPre{ipre};
        for ipost=1:length(cellImPostReg)
            if TTLPre(ipre,1)==TTLPostReg(ipost,1) && TTLPre(ipre,2)==TTLPostReg(ipost,2)
                if abs(CentroidPre(ipre,1)-CentroidPostReg(ipost,1)) <=DCent && abs(CentroidPre(ipre,2)-CentroidPostReg(ipost,2)) <=DCent 
                    c=corrcoef(x,cellImPostReg{ipost});
                    corrMat(ipre,ipost)=c(1,2);
              
                end
            end
        end
        
        
    end
    
    tcorr{pdind}=corrMat;
    corrMax=max(corrMat,[],1);
    [~,orderedIndex]=sort(corrMax,'descend');
    
    for k=orderedIndex
        [CM,IM]=max(corrMat(:,k));
         if CM>corrTh   
             cellIndex(IM,pdind)=k;
             corrMat(IM,:)=0;
             cellImPre{IM}=cellImPre{IM}+cellImPostReg{k};
             cellImages{IM,pdind}=cellImPostReg{k};
         else    
           cellImPre{cellNumber}=cellImPostReg{k};
           TTLPre(cellNumber,:)=TTLPostReg(k,:);
           CentroidPre(cellNumber,:)=CentroidPostReg(k,:);
           cellIndex(cellNumber,pdind)=k; 
           cellImages{cellNumber,pdind}=cellImPostReg{k};
           cellNumber=cellNumber+1;
         end
    end

end



% cellImage=cellImPre;
% TileTopLeft=TTLPre;
% save(strcat(LoadPath{Days(end)},'\registerationUnion'),'cellIndex','cellImage','TileTopLeft','-v7.3');
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Creat Cell Union Traces 
% %  
%  
% load(strcat(LoadPath{AllDays},'\registerationUnion.mat'));
% Traces={};
% sessionLength=zeros(1,length(Days));
%  for ind=1:length(Days)
%      ind
%      A=open(strcat(LoadPath{Days(ind)},'\All_Sessions.mat'));
%      Traces{ind}=A.cellData;
%      sessionLength(ind)=size(A.cellData,2);
%  end
%  progress='data load';
%  AllCells=zeros(size(cellIndex,1),sum(sessionLength));
%  CellsActiveRatio=zeros(size(cellIndex,1),1);
%  for cnum=1:size(cellIndex,1)
%      cnum
%      start=1;
%      for snum=1:size(cellIndex,2)
%          if cellIndex(cnum,snum)>0
%              AllCells(cnum,start:sum(sessionLength(1:snum)))=Traces{snum}(cellIndex(cnum,snum),:);     
%          end
%          start=start+sessionLength(snum);
%      end
%      CellsActiveRatio(cnum)=sum(sessionLength(cellIndex(cnum,:)>0))/sum(sessionLength);
%      
%      Noise=[AllCells(cnum,AllCells(cnum,:)<0),abs(AllCells(cnum,AllCells(cnum,:)<0))];
%      while length(Noise) < max(sessionLength)
%          Noise=[Noise,Noise];
%      end
%      Noise=Noise(randperm(length(Noise)));
%      start=1;
%      for snum=1:size(cellIndex,2)
%          if cellIndex(cnum,snum)==0
%              AllCells(cnum,start:sum(sessionLength(1:snum)))=Noise(randperm(length(Noise),sessionLength(snum)));     
%          end
%          start=start+sessionLength(snum);
%      end
%  end
% save(strcat(LoadPath{AllDays},'\AllCellsUnion'),'AllCells','CellsActiveRatio','-v7.3');
     
%      
     
             
 
 
 





