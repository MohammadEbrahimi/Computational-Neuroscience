AddressSetup;

Days=[33:35,39];
ND=length(Days);


clear tform
clear tcell
clear tcorr

load('C:\Users\Sadegh\Documents\VLMReborn\L367\Data\allDays\registeration.mat')


for i=1:ND
    i
    preIm=open(strcat(LoadPath{Days(i)},'\cellImage.mat'));
    
    
    pdind=i+1;
    if i==ND
        pdind=1;
    end
    
    postIm=open(strcat(LoadPath{Days(pdind)},'\cellImage.mat'));
    
    cellImPre=preIm.cellImage;
    TTLPre=preIm.TileTopLeft;
    
    cellImPost=postIm.cellImage;
    TTLPost=postIm.TileTopLeft;

    

    clear cellImPostReg
    TTLPostReg=zeros(size(TTLPost));
    TTLVal=[1,246,491,736];
    for ipost=1:length(cellImPost) 
        ipost;
        y=zeros(1017,1017);
        indy=TTLPost(ipost,:);
        bufy=cellImPost{ipost};
        sy=size(bufy);
        
        y(indy(1):indy(1)+sy(1)-1,indy(2):indy(2)+sy(2)-1)= ...
            y(indy(1):indy(1)+sy(1)-1,indy(2):indy(2)+sy(2)-1)+bufy;
        
        yReg = imwarp(y,tform{i},'OutputView',imref2d([1017,1017]));
        
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
        

    end
    


    corrMat=tcorr{i};
    for ipre=1:length(cellImPre)
        ipre;
        x=cellImPre{ipre};
        for ipost=1:length(cellImPostReg)

            if TTLPre(ipre,1)==TTLPostReg(ipost,1) && TTLPre(ipre,2)==TTLPostReg(ipost,2)

                %%%%centroid distance 
                [max1,ind1]= max(x,[],1);
                [~,loc2]= max(max1);
                [~,loc1]= max(x(:,loc2));
                [max2,ind2]= max(cellImPostReg{ipost},[],1);
                [~,loc22]= max(max2);
                [~,loc12]= max(cellImPostReg{ipost}(:,loc22));
                if abs(loc2-loc22)>10 || abs(loc1-loc12)>10
                    corrMat(ipre,ipost)=0;
                end
            end
        end
        
        
    end
    
   
    
    [corrMax,indMax]=max(corrMat,[],2);
    tcorr{i}=corrMax;
    indMax(corrMax<0.5)=0;
    tcell{i}=indMax;


end




cellMap=zeros(length(tcell{1}),ND);
cellCounter=0;

for i=1:length(tcell{1})
    start=i;
    for j=1:ND
        vec=tcell{j};
        start=vec(start);
        if start==0
            cellMap(i,:)=0;
            break
        else
        cellMap(i,j)=start;
        end
    end
    
    if start==i
        cellCounter=cellCounter+1;   
    else
        cellMap(i,:)=0;
    end
     
end

save('C:\Users\Sadegh\Documents\VLMReborn\L367\Data\allDays\registeration2','tcell','tcorr','tform','cellMap','-v7.3');




