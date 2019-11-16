TileC=[1,1;1,246;1,491;1,736;246,1;246,246;246,491;246,736;491,...
    1;491,246;491,491;491,736;736,1;736,246;736,491;736,736];

Map123=zeros(length(TTL1),2);


for t=1:16

f1=(TTL1(:,1)==TileC(t,1)) & (TTL1(:,2)==TileC(t,2));
f2=(TTL2(:,1)==TileC(t,1)) & (TTL2(:,2)==TileC(t,2));
f3=(TTL3(:,1)==TileC(t,1)) & (TTL3(:,2)==TileC(t,2));

corrMap12=zeros(length(f1),2);
corrMap23=zeros(length(f2),2);
corrMap31=zeros(length(f3),2);

for i=min(find(f1)):max(find(f1))
   ind=0;
   maxCorr=0;
   for j=min(find(f2)):max(find(f2))
       C=sum(sum(cellImage1{i}.*cellImage2{j}))/( sqrt(sum(sum(cellImage1{i}.*cellImage1{i}))) * ...
           sqrt(sum(sum(cellImage2{j}.*cellImage2{j}))));
       if(C > maxCorr)
           ind=j;
           maxCorr=C;
       end
   end
   corrMap12(i,1)=ind;
   corrMap12(i,2)=maxCorr;
end

for i=min(find(f2)):max(find(f2))
   ind=0;
   maxCorr=0;
   for j=min(find(f3)):max(find(f3))
       C=sum(sum(cellImage2{i}.*cellImage3{j}))/( sqrt(sum(sum(cellImage2{i}.*cellImage2{i}))) * ...
           sqrt(sum(sum(cellImage3{j}.*cellImage3{j}))));
       if(C > maxCorr)
           ind=j;
           maxCorr=C;
       end
   end
   corrMap23(i,1)=ind;
   corrMap23(i,2)=maxCorr;
end

for i=min(find(f3)):max(find(f3))
   ind=0;
   maxCorr=0;
   for j=min(find(f1)):max(find(f1))
       C=sum(sum(cellImage3{i}.*cellImage1{j}))/( sqrt(sum(sum(cellImage3{i}.*cellImage3{i}))) * ...
           sqrt(sum(sum(cellImage1{j}.*cellImage1{j}))));
       if(C > maxCorr)
           ind=j;
           maxCorr=C;
       end
   end
   corrMap31(i,1)=ind;
   corrMap31(i,2)=maxCorr;
end


for i=min(find(f1)):max(find(f1))
    for j=min(find(f1)):max(find(f1))
        if ((i ~= j) && corrMap12(i,1)==corrMap12(j,1) && corrMap12(i,1)>0 )
            if(corrMap12(i,2)>corrMap12(j,2))
                corrMap12(j,:)=0;
            else
                corrMap12(i,:)=0;
            end
        end
    end
end

for i=min(find(f2)):max(find(f2))
    for j=min(find(f2)):max(find(f2))
        if ((i ~= j) && corrMap23(i,1)==corrMap23(j,1) && corrMap23(i,1)>0 )
            if(corrMap23(i,2)>corrMap23(j,2))
                corrMap23(j,:)=0;
            else
                corrMap23(i,:)=0;
            end
        end
    end
end
       
for i=min(find(f3)):max(find(f3))
    for j=min(find(f3)):max(find(f3))
        if ((i ~= j) && corrMap31(i,1)==corrMap31(j,1) && corrMap31(i,1)>0 )
            if(corrMap31(i,2)>corrMap31(j,2))
                corrMap31(j,:)=0;
            else
                corrMap31(i,:)=0;
            end
        end
    end
end

for i=min(find(f1)):max(find(f1))
    if(corrMap12(i,1)>0)
        if(corrMap23(corrMap12(i,1),1) > 0 )
            if(corrMap31(corrMap23(corrMap12(i,1),1))==i)
                Map123(i,1)=corrMap12(i,1);
                Map123(i,2)=corrMap23(corrMap12(i,1),1);
            end
        end
    end
end


end


cn=561;
figure();imagesc(cellImage1{cn});
figure();imagesc(cellImage2{Map123(cn,1)});
figure();imagesc(cellImage3{Map123(cn,2)});





