function U2=normCol(U1)
U2=U1;
for c=1:size(U1,2)
    U2(:,c)=U1(:,c)/norm(U1(:,c));

end

end