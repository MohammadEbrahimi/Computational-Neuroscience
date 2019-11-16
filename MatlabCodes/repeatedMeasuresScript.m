F=zeros(8,2);
P=zeros(8,1);
for a=1:8
filter=(min([hc(1:7,a)';mc(1:7,a)'])>0);
F(a,2)=sum(filter);
F(a,1)=repeatedMeasures([hc(filter,a)';mc(filter,a)']);
[P(a),h(a)]=signrank(hc(filter,a),mc(filter,a));
end
