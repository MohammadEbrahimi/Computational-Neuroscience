function y=randu(p)
y=rand(p,1)*3-1.5;
y=sign(y).*(abs(y)+0.5);
end