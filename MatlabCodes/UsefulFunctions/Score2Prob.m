%%%% Converting scores S= Nby1 of a binary softmax classifier to confidence probablities P=Nby1
function [P]=Score2Prob(S)
P=exp(S) ./ (exp(S)+exp(-1*S));
P=S;
end
