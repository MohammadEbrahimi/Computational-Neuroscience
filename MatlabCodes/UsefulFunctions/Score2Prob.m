function [P]=Score2Prob(S)
P=exp(S) ./ (exp(S)+exp(-1*S));
P=S;
end