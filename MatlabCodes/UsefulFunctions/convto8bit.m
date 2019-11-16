function [B]=convto8bit(A)

B=A-min(min(A));
B=uint8(floor(255*B/max(max(B))));


end