function [B,intercept,error,ev,la]=glmnoreg(Set,group,trialNumber,crossVal,Repeat,doShuff)
ev=0;
la=0;
[B,FitInfo]=lassoglm(Set',group,'binomial','Lambda',0);
intercept=FitInfo.Intercept ;
error=mean( LRClassify(Set' , B , intercept) ~= group);

end


