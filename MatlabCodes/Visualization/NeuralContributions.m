p='E:\Reports2\2018_05_15_NeuralContributions\ROI_Bin\Datasets-5to0_allCells\HitCR_4pixCells\';

 figure();hold on;
 set(0,'defaultlinelinewidth',1.5);
for Mouse=[61,63:67]
    S=load(strcat(p,'Session',num2str(Mouse),'_mode1.mat'));
    D=load(strcat(p,'Session',num2str(Mouse),'_mode2.mat'));
    cellCount=size(S.Info_Val_Fisher,1);
%     SN=(S.Info_Val_Fisher(1:cellCount,:)-repmat(S.Info_Val_Fisher(end,:),[cellCount,1]))...
%         ./repmat(S.Info_Val_Fisher(end,:),[cellCount,1]);
%     DN=(D.Info_Val_Fisher(1:cellCount,:)-repmat(D.Info_Val_Fisher(end,:),[cellCount,1]))...
%         ./repmat(D.Info_Val_Fisher(end,:),[cellCount,1]);
    
    SN=(S.Error_Val_Fisher(1:cellCount,:)-repmat(S.Error_Val_Fisher(end,:),[cellCount,1]));
    DN=(D.Error_Val_Fisher(1:cellCount,:)-repmat(D.Error_Val_Fisher(end,:),[cellCount,1]));
    
    MSN=sort(mean(SN,2),'descend');
    MDN=sort(mean(DN,2),'descend');
    
     th_s=5*sqrt(var([MSN(MSN>0);-1*MSN(MSN>0)]));
     th_d=5*sqrt(var([MDN(MDN>0);-1*MDN(MDN>0)]));
%     th_s=max(MSN);
%     th_d=max(MDN);
    Vari(Mouse-60,1)=var(MSN< -1*th_s);
    Vari(Mouse-60,2)=var(MDN< -1*th_d);
   
    plot(MSN,'b');
    %plot([1,cellCount],-1*[th_s,th_s],'b-.')
    plot(MDN,'r');
    %plot([1,cellCount],-1*[th_d,th_d],'r--')
    
end


