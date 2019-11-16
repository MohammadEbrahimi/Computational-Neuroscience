global StorageArray
objNum=17; 
Day=23;

l=length(objNum);
Stats=zeros(l,5);
GoLick=zeros(l,1500);
NoGoLick=zeros(l,1500);
goCount=0;
nogoCount=0;
clear goTrialLength
clear nogoTrialLength

LORW=zeros(l,2);

for i=1:l
    Angle=StorageArray.Children{1, objNum(i)}.Object.Children{1, 2}.Object.Data;
    Lick=StorageArray.Children{1, objNum(i)}.Object.Children{1, 3}.Object.Data;
    RW=StorageArray.Children{1, objNum(i)}.Object.Children{1, 4}.Object.Data;
    GoTrials=StorageArray.Children{1, objNum(i)}.Object.Children{1, 5}.Object.Data;
    NoGoTrials=StorageArray.Children{1, objNum(i)}.Object.Children{1, 6}.Object.Data;
    AP=StorageArray.Children{1, objNum(i)}.Object.Children{1, 7}.Object.Data;
    WR=StorageArray.Children{1, objNum(i)}.Object.Children{1, 8}.Object.Data;
    T=1:length(GoTrials);
    
    Angle(1)=-1;
    Stats(i,3)=sum(diff(Angle)== -1);%#go
    Stats(i,4)=sum(diff(Angle)== -91);%#nogo
    
    
    GoSE=zeros(Stats(i,3),2);
    NoGoSE=zeros(Stats(i,4),2);
    buf=T(diff(Angle)== 1);
    GoSE=[buf(1:Stats(i,3))',T(diff(Angle)== -1)'];
    buf=T(diff(Angle)== 91);
    NoGoSE=[buf(1:Stats(i,4))',T(diff(Angle)== -91)'];

%%%% only if RW is corrupted    
%     for j=1:Stats(i,3)
%         TimePeriod=T(GoSE(j,1):GoSE(j,2));
%         E=max(TimePeriod(GoTrials(GoSE(j,1):GoSE(j,2))==1));
%         RW(E+30:E+211)=1;
%     end
%     for j=1:Stats(i,4)
%         TimePeriod=T(NoGoSE(j,1):NoGoSE(j,2));
%         E=max(TimePeriod(NoGoTrials(NoGoSE(j,1):NoGoSE(j,2))==1));
%         RW( E+30:min(E+211,NoGoSE(j,2)) )=1;
%     end
    
    
    %
    EOB=max(max([GoSE;NoGoSE]));
    Stats(i,1)=sum(diff(WR(1:EOB))==1)/Stats(i,3);%Hits
    Stats(i,2)=sum(([diff(AP(1:EOB));0].* (RW(1:EOB))==1).*(Angle(1:EOB)==90)) /Stats(i,4);%FA
    Stats(i,5)=EOB;
    

    for j=1:Stats(i,3)
        %%%%Lick Distribution during go
        TimePeriod=T(GoSE(j,1):GoSE(j,2));
        E=max(TimePeriod(GoTrials(GoSE(j,1):GoSE(j,2))==1));
        buf=Lick(GoSE(j,1):E);
        buf=buf(GoTrials(GoSE(j,1):E)==1);
        GoLick(i,1:length(buf))=GoLick(i,1:length(buf))+buf';
        E=min(TimePeriod(RW(GoSE(j,1):GoSE(j,2))==1))-30;
        lo=min(GoSE(j,2)-E+length(buf),1500);
        GoLick(i,length(buf)+1:lo)=GoLick(i,length(buf)+1:lo)+Lick(E+1:min(E+lo-length(buf),GoSE(j,2)))';
        %%%%LORW during go
        if (max(Lick(GoSE(j,1):E))>0)
            LORW(i,1)=LORW(i,1)+(1/Stats(i,3));
        end
        %%%% Pure stimuli before first lick
        firstLick=min(TimePeriod(GoTrials(GoSE(j,1)+1:GoSE(j,2))==0));
        goCount=goCount+1;
        goTrialLength(goCount)=firstLick-GoSE(j,1);
        
    end
    
    for j=1:Stats(i,4)
        %%%%Lick Distribution during nogo
        TimePeriod=T(NoGoSE(j,1):NoGoSE(j,2));
        E=max(TimePeriod(NoGoTrials(NoGoSE(j,1):NoGoSE(j,2))==1));
        buf=Lick(NoGoSE(j,1):E);
        buf=buf(NoGoTrials(NoGoSE(j,1):E)==1);
        NoGoLick(i,1:length(buf))=NoGoLick(i,1:length(buf))+buf';
        E=min(TimePeriod(RW(NoGoSE(j,1):NoGoSE(j,2))==1))-30;
        lo=min(NoGoSE(j,2)-E+length(buf),1500);
        NoGoLick(i,length(buf)+1:lo)=NoGoLick(i,length(buf)+1:lo)+Lick(E+1:min(E+lo-length(buf),NoGoSE(j,2)))';
        %%%%LORW during nogo
        if (max(Lick(NoGoSE(j,1):E))>0)
            LORW(i,2)=LORW(i,2)+(1/Stats(i,4));
        end
        %%%% Pure stimuli before first lick
        firstLick=min(TimePeriod(NoGoTrials(NoGoSE(j,1)+1:NoGoSE(j,2))==0));
        nogoCount=nogoCount+1;
        nogoTrialLength(nogoCount)=firstLick-NoGoSE(j,1);
    end
    
 
    
%     for j=1:T-1
%         if (GoStart(j)==1)
%             for(k=1:min(500,T-j-1))
%                 if (RWStart(j+k)==-1)
%                     GoLick(i,1:k)=GoLick(i,1:k)+Lick(j+1:j+k)';
%                     break;
%                 end
%             end
%         end
%         
%         if (NoGoStart(j)==1)
%             for(k=1:min(500,T-j-1))
%                 if (RWStart(j+k)==-1)
%                     NoGoLick(i,1:k)=NoGoLick(i,1:k)+Lick(j+1:j+k)';
%                     break;
%                 end
%             end
%         end
% 
%     end
%     
%%%%%% Stimuli : 1:120 /  Delay: 121:150 / RW: 150:330
    
end

% LickDist{Day}.go=GoLick;
% LickDist{Day}.nogo=NoGoLick;


%%%%%% plot pure trial length CDF
% goCDF=zeros(120,1);
% nogoCDF=zeros(120,1);
% goTrialLength=goTrialLength+1;
% nogoTrialLength=nogoTrialLength+1;
% 
% for i=1:120
%     goCDF(i)=100*sum(goTrialLength >= i)/length(goTrialLength);
%     nogoCDF(i)=100*sum(nogoTrialLength >= i)/length(nogoTrialLength);
% end
% 
% figure();plot((1:120)/60,goCDF,'b');
% hold on;plot((1:120)/60,nogoCDF,'r--');
% grid on;
% title('Trial Length CDF');
% legend('Go Trials','Nogo Trials');
% xlabel('time (s)')
% ylabel('percentage trials (%)')
% 
% 
% TotalGoTime=sum(goTrialLength)/60
% TotalNoGoTime=sum(nogoTrialLength)/60
% 








    















