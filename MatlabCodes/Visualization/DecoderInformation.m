%%%%%%%%%%%%%%% Information
p='E:\Reports2\2018_03_01_FisherInfo\Raw_RT\Datasets-5to0\HitCR';
Mouse=64;
mode=1;
Sessions=[1:5];
C=load(strcat(p,'\Session',num2str(Mouse),'_mode',num2str(mode),'.mat'));
Decoders=C.PLSRotAll{9}*C.Decoder_All{9};
figure();hold on
plot(mean(C.Info_Val_Fisher(9,:,:),3),'r')
for s=Sessions
    load(strcat(p,'\Session',num2str(Mouse),'_mode',num2str(mode),'_Session',num2str(s),'.mat'));
    plot(mean(Info_Val_Fisher(9,:,:),3),'b')
    Decoders=[Decoders,PLSRotAll{9}*Decoder_All{9}];
end
figure();
imagesc(corrcoef(Decoders));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Projections 
p='E:\Reports2\2018_03_01_FisherInfo\Raw_Bin\Datasets-5to0\HitCR';
AddressSetup;
Mouse=61;
Sessions=[1,4,7,10,13,15];
mode=1;
load(strcat(LoadPath{Mouse-10},'\cellData_ZS.mat'));
load(strcat(LoadPath{Mouse-10},'\Datasets\Datasets-5to0'));
load(strcat(LoadPath{Mouse-10},'\SessionLength.mat'));
load(strcat(p,'\Session',num2str(Mouse),'_mode',num2str(mode),'.mat'));
color1={'b.','r.','g.','b.','m.'};
color2={'bo','ro','go','bo','mo'};
cellData=cellData_Raw_bin;
PLS=PLSRotAll{1,9}(:,1:2);
figure();hold on;
for i=1:length(Sessions)-1
    ts=sum(SessionLength(1:Sessions(i)-1))+1;
    te=sum(SessionLength(1:Sessions(i+1)-1));
    HitIdx=HitDataset{1,mode};
    HitIdx=HitIdx(HitIdx>=ts & HitIdx<=te);
    HitNum=HitTrialNumber{1,mode};
    HitNum=HitNum(HitIdx>=ts & HitIdx<=te);
    dd1=[];
    for trial=unique(HitNum)
        dd1=[dd1,sum(cellData(:,HitIdx(HitNum==trial)),2)];
    end
        
    
    CRIdx=CRDataset{1,mode,1};
    CRIdx=CRIdx(CRIdx>=ts & CRIdx<=te);
    CRNum=CRTrialNumber{1,mode};
    CRNum=CRNum(CRIdx>=ts & CRIdx<=te);
    dd0=[];
    for trial=unique(CRNum)
        dd0=[dd0,sum(cellData(:,CRIdx(CRNum==trial)),2)];
    end
    
%     HitPLSProj=PLS' * cellData(:,HitIdx);
%     CRPLSProj=PLS' * cellData(:,CRIdx);
    HitPLSProj=PLS' * dd1;
    CRPLSProj=PLS' * dd0;
    
    plot(HitPLSProj(1,:),HitPLSProj(2,:),color1{i});
    plot(CRPLSProj(1,:),CRPLSProj(2,:),color2{i});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%Instant Similarity
p='E:\Reports2\2018_03_02_FisherInstant\Raw_Bin\Datasets--5to0\HitCR';
Mouse=66;
rep=1000;
S=load(strcat(p,'\Session',num2str(Mouse),'_mode1.mat'));
D=load(strcat(p,'\Session',num2str(Mouse),'_mode2.mat'));
DecSim=zeros(30,30);

for i=1:30
    i
    for j=1:30
        for r=1:rep
            if i<=25
                x=S.MaxDecoders{9,i,randi(100)};
            else
                x=D.MaxDecoders{9,i-25,randi(100)};
            end
            
            if j<=25
                y=S.MaxDecoders{9,j,randi(100)};
            else
                y=D.MaxDecoders{9,j-25,randi(100)};
            end
            
            C=corrcoef([x,y]);
            DecSim(i,j)=DecSim(i,j)+C(1,2);
        end
    end
end
DecSim=DecSim/rep;


















