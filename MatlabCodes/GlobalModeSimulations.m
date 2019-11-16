%%%%Simulations Outputs
TransMtx={};
AdjMtx={};
CCAModesX={};CCAModesY={};
PCAModes={};
CCAVal={};
CCATrain={};
%%%%%%%%%%%Iterations
Iter=0;

for K=[2,4,6,8,10]
    Iter=Iter+1;
    Iterb=0;
    for beta=[0,0.2,0.4,0.6,0.8,1]
        Iterb=Iterb+1;
        x={};
        xn={};
        W={};
        N=20;
        %K=2;%mean degree (even)
        %beta=0.5;%randomness parameter
        cellPerArea=100;
        x_int=1;
        w_int=1/sqrt(N);
        sOut=0.03;%%% internal random added noise
        TMax=50000;
        xRecording=zeros(N,cellPerArea,TMax);
        %%%% Small-World Network Construction
        A=eye(N);
        for i=1:N
            for j=1:N
                if (min(abs(i-j),N-abs(i-j))<=(K/2) && i~=j)
                    A(i,j)=1;
                end
            end
        end
        for i=1:N
            for j=1:N
                if (A(i,j)==1 && rand(1)<beta  && i~=j)
                    A(i,j)=0;
                    newNode=randperm(N,1);
                    while (newNode==i || newNode==j) newNode=randperm(N,1); end
                    A(i,newNode)=1;
                    
                end
            end
        end
        AdjMtx{Iter,Iterb}=A;
        %%%%Initialization
        for i=1:N
            x{i}=randn(cellPerArea,1);
            x{i}=x{i}/norm(x{i});
            for j=1:N
                if A(i,j)==1
                    W{i,j}=randn(cellPerArea,cellPerArea);
                    for nc=1:cellPerArea
                        W{i,j}(nc,:)=w_int*( W{i,j}(nc,:) / norm(W{i,j}(nc,:)));
                    end
                end
                
            end
        end
        TransMtx{Iter,Iterb}=W;
        %%%%Network Stimulation
        for tc=1:TMax
            for j=1:N
                xRecording(j,:,tc)=x{j};
                noise=randn(cellPerArea,1);
                noise=noise/norm(noise);
                %%xn{i}=W{1,i}*x{1} + W{2,i}*x{2} + W{3,i}*x{3} + (sOut*randn(cellPerArea,1)) ;
                xn{j}=(sOut*randn(cellPerArea,1)) ;
                for i=1:N
                    if A(i,j)==1
                        xn{j}=xn{j}+W{i,j}*x{i};
                    end
                end
                xn{j}=xn{j}/norm(xn{j});
            end
            x=xn;
        end
        
        %%%%Fluctuation Analysis
        Wx={};
        Wy={};
        rTrain=zeros(N,N,20);
        rVal=zeros(N,N,20);
        
        for i=1:N
            i
            for j=1:N
                j
                if i<j
                    Wx{i,j}=zeros(cellPerArea,20);
                    Wy{i,j}=zeros(cellPerArea,20);
                    [Wx{i,j},Wy{i,j},rTrain(i,j,:)]=SparseCCA(squeeze(xRecording(i,:,1:20000))',squeeze(xRecording(j,:,1:20000))',1,1,2,20);
                    for cc=1:20
                        t1=squeeze(xRecording(i,:,20000:40000))'*squeeze(Wx{i,j}(:,cc));
                        t2=squeeze(xRecording(j,:,20000:40000))'*squeeze(Wy{i,j}(:,cc));
                        rVal(i,j,cc)=similarity([t1,t2]);
                    end
                    
                end
            end
        end
        CCAModesX{Iter,Iterb}=Wx;
        CCAModesY{Iter,Iterb}=Wy;
        
        CCAVal{Iter,Iterb}=rVal;
        CCATrain{Iter,Iterb}=rTrain;
        %%%%%%%%%%%%%%%%% PCA Analysis
        
        [COEFF, SCORE, LATENT] = pca(reshape(xRecording,[N*cellPerArea,TMax])','NumComponents',20);
        
        PCAModes{Iter,Iterb}=COEFF;
    end
end
%%%%%%%%%%Save Results
save('E:\Reports2\2019_07_09_GlobalModeSimulations\SmallWorldSimulations','TransMtx','AdjMtx','CCAModesX','CCAModesY','PCAModes','CCAVal','CCATrain','-v7.3');

%%%%%%%%%%%%%%%% Visualization
%%%CCAVal
cc=1;
res=[];
for k=1:5
    buf=[];
    for b=1:6
        buf=[buf;squeeze(CCAVal{k,b}(:,:,cc))];
    end
    res=[res,buf];
end
figure();imagesc(res);
%%%%%%%%%%%%%%
%%%%CCA Modes X,Y
SimMtx=zeros(N,20,6,N,N,5);
for I1=1:5
    I1
    for I2=1:6
        I2
Wx=CCAModesX{I1,I2};
Wy=CCAModesY{I1,I2};

for ccn=1:20
for s=1:N
for d1=1:N
    for d2=1:N
        if s~=d1 && s~=d2 && d1~=d2
            if s<d1
                w1=Wx{s,d1}(:,ccn);
            else
                w1=Wy{d1,s}(:,ccn);
            end
            if s<d2
                w2=Wx{s,d2}(:,ccn);
            else
                w2=Wy{d2,s}(:,ccn);
            end
            SimMtx(d1,ccn,I2,d2,s,I1)=abs(similarity([w1,w2]));
        end
    end
end
end
end
    end
end
imagesc(squeeze(reshape(SimMtx,[N*N*6,N*N*5])))





            




























