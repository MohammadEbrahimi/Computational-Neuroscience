    
function [IcaTraces,IcaFilters]=PCA_ICA(DataIn,nPCs,nICs)


% Do the PCA calculation
% This is to bypass Matlab slow OOP
% We rely on copy on write to avoid unnecessary memory duplication of large
% matrix
DataSize=size(DataIn);
Npixels=DataSize(1)*DataSize(2);
Ntime=DataSize(3);


% Create covariance matrix in the time domain as it is computationaly more
% advantageous than in space and mathematically equivalent.
covmat=zeros(Ntime,Ntime,'single');

Block=200;
IterBlock=1:Block:DataSize(3);
k=0;

% We do a calculation per block to be able to process very large movies
% that do not hold fully in memory. 
for i=IterBlock
    k=k+1;
    EndBlock1=min(i-1+Block,DataSize(3));
    Block1=single(squeeze(DataIn(:,:,i:EndBlock1,1)));
    Block1=reshape(Block1,Npixels,EndBlock1-i+1);
    Block1=bsxfun(@minus,Block1,sum(Block1,1)/Npixels);
    
    covmat(i:EndBlock1,i:EndBlock1)=Block1'*Block1/(Npixels-1);
    for j=((EndBlock1+1):Block:DataSize(3))
        EndBlock2=min(j-1+Block,DataSize(3));
        Block2=single(squeeze(DataIn(:,:,j:EndBlock2,1)));
        Block2=reshape(Block2,Npixels,EndBlock2-j+1);
        Block2=bsxfun(@minus,Block2,sum(Block2,1)/Npixels);
    
        CORR=Block1'*Block2/ (Npixels-1);
        covmat(i:EndBlock1,j:EndBlock2)=CORR;
        covmat(j:EndBlock2,i:EndBlock1)=CORR';
    end
    
    %multiWaitbar('Calculating covariance matrix',i/DataSize(3));
end
covmat=double(covmat);
%multiWaitbar('Calculating covariance matrix','Close');

%multiWaitbar('Performing PCA',1/3);

% Options for the Eigenvectors extraction
opts.issym = 'true';
        
if nPCs<size(covmat,1)
    [PcaTraces, CovEvals] = eigs(covmat, nPCs, 'LM', opts); 
else
    [PcaTraces, CovEvals] = eig(covmat);
    nPCs = size(CovEvals,1);
end

% At this stage PcaTraces is Ntime x nPCs, ie each column is an eigenvector
% CovEvals is a square diagonal matrix with the eigenvalues on the diagonal

% We only keep the diagonal values, ie we no longer have a diagonal matrix
%multiWaitbar('Performing PCA',2/3);

% Throwing out negative eigenvalues
CovEvals=diag(CovEvals);
if nnz(CovEvals<=0)
    nPCs=nPCs - nnz(CovEvals<=0);
    PcaTraces = PcaTraces(:,CovEvals>0);
    CovEvals = CovEvals(CovEvals>0);
end
CovEvals=diag(CovEvals);

% This is because the SVD decomposition is not of the covariance exactly but of
% XX(t). Check wikipedia.
CovEvals=CovEvals*Npixels;

% This is to calculate the PcaFilters. Eigenvalues are the square root.
% Check Wikipedia again. 
Sinv = inv(CovEvals.^(1/2));

% We calculate the corresponding space filters
PcaFilters=zeros(Npixels,nPCs);
for i=IterBlock
    EndBlock=min(i-1+Block,DataSize(3));
    
    DataBlock=double(DataIn(:,:,i:EndBlock,1));
    
    % We substract the spatial mean to each frame to ensure proper
    % calculation of filters
    DataBlock=DataBlock-repmat(mean(mean(DataBlock,1),2),[DataSize(1) DataSize(2)]);
    
    PcaFilters = PcaFilters+reshape(DataBlock,Npixels,1+EndBlock-i)*PcaTraces(i:EndBlock,:)*Sinv;
   % multiWaitbar('Projecting PC filters',i/DataSize(3));
end
%multiWaitbar('Projecting PC filters','Close');

% Now because the Filters are not an EXACT calculation, their variance can
% be slighlty off 1. We make sure it is 1, as this can slightly affect the 
% spatio-temporal ICA that expect normalized eigenvectors.
for i=1:nPCs
    PcaFilters(:,i)=PcaFilters(:,i)/norm(PcaFilters(:,i));
end

%multiWaitbar('Performing PCA','Close');

% Do the ICA calculation
%multiWaitbar('Performing ICA',0);

% We need to take nPCs if some negative PCs appeared and nICs>nPCs


Mu=0.1;
TermTolICs=1e-5;;
MaxRoundsICs=750;

% Seed for the ICs calculation
ica_A_guess = orth(randn(nPCs, nICs));

PcaTraces=PcaTraces';
PcaFilters=PcaFilters';

% Center the data by removing the mean of each PC
% They were already centered after PC but just making sure.
meanTraces = mean(PcaTraces,2);
PcaTraces = PcaTraces - meanTraces * ones(1, size(PcaTraces,2));

meanFilters = mean(PcaFilters,2);
PcaFilters = PcaFilters - meanFilters * ones(1, size(PcaFilters,2));

% Whiten PcaTraces and PcaFilters relative to one another
% ie they were normalized but not of variance 1
% In theory, they are sqrt(N) now each where N is Npixel or Ntime
% Make sure it is variance 1 even after rounding errors
PcaTraces=PcaTraces./repmat(sqrt(var(PcaTraces,0,2)),[1 Ntime]);
PcaFilters=PcaFilters./repmat(sqrt(var(PcaFilters,0,2)),[1 Npixels]);

% Create concatenated data for spatio-temporal ICA
if Mu == 1
    % Pure temporal ICA
    IcaMixed = PcaTraces;
elseif Mu == 0
    % Pure spatial ICA
    IcaMixed = PcaFilters;
else
    % Spatial-temporal ICA
    IcaMixed = [(1-Mu)*PcaFilters, Mu*PcaTraces];
    
    % Again whitening
    IcaMixed=IcaMixed./repmat(sqrt(var(IcaMixed,0,2)),[1 size(IcaMixed,2)]);
end

% Perform ICA
numSamples = size(IcaMixed,2);
ica_A = ica_A_guess;
BOld = zeros(size(ica_A));
numiter = 0;
minAbsCos = 0;

% We preallocate an intermediate matrix
Interm=zeros(size(IcaMixed,2),nICs,'double');

while (numiter<MaxRoundsICs) && ((1-minAbsCos)>TermTolICs)
    numiter = numiter + 1;
    
    % The core of ICA based on Hyvärinen work (1999) on fastICA
    % unmixing matrix is updated based on previous step and using skewness
    % as a contrast function to measure non gaussianity.
    if numiter>1
        Interm=IcaMixed'*ica_A;
        Interm=Interm.^2;
        ica_A = IcaMixed*Interm/numSamples;
    end
    
    % Symmetric orthogonalization.
    ica_A = ica_A * real(inv(ica_A' * ica_A)^(1/2));
    
    % Test for termination condition.
    minAbsCos = min(abs(diag(ica_A' * BOld)));
    
    BOld = ica_A;
    
    %multiWaitbar('Performing ICA',numiter/MaxRoundsICs);
end
%multiWaitbar('Performing ICA','Close');

ica_W = ica_A';

% Add the mean back in.
IcaTraces = ica_W*PcaTraces+ica_W*(meanTraces*ones(1,size(PcaTraces,2)));
IcaFilters = ica_W*PcaFilters+ica_W*(meanFilters*ones(1,size(PcaFilters,2)));

% Sort ICs according to skewness of the temporal component
icskew = Matskewness(IcaTraces');
[~, ICCoord] = sort(icskew, 'descend');
IcaTraces = IcaTraces(ICCoord,:);
IcaFilters = IcaFilters(ICCoord,:);

% Reshape the filter to have a proper image
IcaFilters = reshape(IcaFilters,nICs,DataSize(1),DataSize(2));

end
    
    
  
% This is to get the skewness of a matrix
function s = Matskewness(x)

% Need to tile the output of mean to center X.
tile = ones(1,max(ndims(x),1));
tile(1) = size(x,1);

% Center X, compute its third and second moments, and compute the
% uncorrected skewness.
x0 = x - repmat(mean(x,1), tile);
s2 = mean(x0.^2,1); % this is the biased variance estimator
m3 = mean(x0.^3,1);
s = m3 ./ s2.^(1.5);
end  
    
    