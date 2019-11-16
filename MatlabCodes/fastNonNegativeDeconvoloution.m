clear all
for Mouse=51:57;
    AddressSetup;
    load(strcat(LoadPath{Mouse},'\cellData.mat'));
    
    V.T       = size(cellData_Raw,2); % # of time steps
    V.dt    = 1/10;  % time step size
    V.smc_iter_max = 1;
    
    % initialize params
    P.a     = 1;    % observation scale
    P.b     = 0;    % observation bias
    tau     = 2;    % decay time constant
    P.gam   = 1-V.dt/tau; % C(t) = gam*C(t-1)
    P.lam   =0.05;
    P.sig   = 0.5;  % standard deviation of observation noise
    
    cellData_oopsi=zeros(size(cellData_Raw));
    
    
    for i=1:size(cellData_Raw,1)
        i
        cellData_oopsi(i,:)=(fast_oopsi(double(cellData_Raw(i,:)),V,P));
        
    end
    save(strcat(LoadPath{Mouse},'\cellData_oopsi.mat'),'cellData_oopsi','V','P','-v7.3');
end