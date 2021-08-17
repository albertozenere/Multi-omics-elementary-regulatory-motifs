%This script calculates some of the values in Tab.2

clear all
close all
clc

%% datasets
ALL_FOLDER_FILES = {'data\th1\', 'data\th1_p4\', 'data\Bonneau\hiv\', 'data\Bonneau\mock\'};
ALL_FIG_NAME = {'th1','p4','DC','mock'};

assert(length(ALL_FOLDER_FILES) == length(ALL_FIG_NAME))
n_dataset = length(ALL_FOLDER_FILES);

prc_cohe_FFL = zeros(n_dataset,1);
prc_cohe_null = zeros(n_dataset,1);

%% loop
for n_data = 1:n_dataset

    FOLDER_FILES = ALL_FOLDER_FILES{n_data};
    fprintf('Dataset: %s\n', FOLDER_FILES)
    FIG_NAME = ALL_FIG_NAME{n_data};
    
    %% load
    %data, all replicates
    Xtarget = importdata([FOLDER_FILES, 'Xall.mat']); Ntarget = size(Xtarget,1);
    Xtf = importdata([FOLDER_FILES, 'Xtfall.mat']); Ntf = size(Xtf,1);
    Xatac = importdata([FOLDER_FILES, 'Xatacall.mat']); Natac = size(Xatac,1);
    NTP = size(Xtarget,2);
    %data, mean of replicates, for comparison purposes
    Xtarget_mean = importdata([FOLDER_FILES, 'X.mat']); 
    Xtf_mean = importdata([FOLDER_FILES, 'Xtf.mat']); 
    Xatac_mean = importdata([FOLDER_FILES, 'Xatac.mat']);
    %adjacencies
    peak_tf = importdata([FOLDER_FILES, 'adj_peak_tf.mat']);
    peak_target = importdata([FOLDER_FILES, 'adj_peak_target.mat']);
    tf_target = importdata([FOLDER_FILES, 'adj_tf_target.mat']);
    %motifs
    FFL = importdata([FOLDER_FILES, 'FFL.mat']);
    Rtg = importdata([FOLDER_FILES, 'Rtg.mat']); Rat = importdata([FOLDER_FILES, 'Rat.mat']); Rag = importdata([FOLDER_FILES, 'Rag.mat']);
    ATT = importdata([FOLDER_FILES, 'ATT.mat']);
    Rat2 = importdata([FOLDER_FILES, 'Rat2.mat']); Rat1 = importdata([FOLDER_FILES, 'Rat1.mat']); Rtt = importdata([FOLDER_FILES, 'Rtt.mat']);
    AGG = importdata([FOLDER_FILES, 'AGG.mat']);
    Rag2 = importdata([FOLDER_FILES, 'Rag2.mat']); Rag1 = importdata([FOLDER_FILES, 'Rag1.mat']); Rgg = importdata([FOLDER_FILES, 'Rgg.mat']);
    %genes list
    list = importdata([FOLDER_FILES, 'list.mat']);
    list_tf = importdata([FOLDER_FILES, 'list_tf.mat']);

    if (contains(FOLDER_FILES, 'th1') && NTP == 18) || contains(FOLDER_FILES, 'th1_p4') 
        time = [0 0 0, .5 .5 .5, 1 1 1, 2 2 2, 6 6 6, 24 24 24]; %th1
        ORG = "human";
    elseif contains(FOLDER_FILES, 'hiv') || contains(FOLDER_FILES, 'mock')
        time = [2 2 2, 8 8 8, 24 24 24, 48 48 48]; %hiv
        ORG = "human";
    else
        error()
    end
    
    %% statistical test function
    ntp = NTP/3;
    alpha = 0.05;

    pc_fun = @(pxy,pxz,pyz) (pxy-pxz*pyz)/sqrt((1 - pxz^2)*(1 - pyz^2));    

    t_crit = @(n_control) tinv(.5+alpha/2, ntp-2-n_control);
    t = @(pc,n_control) pc.*sqrt( (ntp-2-n_control)./(1-pc.^2) );

    THS_0 = t_crit(0);
    THS_1 = t_crit(1);
    R_not_zero = @(vec) abs(t(vec,0))>THS_0;
    is_CI = @(vec) abs(t(vec,1))<THS_1; 
    
    %% baseline: random signals
	Niter = 5e4;
	X = randn(Niter,NTP);
	Y = randn(Niter,NTP);
	Z = randn(Niter,NTP);

    Rxy = zeros(Niter,1); Rxz = zeros(Niter,1); Ryz = zeros(Niter,1);
	for iter = 1:Niter
        Rxy(iter) = dyncorr(X(iter,:)', Y(iter,:)', time);
        Rxz(iter) = dyncorr(X(iter,:)', Z(iter,:)', time);
        Ryz(iter) = dyncorr(Y(iter,:)', Z(iter,:)', time);
	end
	to_keep = R_not_zero(Rxy) & R_not_zero(Rxz) & R_not_zero(Ryz);
    Rxy = Rxy(to_keep); Rxz = Rxz(to_keep); Ryz = Ryz(to_keep);
    %bootstrapping
    Niter2 = 1e4;
    Nsample = 1e4;
    Niter = sum(to_keep);
    prc_cohe_vec = zeros(Niter2,1);
    for n = 1:Niter2
        v = randi(Niter,Nsample,1);
        prc_cohe_vec(n) = 100*sum(Rxy(v).*Rxz(v).*Ryz(v)>0)/Nsample;
    end
    
    guass_vec = mean(prc_cohe_vec) + std(prc_cohe_vec)*randn(Niter2,1);
    
    %null distribution
    mu = mean(prc_cohe_vec);
    sigma = std(prc_cohe_vec);
    NULL_DISTRIBUTION = makedist('Normal','mu',mu,'sigma',sigma);
    
	fprintf('\t Null percentage of coherent cycles: %.f%% \n', mu)  
   
	%% FFL
    to_keep = R_not_zero(Rtg) & R_not_zero(Rat) & R_not_zero(Rag);
    FFL = FFL(to_keep,:);
    
    prc_cohe_FFL = 100*sum(FFL.is_cohe)/size(FFL,1);
	pval = cdf(NULL_DISTRIBUTION, prc_cohe_FFL, 'upper');
    
    fprintf('\t Number of FFL: %d\n', size(FFL,1))
    fprintf('\t Percentage of coherent FFLs: %.f%%; p-value: %e; FC: %d\n', prc_cohe_FFL, pval, prc_cohe_FFL/mu)
    
	%% ATT
    to_keep = R_not_zero(Rat1) & R_not_zero(Rat2) & R_not_zero(Rtt);
    ATT = ATT(to_keep,:);
    
    prc_cohe_ATT = 100*sum(ATT.is_cohe)/size(ATT,1);
	pval = cdf(NULL_DISTRIBUTION, prc_cohe_ATT, 'upper');
    
    fprintf('\t Number of ATT: %d\n', size(ATT,1))
    fprintf('\t Percentage of coherent ATTs: %.f%%; p-value: %e; FC: %d\n', prc_cohe_ATT, pval, prc_cohe_ATT/mu)
    
	%% AGG
    to_keep = R_not_zero(Rag1) & R_not_zero(Rag2) & R_not_zero(Rgg);
    AGG = AGG(to_keep,:);
    
    prc_cohe_AGG = 100*sum(AGG.is_cohe)/size(AGG,1);
	pval = cdf(NULL_DISTRIBUTION, prc_cohe_AGG, 'upper');
    
    fprintf('\t Number of AGG: %d\n', size(AGG,1))
    fprintf('\t Percentage of coherent AGGs: %.f%%; p-value: %e; FC: %d\n', prc_cohe_AGG, pval, prc_cohe_AGG/mu)

   
end