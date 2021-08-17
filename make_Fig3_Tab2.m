%This script is used to generate Fig.3-S2-S3-S4 and part Tab.2  

clear all
close all
addpath('bin/');

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

    list = importdata([FOLDER_FILES, 'list.mat']);
    list_tf = importdata([FOLDER_FILES, 'list_tf.mat']);

    if (contains(FOLDER_FILES, 'th1') && NTP == 18) || contains(FOLDER_FILES, 'th1_p4') || contains(FOLDER_FILES, 'control')
        time = [0 0 0, .5 .5 .5, 1 1 1, 2 2 2, 6 6 6, 24 24 24]; %th1
        ORG = "human";
    elseif contains(FOLDER_FILES, 'hiv') || contains(FOLDER_FILES, 'mock')
        time = [2 2 2, 8 8 8, 24 24 24, 48 48 48]; %hiv
        ORG = "human";
    else
        error()
    end
 
    %% statistical test
	pc_fun = @(pxy,pxz,pyz) (pxy-pxz*pyz)/sqrt((1 - pxz^2)*(1 - pyz^2));    
    
    ntp = NTP/3;
    alpha = 0.05;
    t_crit = @(n_control) tinv(.5+alpha/2, ntp-2-n_control);
    t = @(pc,n_control) pc.*sqrt( (ntp-2-n_control)./(1-pc.^2) );
    close_zero_fun = @(pc,n_control) abs( t(pc,n_control) ) < t_crit(n_control);

    THS_0 = t_crit(0);
    THS_1 = t_crit(1);
    R_not_zero = @(vec) abs(t(vec,0))>THS_0;
    is_pc = @(vec) abs(t(vec,1))<THS_1;
    
    pval_fun = @(pc,n_control) tcdf(t(pc,n_control),ntp-2-n_control);
    
    %% indeces of FFLs
    %calculate number of cycles
    n_FFL = 0;
    for n = 1:Natac
       n_tf = sum( peak_tf(n,:) );
       n_target = sum( peak_target(n,:) );

       n_FFL = n_FFL + n_tf*n_target;

    end
    %indeces of all FFLs
    pos_tf = zeros(n_FFL,1);
    pos_peak = zeros(n_FFL,1);
    pos_target = zeros(n_FFL,1);

    [Itf, Jtarget] = find(tf_target);
    Nint = length(Itf);
    c = 1;
    for n = 1:Nint
        Kpeak = find(peak_tf(:,Itf(n)) & peak_target(:,Jtarget(n)));
        n_peak = length(Kpeak);

        pos_tf(c:c+n_peak-1) = Itf(n);
        pos_peak(c:c+n_peak-1) = Kpeak;
        pos_target(c:c+n_peak-1) = Jtarget(n);

        c = c + n_peak;
    end
    assert(c==n_FFL+1)
    
    %% percentage coherent FFLs in data
    Niter = n_FFL;
    v = randperm( n_FFL,Niter );
    Itf = pos_tf(v);
    Jtarget = pos_target(v);
    Kpeak = pos_peak(v);

    is_cohe = false(Niter,1);
    is_not_zero = false(Niter,1);
    Rtg = zeros(Niter,1); 
    Rat = zeros(Niter,1); 
    Rag = zeros(Niter,1); 
    for n = 1:Niter
        Rtg(n) = dyncorr(Xtf(Itf(n),:)', Xtarget(Jtarget(n),:)', time);
        Rat(n) = dyncorr(Xatac(Kpeak(n),:)', Xtf(Itf(n),:)', time);
        Rag(n) = dyncorr(Xatac(Kpeak(n),:)', Xtarget(Jtarget(n),:)', time);

        is_cohe(n) = Rtg(n)*Rag(n)*Rat(n)>0;
        is_not_zero(n) = all(R_not_zero([Rtg(n),Rag(n),Rat(n)]));
    end

    prc_cohe_FFL(n_data) = 100*sum(is_cohe)/Niter;
    prc_cohe_FFL(n_data) = 100*sum(is_cohe & is_not_zero)/sum(is_not_zero);
    fprintf('\t Data: %d/%d=%.f%% FFLs are coherent\n', sum(is_cohe), Niter, prc_cohe_FFL(n_data))
    
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
    %bootstrapping
    Niter2 = 1e4;
    Nsample = 1e4;
    prc_cohe_vec = zeros(Niter2,1);
    for n = 1:Niter2
        v = randi(Niter,Nsample,1);
        prc_cohe_vec(n) = 100*sum(Rxy(v).*Rxz(v).*Ryz(v)>0)/Nsample;
    end
    
    guass_vec = mean(prc_cohe_vec) + std(prc_cohe_vec)*randn(Niter2,1);
    
    mu = mean(prc_cohe_vec);
    sigma = std(prc_cohe_vec);
    pd = makedist('Normal','mu',mu,'sigma',sigma);
    pval = cdf(pd,prc_cohe_FFL(n_data),'upper');
    
    prc_cohe_null(n_data) = mu;
	fprintf('\t Null: %.f%% FFLs are coherent on average\n', mu)
    fprintf('\t Enrichment of coherent FFLs with p-value: %e\n', pval)
   
	%% R(T,G|A), FFL
    Niter = n_FFL;
    
    v = randperm(n_FFL,Niter);
    pos_T = pos_tf(v);
    pos_G = pos_target(v);
    pos_A = pos_peak(v);
    
    pval_tga = zeros(Niter,1);
    Rtg = zeros(Niter,1);
    Rag = zeros(Niter,1);
    Rat = zeros(Niter,1);
    Rtga = zeros(Niter,1);
    is_cohe = false(Niter,1);
    for iter = 1:Niter
        %three dyncorr
        Rat(iter) = dyncorr(Xtf(pos_T(iter),:)', Xatac(pos_A(iter),:)', time);
        Rag(iter) = dyncorr(Xtarget(pos_G(iter),:)', Xatac(pos_A(iter),:)', time);
        Rtg(iter) = dyncorr(Xtf(pos_T(iter),:)', Xtarget(pos_G(iter),:)', time);
        %coherent cycle?
        is_cohe(iter) = Rat(iter)*Rag(iter)*Rtg(iter)>0;
        %partial correlation
        pc = pc_fun(Rtg(iter),Rat(iter),Rag(iter));
        Rtga(iter) = pc;
        assert(isreal(Rtga(iter)))
        %p-value
        pval_tga(iter) = pval_fun(pc,1);
    end
	%remove r close to zero
    r_not_zero = R_not_zero(Rat) & R_not_zero(Rag) & R_not_zero(Rtg);
    is_cohe = is_cohe(r_not_zero);
    Rtga = Rtga(r_not_zero);
    
    fprintf('\t C.i. T->A->G %d/%d=%.2f%% \n', ...
        sum(is_pc(Rtga)), numel(Rtga), 100*sum(is_pc(Rtga))/numel(Rtga))
    
    %plot count
    figure
    print_count({Rtga(is_cohe),Rtga(~is_cohe)},'R(T,G|A)',['IMG/Rtg_a_orig_', FIG_NAME])    
    	
    %plot density
	figure
    print_density({Rtga(is_cohe),Rtga(~is_cohe)},'R(T,G|A)',['IMG/Rtg_a_', FIG_NAME])
    
	%% R(A,G|T), FFL
    Niter = n_FFL;
    
    v = randperm(n_FFL,Niter);
    pos_T = pos_tf(v);
    pos_G = pos_target(v);
    pos_A = pos_peak(v);

    pval_tag = zeros(Niter,1);
    Rag = zeros(Niter,1);
    Rat = zeros(Niter,1);
    Rtg = zeros(Niter,1);
    Ragt = zeros(Niter,1);
    is_cohe = false(Niter,1);
    for iter = 1:Niter
        %three dyncorr
        Rat(iter) = dyncorr(Xtf(pos_T(iter),:)', Xatac(pos_A(iter),:)', time);
        Rtg(iter) = dyncorr(Xtf(pos_T(iter),:)', Xtarget(pos_G(iter),:)', time);
        Rag(iter) = dyncorr(Xtarget(pos_G(iter),:)', Xatac(pos_A(iter),:)', time);
        %coherent cycle?
        is_cohe(iter) = Rat(iter)*Rag(iter)*Rtg(iter)>0;
        %partial correlation
        pc = pc_fun(Rag(iter),Rat(iter),Rtg(iter));
        Ragt(iter) = pc;
        assert(isreal(Ragt(iter)))
        %p-value
        pval_tag(iter) = pval_fun(pc,1);
    end   
    %remove r close to zero
    r_not_zero = R_not_zero(Rtg) & R_not_zero(Rag) & R_not_zero(Rat); 
    Ragt = Ragt(r_not_zero);
    is_cohe = is_cohe(r_not_zero);   
    
	fprintf('\t C.i. A->T->G %d/%d=%.2f%% \n', ...
        sum(is_pc(Ragt)), numel(Ragt), 100*sum(is_pc(Ragt))/numel(Ragt))
    
    %plot count
	figure
    print_count({Ragt(is_cohe),Ragt(~is_cohe)},'R(A,G|T)',['IMG/Rag_t_orig_', FIG_NAME])   
    
    %plot density
	figure
    print_density({Ragt(is_cohe),Ragt(~is_cohe)},'R(A,G|T)',['IMG/Rag_t_', FIG_NAME])

    %% R(T1,T2|A)
    Niter = 1e4;
    
    peak_multiple_tf = find(sum(peak_tf,2)>1);
    pos_A = peak_multiple_tf( randi(length(peak_multiple_tf),Niter,1) );
    
    pos_T1 = zeros(Niter,1);
    pos_T2 = zeros(Niter,1);
    is_cohe = false(Niter,1);
    Rt1t2 = zeros(Niter,1);
    Rt1t2a = zeros(Niter,1);
    not_zero = false(Niter,1);
    for iter = 1:Niter
        %find two tf
        all_T = find(peak_tf(pos_A(iter),:));
        pos_T = all_T( randperm(length(all_T),2) );
        pos_T1(iter) = pos_T(1);
        pos_T2(iter) = pos_T(2);
        %calculate three dyncorr
        Rt1a = dyncorr( Xtf(pos_T1(iter),:)', Xatac(pos_A(iter),:)', time);
        Rt2a = dyncorr( Xtf(pos_T2(iter),:)', Xatac(pos_A(iter),:)', time);
        Rt1t2(iter) = dyncorr( Xtf(pos_T1(iter),:)', Xtf(pos_T2(iter),:)', time);
        %is coherent?
        is_cohe(iter) = Rt1a*Rt2a*Rt1t2(iter)>0;
        %partial correlation
        Rt1t2a(iter) = pc_fun(Rt1t2(iter),Rt1a,Rt2a);
        %is it not zero?
        not_zero(iter) = R_not_zero(Rt1a) && R_not_zero(Rt2a) && R_not_zero(Rt1t2(iter));
    end
    
    Rt1t2a = Rt1t2a(not_zero);
    is_cohe = is_cohe(not_zero);
    
    %plot count
    figure
    print_count({Rt1t2a(is_cohe),Rt1t2a(~is_cohe)},'R(T1,T2|A)',['IMG/Rt1t2_a_orig_',FIG_NAME])     
    
    %plot density
    figure
    print_density({Rt1t2a(is_cohe),Rt1t2a(~is_cohe)},'R(T1,T2|A)',['IMG/Rt1t2_a_',FIG_NAME]) 
 
    fprintf('\t Data: %.f%% T1,T2->A are coherent\n', 100*sum(is_cohe)/numel(is_cohe))
    pval = cdf(pd,100*sum(is_cohe)/numel(is_cohe),'upper');
	fprintf('\t Enrichment of coherent T1,T2->A with p-value: %e\n', pval)
    
    
	%% R(G1,G2|A)
    Niter = 1e4;
    
    peak_multiple_target = find(sum(peak_target,2)>1);
    pos_A = peak_multiple_target( randi(length(peak_multiple_target),Niter,1) );
    
    pos_G1 = zeros(Niter,1);
    pos_G2 = zeros(Niter,1);
    is_cohe = false(Niter,1);
    Rg1g2 = zeros(Niter,1);
    Rg1g2a = zeros(Niter,1);
    not_zero = false(Niter,1);
    for iter = 1:Niter
        %find two tf
        all_G = find(peak_target(pos_A(iter),:));
        pos_G = all_G( randperm(length(all_G),2) );
        pos_G1(iter) = pos_G(1);
        pos_G2(iter) = pos_G(2);
        %calculate three dyncorr
        Rag1 = dyncorr( Xatac(pos_A(iter),:)', Xtarget(pos_G1(iter),:)', time);
        Rag2 = dyncorr( Xatac(pos_A(iter),:)', Xtarget(pos_G2(iter),:)', time);
        Rg1g2(iter) = dyncorr( Xtarget(pos_G1(iter),:)', Xtarget(pos_G2(iter),:)', time);
        %is coherent?
        is_cohe(iter) = Rag1*Rag2*Rg1g2(iter)>0;
        %partial correlation
        Rg1g2a(iter) = pc_fun(Rg1g2(iter),Rag1,Rag2);
        %not zero
        not_zero(iter) = R_not_zero(Rag1) && R_not_zero(Rag2) && R_not_zero(Rg1g2(iter));
    end
    
    Rg1g2a = Rg1g2a(not_zero);
    is_cohe = is_cohe(not_zero);
     
    %plot count
    figure
    print_count({Rg1g2a(is_cohe),Rg1g2a(~is_cohe)},'R(G1,G2|A)',['IMG/Rg1g2_a_orig_',FIG_NAME]) 
    
    %plot density
    figure
    print_density({Rg1g2a(is_cohe),Rg1g2a(~is_cohe)},'R(G1,G2|A)',['IMG/Rg1g2_a_',FIG_NAME]) 

    fprintf('\t Data: %.f%% G1,G2->A are coherent\n', 100*sum(is_cohe)/numel(is_cohe))
    pval = cdf(pd,100*sum(is_cohe)/numel(is_cohe),'upper');
    fprintf('\t Enrichment of coherent G1,G2->A with p-value: %e\n', pval)

end



