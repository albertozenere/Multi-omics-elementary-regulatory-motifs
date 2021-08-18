%This script generates Fig.2

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

sz = [4 3];
all_fig = figure('Position',[50 50 500 300]);
xbin = -1:.02:1;

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
    
	%% statistical test function
    ntp = NTP/3;
    alpha = 0.05;

    pc_fun = @(pxy,pxz,pyz) (pxy-pxz*pyz)/sqrt((1 - pxz^2)*(1 - pyz^2));    

    t_crit = @(n_control) tinv(.5+alpha/2, ntp-2-n_control);
    t = @(pc,n_control) pc.*sqrt( (ntp-2-n_control)./(1-pc.^2) );
    close_zero_fun = @(pc,n_control) abs( t(pc,n_control) ) < t_crit(n_control);

    pval = @(pc,n_control) 2*abs(0.5 - tcdf(t(pc,n_control),ntp-2-n_control));
    
    THS_1 = t_crit(1);
    THS_0 = t_crit(0);
    R_not_zero = @(vec) abs(t(vec,0))>THS_0;
    
    %% functions
    min_fun = @(a,b,c) min(abs([a,b,c]),[],2);
    max_fun = @(a,b,c) max(abs([a,b,c]),[],2);
    geommean_fun = @(a,b,c) geomean(abs([a,b,c]),2);
    
    %% chains and forks
	R1 = [Rtg; Rtt; Rgg];
    R2 = [Rat; Rat1; Rag1];
    R3 = [Rag; Rat2; Rag2];     
    is_cohe = [FFL.is_cohe; ATT.is_cohe; AGG.is_cohe];
    %remove close to zero
    not_zero = abs(R1)>THS_0 & abs(R2)>THS_0 & abs(R3)>THS_0;
    R1 = R1(not_zero); R2 = R2(not_zero); R3 = R3(not_zero); is_cohe = is_cohe(not_zero);
    %calculate functions
    qFFL1 = min_fun(R1,R2,R3);
    qFFL2 = geommean_fun(R1,R2,R3);
    qFFL3 = max_fun(R1,R2,R3);
    %plot
    figure(all_fig)
    subplot(sz(1),sz(2),1+3*(n_data-1))
    hold on
    print_density({qFFL1(is_cohe),qFFL1(~is_cohe)}, [])
    xlim([0 1])
    subplot(sz(1),sz(2),3+3*(n_data-1))
    hold on
    print_density({qFFL2(is_cohe),qFFL2(~is_cohe)}, [])
    xlim([0 1])
    subplot(sz(1),sz(2),2+3*(n_data-1))
    hold on
    print_density({qFFL3(is_cohe),qFFL3(~is_cohe)}, [])
    xlim([0 1])
     
    [h,p] = kstest2(qFFL1(is_cohe), qFFL1(~is_cohe));
    fprintf('\t Are the distributions of the minima of coherent and incoherent cycles similar? P-value:%e\n', p) 
    fprintf('\t\t min mean cohe:%.2f; min mean inco:%.2f\n', mean(qFFL1(is_cohe)), mean(qFFL1(~is_cohe)))    
    [h,p] = kstest2(qFFL2(is_cohe), qFFL2(~is_cohe));
    fprintf('\t Are the distributions of the geometric means of coherent and incoherent cycles similar?% P-value:e\n', p)  
    fprintf('\t\t geomean mean cohe:%.2f; geomean mean inco:%.2f\n', mean(qFFL2(is_cohe)), mean(qFFL2(~is_cohe)))   
    [h,p] = kstest2(qFFL3(is_cohe), qFFL3(~is_cohe));
    fprintf('\t Are the distributions of the maxima of coherent and incoherent cycles similar?% P-value:e\n', p)
    fprintf('\t\t max mean cohe:%.2f; max mean inco:%.2f\n', mean(qFFL3(is_cohe)), mean(qFFL3(~is_cohe))) 
    
end


print('IMG/all_cohe_prop','-depsc') %depsc

