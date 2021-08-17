%This script generates the data in Tab.3-S1 and oart of Tab.2

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
    NTP = size(Xtarget,2)/3;
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
    alpha = 0.05;

    pc_fun = @(pxy,pxz,pyz) (pxy-pxz*pyz)/sqrt((1 - pxz^2)*(1 - pyz^2));    

    t_crit = @(n_control) tinv(.5+alpha/2, NTP-2-n_control);
    t = @(pc,n_control) pc.*sqrt( (NTP-2-n_control)./(1-pc.^2) );

    THS_0 = t_crit(0);
    THS_1 = t_crit(1);
    R_not_zero = @(vec) abs(t(vec,0))>THS_0;
    is_CI = @(vec) abs(t(vec,1))<THS_1;

    %% FFL
    to_keep = R_not_zero(Rag) & R_not_zero(Rtg) & R_not_zero(Rat);
    FFL = FFL(to_keep,:);
    Rtg = Rtg(to_keep); Rat = Rat(to_keep); Rag = Rag(to_keep);

    n_FFL = size(FFL,1);

    %print
    val_agt = is_CI(FFL.Rag_t);
    fprintf('\t A->T->G: TP:%d, FP:%d, FN:%d, TN:%d\n', ...
        sum(val_agt & FFL.is_cohe), sum(~val_agt & FFL.is_cohe), ...
        sum(val_agt & ~FFL.is_cohe), sum(~val_agt & ~FFL.is_cohe))
    fprintf('\t\t Percentage of validated: %.3f%%\n', 100*sum(val_agt)/numel(val_agt))
    fprintf('\t\t Percentage of validated among coherent: %.3f%%\n', 100*sum(val_agt & FFL.is_cohe)/sum(FFL.is_cohe))

    val_tga = abs(t(FFL.Rtg_a,1))<t_crit(1);
    fprintf('\t T->A->G: TP:%d, FP:%d, FN:%d, TN:%d\n', ...
        sum(val_tga & FFL.is_cohe), sum(~val_tga & FFL.is_cohe), ...
        sum(val_tga & ~FFL.is_cohe), sum(~val_tga & ~FFL.is_cohe))
    fprintf('\t\t Percentage of validated: %.3f%%\n', 100*sum(val_tga)/numel(val_tga))
    fprintf('\t\t Percentage of validated among coherent: %.3f%%\n', 100*sum(val_tga & FFL.is_cohe)/sum(FFL.is_cohe))

    F = [sum(val_agt & val_tga),sum(~val_agt & val_tga);
        sum(val_agt & ~val_tga), sum(~val_agt & ~val_tga)];
    p = hygecdf( F(1,1), numel(val_agt), sum(val_agt), F(1,1)+F(1,2));
    fprintf('Percentage of selected ATG=%.4f, percentage of selected ATG when TAG is also selected=%.4f\n', ...
        100*sum(val_agt)/numel(val_agt), 100*F(1,1)/(F(1,1)+F(1,2)))
    fprintf('\t Selected A->T->G vs selected T->A->G: ATG-TAG:%d, !ATG-TAG:%d, ATG-!TAG:%d, !ATG-!TAG:%d\n', ...
        sum(val_agt & val_tga), sum(~val_agt & val_tga), ...
        sum(val_agt & ~val_tga), sum(~val_agt & ~val_tga))
    fprintf('\t Are ATG and TAG disjoint? %d, p-val=%e\n', p<0.05,p)
    
    
    %% ATT
    to_keep = R_not_zero(Rat1) & R_not_zero(Rat2) & R_not_zero(Rtt);
    ATT = ATT(to_keep,:);
    Rtt = Rtt(to_keep); Rat1 = Rat1(to_keep); Rat2 = Rat2(to_keep);

    n_ATT = size(ATT,1);

    %print
    val_tta = is_CI(ATT.Rtt_a);
    fprintf('\t A->T1,T2: TP:%d, FP:%d, FN:%d, TN:%d\n', ...
        sum(val_tta & ATT.is_cohe), sum(~val_tta & ATT.is_cohe), ...
        sum(val_tta & ~ATT.is_cohe), sum(~val_tta & ~ATT.is_cohe))
    fprintf('\t\t Percentage of validated: %.3f%%\n', 100*sum(val_tta)/numel(val_tta))
    fprintf('\t\t Percentage of validated among coherent: %.3f%%\n', 100*sum(val_tta & ATT.is_cohe)/sum(ATT.is_cohe))

	%% AGG
    to_keep = R_not_zero(Rag1) & R_not_zero(Rag2) & R_not_zero(Rgg);
    AGG = AGG(to_keep,:);
    Rgg = Rgg(to_keep); Rag1 = Rag1(to_keep); Rag2 = Rag2(to_keep);

    n_AGG = size(AGG,1);

    %print
    val_gga = is_CI(AGG.Rgg_a);
    fprintf('\t A->G1,G2: TP:%d, FP:%d, FN:%d, TN:%d\n', ...
        sum(val_gga & AGG.is_cohe), sum(~val_gga & AGG.is_cohe), ...
        sum(val_gga & ~AGG.is_cohe), sum(~val_gga & ~AGG.is_cohe))
    fprintf('\t\t Percentage of validated: %.3f%%\n', 100*sum(val_gga)/numel(val_gga))
    fprintf('\t\t Percentage of validated among coherent: %.3f%%\n', 100*sum(val_gga & AGG.is_cohe)/sum(AGG.is_cohe))


    
end
    
  
    
    
    
    
    