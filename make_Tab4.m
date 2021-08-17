%This script generates the data in Tab.4

clear all
close all
clc

%% datasets
FOLDER_TH1 = 'data\th1\';
FOLDER_P4 = 'data\th1_p4\';
time = [0 0 0, .5 .5 .5, 1 1 1, 2 2 2, 6 6 6, 24 24 24];

%genes list
list = importdata([FOLDER_TH1, 'list.mat']);
list_ = importdata([FOLDER_P4, 'list.mat']);
list_tf = importdata([FOLDER_TH1, 'list_tf.mat']);
list_tf_ = importdata([FOLDER_P4, 'list_tf.mat']);
    
assert(isequal(list,list_) && isequal(list_tf,list_tf_))

%% load
%data, all replicates
Xtarget_th1 = importdata([FOLDER_TH1, 'Xall.mat']); Ntarget = size(Xtarget_th1,1);
Xtf_th1 = importdata([FOLDER_TH1, 'Xtfall.mat']); Ntf = size(Xtf_th1,1);
Xatac_th1 = importdata([FOLDER_TH1, 'Xatacall.mat']); Natac = size(Xatac_th1,1);
peaks = importdata([FOLDER_TH1, 'peaks_trimmed.mat']); assert(size(peaks,1)==Natac)

Xtarget_p4 = importdata([FOLDER_P4, 'Xall.mat']); 
Xtf_p4 = importdata([FOLDER_P4, 'Xtfall.mat']); 
Xatac_p4 = importdata([FOLDER_P4, 'Xatacall.mat']); 

assert( size(Xtarget_p4,1)==Ntarget && size(Xtf_p4,1)==Ntf && size(Xatac_p4,1)==Natac )

Xgene_th1 = [Xtarget_th1; Xtf_th1];
Xgene_p4 = [Xtarget_p4; Xtf_p4];
Rgene = arrayfun(@(n) dyncorr(Xgene_th1(n,:)',Xgene_p4(n,:)', time), 1:Ntarget+Ntf)';
Ratac = arrayfun(@(n) dyncorr(Xatac_th1(n,:)',Xatac_p4(n,:)', time), 1:Natac)';

%% adjacencies
peak_tf_th1 = importdata([FOLDER_TH1, 'adj_peak_tf.mat']);
peak_target_th1 = importdata([FOLDER_TH1, 'adj_peak_target.mat']);
tf_target_th1 = importdata([FOLDER_TH1, 'adj_tf_target.mat']);
    
peak_tf_p4 = importdata([FOLDER_P4, 'adj_peak_tf.mat']);
peak_target_p4 = importdata([FOLDER_P4, 'adj_peak_target.mat']);
tf_target_p4 = importdata([FOLDER_P4, 'adj_tf_target.mat']);
    
fprintf('peak-tf: %d in th1, %d in p4, %d in both\n', nnz(peak_tf_th1), nnz(peak_tf_p4), nnz(peak_tf_th1 & peak_tf_p4))
fprintf('peak-target: %d in th1, %d in p4, %d in both\n', nnz(peak_target_th1), nnz(peak_target_p4), nnz(peak_target_th1 & peak_target_p4))
fprintf('tf-target: %d in th1, %d in p4, %d in both\n', nnz(tf_target_th1), nnz(tf_target_p4), nnz(tf_target_th1 & tf_target_p4))

%% statistical test function
NTP = size(Xtarget_th1,2)/3;
alpha = 0.05;

pc_fun = @(pxy,pxz,pyz) (pxy-pxz*pyz)/sqrt((1 - pxz^2)*(1 - pyz^2));    

t_crit = @(n_control) tinv(.5+alpha/2, NTP-2-n_control);
t = @(pc,n_control) pc.*sqrt( (NTP-2-n_control)./(1-pc.^2) );
close_zero_fun = @(pc,n_control) abs( t(pc,n_control) ) < t_crit(n_control);

pval = @(pc,n_control) 2*abs(0.5 - tcdf(t(pc,n_control),NTP-2-n_control));
    
THS_1 = t_crit(1);
THS_0 = t_crit(0);
R_not_zero = @(vec) abs(t(vec,0))>THS_0;
is_CI = @(vec) abs(t(vec,1))<THS_1;

%% compare FFL 
%load th1
FFL_th1 = importdata([FOLDER_TH1, 'FFL.mat']);
Rtg_th1 = importdata([FOLDER_TH1, 'Rtg.mat']); Rat_th1 = importdata([FOLDER_TH1, 'Rat.mat']); Rag_th1 = importdata([FOLDER_TH1, 'Rag.mat']);
%remove corr close to zero
to_keep = R_not_zero(Rtg_th1) & R_not_zero(Rat_th1) & R_not_zero(Rag_th1); 
Rtg_th1 = Rtg_th1(to_keep); Rat_th1 = Rat_th1(to_keep); Rag_th1 = Rag_th1(to_keep);
FFL_th1 = FFL_th1(to_keep,:);
%calculate prior ratio
PRC_PRIOR = 100*sum(Rtg_th1.*Rat_th1.*Rag_th1>0)/numel(Rag_th1);
cohe_prior = sum(Rtg_th1.*Rat_th1.*Rag_th1>0); tot_prior = numel(Rag_th1);
PRC_PRIOR_VAL_ATG = 100*sum(is_CI(FFL_th1.Rag_t))/size(FFL_th1,1);
val_atg_prior = sum(is_CI(FFL_th1.Rag_t)); 
PRC_PRIOR_VAL_TAG = 100*sum(is_CI(FFL_th1.Rtg_a))/size(FFL_th1,1);
val_tag_prior = sum(is_CI(FFL_th1.Rtg_a));  
%load p4
FFL_p4 = importdata([FOLDER_P4, 'FFL.mat']);
Rtg_p4 = importdata([FOLDER_P4, 'Rtg.mat']); Rat_p4 = importdata([FOLDER_P4, 'Rat.mat']); Rag_p4 = importdata([FOLDER_P4, 'Rag.mat']);
%remove corr close to zero
to_keep = R_not_zero(Rtg_p4) & R_not_zero(Rat_p4) & R_not_zero(Rag_p4); 
Rtg_p4 = Rtg_p4(to_keep); Rat_p4 = Rat_p4(to_keep); Rag_p4 = Rag_p4(to_keep);
FFL_p4 = FFL_p4(to_keep,:);
%overview: overlap between datasets
ATG_th1 = [FFL_th1.A,FFL_th1.T,FFL_th1.G];
ATG_p4 = [FFL_p4.A,FFL_p4.T,FFL_p4.G];
ATG_both = intersect(ATG_th1, ATG_p4,'rows');
fprintf('FFL: %d in th1, %d in p4, %d in common \n', size(ATG_th1,1), size(ATG_p4,1), size(ATG_both,1))
%order in the same way
[~, idx_th1] = ismember(ATG_both,ATG_th1,'rows');
[~, idx_p4] = ismember(ATG_both,ATG_p4,'rows');
ATG_th1 = ATG_th1(idx_th1,:);
ATG_p4 = ATG_p4(idx_p4,:);
assert(isequal(ATG_th1,ATG_p4))
FFL_th1 = FFL_th1(idx_th1,:);
FFL_p4 = FFL_p4(idx_p4,:);
assert(isequal(FFL_th1(:,1:3), FFL_p4(:,1:3)))
%compute coherent on overlap
is_cohe_th1 = Rtg_th1(idx_th1).*Rag_th1(idx_th1).*Rat_th1(idx_th1)>0;
is_cohe_p4 = Rtg_p4(idx_p4).*Rag_p4(idx_p4).*Rat_p4(idx_p4)>0;
fprintf('FFL cohe: %d in th1, %d in p4, %d in common \n', sum(is_cohe_th1), sum(is_cohe_p4), sum(is_cohe_th1 & is_cohe_p4))
ATG_cohe = intersect(ATG_th1(is_cohe_th1,:), ATG_p4(is_cohe_p4,:),'rows');
ATG_inco = intersect(ATG_th1(~is_cohe_th1,:), ATG_p4(~is_cohe_p4,:), 'rows');
%compute validated on overlap
is_val_th1_atg = is_CI(FFL_th1.Rag_t);
is_val_p4_atg = is_CI(FFL_p4.Rag_t);
is_val_th1_tag = is_CI(FFL_th1.Rtg_a);
is_val_p4_tag = is_CI(FFL_p4.Rtg_a);
ATG_conf = ATG_th1(is_val_th1_atg & is_val_p4_atg,: );
TAG_conf = ATG_th1(is_val_th1_tag & is_val_p4_tag,: );
%hypergeometric on coherence
p = hygecdf(sum(is_cohe_th1 & is_cohe_p4),tot_prior,cohe_prior,sum(is_cohe_p4),'upper');
fprintf('FFL th1 coherence: prior=%.2f, on overlap=%.2f, on coherent p4=%.2f, p-val=%e\n', ...
    PRC_PRIOR, 100*sum(is_cohe_th1)/numel(is_cohe_th1), 100*sum(is_cohe_th1 & is_cohe_p4)/sum(is_cohe_p4), p)
%calculate hypergeomtric on confirmed ATG
p = hygecdf(sum(is_val_th1_atg & is_val_p4_atg),tot_prior,val_atg_prior,sum(is_val_p4_atg),'upper');
fprintf('ATG th1 confirmed: prior=%.2f, on overlap=%.2f, on confirmed p4=%.2f, p-val=%e\n', ...
    PRC_PRIOR_VAL_ATG, 100*sum(is_val_th1_atg)/numel(is_val_th1_atg), 100*sum(is_val_th1_atg & is_val_p4_atg)/sum(is_val_p4_atg), p)
%calculate hypergeomtric on confirmed TAG
p = hygecdf(sum(is_val_th1_tag & is_val_p4_tag),tot_prior,val_tag_prior,sum(is_val_p4_tag),'upper');
fprintf('TAG th1 confirmed: prior=%.2f, on overlap=%.2f, on confirmed p4=%.2f, p-val=%e\n', ...
    PRC_PRIOR_VAL_TAG, 100*sum(is_val_th1_tag)/numel(is_val_th1_tag), 100*sum(is_val_th1_tag & is_val_p4_tag)/sum(is_val_p4_tag), p)

%% compare ATT 
%load th1
ATT_th1_tab = importdata([FOLDER_TH1, 'ATT.mat']);
Rtt_th1 = importdata([FOLDER_TH1, 'Rtt.mat']); Rat1_th1 = importdata([FOLDER_TH1, 'Rat1.mat']); Rat2_th1 = importdata([FOLDER_TH1, 'Rat2.mat']);
%keep non zero corr
to_keep = R_not_zero(Rtt_th1) & R_not_zero(Rat1_th1) & R_not_zero(Rat2_th1); 
Rtt_th1 = Rtt_th1(to_keep); Rat1_th1 = Rat1_th1(to_keep); Rat2_th1 = Rat2_th1(to_keep);
ATT_th1_tab = ATT_th1_tab(to_keep,:);
%calculate prior ratio
PRC_PRIOR = 100*sum(Rtt_th1.*Rat1_th1.*Rat2_th1>0)/numel(Rat1_th1);
cohe_prior = sum(Rtt_th1.*Rat1_th1.*Rat2_th1>0); tot_prior = numel(Rtt_th1);
PRC_PRIOR_VAL = 100*sum(is_CI(ATT_th1_tab.Rtt_a))/size(ATT_th1_tab,1);
val_tta_prior = sum(is_CI(ATT_th1_tab.Rtt_a)); 
%load p4
ATT_p4_tab = importdata([FOLDER_P4, 'ATT.mat']);
Rtt_p4 = importdata([FOLDER_P4, 'Rtt.mat']); Rat1_p4 = importdata([FOLDER_P4, 'Rat1.mat']); Rat2_p4 = importdata([FOLDER_P4, 'Rat2.mat']);
%keep nonzero corr
to_keep = R_not_zero(Rtt_p4) & R_not_zero(Rat1_p4) & R_not_zero(Rat2_p4); 
Rtt_p4 = Rtt_p4(to_keep); Rat1_p4 = Rat1_p4(to_keep); Rat2_p4 = Rat2_p4(to_keep);
ATT_p4_tab = ATT_p4_tab(to_keep,:);
%keep common only
ATT_th1 = [ATT_th1_tab.A,ATT_th1_tab.T1,ATT_th1_tab.T2];
ATT_p4 = [ATT_p4_tab.A,ATT_p4_tab.T1,ATT_p4_tab.T2];
ATT_both = intersect(ATT_th1,ATT_p4,'rows');
fprintf('ATT: %d in th1, %d in p4, %d in common \n', size(ATT_th1,1), size(ATT_p4,1), size(ATT_both,1))
%order in the same way
[~, idx_th1] = ismember(ATT_both,ATT_th1,'rows');
[~, idx_p4] = ismember(ATT_both,ATT_p4,'rows');
ATT_th1 = ATT_th1(idx_th1,:);
ATT_p4 = ATT_p4(idx_p4,:);
assert(isequal(ATT_th1,ATT_p4))
ATT_th1_tab = ATT_th1_tab(idx_th1,:);
ATT_p4_tab = ATT_p4_tab(idx_p4,:);
assert(isequal(ATT_th1_tab(:,1:3), ATT_p4_tab(:,1:3)))
%compute coherent 
is_cohe_th1 = Rtt_th1(idx_th1).*Rat1_th1(idx_th1).*Rat2_th1(idx_th1)>0 ;
is_cohe_p4 = Rtt_p4(idx_p4).*Rat1_p4(idx_p4).*Rat2_p4(idx_p4)>0;
fprintf('ATT cohe: %d in th1, %d in p4, %d in common \n', sum(is_cohe_th1), sum(is_cohe_p4), sum(is_cohe_th1 & is_cohe_p4))
ATT_cohe = intersect(ATT_th1(is_cohe_th1,:),ATT_p4(is_cohe_p4,:),'rows');
ATT_inco = intersect(ATT_th1(~is_cohe_th1,:),ATT_p4(~is_cohe_p4,:),'rows');
%compute validated
is_val_th1 = is_CI(ATT_th1_tab.Rtt_a);
is_val_p4 = is_CI(ATT_p4_tab.Rtt_a);
ATT_conf = ATT_th1(is_val_th1 & is_val_p4,:);
%hypergeometric on coherence
p = hygecdf(sum(is_cohe_th1 & is_cohe_p4),tot_prior,cohe_prior,sum(is_cohe_p4),'upper');
fprintf('ATT th1 coherence: prior=%.2f, on overlap=%.2f, on coherent p4=%.2f, p-val=%e\n', ...
    PRC_PRIOR, 100*sum(is_cohe_th1)/numel(is_cohe_th1), 100*sum(is_cohe_th1 & is_cohe_p4)/sum(is_cohe_p4), p)
%calculate hypergeomtric on confirmed
p = hygecdf(sum(is_val_th1 & is_val_p4),tot_prior,val_tta_prior,sum(is_val_p4),'upper');
fprintf('ATT th1 confirmed: prior=%.2f, on overlap=%.2f, on confirmed p4=%.2f, p-val=%e\n', ...
    PRC_PRIOR_VAL, 100*sum(is_val_th1)/numel(is_val_th1), 100*sum(is_val_th1 & is_val_p4)/sum(is_val_p4), p)


%% compare AGG 
%load th1
AGG_th1_tab = importdata([FOLDER_TH1, 'AGG.mat']);
Rgg_th1 = importdata([FOLDER_TH1, 'Rgg.mat']); Rag1_th1 = importdata([FOLDER_TH1, 'Rag1.mat']); Rag2_th1 = importdata([FOLDER_TH1, 'Rag2.mat']);
%keep nonzero only
to_keep = R_not_zero(Rgg_th1) & R_not_zero(Rag1_th1) & R_not_zero(Rag2_th1); 
Rgg_th1 = Rgg_th1(to_keep); Rag1_th1 = Rag1_th1(to_keep); Rag2_th1 = Rag2_th1(to_keep);
AGG_th1_tab = AGG_th1_tab(to_keep,:);
%prior ratio
PRC_PRIOR = 100*sum(Rgg_th1.*Rag1_th1.*Rag2_th1>0)/numel(Rag2_th1);
cohe_prior = sum(Rgg_th1.*Rag1_th1.*Rag2_th1>0); tot_prior = numel(Rgg_th1);
PRC_PRIOR_VAL = 100*sum(is_CI(AGG_th1_tab.Rgg_a))/size(AGG_th1_tab,1);
val_gga_prior = sum(is_CI(AGG_th1_tab.Rgg_a)); 
%load p4
AGG_p4_tab = importdata([FOLDER_P4, 'AGG.mat']);
Rgg_p4 = importdata([FOLDER_P4, 'Rgg.mat']); Rag1_p4 = importdata([FOLDER_P4, 'Rag1.mat']); Rag2_p4 = importdata([FOLDER_P4, 'Rag2.mat']);
%keep nonzero only
to_keep = R_not_zero(Rgg_p4) & R_not_zero(Rag1_p4) & R_not_zero(Rag2_p4); 
Rgg_p4 = Rgg_p4(to_keep); Rag1_p4 = Rag1_p4(to_keep); Rag2_p4 = Rag2_p4(to_keep);
AGG_p4_tab = AGG_p4_tab(to_keep,:);
%keep common only
AGG_th1 = [AGG_th1_tab.A,AGG_th1_tab.G1,AGG_th1_tab.G2];
AGG_p4 = [AGG_p4_tab.A,AGG_p4_tab.G1,AGG_p4_tab.G2];
AGG_both = intersect(AGG_th1, AGG_p4,'rows');
fprintf('AGG: %d in th1, %d in p4, %d in common \n', size(AGG_th1,1), size(AGG_p4,1), size(AGG_both,1))
%order in the same way
[~, idx_th1] = ismember(AGG_both,AGG_th1,'rows');
[~, idx_p4] = ismember(AGG_both,AGG_p4,'rows');
AGG_th1 = AGG_th1(idx_th1,:);
AGG_p4 = AGG_p4(idx_p4,:);
assert(isequal(AGG_th1,AGG_p4))
AGG_th1_tab = AGG_th1_tab(idx_th1,:);
AGG_p4_tab = AGG_p4_tab(idx_p4,:);
assert(isequal(AGG_th1_tab(:,1:3), AGG_p4_tab(:,1:3)))
%compute coherent
is_cohe_th1 = Rgg_th1(idx_th1).*Rag1_th1(idx_th1).*Rag2_th1(idx_th1)>0;
is_cohe_p4 = Rgg_p4(idx_p4).*Rag1_p4(idx_p4).*Rag2_p4(idx_p4)>0;
fprintf('AGG cohe: %d in th1, %d in p4, %d in common \n', sum(is_cohe_th1), sum(is_cohe_p4), sum(is_cohe_th1 & is_cohe_p4))
AGG_cohe = intersect(AGG_th1(is_cohe_th1,:), AGG_p4(is_cohe_p4,:),'rows');
AGG_inco = intersect(AGG_th1(~is_cohe_th1,:), AGG_p4(~is_cohe_p4,:),'rows');
%compute validated
is_val_th1 = is_CI(AGG_th1_tab.Rgg_a);
is_val_p4 = is_CI(AGG_p4_tab.Rgg_a);
AGG_conf = AGG_th1(is_val_th1 & is_val_p4,:);
%calculate hypergeometric on coherent
p = hygecdf(sum(is_cohe_th1 & is_cohe_p4),tot_prior,cohe_prior,sum(is_cohe_p4),'upper');
fprintf('AGG th1 coherence: prior=%.2f, on overlap=%.2f, on coherent p4=%.2f, p-val=%e\n', ...
    PRC_PRIOR, 100*sum(is_cohe_th1)/numel(is_cohe_th1), 100*sum(is_cohe_th1 & is_cohe_p4)/sum(is_cohe_p4), p)
%calculate hypergeomtric on confirmed
p = hygecdf(sum(is_val_th1 & is_val_p4),tot_prior,val_gga_prior,sum(is_val_p4),'upper');
fprintf('AGG th1 confirmed: prior=%.2f, on overlap=%.2f, on confirmed p4=%.2f, p-val=%e\n', ...
    PRC_PRIOR_VAL, 100*sum(is_val_th1)/numel(is_val_th1), 100*sum(is_val_th1 & is_val_p4)/sum(is_val_p4), p)
