%This script generates the data in Tab.S2

clear all
close all
clc

%% datasets
FOLDER_HIV = 'data\Bonneau\hiv\';
FOLDER_MOCK = 'data\Bonneau\mock\';
time = [2 2 2, 8 8 8, 24 24 24, 48 48 48];

%genes list
list = importdata([FOLDER_HIV, 'list.mat']);
list_ = importdata([FOLDER_MOCK, 'list.mat']);
list_tf = importdata([FOLDER_HIV, 'list_tf.mat']);
list_tf_ = importdata([FOLDER_MOCK, 'list_tf.mat']);
    
assert(isequal(list,list_) && isequal(list_tf,list_tf_))

%% load
%data, all replicates
Xtarget_hiv = importdata([FOLDER_HIV, 'Xall.mat']); Ntarget = size(Xtarget_hiv,1);
Xtf_hiv = importdata([FOLDER_HIV, 'Xtfall.mat']); Ntf = size(Xtf_hiv,1);
Xatac_hiv = importdata([FOLDER_HIV, 'Xatacall.mat']); Natac = size(Xatac_hiv,1);
peaks = importdata([FOLDER_HIV, 'peaks_trimmed.mat']); assert(size(peaks,1)==Natac)

Xtarget_mock = importdata([FOLDER_MOCK, 'Xall.mat']); 
Xtf_mock = importdata([FOLDER_MOCK, 'Xtfall.mat']); 
Xatac_mock = importdata([FOLDER_MOCK, 'Xatacall.mat']); 

assert( size(Xtarget_mock,1)==Ntarget && size(Xtf_mock,1)==Ntf && size(Xatac_mock,1)==Natac )

Xgene_hiv = [Xtarget_hiv; Xtf_hiv];
Xgene_mock = [Xtarget_mock; Xtf_mock];
Rgene = arrayfun(@(n) dyncorr(Xgene_hiv(n,:)',Xgene_mock(n,:)', time), 1:Ntarget+Ntf)';
Ratac = arrayfun(@(n) dyncorr(Xatac_hiv(n,:)',Xatac_mock(n,:)', time), 1:Natac)';

%% adjacencies
peak_tf_hiv = importdata([FOLDER_HIV, 'adj_peak_tf.mat']);
peak_target_hiv = importdata([FOLDER_HIV, 'adj_peak_target.mat']);
tf_target_hiv = importdata([FOLDER_HIV, 'adj_tf_target.mat']);
    
peak_tf_mock = importdata([FOLDER_MOCK, 'adj_peak_tf.mat']);
peak_target_mock = importdata([FOLDER_MOCK, 'adj_peak_target.mat']);
tf_target_mock = importdata([FOLDER_MOCK, 'adj_tf_target.mat']);
    
fprintf('peak-tf: %d in hiv, %d in mock, %d in both\n', nnz(peak_tf_hiv), nnz(peak_tf_mock), nnz(peak_tf_hiv & peak_tf_mock))
fprintf('peak-target: %d in hiv, %d in mock, %d in both\n', nnz(peak_target_hiv), nnz(peak_target_mock), nnz(peak_target_hiv & peak_target_mock))
fprintf('tf-target: %d in hiv, %d in mock, %d in both\n', nnz(tf_target_hiv), nnz(tf_target_mock), nnz(tf_target_hiv & tf_target_mock))

%% statistical test function
NTP = size(Xtarget_hiv,2)/3;
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
%load hiv
FFL_hiv = importdata([FOLDER_HIV, 'FFL.mat']);
Rtg_hiv = importdata([FOLDER_HIV, 'Rtg.mat']); Rat_hiv = importdata([FOLDER_HIV, 'Rat.mat']); Rag_hiv = importdata([FOLDER_HIV, 'Rag.mat']);
%remove corr close to zero
to_keep = R_not_zero(Rtg_hiv) & R_not_zero(Rat_hiv) & R_not_zero(Rag_hiv); 
Rtg_hiv = Rtg_hiv(to_keep); Rat_hiv = Rat_hiv(to_keep); Rag_hiv = Rag_hiv(to_keep);
FFL_hiv = FFL_hiv(to_keep,:);
%calculate prior ratio
PRC_PRIOR = 100*sum(Rtg_hiv.*Rat_hiv.*Rag_hiv>0)/numel(Rag_hiv);
cohe_prior = sum(Rtg_hiv.*Rat_hiv.*Rag_hiv>0); tot_prior = numel(Rag_hiv);
PRC_PRIOR_VAL_ATG = 100*sum(is_CI(FFL_hiv.Rag_t))/size(FFL_hiv,1);
val_atg_prior = sum(is_CI(FFL_hiv.Rag_t)); 
PRC_PRIOR_VAL_TAG = 100*sum(is_CI(FFL_hiv.Rtg_a))/size(FFL_hiv,1);
val_tag_prior = sum(is_CI(FFL_hiv.Rtg_a));  
%load mock
FFL_mock = importdata([FOLDER_MOCK, 'FFL.mat']);
Rtg_mock = importdata([FOLDER_MOCK, 'Rtg.mat']); Rat_mock = importdata([FOLDER_MOCK, 'Rat.mat']); Rag_mock = importdata([FOLDER_MOCK, 'Rag.mat']);
%remove corr close to zero
to_keep = R_not_zero(Rtg_mock) & R_not_zero(Rat_mock) & R_not_zero(Rag_mock); 
Rtg_mock = Rtg_mock(to_keep); Rat_mock = Rat_mock(to_keep); Rag_mock = Rag_mock(to_keep);
FFL_mock = FFL_mock(to_keep,:);
%overview: overlap between datasets
ATG_hiv = [FFL_hiv.A,FFL_hiv.T,FFL_hiv.G];
ATG_mock = [FFL_mock.A,FFL_mock.T,FFL_mock.G];
ATG_both = intersect(ATG_hiv, ATG_mock,'rows');
fprintf('FFL: %d in hiv, %d in mock, %d in common \n', size(ATG_hiv,1), size(ATG_mock,1), size(ATG_both,1))
%order in the same way
[~, idx_hiv] = ismember(ATG_both,ATG_hiv,'rows');
[~, idx_mock] = ismember(ATG_both,ATG_mock,'rows');
ATG_hiv = ATG_hiv(idx_hiv,:);
ATG_mock = ATG_mock(idx_mock,:);
assert(isequal(ATG_hiv,ATG_mock))
FFL_hiv = FFL_hiv(idx_hiv,:);
FFL_mock = FFL_mock(idx_mock,:);
assert(isequal(FFL_hiv(:,1:3), FFL_mock(:,1:3)))
%compute coherent on overlap
is_cohe_hiv = Rtg_hiv(idx_hiv).*Rag_hiv(idx_hiv).*Rat_hiv(idx_hiv)>0;
is_cohe_mock = Rtg_mock(idx_mock).*Rag_mock(idx_mock).*Rat_mock(idx_mock)>0;
fprintf('FFL cohe: %d in hiv, %d in mock, %d in common \n', sum(is_cohe_hiv), sum(is_cohe_mock), sum(is_cohe_hiv & is_cohe_mock))
ATG_cohe = intersect(ATG_hiv(is_cohe_hiv,:), ATG_mock(is_cohe_mock,:),'rows');
ATG_inco = intersect(ATG_hiv(~is_cohe_hiv,:), ATG_mock(~is_cohe_mock,:), 'rows');
%compute validated on overlap
is_val_hiv_atg = is_CI(FFL_hiv.Rag_t);
is_val_mock_atg = is_CI(FFL_mock.Rag_t);
is_val_hiv_tag = is_CI(FFL_hiv.Rtg_a);
is_val_mock_tag = is_CI(FFL_mock.Rtg_a);
ATG_conf = ATG_hiv(is_val_hiv_atg & is_val_mock_atg,: );
TAG_conf = ATG_hiv(is_val_hiv_tag & is_val_mock_tag,: );
%hypergeometric on coherence
p = hygecdf(sum(is_cohe_hiv & is_cohe_mock),tot_prior,cohe_prior,sum(is_cohe_mock),'upper');
fprintf('FFL hiv coherence: prior=%.2f, on overlap=%.2f, on coherent mock=%.2f, p-val=%e\n', ...
    PRC_PRIOR, 100*sum(is_cohe_hiv)/numel(is_cohe_hiv), 100*sum(is_cohe_hiv & is_cohe_mock)/sum(is_cohe_mock), p)
%calculate hypergeomtric on confirmed ATG
p = hygecdf(sum(is_val_hiv_atg & is_val_mock_atg),tot_prior,val_atg_prior,sum(is_val_mock_atg),'upper');
fprintf('ATG hiv confirmed: prior=%.2f, on overlap=%.2f, on confirmed mock=%.2f, p-val=%e\n', ...
    PRC_PRIOR_VAL_ATG, 100*sum(is_val_hiv_atg)/numel(is_val_hiv_atg), 100*sum(is_val_hiv_atg & is_val_mock_atg)/sum(is_val_mock_atg), p)
%calculate hypergeomtric on confirmed TAG
p = hygecdf(sum(is_val_hiv_tag & is_val_mock_tag),tot_prior,val_tag_prior,sum(is_val_mock_tag),'upper');
fprintf('TAG hiv confirmed: prior=%.2f, on overlap=%.2f, on confirmed mock=%.2f, p-val=%e\n', ...
    PRC_PRIOR_VAL_TAG, 100*sum(is_val_hiv_tag)/numel(is_val_hiv_tag), 100*sum(is_val_hiv_tag & is_val_mock_tag)/sum(is_val_mock_tag), p)

%% compare ATT 
%load hiv
ATT_hiv_tab = importdata([FOLDER_HIV, 'ATT.mat']);
Rtt_hiv = importdata([FOLDER_HIV, 'Rtt.mat']); Rat1_hiv = importdata([FOLDER_HIV, 'Rat1.mat']); Rat2_hiv = importdata([FOLDER_HIV, 'Rat2.mat']);
%keep non zero corr
to_keep = R_not_zero(Rtt_hiv) & R_not_zero(Rat1_hiv) & R_not_zero(Rat2_hiv); 
Rtt_hiv = Rtt_hiv(to_keep); Rat1_hiv = Rat1_hiv(to_keep); Rat2_hiv = Rat2_hiv(to_keep);
ATT_hiv_tab = ATT_hiv_tab(to_keep,:);
%calculate prior ratio
PRC_PRIOR = 100*sum(Rtt_hiv.*Rat1_hiv.*Rat2_hiv>0)/numel(Rat1_hiv);
cohe_prior = sum(Rtt_hiv.*Rat1_hiv.*Rat2_hiv>0); tot_prior = numel(Rtt_hiv);
PRC_PRIOR_VAL = 100*sum(is_CI(ATT_hiv_tab.Rtt_a))/size(ATT_hiv_tab,1);
val_tta_prior = sum(is_CI(ATT_hiv_tab.Rtt_a)); 
%load mock
ATT_mock_tab = importdata([FOLDER_MOCK, 'ATT.mat']);
Rtt_mock = importdata([FOLDER_MOCK, 'Rtt.mat']); Rat1_mock = importdata([FOLDER_MOCK, 'Rat1.mat']); Rat2_mock = importdata([FOLDER_MOCK, 'Rat2.mat']);
%keep nonzero corr
to_keep = R_not_zero(Rtt_mock) & R_not_zero(Rat1_mock) & R_not_zero(Rat2_mock); 
Rtt_mock = Rtt_mock(to_keep); Rat1_mock = Rat1_mock(to_keep); Rat2_mock = Rat2_mock(to_keep);
ATT_mock_tab = ATT_mock_tab(to_keep,:);
%keep common only
ATT_hiv = [ATT_hiv_tab.A,ATT_hiv_tab.T1,ATT_hiv_tab.T2];
ATT_mock = [ATT_mock_tab.A,ATT_mock_tab.T1,ATT_mock_tab.T2];
ATT_both = intersect(ATT_hiv,ATT_mock,'rows');
fprintf('ATT: %d in hiv, %d in mock, %d in common \n', size(ATT_hiv,1), size(ATT_mock,1), size(ATT_both,1))
%order in the same way
[~, idx_hiv] = ismember(ATT_both,ATT_hiv,'rows');
[~, idx_mock] = ismember(ATT_both,ATT_mock,'rows');
ATT_hiv = ATT_hiv(idx_hiv,:);
ATT_mock = ATT_mock(idx_mock,:);
assert(isequal(ATT_hiv,ATT_mock))
ATT_hiv_tab = ATT_hiv_tab(idx_hiv,:);
ATT_mock_tab = ATT_mock_tab(idx_mock,:);
assert(isequal(ATT_hiv_tab(:,1:3), ATT_mock_tab(:,1:3)))
%compute coherent 
is_cohe_hiv = Rtt_hiv(idx_hiv).*Rat1_hiv(idx_hiv).*Rat2_hiv(idx_hiv)>0 ;
is_cohe_mock = Rtt_mock(idx_mock).*Rat1_mock(idx_mock).*Rat2_mock(idx_mock)>0;
fprintf('ATT cohe: %d in hiv, %d in mock, %d in common \n', sum(is_cohe_hiv), sum(is_cohe_mock), sum(is_cohe_hiv & is_cohe_mock))
ATT_cohe = intersect(ATT_hiv(is_cohe_hiv,:),ATT_mock(is_cohe_mock,:),'rows');
ATT_inco = intersect(ATT_hiv(~is_cohe_hiv,:),ATT_mock(~is_cohe_mock,:),'rows');
%compute validated
is_val_hiv = is_CI(ATT_hiv_tab.Rtt_a);
is_val_mock = is_CI(ATT_mock_tab.Rtt_a);
ATT_conf = ATT_hiv(is_val_hiv & is_val_mock,:);
%hypergeometric on coherence
p = hygecdf(sum(is_cohe_hiv & is_cohe_mock),tot_prior,cohe_prior,sum(is_cohe_mock),'upper');
fprintf('ATT hiv coherence: prior=%.2f, on overlap=%.2f, on coherent mock=%.2f, p-val=%e\n', ...
    PRC_PRIOR, 100*sum(is_cohe_hiv)/numel(is_cohe_hiv), 100*sum(is_cohe_hiv & is_cohe_mock)/sum(is_cohe_mock), p)
%calculate hypergeomtric on confirmed
p = hygecdf(sum(is_val_hiv & is_val_mock),tot_prior,val_tta_prior,sum(is_val_mock),'upper');
fprintf('ATT hiv confirmed: prior=%.2f, on overlap=%.2f, on confirmed mock=%.2f, p-val=%e\n', ...
    PRC_PRIOR_VAL, 100*sum(is_val_hiv)/numel(is_val_hiv), 100*sum(is_val_hiv & is_val_mock)/sum(is_val_mock), p)


%% compare AGG 
%load hiv
AGG_hiv_tab = importdata([FOLDER_HIV, 'AGG.mat']);
Rgg_hiv = importdata([FOLDER_HIV, 'Rgg.mat']); Rag1_hiv = importdata([FOLDER_HIV, 'Rag1.mat']); Rag2_hiv = importdata([FOLDER_HIV, 'Rag2.mat']);
%keep nonzero only
to_keep = R_not_zero(Rgg_hiv) & R_not_zero(Rag1_hiv) & R_not_zero(Rag2_hiv); 
Rgg_hiv = Rgg_hiv(to_keep); Rag1_hiv = Rag1_hiv(to_keep); Rag2_hiv = Rag2_hiv(to_keep);
AGG_hiv_tab = AGG_hiv_tab(to_keep,:);
%prior ratio
PRC_PRIOR = 100*sum(Rgg_hiv.*Rag1_hiv.*Rag2_hiv>0)/numel(Rag2_hiv);
cohe_prior = sum(Rgg_hiv.*Rag1_hiv.*Rag2_hiv>0); tot_prior = numel(Rgg_hiv);
PRC_PRIOR_VAL = 100*sum(is_CI(AGG_hiv_tab.Rgg_a))/size(AGG_hiv_tab,1);
val_gga_prior = sum(is_CI(AGG_hiv_tab.Rgg_a)); 
%load mock
AGG_mock_tab = importdata([FOLDER_MOCK, 'AGG.mat']);
Rgg_mock = importdata([FOLDER_MOCK, 'Rgg.mat']); Rag1_mock = importdata([FOLDER_MOCK, 'Rag1.mat']); Rag2_mock = importdata([FOLDER_MOCK, 'Rag2.mat']);
%keep nonzero only
to_keep = R_not_zero(Rgg_mock) & R_not_zero(Rag1_mock) & R_not_zero(Rag2_mock); 
Rgg_mock = Rgg_mock(to_keep); Rag1_mock = Rag1_mock(to_keep); Rag2_mock = Rag2_mock(to_keep);
AGG_mock_tab = AGG_mock_tab(to_keep,:);
%keep common only
AGG_hiv = [AGG_hiv_tab.A,AGG_hiv_tab.G1,AGG_hiv_tab.G2];
AGG_mock = [AGG_mock_tab.A,AGG_mock_tab.G1,AGG_mock_tab.G2];
AGG_both = intersect(AGG_hiv, AGG_mock,'rows');
fprintf('AGG: %d in hiv, %d in mock, %d in common \n', size(AGG_hiv,1), size(AGG_mock,1), size(AGG_both,1))
%order in the same way
[~, idx_hiv] = ismember(AGG_both,AGG_hiv,'rows');
[~, idx_mock] = ismember(AGG_both,AGG_mock,'rows');
AGG_hiv = AGG_hiv(idx_hiv,:);
AGG_mock = AGG_mock(idx_mock,:);
assert(isequal(AGG_hiv,AGG_mock))
AGG_hiv_tab = AGG_hiv_tab(idx_hiv,:);
AGG_mock_tab = AGG_mock_tab(idx_mock,:);
assert(isequal(AGG_hiv_tab(:,1:3), AGG_mock_tab(:,1:3)))
%compute coherent
is_cohe_hiv = Rgg_hiv(idx_hiv).*Rag1_hiv(idx_hiv).*Rag2_hiv(idx_hiv)>0;
is_cohe_mock = Rgg_mock(idx_mock).*Rag1_mock(idx_mock).*Rag2_mock(idx_mock)>0;
fprintf('AGG cohe: %d in hiv, %d in mock, %d in common \n', sum(is_cohe_hiv), sum(is_cohe_mock), sum(is_cohe_hiv & is_cohe_mock))
AGG_cohe = intersect(AGG_hiv(is_cohe_hiv,:), AGG_mock(is_cohe_mock,:),'rows');
AGG_inco = intersect(AGG_hiv(~is_cohe_hiv,:), AGG_mock(~is_cohe_mock,:),'rows');
%compute validated
is_val_hiv = is_CI(AGG_hiv_tab.Rgg_a);
is_val_mock = is_CI(AGG_mock_tab.Rgg_a);
AGG_conf = AGG_hiv(is_val_hiv & is_val_mock,:);
%calculate hypergeometric on coherent
p = hygecdf(sum(is_cohe_hiv & is_cohe_mock),tot_prior,cohe_prior,sum(is_cohe_mock),'upper');
fprintf('AGG hiv coherence: prior=%.2f, on overlap=%.2f, on coherent mock=%.2f, p-val=%e\n', ...
    PRC_PRIOR, 100*sum(is_cohe_hiv)/numel(is_cohe_hiv), 100*sum(is_cohe_hiv & is_cohe_mock)/sum(is_cohe_mock), p)
%calculate hypergeomtric on confirmed
p = hygecdf(sum(is_val_hiv & is_val_mock),tot_prior,val_gga_prior,sum(is_val_mock),'upper');
fprintf('AGG hiv confirmed: prior=%.2f, on overlap=%.2f, on confirmed mock=%.2f, p-val=%e\n', ...
    PRC_PRIOR_VAL, 100*sum(is_val_hiv)/numel(is_val_hiv), 100*sum(is_val_hiv & is_val_mock)/sum(is_val_mock), p)
