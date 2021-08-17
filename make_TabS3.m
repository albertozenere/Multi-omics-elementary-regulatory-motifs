%Script to generate the data of Fig.S3

clear all
close all
clc

addpath('P:\server\P1')

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
    
    peaks = importdata([FOLDER_FILES, 'peaks_trimmed.mat']);
	assert(size(peaks,1)==size(Xatac,1))
    
    NTP = size(Xtarget,2);
    %data, mean of replicates, for comparison purposes
    list = importdata([FOLDER_FILES, 'list.mat']); 
    list_tf = importdata([FOLDER_FILES, 'list_tf.mat']); 
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

	if (contains(FOLDER_FILES, 'th1') && NTP == 18) || contains(FOLDER_FILES, 'th1_p4') || contains(FOLDER_FILES, 'control')
        time = [0 0 0, .5 .5 .5, 1 1 1, 2 2 2, 6 6 6, 24 24 24]; %th1
        ORG = "human";
        
        dist_TSS = importdata('P:\server\P1\atacseq\th1_P4\FinalConsensuspeaksannotatedtotranscripts.mat');

	elseif contains(FOLDER_FILES, 'hiv') || contains(FOLDER_FILES, 'mock')
        time = [2 2 2, 8 8 8, 24 24 24, 48 48 48]; %hiv
        ORG = "human";
     
        dist_TSS = importdata('P:\server\P1\sign_analysis\external_data\Bonneau\atacseq\FinalConsensuspeaksannotatedtotranscripts.mat');
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
    
    %% c.i. 
    %chains
	idx_is_pc = (is_pc(FFL.Rag_t) | is_pc(FFL.Rtg_a)) & ...
        (R_not_zero(Rtg) & R_not_zero(Rag) & R_not_zero(Rat) );
  
    %% filter dist_TSS     
    [i,o] = ismember( [peaks.Start-1 peaks.End], [dist_TSS.start, dist_TSS.end1] , 'rows');
    assert(all(i))
    
    dist_TSS = dist_TSS(o, :);
       
    %% chains
    %selected pairs
    AG_pairs = [FFL.A(idx_is_pc), FFL.G(idx_is_pc)];
   
    AG_pairs = unique(AG_pairs, 'rows');
    Niter = size(AG_pairs,1);
    
    tit = 'selected';
    
    distAG = nan(Niter,1);
    for n = 1:Niter
        posA = AG_pairs(n,1); 
        posG = AG_pairs(n,2); 
        g = list(posG);
       
        %check
        cond = contains(dist_TSS.geneIds(posA), g) ...
           | contains(dist_TSS.flank_geneIds(posA), g);  
        assert(cond)
       
        if contains(dist_TSS.geneIds(posA), g)         
            str = string(dist_TSS.geneIds(posA));
            str_vec = strsplit(str,";");
                     
            dist = string(dist_TSS.distancetoTSS(posA));
            dist_vec = str2double((strsplit(dist,";")));
           
            assert(length(str_vec)==length(dist_vec))
            
            i = find(ismember(str_vec,g),1);
            distAG(n) = dist_vec(i);
            
        elseif contains(dist_TSS.flank_geneIds(posA), g)
            str = string(dist_TSS.flank_geneIds(posA));
            str_vec = strsplit(str,";");
                     
            dist = string(dist_TSS.flank_distancetoTSS(posA));
            dist_vec = str2double((strsplit(dist,";")));
           
            assert(length(str_vec)==length(dist_vec))
            
            i = find(ismember(str_vec,g),1);
            distAG(n) = dist_vec(i);
        end   
    end
    distAG(distAG>1e4) = []; 
    assert(~any(isnan(distAG)))
    
    fprintf('Selected ATG/TAG: Percentage of TSS inside the peak: %.2f%%\n', 100*sum(distAG==0)/numel(distAG))
    fprintf('Selected ATG/TAG: Mean distance between TSS and peak: %.f bp \n', mean(abs(distAG)))
    
	%% distance to TSS
    %remaining pairs
    [I,J] = find(peak_target);
    AG_rand = [I, J];
	AG_rand(ismember(AG_rand,AG_pairs,'rows'),:) = [];
    AG_pairs = AG_rand;
    
    Niter = size(AG_rand,1);
    tit = 'rand';
    
    distAG = nan(Niter,1);
    for n = 1:Niter
        posA = AG_pairs(n,1); 
        posG = AG_pairs(n,2); 
        g = list(posG);
       
        %check
        cond = contains(dist_TSS.geneIds(posA), g) ...
           | contains(dist_TSS.flank_geneIds(posA), g);  
        assert(cond)
       
        if contains(dist_TSS.geneIds(posA), g)         
            str = string(dist_TSS.geneIds(posA));
            str_vec = strsplit(str,";");
                     
            dist = string(dist_TSS.distancetoTSS(posA));
            dist_vec = str2double((strsplit(dist,";")));
           
            assert(length(str_vec)==length(dist_vec))
            
            i = find(ismember(str_vec,g),1);
            distAG(n) = dist_vec(i);
            
        elseif contains(dist_TSS.flank_geneIds(posA), g)
            str = string(dist_TSS.flank_geneIds(posA));
            str_vec = strsplit(str,";");
                     
            dist = string(dist_TSS.flank_distancetoTSS(posA));
            dist_vec = str2double((strsplit(dist,";")));
           
            assert(length(str_vec)==length(dist_vec))
            
            i = find(ismember(str_vec,g),1);
            distAG(n) = dist_vec(i);
        end   
    end
    distAG(distAG>1e4) = [];
    assert(~any(isnan(distAG)))
    
       
    fprintf('Random AG pairs: Percentage of TSS inside the peak: %.2f%%\n', 100*sum(distAG==0)/numel(distAG))
    fprintf('Random AG pairs: Mean distance between TSS and peak: %.f bp \n', mean(abs(distAG)))
  
end