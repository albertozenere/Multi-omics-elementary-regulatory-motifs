%This is the first script that needs to be run, in order to create
%necessary .mat files

clear all
close all
clc

%% datasets
ALL_FOLDER_FILES = {'data\th1\', 'data\th1_p4\', 'data\Bonneau\hiv\', 'data\Bonneau\mock\'};
n_dataset = length(ALL_FOLDER_FILES);

%% loop
for n_data = 1:n_dataset

    FOLDER_FILES = ALL_FOLDER_FILES{n_data};
    fprintf('Dataset: %s\n', FOLDER_FILES)
    
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

    list = importdata([FOLDER_FILES, 'list.mat']);
    list_tf = importdata([FOLDER_FILES, 'list_tf.mat']);

    if contains(FOLDER_FILES, 'th1')
        time = [0 0 0, .5 .5 .5, 1 1 1, 2 2 2, 6 6 6, 24 24 24]; %th1
        ORG = "human";
    elseif contains(FOLDER_FILES, 'hiv') || contains(FOLDER_FILES, 'mock')
        time = [2 2 2, 8 8 8, 24 24 24, 48 48 48]; %hiv
        ORG = "human";
    else
        error()
    end

    %% indeces of FFLs
    fprintf('\t Started computing Chain files..\n')
    %calculate number of cycles
    n_FFL = 0;
    for n = 1:Natac
       n_tf = sum( peak_tf(n,:) );
       n_target = sum( peak_target(n,:) );

       n_FFL = n_FFL + n_tf*n_target;
    end
    
    %indeces of all FFLs
    T = zeros(n_FFL,1);
    A = zeros(n_FFL,1);
    G = zeros(n_FFL,1);

    [Itf, Jtarget] = find(tf_target);
    Nint = length(Itf);
    c = 1;
    for n = 1:Nint
        Kpeak = find(peak_tf(:,Itf(n)) & peak_target(:,Jtarget(n)));
        n_peak = length(Kpeak);

        T(c:c+n_peak-1) = Itf(n);
        A(c:c+n_peak-1) = Kpeak;
        G(c:c+n_peak-1) = Jtarget(n);

        c = c + n_peak;
    end
    assert(c==n_FFL+1)

    % calculate PC   
    pc_fun = @(pxy,pxz,pyz) (pxy-pxz*pyz)/sqrt((1 - pxz^2)*(1 - pyz^2)); %parcorr function, R(x,y|z)   
    
    is_cohe = false(n_FFL,1);
    Rtg_a = zeros(n_FFL,1); 
    Rag_t = zeros(n_FFL,1);
    Rtg = zeros(n_FFL,1);
    Rag = zeros(n_FFL,1);
    Rat = zeros(n_FFL,1);
    
    pool = parpool([2,10]);
    parfor n = 1:n_FFL
        Rtg_ = dyncorr(Xtf(T(n),:)', Xtarget(G(n),:)', time);
        Rag_ = dyncorr(Xatac(A(n),:)', Xtarget(G(n),:)', time);
        Rat_ = dyncorr(Xtf(T(n),:)', Xatac(A(n),:)', time);
        
        is_cohe(n) = Rtg_*Rag_*Rat_>0;
        Rtg_a(n) = pc_fun(Rtg_,Rag_,Rat_);
        Rag_t(n) = pc_fun(Rag_,Rtg_,Rat_);
        
        Rtg(n) = Rtg_;
        Rat(n) = Rat_;
        Rag(n) = Rag_;
    end
    delete(pool)
    
    FFL = table(A,T,G,Rtg_a,Rag_t,is_cohe);
    save([FOLDER_FILES,'FFL.mat'],'FFL')
    save([FOLDER_FILES,'Rtg.mat'],'Rtg')
    save([FOLDER_FILES,'Rat.mat'],'Rat')
    save([FOLDER_FILES,'Rag.mat'],'Rag')

    clear A T G
    fprintf('\t Finished computing Chain files..\n')
    
    %% indeces of A->T1,T2
    fprintf('\t Started computing T1<-A->T2 files..\n')
    %calculate number of cycles
    n_ATT = 0;
    for n = 1:Natac
        n_tf = sum(peak_tf(n,:));
        if n_tf>1
            n_ATT = n_ATT + size( combnk(1:n_tf,2),1 );
        end
    end
    
    %indeces of all A->T1,T2
    A = zeros(n_ATT,1);
    T1 = zeros(n_ATT,1);
    T2 = zeros(n_ATT,1);
    
    c = 1;
	for n = 1:Natac
        n_tf = sum(peak_tf(n,:));
        if n_tf>1
            pos_tf = find(peak_tf(n,:));
            pos_pair = combnk(pos_tf,2);
            n_pair = size(pos_pair,1);
            
            T1(c:c+n_pair-1) = pos_pair(:,1);
            T2(c:c+n_pair-1) = pos_pair(:,2);
            A(c:c+n_pair-1) = n;
            
            c = c + n_pair;
        end
 	end
    
    assert(c==n_ATT+1)

    % calculate
    is_cohe = false(n_ATT,1);
    Rtt_a = zeros(n_ATT,1); 
    Rtt = zeros(n_ATT,1);
    Rat1 = zeros(n_ATT,1);
    Rat2 = zeros(n_ATT,1);
    
    pool = parpool([2,10]);
    parfor n = 1:n_ATT
        Rtt_ = dyncorr(Xtf(T1(n),:)', Xtf(T2(n),:)', time);
        Rat1_ = dyncorr(Xtf(T1(n),:)', Xatac(A(n),:)', time);
        Rat2_ = dyncorr(Xtf(T2(n),:)', Xatac(A(n),:)', time);
        
        is_cohe(n) = Rtt_*Rat1_*Rat2_>0;
        Rtt_a(n) = pc_fun(Rtt_,Rat1_,Rat2_);
        
        Rtt(n) = Rtt_;
        Rat1(n) = Rat1_
        Rat2(n) = Rat2_;
    end    
    delete(pool) 
    
    ATT = table(A,T1,T2,Rtt_a,is_cohe);
    save([FOLDER_FILES,'ATT.mat'],'ATT') 
    save([FOLDER_FILES,'Rtt.mat'],'Rtt') 
    save([FOLDER_FILES,'Rat1.mat'],'Rat1') 
    save([FOLDER_FILES,'Rat2.mat'],'Rat2')    

    clear A T1 T2 n_tf
    fprintf('\t Finished computing T1<-A->T2 files..\n')
    
    %% indeces of A->G1,G2
    fprintf('\t Started computing G1<-A->G2 files..\n')
    %calculate number of cycles
    n_AGG = 0;
    for n = 1:Natac
        n_target = sum(peak_target(n,:));
        if n_target>1
            n_AGG = n_AGG + size( combnk(1:n_target,2),1 );
        end
    end
    
    %indeces of all A->T1,T2
    A = zeros(n_AGG,1);
    G1 = zeros(n_AGG,1);
    G2 = zeros(n_AGG,1);
    
    c = 1;
	for n = 1:Natac
        n_target = sum(peak_target(n,:));
        if n_target>1
            pos_target = find(peak_target(n,:));
            pos_pair = combnk(pos_target,2);
            n_pair = size(pos_pair,1);
            
            G1(c:c+n_pair-1) = pos_pair(:,1);
            G2(c:c+n_pair-1) = pos_pair(:,2);
            A(c:c+n_pair-1) = n;
            
            c = c + n_pair;
        end
 	end
    
    assert(c==n_AGG+1)

    % calculate
    is_cohe = false(n_AGG,1);
    Rgg_a = zeros(n_AGG,1); 
    Rgg = zeros(n_AGG,1);
    Rag1 = zeros(n_AGG,1);
    Rag2 = zeros(n_AGG,1);
    
    fprintf('\t AGG started..\n')
    pool = parpool([2,10]);
    parfor n = 1:n_AGG
        Rgg_ = dyncorr(Xtarget(G1(n),:)', Xtarget(G2(n),:)', time);
        Rag1_ = dyncorr(Xtarget(G1(n),:)', Xatac(A(n),:)', time);
        Rag2_ = dyncorr(Xtarget(G2(n),:)', Xatac(A(n),:)', time);
        
        is_cohe(n) = Rgg_*Rag1_*Rag2_>0;
        Rgg_a(n) = pc_fun(Rgg_,Rag1_,Rag2_);
        
        Rgg(n) = Rgg_;
        Rag1(n) = Rag1_
        Rag2(n) = Rag2_;
    end    
    delete(pool) 
    
    AGG = table(A,G1,G2,Rgg_a,is_cohe);
    save([FOLDER_FILES,'AGG.mat'],'AGG') 
    save([FOLDER_FILES,'Rgg.mat'],'Rgg') 
    save([FOLDER_FILES,'Rag1.mat'],'Rag1') 
    save([FOLDER_FILES,'Rag2.mat'],'Rag2')    

    clear A G1 G2 n_target
    fprintf('\t Finished computing G1<-A->G2 files..\n')    
    
 

end
