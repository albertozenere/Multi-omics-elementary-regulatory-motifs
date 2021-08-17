%Creates the list of genes in Tab.S8, look at bottom of script

clear all
close all
clc

ALL_FOLDER_FILES = {'data\th1\', 'data\th1_p4\', 'data\Bonneau\hiv\', 'data\Bonneau\mock\'};
ALL_FIG_NAME = {'th1','p4','DC','mock'};

stat1_tg = cell(2,1);
stat4_tg = cell(2,1);
tbet_tg = cell(2,1);

irf3_tg = cell(2,1);
nfkb_tg = cell(2,1);

for n_data = 1:length(ALL_FOLDER_FILES)
    
    %% load
    FOLDER_FILES = ALL_FOLDER_FILES{n_data};
    
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

    if contains(FOLDER_FILES, "th1")
        time = [0 0 0, .5 .5 .5, 1 1 1, 2 2 2, 6 6 6, 24 24 24]; %th1
    else
        time = [2 2 2, 8 8 8, 24 24 24, 48 48 48]; %hiv
    end
    ORG = "human";

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

    pval_fun = @(pc,n_control) tcdf(t(pc,n_control),ntp-2-n_control);


    if contains(FOLDER_FILES, "th1")

        %% TBX21
        tbet = "ENSG00000073861";
        pos_tbet = find(ismember(list_tf,tbet));
        pos_chr_tbet = find(peak_tf(:,pos_tbet)); 
        pos_tg_tbet = find(tf_target(pos_tbet,:))';

        list_tg = make_TabS8_sub(pos_tbet,pos_chr_tbet,peak_target,Xtarget,Xtf,Xatac,time);
        list_tg_name = convert_id(list(list_tg),'ensg','name');
        list_tg_name = string(list_tg_name(arrayfun(@(n) ~isempty(list_tg_name{n}), 1:length(list_tg_name))));

        tbet_tg{n_data} = list_tg_name;
        
        %% stat1
        stat1 = "ENSG00000115415";
        pos_stat1 = find(ismember(list_tf,stat1));
        pos_chr_stat1 = find(peak_tf(:,pos_stat1)); 
        pos_tg_stat1 = find(tf_target(pos_stat1,:));

        list_tg = make_TabS8_sub(pos_stat1,pos_chr_stat1,peak_target,Xtarget,Xtf,Xatac,time);
        list_tg_name = convert_id(list(list_tg),'ensg','name');
        list_tg_name = string( list_tg_name(arrayfun(@(n) ~isempty( list_tg_name{n} ), 1:length(list_tg_name) )));

        stat1_tg{n_data} = list_tg_name;
        
        %% stat4
        stat4 = "ENSG00000138378";
        pos_stat4 = find(ismember(list_tf,stat4));
        pos_chr_stat4 = find(peak_tf(:,pos_stat4)); 
        pos_tg_stat4 = find(tf_target(pos_stat4,:));

        list_tg = make_TabS8_sub(pos_stat4,pos_chr_stat4,peak_target,Xtarget,Xtf,Xatac,time);
        list_tg_name = convert_id(list(list_tg),'ensg','name');
        list_tg_name = string(list_tg_name(arrayfun(@(n) ~isempty(list_tg_name{n}), 1:length(list_tg_name))));

        stat4_tg{n_data} = list_tg_name;
        
    else
        
        %% NFKB1 
        nfkb1 = "ENSG00000170345";
        pos_nfkb1 = find(ismember(list_tf,nfkb1));
        pos_chr_nfkb1 = find(peak_tf(:,pos_nfkb1));
        pos_tg_nfkb1 = find(tf_target(pos_nfkb1,:)); 

        list_tg = make_TabS8_sub(pos_nfkb1,pos_chr_nfkb1,peak_target,Xtarget,Xtf,Xatac,time);
        list_tg_name = convert_id(list(list_tg),'ensg','name');
        list_tg_name = string(list_tg_name(arrayfun(@(n) ~isempty(list_tg_name{n}), 1:length(list_tg_name))));

        nfkb_tg{n_data-2} = list_tg_name;

        %% IRF3
        irf3 = "ENSG00000126456";
        pos_irf3 = find(ismember(list_tf,irf3));
        pos_chr_irf3 = find(peak_tf(:,pos_irf3)); 
        pos_tg_irf3 = find(tf_target(pos_irf3,:)); 

        list_tg = make_TabS8_sub(pos_irf3,pos_chr_irf3,peak_target,Xtarget,Xtf,Xatac,time);
        list_tg_name = convert_id(list(list_tg),'ensg','name');
        list_tg_name = string( list_tg_name(arrayfun(@(n) ~isempty( list_tg_name{n} ), 1:length(list_tg_name) )));
        
        irf3_tg{n_data-2} = list_tg_name;
    
    end
    
end

%% Results

%% dataset A
%selected targets of STAT1 in dataset A 
stat1_A = sort(stat1_tg{1});
%selected targets of STAT4 in dataset A 
stat4_A = sort(stat4_tg{1});
%selected targets of t-bet in dataset A 
tbet_A = sort(tbet_tg{1});

%% dataset B
%selected targets of STAT1 in dataset B 
stat1_B = sort(stat1_tg{2});
%selected targets of STAT4 in dataset B 
stat4_B = sort(stat4_tg{2});
%selected targets of t-bet in dataset B 
tbet_B = sort(tbet_tg{2});

%% dataset C
%selected targets of IRF3 in dataset C
irf3_C = sort(irf3_tg{1});
%selected targets of nfkb in dataset C
nfkb_C = sort(nfkb_tg{1});

%% dataset D
%selected targets of IRF3 in dataset D
irf3_D = sort(irf3_tg{2});
%selected targets of nfkb in dataset D
nfkb_D = sort(nfkb_tg{2});




