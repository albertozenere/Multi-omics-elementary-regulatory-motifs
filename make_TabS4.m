%Generate the lists of genes in Tab.S4-S5-S6-S7, in a format ready to be copy-pasted in
%Latex

clear all
close all
clc

ALL_FOLDER_FILES = {'data\th1\', 'data\th1_p4\', 'data\Bonneau\hiv\', 'data\Bonneau\mock\'};
ALL_FIG_NAME = {'th1','p4','DC','mock'};

%% KEGG pathways
pathway_gene = importdata('KEGG_pathways/pathway_gene.mat');

keyword_th1 = ["th1_"; "t_cell"];
idx = contains(pathway_gene(:,1),keyword_th1,'IgnoreCase',true);
pathways_th1 = pathway_gene(idx, 2)';
pathways_th1 = unique([pathways_th1{:}])';

keyword_dc = ["antigen"; "toll"; "kappa"];
idx = contains(pathway_gene(:,1),keyword_dc,'IgnoreCase',true);
pathways_dc = pathway_gene(idx, 2)';
pathways_dc = unique([pathways_dc{:}])';

%% do
tf_th1 = strings(length(pathways_th1),1);
tg_th1 = cell(length(pathways_th1),2);
tf_dc = strings(length(pathways_dc),1);
tg_dc = cell(length(pathways_dc),2);

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

    %% do
    if contains(FOLDER_FILES, "th1")

        for n_gene = 1:length(pathways_th1)
            g = pathways_th1(n_gene);
            if ismember(g,list_tf)
                tf = g;
                tf_th1(n_gene) = g;
                pos_tf = find(ismember(list_tf,tf));
                pos_chr = find(peak_tf(:,pos_tf));
                pos_tg = find(tf_target(pos_tf,:))';

                list_tg = example_sub(pos_tf,pos_chr,peak_target,Xtarget,Xtf,Xatac,time,pc_fun, THS_0,THS_1);
                list_tg_name = convert_id(list(list_tg),'ensg','name');
                list_tg_name = string(list_tg_name(arrayfun(@(n) ~isempty(list_tg_name{n}), 1:length(list_tg_name))));

                tg_th1{n_gene,n_data} = unique(list_tg_name);
            end
        end
       
    else
        
        for n_gene = 1:length(pathways_dc)
            g = pathways_dc(n_gene);
            tf_dc(n_gene) = g;
            if ismember(g,list_tf)
                tf = g;
                pos_tf = find(ismember(list_tf,tf));
                pos_chr = find(peak_tf(:,pos_tf));
                pos_tg = find(tf_target(pos_tf,:))';

                list_tg = example_sub(pos_tf,pos_chr,peak_target,Xtarget,Xtf,Xatac,time,pc_fun, THS_0,THS_1);
                list_tg_name = convert_id(list(list_tg),'ensg','name');
                list_tg_name = string(list_tg_name(arrayfun(@(n) ~isempty(list_tg_name{n}), 1:length(list_tg_name))));

                tg_dc{n_gene,n_data-2} = unique(list_tg_name);
            end
        end
    
    end
    
end


%% remove empty
i = arrayfun(@(n) isempty(tg_th1{n,1}), 1:size(tg_th1,1) );
tg_th1(i,:) = [];
tf_th1(i) = [];

i = arrayfun(@(n) isempty(tg_dc{n,1}), 1:size(tg_dc,1) );
tg_dc(i,:) = [];
tf_dc(i) = [];

%% overlap
overlap_th1 = zeros(size(tg_th1,1),1);
for n = 1:length(overlap_th1)
    overlap_th1(n) = length(intersect(tg_th1{n,1},tg_th1{n,2}))/length(tg_th1{n,1});
end

overlap_dc = zeros(size(tg_dc,1),1);
for n = 1:length(overlap_dc)
    overlap_dc(n) = length(intersect(tg_dc{n,1},tg_dc{n,2}))/length(tg_dc{n,1});
end

%% convert id
tf_th1 = string(convert_id(tf_th1, 'ensg', 'name'));
tf_dc = string(convert_id(tf_dc, 'ensg', 'name'));

assert(length(tf_th1)==size(tg_th1,1))
assert(length(tf_dc)==size(tg_dc,1))

[~,s] = sort(tf_th1);
tf_th1 = tf_th1(s);
tg_th1 = tg_th1(s,:);

[~,s] = sort(tf_dc);
tf_dc = tf_dc(s);
tg_dc = tg_dc(s,:);

%% write on text file
fileID = fopen('tf_tg_th1.txt','w');
for n = 1:length(tf_th1)
    fprintf(fileID, '%s & \t', tf_th1(n));
    fprintf(fileID, '%s, \t', tg_th1{n,1}(1:end-1));
    fprintf(fileID, '%s \\\\ \n', tg_th1{n,1}(end));
end
fclose(fileID);

fileID = fopen('tf_tg_p4.txt','w');
for n = 1:length(tf_th1)
    fprintf(fileID, '%s & \t', tf_th1(n));
    fprintf(fileID, '%s, \t', tg_th1{n,2}(1:end-1));
    fprintf(fileID, '%s \\\\ \n', tg_th1{n,2}(end));
end
fclose(fileID);

fileID = fopen('tf_tg_dc.txt','w');
for n = 1:length(tf_dc)
    fprintf(fileID, '%s & \t', tf_dc(n));
    fprintf(fileID, '%s, \t', tg_dc{n,1}(1:end-1));
    fprintf(fileID, '%s \\\\ \n', tg_dc{n,1}(end));
end
fclose(fileID);

fileID = fopen('tf_tg_mock.txt','w');
for n = 1:length(tf_dc)
    fprintf(fileID, '%s & \t', tf_dc(n));
    fprintf(fileID, '%s, \t', tg_dc{n,2}(1:end-1));
    fprintf(fileID, '%s \\\\ \n', tg_dc{n,2}(end));
end
fclose(fileID);
