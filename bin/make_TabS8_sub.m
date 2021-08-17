%calculate the percentage of coherent/c.i. triplets for a given TF

function list_tg = make_TabS8_sub(pos_tf,pos_chr,peak_target,Xtarget,Xtf,Xatac,time)

assert(pos_tf<=size(Xtf,1), 'this gene is not a TF')

%% statistical test
pc_fun = @(pxy,pxz,pyz) (pxy-pxz*pyz)/sqrt((1 - pxz^2)*(1 - pyz^2));    

NR = 3; %number replicates (same for all datasets considered)
ntp = size(Xatac,2)/NR;
alpha = 0.05;
t_crit = @(n_control) tinv(.5+alpha/2, ntp-2-n_control);
t = @(pc,n_control) pc.*sqrt( (ntp-2-n_control)./(1-pc.^2) );

THS_0 = t_crit(0);
THS_1 = t_crit(1);
R_not_zero = @(vec) abs(t(vec,0))>THS_0;
is_pc = @(vec) abs(t(vec,1))<THS_1;

%% do
N = length(pos_chr);
Rtg_cell = cell(1,N); Rag_cell = cell(1,N); Rat_cell = cell(1,N);
pc_a = cell(1,N); pc_t = cell(1,N);

pos_tg_vec = cell(1,N);

for n = 1:N %loop each chr
    pos_tg = find(peak_target(pos_chr(n),:));
    M = length(pos_tg);
    Rtg = zeros(M,1); Rag = zeros(M,1); Rat = zeros(M,1);
    Rtg_a = zeros(M,1); Rag_t = zeros(M,1);
    
    pos_tg_vec{n} = pos_tg;
    for k = 1:M %loop each target
        Rtg(k) = dyncorr( Xtarget(pos_tg(k),:)', Xtf(pos_tf,:)', time);
        Rag(k) = dyncorr( Xtarget(pos_tg(k),:)', Xatac(pos_chr(n),:)', time);
        Rat(k) = dyncorr( Xtf(pos_tf,:)', Xatac(pos_chr(n),:)', time);

        Rtg_a(k) = pc_fun(Rtg(k), Rat(k), Rag(k));
        Rag_t(k) = pc_fun(Rag(k), Rat(k), Rtg(k));
        
    end
    Rtg_cell{n} = Rtg'; Rag_cell{n} = Rag'; Rat_cell{n} = Rat';
    pc_a{n} = Rtg_a';
    pc_t{n} = Rag_t';
    
end
Rtg_cell = [Rtg_cell{:}]'; Rag_cell = [Rag_cell{:}]'; Rat_cell = [Rat_cell{:}]';
pos_tg_vec = [pos_tg_vec{:}]';

pc_a = [pc_a{:}]';
pc_t = [pc_t{:}]';

to_keep = R_not_zero(Rtg_cell) & R_not_zero(Rag_cell) & R_not_zero(Rat_cell); 

pc_a = pc_a(to_keep);
pc_t = pc_t(to_keep);
pos_tg_vec = pos_tg_vec(to_keep);

assert(length(pos_tg_vec)==length(pc_a))

list_tg = pos_tg_vec(is_pc(pc_a));



