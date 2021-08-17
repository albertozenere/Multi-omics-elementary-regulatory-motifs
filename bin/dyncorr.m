%From Rainer Opgen-Rhein and Korbinian Strimmer. REVSTAT – Statistical Journal Volume 4, Number 1, March 2006, 53–65

function dynR = dyncorr(Xmat,varargin)
%   dynR = dyncor(Xmat,Ymat,T) calculates the dynamical correlation between
%   every column of Xmat and every column of Ymat
%   dynR = dyncor(Xmat,T) calculates the dynamical correlation between
%   columns of Xmat
% Xmat,Ymat: n rows (measurements), m columns (variables). T: vector of length n
% containing the timepoint of each measurements. The function return a m x m
% matrix containing the dynamical correlations. Given how the dynamical
% correlation is defined in their paper, T must be the same for X and Y.
% Also, every timepoint must have the same number of replicates.

%%
if nargin == 3
    Ymat = varargin{1};
    T = varargin{2};
elseif nargin == 2
    T = varargin{1};
    Ymat = Xmat;
else
    error('wrong number of inputs')
end
if size(T,1)==1
    T = T';
end

%% checks
assert(min(size(T))==1)
assert(size(Xmat,1)==size(Ymat,1))
assert(size(Xmat,1)==size(T,1))

Tu = unique(T);
replicates_per_timepoint = arrayfun(@(t) sum(T==t), Tu);
assert(~any(replicates_per_timepoint-replicates_per_timepoint(1)),'every timepoint must have same number of replicates')
% assert(replicates_per_timepoint(1)>1, 'The function needs at least 2 replicates per timepoint to calculate Eq.2.4-2.5')

%% do
N = size(Xmat,1);
MX = size(Xmat,2);
MY = size(Ymat,2);
dynR = zeros(MX,MY);
for i = 1:MX
    for j = 1:MY
        dynR(i,j) = dyncor_ij(Xmat(:,i),Ymat(:,j),T);
    end
end

%% dyncorr code
function dynRij = dyncor_ij(X,Y,T)
    Tu = sort(unique(T),'ascend');

    %% group samples in X and Y based on T
    Xg = cell2mat(arrayfun(@(t) X(T==t)', Tu, 'uni', false));
    Yg = cell2mat(arrayfun(@(t) Y(T==t)', Tu, 'uni', false));
    assert(size(Xg,1)==length(Tu) && size(Xg,2)== replicates_per_timepoint(1))
    assert(all(size(Xg)==size(Yg)))
    dt = [Tu; Tu(end)] - [Tu(1); Tu]; %Eq 2.2
    wt =  (dt(1:end-1)+dt(2:end))/(2*(Tu(end)-Tu(1))) ; %Eq 2.2
    if ~any(Tu-Tu(1))
        wt = 1;
    end
    
    %% center each replicate with mean of mean
    Xm = mean(Xg,2);
    Xc = Xg - Xm'*wt;
    Ym = mean(Yg,2);
    Yc = Yg - Ym'*wt;

    %% calculate variance
    nr = size(Xg,2);
    assert(nr == size(Yg,2));
    varX = 0;
    for n = 1:nr %Eq.2.4
        varX = varX + Xc(:,n)' * (Xc(:,n).*wt);
    end
	varX = varX/( max(2,nr)-1); %to account for the case nr=1
    
    varY = 0;
    for n = 1:nr %Eq.2.4
        varY = varY + Yc(:,n)' * (Yc(:,n).*wt);
    end
	varY = varY/( max(2,nr)-1);
    
    Xs = Xc/sqrt(varX);
    Ys = Yc/sqrt(varY);

    %% dynamic correlation
    dynRij = 0;
    for n = 1:nr
        dynRij = dynRij + Xs(:,n)'*(Ys(:,n).*wt);
    end
	dynRij = dynRij/( max(2,nr)-1);
    
end
end