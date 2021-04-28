% AUROC: [auc, pval, shuff_auc] = auROC(x, y, N, NUM_ITER)
%  auROC computes the area under the ROC given two distributions x and y, 
%   and optionally does significance tests.  Significance tests (by 
%   permutation tests) are only performed if one of the optional 2nd or 3rd
%   outputs are saved in the call.
%
%  Inputs (required):
%   - x and y: two DATA SERIES (NOT DISTRIBUTIONS!). Need not be the same size or dimensions.
% 
%  Inputs (optional):
%   - N: number of threshold values used to construct ROC. [def=100]
%   - NUM_ITER: # of shufflings to perform for any signif tests [def=5000]
%  Outputs:
%   - auc: area under ROC curve, such that 1 indicates x > y, reliably.
%   - pval: fraction of shufflings which are further from 0.5 than true auc
%   - shuff_auc: the vector of shuffled auc's.  This is useful for
%                constructing confidence intervals.

% Adapted from code of Vinod Rao by A.M.

function [auc, pval, shuff_auc] = auROC_withStates_cxf(timeseries,behavioralstate,zerostate, N, NUM_ITER)

%% Check variables
if nargin<3 || isempty(N)
    N = 100; %default number of threshold values used to construct ROC.
end
if nargin<4
    NUM_ITER = 15000;
end
if nargout>=2
    testflag = 1; %flag to do the simulations
else
    testflag = 0; 
end

%% check
nsample = length(timeseries);
if islogical(behavioralstate)
    assert(length(behavioralstate)==nsample);
    behavioralstate = find(behavioralstate); 
end
if islogical(zerostate)
    assert(length(zerostate)==nsample);
    zerostate = find(zerostate); 
end
timeseries = timeseries(:);  %as column vector
behavioralstate = behavioralstate(:);
zerostate = zerostate(:);
nbstate = length(behavioralstate);
nzstate = length(zerostate);

%% DON'T balence the number of behavioralstate v.s. zerostate
if nbstate<nzstate
    % nzstate large
    nrepeat = floor(nzstate / nbstate);
    behavioralstate   = repmat(behavioralstate, nrepeat, 1);
else  
    % nzstate large
    nrepeat = floor(nbstate / nzstate);
    zerostate   = repmat(zerostate, nrepeat, 1);
end
nbstate = length(behavioralstate);
nzstate = length(zerostate);

%% bootstrap in 'behavioralstate'
x0=timeseries(behavioralstate); %n_by_bootstrap
y =timeseries(zerostate);
rand_bstate = zeros(nbstate, NUM_ITER+1);
rand_bstate(:,1) = (1:nbstate)';  %  帧数
reservoir = [x0; y];  % 每一帧的电信号
% reservoir = timeseries;
nreservoir = length(reservoir);
for i=1:NUM_ITER
    rand_bstate(:,i+1) = randperm(nreservoir, nbstate);
end

%% Compute ROC and calculate the area under the curve
x=reservoir(rand_bstate); %n_by_bootstrap
y=timeseries(zerostate);  %n_by_1

threshlo = min([x(:); y]);
threshhi = max([x(:); y]);
thresh = linspace(threshlo,threshhi,N); %thresh = linspace(threshlo,threshhi,N*5);
fa = cumsum(flip(histc(y, thresh, 1) / length(y))); % flip of cdf
hit = cumsum(flip(histc(x, thresh, 1) / nbstate));
aucs = trapz(fa,hit); 
auc = aucs(1);


%% if necessary, bootstrap the ROC for significance testing
% The idea is to randomly assign the same values into two groups and
% recompute ROCs. But 
shuff_auc = aucs(2:end);
% shuff_auc2 = abs(shuff_auc-0.5);

pval = sum( abs(shuff_auc-0.5) > abs(auc-0.5) ) / NUM_ITER; 
    
if nargout < 2
    varargout = [];
elseif nargout == 2
    varargout{1} = pval;
elseif nargout == 3
    varargout{1} = pval;
    varargout{2} = shuff_auc;
end

end
