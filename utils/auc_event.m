function event_auc_result = auc_event(data,neuron,auc_thresh,filefolder,filelist)
% Calculate auc by each cue.

Fs = neuron.Fs;
nframe = neuron.nframe;
num_neuron = neuron.num_neuron;
events = neuron.events;
action_label = neuron.action_label;
%% check the file
fname = sprintf('auc_event_%.02f.mat',auc_thresh);

for i =1:length(filelist)
if contains(filelist{i},fname)
    load([filefolder,filesep,filelist{i}]);
    return
end
end
%% select event state
tlen = nframe / Fs;  % 总时长
% zerostates = zeros(1, nframe);
intruder = events{neuron.intruder_idx};
nintruder = size(intruder,1);  % intruder总数
nstate = size(action_label,2);  % action总数
zerostates = zeros(nintruder, nframe);
bstates = zeros(1, nframe);
%% set zerostate and bstate
% Extrace the waveform of 80-2 seconds before the event as zerostate
base_begin = -100;
base_end = -2;
for i=1:nintruder
    zerostate = zeros(1, nframe);
    tRise = intruder(i,1)+base_begin; tDur = base_end-base_begin;
    if i == 1 && tRise < 0
        tRise = 1/Fs;
    end
    if i>1 && tRise<intruder(i-1,2)
        zerostates(i,:) = zerostates(i-1,:);  %  baseline less than 80sec, extrace the baseline of the previous cue
    end
    flag = revertTTL2bin(tRise, tDur, Fs, tlen);
    zerostate(flag) = 1;
    zerostates(i,:) = zerostate(1:nframe);
   
    bRise = intruder(i,1); bDur = intruder(i,2) - intruder(i,1);
    bflag = revertTTL2bin(bRise, bDur, Fs, tlen);
    bstates(bflag) = i;
end
zerostates = logical(zerostates);
bstates = bstates(1:nframe);

%% calculate auROC
tic
aucs = zeros(num_neuron, nstate);
pvals = zeros(num_neuron, nstate);
rankps = zeros(num_neuron, nstate);
for iintruder = 1:nintruder
%     zerostate = zerostates==istate;
%     zerostate = zerostates>0;
    behavioralstate = (bstates==iintruder);
    for icell = 1:num_neuron
        timeseries = data(icell,:);
        [aucs(icell, iintruder), pvals(icell, iintruder), ~] = auROC_withStates_cxf(timeseries,behavioralstate,zerostates(iintruder,:),100,5000);
        rankps(icell, iintruder) = ranksum(timeseries(behavioralstate), timeseries(zerostates(iintruder,:)));
    end
end
toc

if auc_thresh >= 0.5
    out_th=(aucs > auc_thresh) & aucs~=0;
else
    out_th=(aucs < auc_thresh) & aucs~=0;
end
h_signifi = pvals<0.05 & out_th; %h==1, significant; h==0, non-significant. 1是真的

event_auc_result.bstates = bstates;
event_auc_result.aucs = aucs;
event_auc_result.pvals = pvals;
event_auc_result.rankps = rankps;
event_auc_result.h_signifi = h_signifi;
event_auc_result.auc_thresh = auc_thresh;

eval(sprintf('save %s%s%s_auc_event_%.02f.mat event_auc_result', filefolder, filesep, neuron.name, auc_thresh));

% %% static
% fprintf('\n');
% fprintf('Neuron all = %.0f \n', num_neuron);
% fprintf('Neuron to event 1 = %.0f \n', sum(h_signifi(:,1)));
% fprintf('                2 = %.0f \n', sum(h_signifi(:,2)));
% fprintf('                3 = %.0f \n', sum(h_signifi(:,3)));
% fprintf(' (None response)  = %.0f \n', sum(sum(h_signifi,2) == 0));
% fprintf('\n');
% fprintf('Neuron to event 1&2 = %.0f \n', sum(h_signifi(:,1)&h_signifi(:,2)));
% fprintf('                1&3 = %.0f \n', sum(h_signifi(:,1)&h_signifi(:,3)));
% fprintf('                2&3 = %.0f \n', sum(h_signifi(:,2)&h_signifi(:,3)));
% fprintf('              1&2&3 = %.0f \n', sum(h_signifi(:,1)&h_signifi(:,2)&h_signifi(:,3)));
% fprintf('\n');