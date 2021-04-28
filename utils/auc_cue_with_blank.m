function [auc_result, neuron] = auc_cue_with_blank(data,neuron,auc_thresh,filefolder,filelist,varargin)
% Calculate auc by each action in each cue.

Fs = neuron.Fs;
fname = sprintf('auc+blank_%.02f.mat',auc_thresh);
events = neuron.events;
action_label = neuron.action_label;

if ~isempty(varargin)
    if strcmp(varargin{1},'sex')
        fname = sprintf('auc_sex_%.02f.mat',auc_thresh);
        delete_event = neuron.delete_event;
        for i = 1:length(delete_event)
            del_id = delete_event(i);
            events{neuron.intruder_idx}(del_id,:) = [];
            neuron.event_label(del_id) = [];
                delete_event(i:end) = delete_event(i:end) - 1;
        end
    else
        disp('Type unmatch, calculate auc for all events.');
    end
end
        
%% check the file
for i =1:length(filelist)
if contains(filelist{i},fname)
    load([filefolder,filesep,filelist{i}]);
    return
end
end
%% select event state
[num_neuron, nframe] = size(data);
tlen = nframe / Fs;  % total time
intruder = events{neuron.intruder_idx};
nintruder = size(intruder,1);  % intruder总数
nstate = size(action_label,2);  % action总数
zerostates = zeros(nintruder, nframe);
eventstates = zeros(nintruder,nframe);
%% set baseline, intruders and actions
base_begin = -100;
base_end = -2;
for i=1:nintruder
    eventstate = zeros(1,nframe);
    tRise = intruder(i,1); tDur = intruder(i,2) - intruder(i,1);
    flag = revertTTL2bin(tRise, tDur, Fs, tlen);
    eventstate(flag) = 1;
    eventstates(i,:) = i*10*eventstate;
    
    zerostate = zeros(1, nframe);
    tRise = intruder(i,1)+base_begin; tDur = base_end-base_begin;
    if i == 1 && tRise < 0
        tRise = 1/Fs;
    end
    if i>1 && tRise<intruder(i-1,2)
        zerostates(i,:) = zerostates(i-1,:);  %  如果baseline不足100秒，则取前一个cue的baseline
        continue
    end
    flag = revertTTL2bin(tRise, tDur, Fs, tlen);
    zerostate(flag) = 1;
    zerostates(i,:) = zerostate(1:nframe);
end

zerostates = logical(zerostates);

none_action = [];
bstates = zeros(1, nframe);

for i = 1:nstate
    event_now = events{i+2};
    tRise = event_now(:,1); tDur = event_now(:,2) - tRise;
    flag = revertTTL2bin(tRise, tDur, Fs, tlen);
    for j=1:nintruder
        cue_now = intruder(j,:);  % 取第j行的开始和结束时间。
        cue_dur = cue_now(2) - cue_now(1);
        flag_cue = revertTTL2bin(cue_now(1), cue_dur, Fs, tlen);
        tot_flag = flag & flag_cue;
        bstates(tot_flag) = i + j*10;
        if sum(tot_flag)==0
            none_action = [none_action; i j];
        end
%     serial = (1:nframe)'/Fs; %'
%     flagzero = serial>tseg(i) & serial<tseg(i+1) & ~flag;
%     zerostates(flagzero) = i; 
    end
    
end
bstates = bstates(1:nframe);

for i = 1:nintruder
    eventstate = eventstates(i,:) - bstates;
    eventstates(i,:) = eventstate>0;
    zerostates(i,:) = zerostates(i,:)+eventstates(i,:);
end

%% calculate auROC
tic;
aucs = zeros(num_neuron, nstate, nintruder);
pvals = zeros(num_neuron, nstate, nintruder);
rankps = zeros(num_neuron, nstate, nintruder);
for iintruder = 1:nintruder
    for istate = 1:nstate
    %     zerostate = zerostates==istate;
    %     zerostate = zerostates>0;
        if sum(ismember(none_action, [istate, iintruder],'rows'))~=0  % 如若此行为在此intruder中没有出现，则不用计算
            aucs(:, istate, iintruder) = 0;
            pvals(:, istate, iintruder) = 0;
            rankps(:, istate, iintruder) = 0;
            continue
        end
        flag = istate + iintruder*10;
        behavioralstate = (bstates==flag);
        parfor icell = 1:num_neuron
            timeseries = data(icell,:);
            [aucs(icell, istate, iintruder), pvals(icell, istate, iintruder), ~] = auROC_withStates_cxf(timeseries,behavioralstate,zerostates(iintruder,:),100,5000);
            rankps(icell, istate, iintruder) = ranksum(timeseries(behavioralstate), timeseries(zerostates(iintruder,:)));
        end
    end
end
toc

for i = 1:size(none_action,1)
    aucs(:,none_action(i,1),none_action(i,2))=NaN;
    pvals(:,none_action(i,1),none_action(i,2))=NaN;
    rankps(:,none_action(i,1),none_action(i,2))=NaN;
end

if auc_thresh >= 0.5
    out_th=(aucs > auc_thresh) & ~isnan(aucs);
else
    out_th=(aucs < auc_thresh) & ~isnan(aucs);
end
h_signifi = pvals<0.01 & out_th;

auc_result.bstates = bstates;
auc_result.aucs = aucs;
auc_result.pvals = pvals;
auc_result.rankps = rankps;
auc_result.h_signifi = h_signifi;
auc_result.auc_thresh = auc_thresh;

eval(sprintf('save %s%s%s_%s auc_result', filefolder, filesep, neuron.name, fname))
% %% static
% fprintf('\n');
% fprintf('Neuron all = %.0f \n', num_neuron);
% fprintf('Neuron to event 1 = %.0f \n', sum(h_signifi(:,1))/num_neuron);
% fprintf('                2 = %.0f \n', sum(h_signifi(:,2)));
% fprintf('                3 = %.0f \n', sum(h_signifi(:,3)));
% fprintf(' (None response)  = %.0f \n', sum(sum(h_signifi,2) == 0));
% fprintf('\n');
% fprintf('Neuron to event 1&2 = %.0f \n', sum(h_signifi(:,1)&h_signifi(:,2)));
% fprintf('                1&3 = %.0f \n', sum(h_signifi(:,1)&h_signifi(:,3)));
% fprintf('                2&3 = %.0f \n', sum(h_signifi(:,2)&h_signifi(:,3)));
% fprintf('              1&2&3 = %.0f \n', sum(h_signifi(:,1)&h_signifi(:,2)&h_signifi(:,3)));
% fprintf('\n');