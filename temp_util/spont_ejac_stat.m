clc; clear;
addpath('E:\wupeixuan\CalImgProcess\CalciumPlot\utils');
trace_dir =  'E:\wupeixuan\auc_plot\data\Esr29m\Esr290218\Esr290218_7_prime_combined.mat';
result_folder = 'E:\wupeixuan\auc_plot\data\Esr29m\Esr290218\spont_ejact_result';
filelist = dir(result_folder);
filelist = {filelist.name};
%% load the behavioral annotation
slx_path = 'E:\wupeixuan\auc_plot\data\Esr29m\Esr290218\Esr2-9+20210218+analysis2.xlsx';
[neuron, form_label, events] = behav_read(slx_path);
neuron.Fs = 10; %sample rate of final calcium data.
neuron.event_label = {'pups','R female','Another R female'};
neuron.action_label = {'positive','passive','mount','intromission','ejaculation','attack','biting'};
neuron.color = {[1 0 0],[1 0.63 0],[0 1 0],[0.4940 0.1840 0.5560],[0 0 1],[0 1 1], [1 1 0]};
%% Read the traces
[rpath, rname, rtype] = fileparts(trace_dir);
rname_ind = strfind(rname,'_');
neuron.name = rname(1:rname_ind(1)-1);
load(trace_dir);
num_neuron = size(neuron_id,2);
data = trace(2:num_neuron+1,:);
neuron.neuron_id = neuron_id;
[neuron.num_neuron, neuron.nframe] = size(data);
[dec_data, sig_data] = calculate_deconv(data,neuron, result_folder, filelist,1);
% dec_data = dec_data(1,:);
% sig_data = sig_data(1,:);

%% Refine the behavioral data
neuron = waive_overlap_behav(neuron, form_label, 4);

tlen = neuron.nframe/neuron.Fs;

%% Find the neuron that spikes during ejaculation
ejact_idx = find(ismember(neuron.action_label,'ejaculation'));
tRise = neuron.events{ejact_idx+2}(1)-2;
tDur = neuron.events{ejact_idx+2}(2)-tRise+2;
ejact_time = [tRise, tDur];
ejact_state = revertTTL2bin(tRise, tDur, neuron.Fs, tlen);
ejact_state = ejact_state';
ejact_thresh = 0.4;
norm_data = normalize(dec_data,2,'range');
ejact_data = norm_data(:,ejact_state);
neuron_ejact = zeros(num_neuron,1);
for i=1:num_neuron
    response_ind = find(ejact_data(i,:)>ejact_thresh,1);
    if ~isempty(response_ind)
        neuron_ejact(i) = 1;
    end
end

%% Find the neuron that spikes spontaneous before and after the tests
% waive the largest peak and normalize the trace, to filter some small spikes
sig_waive_peak = zeros(num_neuron,neuron.nframe);
for i=1:neuron.num_neuron
    large_idx = find(norm_data>0.99,1);
    large_time = large_idx/neuron.Fs;
    tRise = large_time-5; tDur=10;
    large_time = revertTTL2bin(tRise, tDur, neuron.Fs, tlen);
    revert_large = (~large_time)';
    sig_waive_peak(i,:) = normalize(sig_data(i,:).*revert_large,'range');
    sig_waive_peak = (sig_waive_peak>0.1);
end
sig_data0 = sig_data.*sig_waive_peak;

% First spontaneous
spont_end1 = neuron.events{neuron.intruder_idx}(1);
tRise = 1/neuron.Fs; tDur = spont_end1-2;
spont_time1 = [tRise, tDur];
spont1 = revertTTL2bin(tRise, tDur, neuron.Fs, tlen);
spont1 = spont1';
spont_trace1 = dec_data(:,spont1);
spont_spike01 = sig_data0(:,spont1);

% Last spontaneous
tRise = neuron.events{neuron.intruder_idx}(end);
tDur = neuron.nframe/neuron.Fs - tRise;
spont_time2 = [tRise, tDur];
spont2 = revertTTL2bin(tRise, tDur, neuron.Fs, tlen);
spont2 = spont2';
spont_trace2 = dec_data(:,spont2);
spont_spike02 = sig_data0(:,spont2);
spont_spike0 = [spont_spike01, spont_spike02];

% a. normalize two pieces as a whole
% spont_trace = normalize([spont_trace1,spont_trace2],'range');
% spont_spike = normalize([spont_spike1,spont_spike2],'range');

% b. normalize two pieces seperately
spont_trace_norm1 = normalize(spont_trace1,2,'range');
spont_spike_norm1 = normalize(spont_spike01,2,'range');
spont_trace_norm2 = normalize(spont_trace2,2,'range');
spont_spike_norm2 = normalize(spont_spike02,2,'range');
spont_trace = [spont_trace_norm1,spont_trace_norm2];

threshold = 0.3;
spont_spike1 = (spont_spike_norm1>threshold);
spont_spike2 = (spont_spike_norm2>threshold);
spont_spike = [spont_spike1, spont_spike2];

%% plot the data as example
neuron_num = 5;
figure('name',sprintf('Raster of spontaneous,neuron %d',neuron.neuron_id(neuron_num))); 
subplot(3,1,1);
plotHz(neuron.Fs, spont_trace(neuron_num,:));hold on;
barline(spont_end1,[0,1], 'k', 'linewidth', 0.05);
subplot(3,1,2);
plotHz(neuron.Fs, spont_spike0(neuron_num,:));hold on;
barline(spont_end1,[0,1], 'k', 'linewidth', 0.05);
subplot(3,1,3);
plotHz(neuron.Fs, spont_spike(neuron_num,:));hold on;
barline(spont_end1,[0,1], 'k', 'linewidth', 0.05);
last_num = ceil(length(spont_spike)/neuron.Fs/100)*100;
xticks([0:75:last_num]);xlim([0 last_num]); 
yticks([0 0.5 1]);

%% Calculate the spike count
spont_count1 = sum(spont_spike1,2);
spont_count2 = sum(spont_spike2,2);
%% Plot the histogram of the spike count
edge = max([spont_count1;spont_count2]);
space = 15;
bin = [0:space:ceil(edge/space)*space];
space2 = 30;
bin2 = [0:space2:ceil(edge/space)*space];
figure;
subplot(1,2,1);
plot(histcounts(spont_count1,bin));
yticks([0:5:60]);
xticks([1:2:length(bin)]);
xticklabels(bin2);
title('Spontaneous count1');
subplot(1,2,2);
plot(histcounts(spont_count2,bin));
yticks([0:5:60]);
xticks([1:2:length(bin)]);
xticklabels(bin2);
title('Spontaneous count2');

cut = 30;  % Most of the spike_count concentrates in [0, 30]
% find the neuron_id with spike above the cut
neuron_above1 = spont_count1>cut;
neuron_above2 = spont_count2>cut;

%% Group the neurons
spont_ejact_state = [neuron_above1, neuron_ejact, neuron_above2];
neuron123 = neuron_above1&neuron_above2&neuron_ejact;
neuron12 = neuron_above1&neuron_ejact-neuron123;
neuron13 = neuron_above1&neuron_above2-neuron123;
neuron23 = neuron_ejact&neuron_above2-neuron123;
neuron1 = neuron_above1-neuron123-neuron12-neuron13;
neuron2 = neuron_ejact-neuron123-neuron12-neuron23;
neuron3 = neuron_above2-neuron123-neuron23-neuron13;
neuron_single = neuron1&neuron2&neuron3;
rest_neuron = ~(sum(spont_ejact_state,2)>0);

neuron_order = (1:num_neuron)';
seg_ind = {neuron_order(neuron123),neuron_order(neuron12),neuron_order(neuron13),neuron_order(neuron23),neuron_order(neuron_single),neuron_order(rest_neuron)};
fig_ttl = {'spontaneous+ejaculation+spontaneous';'spontaneous+ejaculation';'spontaneous+spontaneous';'ejaculation+spontaneous';'single event';'other neurons'};
ind = ones(1,6);
for i=1:6
    if isempty(seg_ind{i}) == 1
        ind(i)=0;
    end
end
seg_ind=seg_ind(logical(ind));
fig_ttl=fig_ttl(logical(ind));
boxes = {spont_ejact_state,[spont_time1;ejact_time;spont_time2]};
plot_data(80*norm_data,neuron,seg_ind,'box',boxes,'ttl',{fig_ttl,'spont_ejact'});