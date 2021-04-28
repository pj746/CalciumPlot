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

%% Group the neurons
neuron_order = (1:num_neuron)';
seg_ind = {neuron_order(logical(neuron_ejact)),neuron_order(logical(~neuron_ejact))};
fig_ttl = {'ejaculation';'doesnt ejact'};
plot_data(dec_data,neuron,seg_ind,'deconv','ttl',{fig_ttl,'spont_ejact'});