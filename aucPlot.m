function [processed_trace] = aucPlot(result_folder, Fs, thresh,event_label,trace_per_plot,first_event,second_event, auc_thresh)
%% Process the raw trace(auc,raster) and plot the result.

%% inputs:
% folder: the path of the folder, where there contains raw traces and behavior.
% options:
% 'raw': plot the raw trace in original order
% 'specific': plot specific traces. Additional input: n x 1 vector, each element is the id of the neuron.
% 'deconv'�? calculate the deconvolution of the traces
% 'auc_event': calculate the auc by events. Additional input: 0:inhibition 1: excitation
% 'auc_cues': calculate the auc by cues in each events. Additional input: 0: inhibition, 1: excitation
% 'auc_post': calculate the auc: no cue / cue

addpath('E:\wupeixuan\CalImgProcess\CalciumPlot\utils');

%% Get the current list of data
cd(result_folder);
filelist = dir(result_folder);
filelist = {filelist.name};

%% check the file
for i =1:length(filelist)
    if contains(filelist{i},'analysis.xls')
        slx_path = filelist{i};
        continue
    end
    if contains(filelist{i},'_deconv.mat')
        load([result_folder,filesep,filelist{i}]);
        continue
    end
    if contains(filelist{i},'_trace.mat')
        trace_dir = filelist{i};
        continue
    end
end
if ~exist('trace_dir','var')
    trace_dir = input('Please specify the path of the trace:');
end

%% Read the behavior tabel
[~,event_label0,~]=xlsread(slx_path,3,'A:A');
event_label0 = rmmissing(event_label0);
form_label = lower(event_label0');
num_event = length(form_label);
events = cell(num_event,1);
[num,~,~]=xlsread(slx_path,3);
u = isnan(num);
j = 1;
ini_ind = 1;
for i= 1:length(~u)
    if j == num_event
        events{j} = num(ini_ind:length(u),:);
    end
    if ~u(i,1)==0
        events{j}=num(ini_ind:i-1,:);
        ini_ind = i+1;
        j = j+1;
    end
end
neuron.intruder_idx = find(ismember(form_label,'intruder'));
neuron.events = events;
neuron.event_label = event_label;
neuron.action_label = {'positive','passive','mount','intromission','ejaculation','attack','biting','retrieve','beaten','sniff ag','sniff other','sniff pups at nest', 'feeding'};
neuron.color = {[1 0 0],[1 0.63 0],[0 1 0],[0.4940 0.1840 0.5560],[0 0 1],[0 1 1], [1 1 0],[1 0.62 0.83], [1 0.5 0.5], [0.46 0.078 0.57], [1 0.155 0.155], [0.15 1 0.15], [0.82 0.41 0.12]};
neuron.Fs = Fs; %sample rate of final calcium data.
%% Read the traces
[rpath, rname, rtype] = fileparts(trace_dir);
rname_ind = strfind(rname,'_');
neuron.name = rname(1:rname_ind(1)-1);
load(trace_dir);
num_neuron = size(neuron_id,2);
data = trace(2:num_neuron+1,:);
neuron.neuron_id = neuron_id;
[neuron.num_neuron, neuron.nframe] = size(data);

%% Refine the behavioral data
neuron = waive_overlap_behav(neuron, form_label, 4);

%% calculate deconvolution
if ~exist('dec_data','var')
    [dec_data, sig_data] = calculate_deconv(data,neuron, result_folder, filelist,1);
end

%% plot deconvolution data
[dec_ind, dec_seg_ind] = thresh_order(dec_data,neuron,thresh,trace_per_plot);
% plot_data(dec_data,neuron,dec_seg_ind,dec_num_plot,'deconv','plot_thresh',thresh)

% [dec_ind2, dec_seg_ind2, dec_num_plot2] = double_sort(dec_data,neuron,thresh,first_event, second_event);
% plot_data(dec_data,neuron,dec_seg_ind2,dec_num_plot2,'double_sort',thresh,[first_event,second_event],'zeroline','save_path', result_folder)

%% calculate auc
auc_thresh1 = auc_thresh(1);
auc_result = auc_cue_with_blank(dec_data,neuron,auc_thresh1,result_folder,filelist);
plot_data(dec_data,neuron,dec_seg_ind,'auc',auc_result,'zeroline','zero_blank','save_path', result_folder);
% auc_result = auc_cue(dec_data,neuron,auc_thresh1,result_folder,filelist);
% plot_data(dec_data,neuron,dec_seg_ind,'auc',auc_result,'zeroline','save_path', result_folder);
% 
% auc_thresh2 = auc_thresh(2);
% auc_result = auc_cue_with_blank(dec_data,neuron,auc_thresh2,result_folder,filelist);
% plot_data(dec_data,neuron,dec_seg_ind,'auc',auc_result,'zeroline','zero_blank','save_path', result_folder);
% auc_result = auc_cue(dec_data,neuron,auc_thresh2,result_folder,filelist);
% plot_data(dec_data,neuron,dec_seg_ind,'auc',auc_result,'zeroline','save_path', result_folder);
% 
% eval(sprintf('save %s%s%s_neuron.mat neuron', result_folder, filesep, neuron.name))