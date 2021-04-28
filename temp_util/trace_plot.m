result_folder = 'E:\wupeixuan\auc_plot\data\St187f\St1870426_spontaneous';
%% Get the current list of data
cd(result_folder);
filelist = dir(result_folder);
filelist = {filelist.name};

%% check the file
for i =1:length(filelist)
    if contains(filelist{i},'_deconv.mat')
        load([result_folder,filesep,filelist{i}]);
        continue
    end
    if contains(filelist{i},'_trace.mat')
        trace_dir = filelist{i};
        continue
    end
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

%% calculate deconvolution
if ~exist('dec_data','var')
    [dec_data, sig_data] = calculate_deconv(data,neuron, result_folder, filelist,1);
end

[dec_ind, dec_seg_ind] = thresh_order(dec_data,neuron,thresh,trace_per_plot);
plot_data(dec_data,neuron,dec_seg_ind,dec_num_plot,'deconv','plot_thresh',thresh)

