clc; clear;
addpath('E:\wupeixuan\CalImgProcess\CalciumPlot\utils');
%% load the behavioral annotation
slx_path = 'E:\wupeixuan\auc_plot\data\Esr29m\Esr290218\Esr2-9+20210218+updated_analysis.xlsx';
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
for i = 1:length(events)
    if isnan(events{i}(1,2))
        events{i}(:,2) = events{i}(:,1)+0.01;
    end
end
neuron.Fs = 10; %sample rate of final calcium data.
neuron.intruder_idx = find(ismember(form_label,'intruder'));
neuron.events = events;
neuron.event_label = {'pups','R female','Another R female'};
action_label = {'positive','negative','mount','intromission','ejaculation','attack','biting','thrust','penile groom'};
color = {[1 0 0],[1 0.63 0],[0 1 0],[0.4940 0.1840 0.5560],[0 0 1],[0 1 1], [1 1 0], [0.12 0.56 1],[0.82 0.41 0.12]};
%% refine the action_label and reorganize the events
u = zeros(1,length(action_label));
for i = 3:length(form_label)
    contain = strncmp(action_label,form_label{i},4);
    u = u+contain;
end
u = logical(u);
neuron.action_label = action_label(u);
neuron.color = color(u);
for i = 3:length(form_label)
    contain = strncmp(neuron.action_label,form_label{i},4);
    idx = find(contain,1);    
    neuron.events{idx+2} = events{i};
end

trace_dir =  'E:\wupeixuan\auc_plot\data\Esr29m\Esr290218\Esr290218_7_prime_combined.mat';
result_folder = 'E:\wupeixuan\auc_plot\temp_util';
[rpath, rname, rtype] = fileparts(trace_dir);
rname_ind = strfind(rname,'_');
neuron.name = rname(1:rname_ind(1)-1);
load(trace_dir);
num_neuron = size(neuron_id,2);
data = trace(2:num_neuron+1,:);
neuron.neuron_id = neuron_id;
[neuron.num_neuron, neuron.nframe] = size(data);

filelist = dir(result_folder);
filelist = {filelist.name};

%% calculate deconvolution
[dec_data, sig_data] = calculate_deconv(data,neuron, result_folder, filelist,1);

seg_ind = {[22 24 51 56 65 67 84]};
num_plot = 1;

data1 = sig_data;
data2 = dec_data;
for k=1:num_plot
%     data_plot = cell2mat(seg_trace(k));
    ind_plot = cell2mat(seg_ind(k));
    figure; 
    ind_choose = ind_plot(:);
    n_neu = length(ind_choose);
    hold on;
    for i=1:n_neu
        ii = ind_choose(i);
        data_now = 5*data1(ii, :)-(i-1)*10*10.5;%  scale ����trace
        h = plotHz(neuron.Fs, data_now, 'k');
        set(h, 'xdata', get(h, 'xdata'));
        MakeDrag(h, 'y');
        data_now = 0.75*data2(ii, :)-(i-0.5)*10*10.5;%  scale ����trace
        h = plotHz(neuron.Fs, data_now, 'k');
        set(h, 'xdata', get(h, 'xdata'));
        MakeDrag(h, 'y');
    end

h8 =  barline(neuron.events{neuron.intruder_idx}(:,1), [-19700 40], 'k', 'linewidth', 0.05);
h9 = barline(neuron.events{neuron.intruder_idx}(:,2), [-19700 40], '--k', 'linewidth', 0.05);

for i=1:size(neuron.event_label,2)
    text(neuron.events{neuron.intruder_idx}(i,1),80,char(neuron.event_label(i)));
end
for i=3:size(neuron.events,1)
    if neuron.events{i}(1,1) == 0 continue; end
    h = barpatch(neuron.events{i}(:,1), neuron.events{i}(:,2)-neuron.events{i}(:,1), [-19700 40],cell2mat(neuron.color(i-2)));
    set(h, 'facealpha', 0.5);
    text(neuron.events{i}(1,2),65,char(neuron.action_label(i-2)),'FontSize',6);
end
ttl = sprintf('%s    Raster',neuron.name);
title(ttl);
xlabel('Time (sec)');
last_num = ceil(length(data)/neuron.Fs/100)*100;
% xticks([0:75:last_num]);xlim([0 last_num]);
xticks([825:50:1800]);xlim([825 1800]);
yticklabel = arrayfun(@num2str, neuron.neuron_id(ind_choose(n_neu:-1:1)), 'UniformOutput', false); 
set(gca, 'ytick',(-n_neu+1:0)*10*10.5,'YTickLabel',yticklabel,'Position',[0.03,0.06,0.96,0.92]);
set(gca, 'tickdir', 'out');
ylim([min(data_now)- 5, 100]);
ylabel('Cell (#)')
fig_ttl = sprintf('./raster_intromission_%d.fig',k);
savefig(fig_ttl);
end