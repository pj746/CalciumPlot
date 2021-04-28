function plot_general(trace_dir,result_folder,neuron,thresh,trace_per_plot,first_event,second_event)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%% load data
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

thresh = 0.99;
trace_per_plot = 19;
%% plot raw data
[raw_ind, raw_seg_ind, raw_num_plot] = thresh_order(data, neuron, thresh, trace_per_plot);
plot_data(data,neuron,raw_seg_ind,raw_num_plot,'raw','plot_thresh','zero_line','save_path', result_folder);

%% calculate deconvolution
[dec_data, sig_data] = calculate_deconv(data,neuron, result_folder, filelist,1);

%% plot deconvolution data
[dec_ind, dec_seg_ind, dec_num_plot] = thresh_order(dec_data,neuron,thresh,trace_per_plot);
plot_data(dec_data,neuron,dec_seg_ind,dec_num_plot,'deconv','plot_thresh',thresh,'save_path', result_folder)

[dec_ind2, dec_seg_ind2, dec_num_plot2] = double_sort(dec_data,neuron,thresh,first_event, second_event);
plot_data(dec_data,neuron,dec_seg_ind2,dec_num_plot2,'double_sort',thresh,[first_event,second_event],'save_path', result_folder)

%% calculate auc
auc_thresh = 0.65;
auc_result = auc_cue(dec_data,neuron,auc_thresh,result_folder,filelist);
plot_data(data,neuron,dec_seg_ind,dec_num_plot,'auc',auc_result);

auc_result = auc_event(dec_data,neuron,auc_thresh,result_folder,filelist);
plot_data(data,neuron,dec_seg_ind,dec_num_plot,'auc',auc_result);

auc_thresh = 0.35;
auc_result = auc_cue(dec_data,neuron,auc_thresh,result_folder,filelist);
plot_data(data,neuron,dec_seg_ind,dec_num_plot,'auc',auc_result);

auc_result = auc_event(dec_data,neuron,auc_thresh,result_folder,filelist);
plot_data(data,neuron,dec_seg_ind,dec_num_plot,'auc',auc_result);
end

