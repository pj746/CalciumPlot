addpath('E:\wupeixuan\CalImgProcess\CalciumPlot\utils');

%  a. ����Inscopix���ҵ����ʵ���Ԫ��������С��һ
%  b. neuron_id��Ӧ��deconv��raw�Ĵ���
%  b. д��raw��deconv
%  c. ʱ���߸�raw��ʱ���߶��룺��ͷ�ӵ�β����հף����Ϊraw�ĵ�0�͵�1��ʱ��ļ��
%  Esr290222 neuron 061 136
Fs = 10;
neuron = [61, 136];
result_folder = 'E:\wupeixuan\auc_plot\data\Esr29m\Esr290222';
filelist = dir(result_folder);
filelist = {filelist.name};

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
%% Read the traces
[rpath, rname, rtype] = fileparts(trace_dir);
rname_ind = strfind(rname,'_');
neuron.name = rname(1:rname_ind(1)-1);
load(trace_dir);
num_neuron = size(neuron_id,2);
data = trace(2:num_neuron+1,:);
neuron.neuron_id = neuron_id;
[neuron.num_neuron, neuron.nframe] = size(data);

various={'Time','61_raw','61_deconv','136_raw','136_deconv'};
time0 = trace(1,:)';
trace_table=table(time0,trace1',c_cvx0',c_cvx,dec_trace6,dec_trace7,dec_trace1_0,dec_trace1,dec_trace2,dec_trace8,dec_trace9,'VariableNames',various);
writetable(trace_table, '\\liying.cibr.ac.cn\public\Code\auROC\wpx_code\traces_plot_with_events\Esr29\Esr290201(no behave annnot)\neuron136_deconv.csv');
