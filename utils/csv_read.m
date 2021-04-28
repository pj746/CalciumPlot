function [neuron, trace] = csv_read(fname)
% 用于读取由inscopix识别神经元活动产生的csv文件
% Only the neuron that annotated as 'accepted' would be counted in the
% trace.
% neuron; a list of id that annotated as 'accepted'
% trace = [N+1, I], where N = num(neuron), I is the intensity of the
% neuronal activity. First row is the recording time.
% 注意：转完之后第一行是时间，后面为神经元的trace
T = readtable(fname);
vname = char(T.Properties.VariableNames);
trace = str2num(char(T.(vname(1,:))(2:length(T.(vname(1,:))))))';  % time, x axis
neuron = [0];
for i=2:size(T.Properties.VariableNames,2)
    if strcmp(char(T.(strtrim(vname(i,:)))(1)),'accepted')
        trace = [trace;str2num(char(T.(strtrim(vname(i,:)))(2:length(T.(strtrim(vname(i,:)))))))'];
        % 另外还要存储神经元的编号
        neuron = [neuron,str2num(strtrim(vname(i,2:size(vname(i,:),2))))];
    end
num_neuron = size(neuron, 2);
end
neuron = neuron(2:num_neuron);