function [ind, seg_ind] = thresh_order(data, neuron, thresh, trace_per_plot)
% Sort the trace at the time when it's normalized amplitude first reaches thresh.
% ind: n x 2 matrix. First row: time that first reached the thresh; second row: neuron id
% order from smallest to biggest
data_norm = normalize(data,2,'range');  % normalize the data to [0,1]
ind = zeros(neuron.num_neuron,2);
for i=1:neuron.num_neuron
    ind(i,1) = find(data_norm(i,:)>thresh,1);
    ind(i,2)=i;
end
ind = sortrows(ind);
num_plot0 = floor(neuron.num_neuron/trace_per_plot);
last_num = mod(neuron.num_neuron,trace_per_plot);
if last_num==0 last_num = []; end
mtr_seg = ones(1,num_plot0)*trace_per_plot;
mtr_seg = [mtr_seg, last_num];
seg_ind = mat2cell(ind,mtr_seg);