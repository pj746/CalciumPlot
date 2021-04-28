function [ind, seg_ind, num_plot] = double_sort(data,neuron,thresh, first_event, second_event)
    % Sort the trace  by two event. a. Find the traces that reaches the thresh in the first event and sort them;
    % b. Find traces that reaches the thresh in the second event among the rest traces and sort them.
    % c. Sort the rest traces.
    % inputs:
    % data: n x T matrix
    % thresh: double
    % first_event & second event: integers that indicate the event.

num_neuron = neuron.num_neuron;
data_norm = normalize(data,2,'range');  % normalize the data to [0,1]
intruder = cell2mat(neuron.events(neuron.intruder_idx));
ind1 = [0,0];
ind2 = [0,0];
ind3 = [0,0];
for i=1:num_neuron
%     ind(i,1) = find(data_norm(i,:)>thresh,1);
    frame_idx = find(data_norm(i,:)>thresh);
    time_idx = frame_idx/neuron.Fs;
    first_ind = time_idx>intruder(first_event,1) & time_idx<intruder(first_event,2);
    second_ind = time_idx>intruder(second_event,1) & time_idx<intruder(second_event,2);
    
    if sum(first_ind)>0
        ind1 = [ind1;frame_idx(find(first_ind,1)),i];
    elseif sum(second_ind)>0
        ind2 = [ind2;frame_idx(find(second_ind,1)),i];
    else
        ind3 = [ind3;find(data_norm(i,:)>thresh,1),i];
    end
end
try
    ind1 = ind1(2:size(ind1,1),:);
    ind1 = sortrows(ind1);
catch
    ind1 = [];
end
try
    ind2 = ind2(2:size(ind2,1),:);
    ind2 = sortrows(ind2);
catch
    ind2 = [];
end
try
    ind3 = ind3(2:size(ind3,1),:);
    ind3 = sortrows(ind3);
catch
    ind3 = [];
end
% ind2 = ind2(2,2),:);
% ind3 = ind3(2:size(ind3,2),:);
% ind2 = sortrows(ind2);
% ind3 = sortrows(ind3);
ind = [ind1;ind2;ind3];
seg_ind = {ind1,ind2,ind3};