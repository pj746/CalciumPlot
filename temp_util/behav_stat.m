St1810108_init;
% TODO 打印refined_action和intruder,手动输入cue和行为
load('St1810108_auc+blank_0.70.mat');
h_signifi = auc_result.h_signifi;
% idx = [1,3;4,9];  % 第一列: action，第二列: intruder
idx = [4,9];
nstat = size(idx,1);
sum_piece = 'sum(h_signifi(:,idx(1,1),idx(1,2))';
ans_piece = sprintf('Neurons to %s of %s',neuron.action_label{idx(1,1)},neuron.event_label{idx(1,2)});
if nstat>1
    for i=2:nstat
        sum_piece =[sum_piece,sprintf('&h_signifi(:,idx(%d,1),idx(%d,2))',i,i)];
        ans_piece = [ans_piece, sprintf(' and %s of %s',neuron.action_label{idx(i,1)}),neuron.event_label{idx(i,2)}];
    end
end
sum_piece =[sum_piece,')'];
ans_piece = [ans_piece, ' = %.2f%%'];
stat_neuron_sum = eval(sum_piece);
sprintf(ans_piece,100*stat_neuron_sum/num_neuron)