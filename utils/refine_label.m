function neuron = refine_label(neuron, action_label, form_label, color, events, match_num)
u = zeros(1,length(action_label));
for i = 3:length(form_label)
    contain = strncmp(action_label,form_label{i},match_num);
    u = u+contain;
end
u = logical(u);
neuron.action_label = action_label(u);
neuron.color = color(u);
for i = 3:length(form_label)
    contain = strncmp(neuron.action_label,form_label{i},match_num);
    idx = find(contain,1);    
    neuron.events{idx+2} = events{i};
end
end