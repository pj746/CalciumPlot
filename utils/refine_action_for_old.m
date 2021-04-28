function neuron = refine_action_for_old(neuron)
% The function is used for refining the cues for the old data format.
% Originally, the behaviors are copy-pasted from the .xlsx to the .m code
% to plot the traces. For the behaviors that doesn't appear in this day,
% researcher uses [0 0.1] to mark as omitting action_label. This file is
% specifically used to adapt the old formate to the newer, tidier way to plot
% the data.

action_label = neuron.action_label;
color = neuron.color;
events = neuron.events;
nstate = size(action_label,2);
minu=0;
for i = 1:nstate
    if events{i+2-minu}(1) == 0
        action_label(i-minu) = [];
        color(i-minu) = [];
        events(i+2-minu) = [];
        minu=minu+1;
        nstate = nstate-1;
        continue
    end
end
neuron.action_label = action_label;
neuron.color = color;
neuron.events = events;
end