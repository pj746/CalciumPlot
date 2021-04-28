function neuron = waive_overlap_behav(neuron, form_label, match_num)
% The function is to used for refining the wrong annotations that makes the
% behaviors overlap.
Fs = neuron.Fs;
action_label = neuron.action_label;
events = neuron.events;
nframe = neuron.nframe;
color = neuron.color;

neuron = refine_label(neuron, action_label, form_label, color, events, match_num);
assert(length(neuron.action_label)+2 == length(neuron.events),'Failed to refine label and events!');
nstate = length(neuron.action_label);
bstates = zeros(1,nframe);
tlen = nframe / Fs;  % total time
for i = 1:nstate
    event_now = neuron.events{i+2};
    tRise = event_now(:,1); tDur = event_now(:,2) - tRise;
    flag = revertTTL2bin(tRise,tDur, Fs, tlen);
    bstates(flag) = i;
end
events{1} = neuron.events{1};
events{2} = neuron.events{2};
for i = 1:nstate
    behav_state = (bstates==i);
    [tRise, tDur] = detectTTL(behav_state);
    tEnd = tRise + tDur;
    event = tRise/Fs;
    event(:,2) = tEnd/Fs;
    events{i+2} = event;
end
neuron.events = events;
end