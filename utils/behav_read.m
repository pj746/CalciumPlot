function [neuron,form_label, events] = behav_read(slx_path)
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
neuron.intruder_idx = find(ismember(form_label,'intruder'));
neuron.events = events;
end