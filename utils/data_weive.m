function [data, neuron] = data_weive(data, neuron)
events = neuron.events;
event_label = neuron.event_label;
intruder_idx = neuron.intruder_idx;
Fs = neuron.Fs;
assert(length(events)==length(event_label), "Event label doesn't fit the annotation of intruder!");
delete_event = [];
for i=1:length(neuron.event_label)
    % find the event that is not relate to mating
    if contains(neuron.event_label{i},["bedding","pups"])
%         event_label(i) = [];
        blank_dur = events{intruder_idx}(i+1,1) - events{intruder_idx}(i,2);
        start_time = events{intruder_idx}(i,1);
        start_waive = round(start_time*Fs);
        % delete these event
        if blank_dur > 180
            end_time = events{intruder_idx}(i+1,1)-180;
            end_waive = round(end_time*Fs);
            data(:,start_waive:end_waive) = [];
        else
            end_time = events{intruder_idx}(i,2);
            end_waive = round(end_time*Fs);
            data(:,start_waive:end_waive) = [];
        end
        delete_event = [delete_event, i];
        delete_time = end_time - start_time;
        for k = 1:length(events)
                inside_idx = (events{k}(:,1) > start_time) & (events{k}(:,1) < end_time);
                beyond_idx = events{k}(:,1) > end_time;
                events{k}(beyond_idx,:) = events{k}(beyond_idx,:) - delete_time;
                events{k}(inside_idx,:) = [];
        end
        events{intruder_idx}(i,2) = events{intruder_idx}(i,1);
    end
end

neuron.events = events;
% neuron.event_label = event_label;
neuron.delete_event = delete_event;
neuron.nframe = size(data,2);
end

%         tRise = events{intruder_idx}(i,1);
%     
%         if blank_dur > 180
%             tDur = events{intruder_idx}(i+1,1) - 180 - tRise;
%             flag = revertTTL2bin(tRise, tDur, neuron.Fs, tlen);
%             data(:,~flag) = [];
%         else
%             tDur = events{intruder_idx}(i,2) - tRise;
%             flag = revertTTL2bin(tRise, tDur, neuron.Fs, tlen);
%             data(:,~flag) = [];
%         end
%         delete_event = [delete_event, i];
%         for k = 1:length(events)
%                 inside_idx = (events{k}(:,1) > tRise) & (events{k}(:,1) < tRise+tDur);
%                 beyond_idx = events{k}(:,1) > tRise+tDur;