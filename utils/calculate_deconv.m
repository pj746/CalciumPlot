function [dec_data, sig_data] = calculate_deconv(data,neuron,filefolder,filelist,save_or_not)
num_neuron = neuron.num_neuron;
nframe = neuron.nframe;
for i =1:length(filelist)
if contains(filelist{i},'_deconv.mat')
    load([filefolder,filesep,filelist{i}]);
    return
end
end
dec_data = zeros(num_neuron,nframe);
sig_data = zeros(num_neuron,nframe);
lambda = 1.5;
for i=1:num_neuron
    try
        [dec_trace, signal, ~] = deconvolveCa(data(i,:)', 'ar1', 'constrained', 'optimize_b', 'optimize_pars');
        dec_data(i,:) = dec_trace';
        sig_data(i,:) = signal';
    catch
        [dec_trace, signal, ~] = deconvolveCa(data(i,:)', 'ar1', 'foopsi', 'lambda', lambda, 'optimize_pars');
        dec_data(i,:) = dec_trace';
        sig_data(i,:) = signal';
    end
end
if save_or_not
    eval(sprintf('save %s%s%s_deconv.mat dec_data sig_data', filefolder, filesep, neuron.name));
end
end