% 来自：demo_endoscope

neuron_full.options.deconv_flag = true; 
neuron_full.options.deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'foopsi', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'smin', -5, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', true, ...% optimize the baseline);
    'max_tau', 100);    % maximum decay time (unit: frame);

% 来自 updateTemporal_endoscope
% deconvolution
        if obj.options.deconv_flag
            [ck, sk, deconv_options]= deconvolveCa(temp, deconv_options_0, 'maxIter', 2, 'sn', tmp_sn);
            smin(k) = deconv_options.smin;
            kernel_pars{k} = reshape(deconv_options.pars, 1, []);
            temp = temp - deconv_options.b; 
        else
            ck = max(0, temp);
        end