% Wrapper script to post-process EPI data using physiological
% signals (e.g. WM/CSF signal, RETROICOR, COMPCORR) and head motion
% parameteres.
% Optional:
%  - compute RSFA maps 
%  - functional connectivity (for parcellated data only)
%  - process EPI data

parfor isub = 1:size(T,1)
    try
        funcdir     = fullfile(S.paths.genfi,'func/volumetric',T.SubID{isub});% Directory path to EPI data
        X           = [];
        X           = kat_fmri_postprocessing_GLM_params_template_fc(cfg);
        X.dataset   = 'rsfa_aroma';
        X.SubID     = T.SubID{isub}; % Subject ID
        X.f_in      = fullfile(funcdir,'func_mnispace.nii.gz'); % Filepath to pre-processed EPI data 
        X.f_rp      = fullfile(funcdir,'func.1D'); % Filepath to realignment parameters
        X.f_out     = fullfile(S.paths.genfi,'func/rsfa', T.SubID{isub},[X.name '_func_mnispace.nii']);% Directory where to save output 
        X.f_mask    = '/imaging/kt03/templates/masks/brain/brainmask_SPM_PT50_91x109x91.nii'; % Filepath to Brain mask image
        X.TR        = T.repetition_time(isub); % Repetition time in seconds
        X.smooth    = [8]; % Smoothing paramter in mm
        X.wm        = readmatrix(fullfile(funcdir,'func_WM_mean.txt')); % Filepath to file with WM timecourse. Otherwise, give in X.f_wm the fullpath to WM mask in the same resolution as EPI data 
        X.csf       = readmatrix(fullfile(funcdir,'func_CSF_mean.txt'));% Filepath to file with CSF timecourse. Otherwise, give in X.f_csf the fullpath to CSF mask in the same resolution as EPI data 
        X.compcor   = [readmatrix(fullfile(funcdir,'func_WM_pca.txt')) readmatrix(fullfile(funcdir,'func_CSF_pca.txt'))]; % Filepath to WM and CSF pca. if using X.f_wm and X.f_csf, leave blank
%             if ~exist(regexprep(X.f_out,'wcm2f11','cov_wcm2f11'),'file')
                kat_fmri_postprocessing_GLM(X);
%             end
        X=[];
    catch
        errmsg(isub) = 1;
    end
end
