function kat_fmri_postprocessing_GLM_relCVR(SS,aY,V)
% Auxiliary function to fmri_processing_GLM to generate relative CVR maps
% Compute and save relative CVR maps as in Liu et al 2017 Neuroimage
% 

Ns = SS.Ns;
idxMaskBrain = logical(spm_read_vols(spm_vol(SS.f_mask))); % Read grey matter reference mask to derive reference signal time course
idxMaskGM = logical(spm_read_vols(spm_vol(SS.f_maskGM))); % Read grey matter reference mask to derive reference signal time course
dims      = size(idxMaskBrain);
idxMaskGMvec = logical(reshape(idxMaskBrain,prod(dims),1));
yBrain = aY(:,idxMaskGMvec);
yGM    = aY(:,idxMaskGM);
[~,scoreGM, ~, ~, explainedGM] = pca(yGM);
[Tval,df,Bval]   = spm_ancova(zscore(scoreGM(:,1)),speye(Ns,Ns),zscore(yBrain),[1]);
Tval3d = nan(dims);
Bval3d = nan(dims);
Tval3d(idxMaskBrain) = Tval;
Bval3d(idxMaskBrain) = atanh(Bval);


% Use whole-brain mask images
%     idxMaskBrain = spm_read_vols(spm_vol(SS.f_mask)); % Read grey matter reference mask to derive reference signal time course
%     Tval = Tval.*idxMaskGM;
%     
[d_out,f_out] = fileparts(SS.f_out);
V(1).fname = fullfile(d_out,['aCVR_' f_out '_tvals.nii']);
spm_write_vol(V(1),Tval3d); 
smoothFilt = SS.smooth;
for iSmooth = 1:numel(smoothFilt)
    sfilt = smoothFilt(iSmooth);
    f_new = regexprep(V(1).fname,'.nii',sprintf('_s%d.nii',sfilt));
    spm_smooth(V(1).fname,f_new,[sfilt sfilt sfilt]);
    % Mask again
    Vtemp = spm_vol(f_new);
    Y = spm_read_vols(Vtemp);
    Ymasked = Y.*idxMaskBrain;
    spm_write_vol(Vtemp,Ymasked);
end

V(1).fname = fullfile(d_out,['aCVR_' f_out '_betas.nii']);
spm_write_vol(V(1),Bval3d); 
for iSmooth = 1:numel(smoothFilt)
    sfilt = smoothFilt(iSmooth);
    f_new = regexprep(V(1).fname,'.nii',sprintf('_s%d.nii',sfilt));
    spm_smooth(V(1).fname,f_new,[sfilt sfilt sfilt]);
    % Mask again
    Vtemp = spm_vol(f_new);
    Y = spm_read_vols(Vtemp);
    Ymasked = Y.*idxMaskBrain;
    spm_write_vol(Vtemp,Ymasked);
end

% -----------------
% Relative CVR
% -----------------
Tval3d = nan(dims);
Bval3d = nan(dims);
Tval3d(idxMaskBrain) = Tval/mean(Tval);
Bval3d(idxMaskBrain) = Bval/mean(Bval);

V(1).fname = fullfile(d_out,['rCVR_' f_out '_tvals.nii']);
spm_write_vol(V(1),Tval3d); 
for iSmooth = 1:numel(smoothFilt)
    sfilt = smoothFilt(iSmooth);
    f_new = regexprep(V(1).fname,'.nii',sprintf('_s%d.nii',sfilt));
    spm_smooth(V(1).fname,f_new,[sfilt sfilt sfilt]);
    % Mask again
    Vtemp = spm_vol(f_new);
    Y = spm_read_vols(Vtemp);
    Ymasked = Y.*idxMaskBrain;
    spm_write_vol(Vtemp,Ymasked);
end

V(1).fname = fullfile(d_out,['rCVR_' f_out '_betas.nii']);
spm_write_vol(V(1),Bval3d); 
for iSmooth = 1:numel(smoothFilt)
    sfilt = smoothFilt(iSmooth);
    f_new = regexprep(V(1).fname,'.nii',sprintf('_s%d.nii',sfilt));
    spm_smooth(V(1).fname,f_new,[sfilt sfilt sfilt]);
    % Mask again
    Vtemp = spm_vol(f_new);
    Y = spm_read_vols(Vtemp);
    Ymasked = Y.*idxMaskBrain;
    spm_write_vol(Vtemp,Ymasked);
end