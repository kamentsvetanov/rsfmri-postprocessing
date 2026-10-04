function kat_fmri_postprocessing_GLM_rsfa(SS, aY,V)
% Auxiliary function to fmri_processing_GLM to generate RSFA maps
%
%
% 

[dirOut,fn] = fileparts(SS.f_out); 
mkdir(dirOut);
Ym = spm_read_vols(spm_vol(SS.f_mask));
Ysd = std(aY);
Ysd3d = reshape(Ysd,V(1).dim);
Ysd3d = (Ysd3d.*Ym);
V(1).fname = regexprep(SS.f_out,fn,['sd_' fn]);
V(1).dt = [16 0];
spm_write_vol(V(1),Ysd3d);

% Apply smoothing
try smoothFilt = SS.smooth; catch smoothFilt = []; end
for iSmooth = 1:numel(smoothFilt)
    sfilt = smoothFilt(iSmooth);
    f_new = regexprep(V(1).fname,'.nii',sprintf('_s%d.nii',sfilt));
    spm_smooth(V(1).fname,f_new,[sfilt sfilt sfilt]);
    % Mask again
    Vtemp = spm_vol(f_new);
    Y = spm_read_vols(Vtemp);
    Ymasked = Y.*Ym;
    spm_write_vol(Vtemp,Ymasked);
end

% -------------------------------------------
% Normalised SD, i.e. Coefficient of Variation
% -------------------------------------------
SDglobal = mean(Ysd3d(logical(Ym)));
Ysd3Dnorm = Ysd3d./SDglobal;
V(1).fname = regexprep(SS.f_out,fn,['cov_' fn]);
spm_write_vol(V(1),Ysd3Dnorm);

for iSmooth = 1:numel(smoothFilt)
    sfilt = smoothFilt(iSmooth);
    f_new = regexprep(V(1).fname,'.nii',sprintf('_s%d.nii',sfilt));
    spm_smooth(V(1).fname,f_new,[sfilt sfilt sfilt]);
    % Mask again
    Vtemp = spm_vol(f_new);
    Y = spm_read_vols(Vtemp);
    Ymasked = Y.*Ym;
    spm_write_vol(Vtemp,Ymasked); 
end


% % ---------------------------------------------------
% % save mean signal (not SD) for diferent tissue types
% % ---------------------------------------------------
% Ymean = mean(aY);
% 
% threshold   = .8;    
% dir_tissues = '/imaging/camcan/sandbox/kt03/templates/masks/brain/spm_fieldmap';
% f_csf       = fullfile(dir_tissues,'csf_61x73x61.nii');
% f_wm        = fullfile(dir_tissues,'white_61x73x61.nii');
% f_gm        = fullfile(dir_tissues,'grey_61x73x61.nii');
% numComp     = 5; % Number of PCs to extract
% Ref         = { 'gm'  f_gm;...
%                 'wm'  f_wm;...
%                 'csf' f_csf};
% compSign = nan(numComp+1,size(Ref,1));
% for iDat = 1:length(Ref)
%     mask = find(spm_read_vols(spm_vol(Ref{iDat,2}))>threshold);
%     dat  = aY(:,mask);
%     [coeff, score, latent, tsquared, explained] = pca(dat');
% 
%     compSign(1:numComp,iDat) = mean(score(:,1:numComp)); 
%     compSign(end,iDat) = mean(dat(:));
% end
% compSign = array2table(compSign);
% compSign.Properties.VariableNames = Ref(:,1);
% compSign.Properties.RowNames(1:numComp)= cellstr(strcat('PC', num2str([1:numComp]')));
% compSign.Properties.RowNames(end)      = {'mean'};
% f_out = spm_file(SS.f_out,'prefix','compSign_','ext','mat');
% save(f_out,'compSign');