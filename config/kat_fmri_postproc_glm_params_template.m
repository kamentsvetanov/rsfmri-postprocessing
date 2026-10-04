function Pout = kat_fmri_postprocessing_GLM_params(cfg)
% Template for parameters P to be used in
% fmri_postprocessing_GLM_params_template.m



% try smooth_filter   = cfg.smooth_filter;    catch smooth_filter = [];   end
% try f_brainmask     = cfg.f_brainmask;      catch f_brainmask = '';     end
% try TR              = cfg.TR;               catch error('Specify TR in Seconds'); end
%     

% First column defines the name of the option. Following columns take the
% values for these options. Every column is a setup for a new processing
% pipeline
P1    = {'Name','m2f21_10_128' 'test';...  1. Name for processing combination
         'detrend'          ,1 ,1 ;...  2. Add detrend
         'regressCovs'      ,1 ,1 ;...  3. Regress Covariates in 4D
         'regressGlobal'    ,0 ,1 ;...  4. Regress Global signal
         'regressWM'        ,1 ,1 ;...  5. Regress WM signal
         'regressCSF'       ,1 ,1 ;...  6. Regress CSF signal
         'regressCompCor'   ,1 ,1 ;...  7. Regress CompCorr Components (Behzadi et al 2007)
         'derivatives'      ,1 ,1 ;...  8. Add derivatives
         'squaredTerms'     ,1 ,1 ;...  9. Add squared terms
         'motionOpt'        ,2 ,1 ;... 10. Motion option 0 - No motion;  1- Power et al'13; 2 - Satthertwaite et al '12 (24 parameters in GLM); 3 - 6 RPs
         'filterMethod'     ,1 ,1 ;... 11. Add filtering 0 - No filtering;  1 - Do filtering part of regression (Hallquist et al 2013) | 2 - Butterworth filter, sequential regression
         'filterType'       ,1 ,1 ;... 12. 1 - Band-pass; 2 - High-pass; 3 - Low-pass 
         'PreWhiten'        ,0 ,1 ;... 13. PreWhiten
         'PartialFC'        ,0 ,1 ;... 14. Do Partial FC, i.e. regress all other node timeseries
         'zscoreY'          ,1 ,1 ;... 15. Zscore input EPI timeseries before processing
         'zscoreYr'         ,0 ,1 ;... 16. Zscore output EPI timeseries after processing
         'zsoreCovs'        ,1 ,1 ;... 17. Zscore covariates/all regressors
         'wantZval'         ,0 ,1 ;... 18. Zval connectivity
         'trimVolumes'      ,5 ,1 ;... 19. Trim N initial volumesn
         'addConstant'      ,0 ,1 ;... 20. Add constatnt value to all voxels to hack 1st level SPM
         'HPC'          ,.0078 ,1 ;... 21. High-pass filter
         'LPC'             ,.1 ,1 ;... 22. Low-pass filter
         'saveSD'           ,0 ,1 ;... 23. Compute temporal SD of timeseries (e.g. for QC) or RSFA maps
         'save4D'           ,0 ,1 ;... 24. Save processed 4D maps -timeseries (or rather not if interested only in RSFA maps)
         'saveCVR'          ,0 ,1 ;... 25. Compute rCVR as per Liu et al 2017 Neuroimage
         'saveFCnode2node'  ,0 ,1 ;... 26. Compute Functional connectivity between ROIs
         'saveFCnode2voxel' ,1 ,1 ;... 27. Compute Functional connectivity b/n ROIs and all voxels in the brain
         };

P = cell2table(P1(:,2:end)','VariableNames',P1(:,1));

% Exclude 'test' settings
P(ismember(P.Name,'test'),:) = [];

for pp = 1:size(P,1)%numel(postproc_name)
    Pout(pp).name               = P.Name{pp};
    Pout(pp).do                 = table2struct(P(pp,:));
    Pout(pp).LPC                = 1/P.LPC(pp);
    Pout(pp).HPC                = 1/P.HPC(pp);
%     S.process.postproc(pp).f_mask             = f_brainmask ;% Optional brain Mask used for 4D data 
%     S.process.postproc(pp).TR                 = TR; % Repetition time in seconds
%     S.process.postproc(pp).smooth             = smooth_filter;% Smoothing filters for 4D data
end

