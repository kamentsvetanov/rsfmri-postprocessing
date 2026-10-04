function Pout = kat_fmri_postprocessing_GLM_params(cfg)


% try smooth_filter   = cfg.smooth_filter;    catch smooth_filter = [];   end
% try f_brainmask     = cfg.f_brainmask;      catch f_brainmask = '';     end
% try TR              = cfg.TR;               catch error('Specify TR in Seconds'); end
%     
    
P_name = {'Name'        ...  1. Name for processing combination
         'detrend'      ...  2. Add detrend
         'regressCovs'  ...  3. Regress Covariates
         'regressGlobal'...  4. Regress Global signal
         'regressWM'    ...  5. Regress WM signal
         'regressCSF'   ...  6. Regress CSF signal
         'CompCor'      ...  7. Regress CompCorr Components (Behzadi et al 2007)
         'derivatives'  ...  8. Add derivatives
         'squaredTerms' ...  9. Add squared terms
         'motionOpt',   ... 10. Motion option 0 - No motion;  1- Power et al'13; 2 - Satthertwaite et al '12 (24 parameters in GLM); 3 - 6 RPs
         'filterMethod' ... 11. Add filtering 0 - No filtering;  1 - Do filtering part of regression (Hallquist et al 2013) | 2 - Butterworth filter, sequential regression
         'filterType'   ... 12. 1 - Band-pass; 2 - High-pass; 3 - Low-pass 
         'PreWhiten'    ... 13. PreWhiten
         'PartialFC'    ... 14. Do Partial FC, i.e. regress all other node timeseries
         'zscoreY'      ... 15. Zscore input EPI timeseries before processing
         'zscoreYr'     ... 16. Zscore output EPI timeseries after processing
         'zsoreCovs'    ... 17. Zscore covariates/all regressors
         'wantZval'     ... 18. Zval connectivity
         'trimVolumes'  ... 19. Trim N initial volumesn
         'addConstant'  ... 20. Add constatnt value to all voxels to hack 1st level SPM
         'saveSD'       ... 21. Compute temporal SD of timeseries (e.g. for QC) or RSFA maps
         'saveTS'       ... 22. Save processed timeseries (or rather not if interested only in RSFA maps)
         'saveFC'       ... 23. Compute Functional connectivity between ROIs
         'saveCVR'      ... 24. Compute rCVR as per Liu et al 2017 Neuroimage
         'HPC'          ... 25. High-pass filter
         'LPC'          ... 26. Low-pass filter
         };


%    1.              2     3   4   5   6    7    8    9    10       11       12   13   14      15       16      17         18       19      20      21     22       23   24   25    26                    
%   Name           Detr regCov GS WM  CSF  Cc  Deriv ^2  movOpt FiltMeth  FiltTyp PW  pFC   zscoreY zscoreYr zsoreCovs wantZval trimVol addConst saveSD saveTS saveFC saveCVR HPC  LPC   
P ={
 ...'m2f21_10_128'   1     1    0  0   0   0     1    1    2         2       1     0    1       0       0        1          0       5        0       1      0      0     1   .0078 .1 ;   ...5. High-pass 1 + Movement 2 + CSF + WM
 ...'m2f21_50_100'   1     1    0  0   0   0     1    1    2         2       1     0    1       0       0        1          0       5        0       1      0      0     1    .01  .02;   ...5. High-pass 1 + Movement 2 + CSF + WM
 ...'m2f21_25_50'    1     1    0  0   0   0     1    1    2         2       1     0    1       0       0        1          0       5        0       1      0      0     1    .02  .04;   ...5. High-pass 1 + Movement 2 + CSF + WM
    'wcm2f11_10_128' 1     1    0  1   1   0     1    1    2         1       1     0    1       0       0        1          0       5        0       0      0      0     0   .0078 .1 ;   ...5. High-pass 1 + Movement 2 + CSF + WM
 ...'wcm2f21_10_128' 1     1    0  1   1   0     1    1    2         2       1     0    1       0       0        1          0       5        0       1      0      0     1   .0078 .1 ;   ...5. High-pass 1 + Movement 2 + CSF + WM
 ...'wcm2f21_25_50'  1     1    0  1   1   0     1    1    2         2       1     0    1       0       0        1          0       5        0       1      0      0     1    .02  .04;   ...5. High-pass 1 + Movement 2 + CSF + WM
 ...'m2f21_12_25'    1     1    0  0   0   0     1    1    2         2       1     0    1       0       0        1          0       5        0       1      0      0     1    .04  .08;   ...5. High-pass 1 + Movement 2 + CSF + WM
 ...'m2f21_5_12'     1     1    0  0   0   0     1    1    2         2       1     0    1       0       0        1          0       5        0       1      0      0     1    .08  .2 ;   ...5. High-pass 1 + Movement 2 + CSF + WM
    };

P = cell2table(P,'VariableNames',P_name);

for pp = 1:size(P,1)%numel(postproc_name)
    Pout(pp).name               = P.Name{pp};
    Pout(pp).do                 = table2struct(P(pp,:));
    Pout(pp).LPC                = 1/P.LPC(pp);
    Pout(pp).HPC                = 1/P.HPC(pp);
%     S.process.postproc(pp).f_mask             = f_brainmask ;% Optional brain Mask used for 4D data 
%     S.process.postproc(pp).TR                 = TR; % Repetition time in seconds
%     S.process.postproc(pp).smooth             = smooth_filter;% Smoothing filters for 4D data
end

