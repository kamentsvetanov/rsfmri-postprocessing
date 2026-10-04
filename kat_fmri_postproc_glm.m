function [SS]  = kat_fmri_postprocessing_GLM(SS)

% function [SS] = kat_fmri_postprocessing_GLM(SS);
% % Originally based on GLM code by Rik Henson (MRC CBU):
% https://github.com/MRC-CBU/riksneurotools/blob/master/GLM/glm.m
%
% Subsequently extended to more flexible nuisance regression, temporal filtering, RSFA and functional
% connectivity analyses.
%
% Function (using SPM8 functions) for estimating linear regressions
% between fMRI timeseries in each pair of Nr ROIs, adjusting for bandpass
% filter, confounding timeseries (eg CSF) and (SVD of) various expansions of
% movement parameters, and properly modelling dfs based on comprehensive
% model of error autocorrelation.

%
% S.Y = [Ns x Nr] data matrix, where Ns = number of scans (ie resting-state fMRI timeseries) and Nr = number of ROIs
% S.M = [Ns x 6] matrix of 6 movement parameters from realignment (x,y,z,pitch,roll,yaw)
% S.C = [Ns x Nc] matrix of confounding timeseries, Nc = number of confounds (eg extracted from WM, CSF or Global masks)
% S.G = [0/1] - whether to include global over ROIs (ie given current data) (default = 0)
% S.TR = time between volumes (TR), in seconds
%
% (S.HPC (default = 100) = highpass cut-off (in seconds)
% (S.LPC (default = 10) = lowpass cut-off (in seconds)
% (S.CY (default = {}) = precomputed covariance over pooled voxels (optional)
% (S.pflag (default=0) is whether to calculate partial regressions too (takes longer))
% (S.svd_thr (default=.99) is threshold for SVD of confounds)
% (S.SpikeMovAbsThr (default='', ie none) = absolute threshold for outliers based on RMS of Difference of Translations
% (S.SpikeMovRelThr (default='', ie none) = relative threshold (in SDs) for outliers based on RMS of Difference of Translations or Rotations
% (S.SpikeDatRelThr (default='', ie none) = relative threshold (in SDs) for mean of Y over voxels (outliers, though arbitrary?)
% (S.SpikeLag (default=1) = how many TRs after a spike are modelled out as separate regressors
% (S.StandardiseY (default = 0) = whether to Z-score each ROI's timeseries to get standardised Betas)
% (S.PreWhiten (default = 1) = whether to estimate autocorrelation of error and prewhiten (slower, but better Z-values (less important for Betas?))
% (S.optmotion (default = 1) = whether to use the Satterthwaite motion confounds (optmotion==2) 
%                              or Rik's Volterra expansion (optmotion==1) 
% (S.WantZval = whether Z-values need to be computed

% Zmat  = [Nr x Nr] matrix of Z-statistics for linear regression from seed (row) to target (column) ROI
% Bmat  = [Nr x Nr] matrix of betas for linear regression from seed (row) to target (column) ROI
% pZmat = [Nr x Nr] matrix of Z-statistics for partial linear regression from seed (row) to target (column) ROI
% pZmat = [Nr x Nr] matrix of betas for partial linear regression from seed (row) to target (column) ROI
% aY    = data adjusted for confounds
% X0r   = confounds (filter, confound ROIs, movement expansion)
%
% Note:
%    some Inf Z-values can be returned in Zmat (when p-value so low than inv_Ncdf is infinite)
%       (these could be replaced by the maximum Z-value in matrix?)
%
% Potential improvements:
%   regularised regression (particularly for pZmat), eg L1 with LASSO - or L2 implementable with spm_reml_sc?

% Script edit log (kat)
% 09/10/2014 - can take 4D images as input (paths to nii)
% 20/06/2015 - converted into function with a single structure as input
% 23/11/2016 - Option to remove initial/dummy volumes
% 12/03/2017 - Additional filtering option for Low/High/Band-pass filtering 
% 15/05/2018 - Added CompCor
% 23/10/2019 - Integration with estimation of RSFA and relCVR maps
% 25/01/2021 - Moved FC computaion to a separate function
% 15/05/2024 - Added Seed-based connectivity analysis with all voxels
% 
% STRUCTURE TEMPLATE
% P     = {'Name' ,'m2f21_10_128';...  1. Name for processing combination
%          'detrend'          ,1 ;...  2. Add detrend
%          'regressCovs'      ,2 ;...  3. Regress Covariates
%          'regressGlobal'    ,1 ;...  4. Regress Global signal
%          'regressWM'        ,1 ;...  5. Regress WM signal
%          'regressCSF'       ,1 ;...  6. Regress CSF signal
%          'regressCompCor'   ,1 ;...  7. Regress CompCorr Components (Behzadi et al 2007)
%          'derivatives'      ,1 ;...  8. Add derivatives
%          'squaredTerms'     ,1 ;...  9. Add squared terms
%          'motionOpt'        ,2 ;... 10. Motion option 0 - No motion;  1- Power et al'13; 2 - Satthertwaite et al '12 (24 parameters in GLM); 3 - 6 RPs
%          'filterMethod'     ,1 ;... 11. Add filtering 0 - No filtering;  1 - Do filtering part of regression (Hallquist et al 2013) | 2 - Butterworth filter, sequential regression
%          'filterType'       ,1 ;... 12. 1 - Band-pass; 2 - High-pass; 3 - Low-pass 
%          'PreWhiten'        ,0 ;... 13. PreWhiten
%          'PartialFC'        ,0 ;... 14. Do Partial FC, i.e. regress all other node timeseries
%          'zscoreY'          ,1 ;... 15. Zscore input EPI timeseries before processing
%          'zscoreYr'         ,0 ;... 16. Zscore output EPI timeseries after processing
%          'zsoreCovs'        ,1 ;... 17. Zscore covariates/all regressors
%          'wantZval'         ,0 ;... 18. Zval connectivity
%          'trimVolumes'      ,5 ;... 19. Trim N initial volumesn
%          'addConstant'      ,0 ;... 20. Add constatnt value to all voxels to hack 1st level SPM
%          'HPC'          ,.0078 ;... 21. High-pass filter
%          'LPC'             ,.1 ;... 22. Low-pass filter
%          'saveSD'           ,0 ;... 23. Compute temporal SD of timeseries (e.g. for QC) or RSFA maps
%          'save4D'           ,0 ;... 24. Save processed 4D timeseries (or rather not if interested only in RSFA maps)
%          'saveCVR'          ,0 ;... 25. Compute rCVR as per Liu et al 2017 Neuroimage
%          'saveFCnode2node'  ,0 ;... 26. Compute Functional connectivity between ROIs
%          'saveFCnode2voxel' ,1 ;... 27. Compute Functional connectivity b/n ROIs and all voxels in the brain
%          };
%
% --------------
% TO DO:
% - introduce a variable keeping track of the names for regressors in X0r




try fprintf('Analysis name: %s  \n',SS.name); end

try Y                   = SS.Y;                 catch Y                     = [];end
try CY                  = SS.CY;                catch CY                    = [];end % Covariance over ROIs 
try M                   = SS.M;                 catch M                     = [];end % Read movement parameters from rp.txt file
try C                   = SS.C;                 catch C                     = []; warning('No additional covariates considered.'); end
try doZval              = SS.do.Zval;           catch doZval                = 0; end
try doZscoreY           = SS.do.zscoreY;        catch doZscoreY             = 0; end
try doZscoreYr          = SS.do.zscoreYr;       catch doZscoreYr            = 0; end

try doRegressCovariates = SS.do.regressCovs;    catch doRegressCovariates   = 1; end
try doRegressGlobal     = SS.do.regressGlobal;  catch doRegressGlobal       = 0; end
try doRegressWM         = SS.do.regressWM;      catch doRegressWM           = 0; end 
try doRegressCSF        = SS.do.regressCSF;     catch doRegressCSF          = 0; end
try doRegressCompCor    = SS.do.regressCompCor; catch doRegressCompCor      = 0; end

try doDerivatives       = SS.do.derivatives;    catch doDerivatives         = 0; end
try doSquaredTerms      = SS.do.squaredTerms;   catch doSquaredTerms        = 0; end
try doDetrend           = SS.do.detrend;        catch doDetrend             = 0; end% Detrend linear and quadratic terms
try doZsoreCovariates   = SS.do.zsoreCovs;      catch doZsoreCovariates     = 0; end

try doMotionOpt         = SS.do.motionOpt;      catch doMotionOpt           = 0; end
try SpikeMovRelThr      = SS.SpikeMovRelThr;    catch SpikeMovRelThr        = '';end % 5 SDs of mean?
try SpikeDatRelThr      = SS.SpikeDatRelThr;    catch SpikeDatRelThr        = '';end
try SpikeMovAbsThr      = SS.SpikeMovAbsThr;    catch SpikeMovAbsThr        = '';end % SpikeMovAbsThr = 0.25;  % mm from Satterthwaite et al 2013
try SpikeLag            = SS.SpikeLag;          catch SpikeLag              = 1; end % SpikeLag = 5;     % 5 TRs after spike?

try doFilterMethod      = SS.do.filterMethod;   catch doFilterMethod        = 0; end % 1 - Power et al; 2 - Satterthwaite et al
try doFilterType        = SS.do.filterType;     catch doFilterType          = 1; end % 1 - Band-pass; 2 - High-pass 
try HPC                 = SS.HPC;               catch HPC                   = 1/0.01; warning('Assuming highpass cut-off of %d',HPC); end
try LPC                 = SS.LPC;               catch LPC                   = 1/0.2;  warning('Assuming lowpass cut-off of %d',LPC); end

try doSaveSD            = SS.do.saveSD;         catch doSaveSD              = 0; end
try doSave4D            = SS.do.save4D;         catch doSave4D              = 0; end
try doSaveFCnode2node   = SS.do.saveFCnode2node;catch doSaveFCnode2node     = 0; end
try doSaveFCnode2voxel  = SS.do.saveFCnode2voxel;catch doSaveFCnode2voxel   = 0; end
try doSaveCVR           = SS.do.saveCVR;        catch doSaveCVR             = 0; end

try TR                  = SS.TR;                catch error('N.B. Need TR (in seconds).'); end % TR in seconds
try doTrimVolumes       = SS.do.trimVolumes;    catch doTrimVolumes         = 0; end
try doAddConstant       = SS.do.addConstant;    catch addConstant           = 0; end

try SS.svd_thr;           catch SS.svd_thr              =.99;end
% try SS.do.PreWhiten;      catch SS.do.PreWhiten         = 0; end
% try SS.do.wantZval;       catch SS.do.wantZval          = 1; end
% try SS.do.PartialFC;      catch SS.do.PartialFC         = 0; end

% Generate variable with names of the Regressors Types . Starting with
% additional Covariates, e.g. RETROICOR regressors
if ~isempty(C)
    nameReg =  string(strcat('AdditionalCov', num2str([1:size(C,2)]','%03d')));
else
    nameReg = string();
end
% ------------------------------------------
% Input timeseries and transforation options
% ------------------------------------------
if isempty(Y) || ischar(Y)
        % Read in image file, concatenating the 4D image into a compressed 2D voxel x time matrix
%     [Y1, Info] = ParseInNii(SS.f_in, 'compress', 0); 
    V = spm_vol(SS.f_in);
    Y = spm_read_vols(V);
    Y = permute(Y,[4 1 2 3]);
    Y = Y(:,:);
end

% ----------------------
% Motion and its options
% ----------------------
if isempty(M) 
    try
        M = load(SS.f_rp); 
    catch
        disp('N.B. Need Nscans x 6 movement parameter matrix (Ignore of not motion regression is required)'); 
        M = [];
    end
end

if doMotionOpt == 0
    doRegressMotion = 0;
else
    doRegressMotion = 1;
end

% ---------------------------
% Covariates and their orders
% ---------------------------

% -------------------------------------------------------------------------
% If Global Signal timecourse is not provided try to load it from a file
% -------------------------------------------------------------------------
if ~isfield(SS,'global') && doRegressGlobal
    % from ROI file containing the signal
    Global   = load(SS.f_global); Global = Global.ROI;
    SS.global = Global(1).mean;
    
    % From a mask
%     mask = logical(ParseInNii(SS.f_mask, 'compress', 0));
    Ym = logical(spm_read_vols(spm_vol(SS.f_mask)));
    SS.global = mean(Y(:,Ym(:)),2);
end

% -------------------------------------------------------------------------
% if WM signal is not provided as a vector in SS.wm, try
%  - to load wm mask and extract signal from EPI image
%  or
%  - to load filepath to file containing wm signal (accoridng to AA's output)
% -------------------------------------------------------------------------
if ~isfield(SS,'wm') && doRegressWM
    [~,~,ext] = fileparts(SS.f_wm); % Get extension to determine where 
    if strcmp(ext,'.nii')
        WMmask = logical(spm_read_vols(spm_vol(SS.f_wm)));
        WM(1).rawdata = Y(:,WMmask(:));
        SS.wm = mean(WM(1).rawdata,2); 
    else 
        WM   = load(SS.f_wm); WM = WM.ROI;
        SS.wm = WM(1).mean;
    end
end


% if CSF signal is not provided as a vector in SS.csf, try to
%  - to load CSF mask and extract signal from EPI image
%  or
%  - to load filepath to file containing CSF signal (accoridng to AA's output)
if ~isfield(SS,'csf') && doRegressCSF
    [~,~,ext] = fileparts(SS.f_csf); % Get extension to determine where 
    if strcmp(ext,'.nii')
        CSFmask = logical(spm_read_vols(spm_vol(SS.f_csf)));
        CSF(1).rawdata = Y(:,CSFmask(:));
        SS.csf = mean(CSF(1).rawdata,2);
    else 
        CSF   = load(SS.f_csf); CSF = CSF.ROI;
        SS.csf = CSF(1).mean;
    end
end

% if CompCor signal is not provided as a vector in SS.compcor, try to load
% filepath to file containing CSF signal (accoridng to AA's output)
if ~isfield(SS,'compcor') && doRegressCompCor
    try % Using data from filepaths to ROIs
        [~,WMscore, ~, ~, explainedWM] = pca(WM(1).rawdata);
        [~,CSFscore, ~, ~, explainedCSF] = pca(CSF(1).rawdata);
    catch % 
        [~,WMscore, ~, ~, explainedWM] = pca(SS.wmRawdata);
        [~,CSFscore, ~, ~, explainedCSF] = pca(SS.csfRawdata);
    end
    CompCor = [WMscore(:,1:5) CSFscore(:,1:5)];
elseif isfield(SS,'compcor') && doRegressCompCor
    CompCor = SS.compcor; 
end
clear WM CSF;

% ------------------------
% Filtering options
% ------------------------
if doFilterMethod == 0   
    doFilter = 0;
else
    doFilter = 1;
end

%% If want to try Matlab's LASSO (takes ages though)
% lassoflag = 0;
% if lassoflag
%     matlabpool open
%     opts = statset('UseParallel','always');
% end

if doTrimVolumes 
    trim_initvols = doTrimVolumes;
    Y = Y(trim_initvols+1:end,:);
    try M           = M(trim_initvols+1:end,:);         catch warning('Attempted to trim RPs'); end
    try C           = C(trim_initvols+1:end,:);         catch warning('Attempted to trim Covs');end
    try SS.wm       = SS.wm(trim_initvols+1:end,:);     catch warning('Attempted to trim WM');end
    try SS.csf      = SS.csf(trim_initvols+1:end,:);    catch warning('Attempted to trim CSF'); end
    try SS.gm       = SS.gm(trim_initvols+1:end,:);     catch warning('Attempted to trim GM');end
    try SS.global   = SS.global(trim_initvols+1:end,:); catch warning('Attempted to trim GM');end
    try CompCor     = CompCor(trim_initvols+1:end,:);   catch warning('Attempted to trim CompCor');end
end

Ns = size(Y,1);
Nr = size(Y,2);

% In case Y's first volumes have already been trimmed, but design matrix
% and other regressors not
if Ns~=size(M,1)
    disp('NB: Trimming first datapoints in covariates to match Y');
    trim_initvols = size(M,1)-Ns;
    try M           = M(trim_initvols+1:end,:);         catch warning('Attempted to trim RPs');end
    try C           = C(trim_initvols+1:end,:);         catch warning('Attempted to trim Covs'); end
    try SS.wm       = SS.wm(trim_initvols+1:end);       catch warning('Attempted to trim WM'); end
    try SS.csf      = SS.csf(trim_initvols+1:end,:);    catch warning('Attempted to trim CSF');end
    try SS.gm       = SS.gm(trim_initvols+1:end,:);     catch warning('Attempted to trim GM');end
    try SS.global   = SS.global(trim_initvols+1:end,:); catch warning('Attempted to trim GM');end
    try CompCor     = CompCor(trim_initvols+1:end,:);   catch warning('Attempted to trim CompCor');end
end


%% Create a DCT bandpass filter (so filtering part of model, countering Hallquist et al 2013 Neuroimage)
if doFilter
    if doFilterMethod == 1
        K   = spm_dctmtx(Ns,Ns);
        nHP = fix(2*(Ns*TR)/HPC + 1);
        nLP = fix(2*(Ns*TR)/LPC + 1);
        
        if doFilterType == 1 % Band-pass filter
            K   = K(:,[2:nHP nLP:Ns]);      % Remove initial constant
        elseif doFilterType == 2 
            K   = K(:,[2:nHP]);   %Only high-pass filter [2:nHP ]);      % Remove initial constant
        elseif doFilterType == 3
            K   = K(:,[nLP:Ns]);  %Only low pass filter [2:nHP ]);      % Remove initial constant
        end
        Nk  = size(K,2);
   
    elseif doFilterMethod == 2
        % or use a first order Butterworth filter 
        HPChz = 1/HPC;
        LPChz = 1/LPC;
        TRhz  = 1/TR;
        ny    = (1/2)*TRhz; % Nyqust Frequency
        Wn    = [HPChz LPChz]/ny;
        n = 2;
        % Design a bandpass filter with a passband from HPC to LPC Hz with 
        % at most 3 dB of passband ripple and at least 40 dB attenuation in
        % the stopbands. Specify a sample rate of 1/TR Hz. Set the stopband 
        % width to half freqyency Hz on both sides of the passband. 
        % Find the filter order and cutoff frequencies.
%         Wp = [HPChz LPChz]/ny;
%         Ws = [HPChz*.5 LPChz*1.5]/ny;
%         Rp = 3;
%         Rs = 40;
%         [n,Wn] = buttord(Wp,Ws,Rp,Rs);
%         
%         [a,b]= butter(n,Wn);
%         [z,p,k] = butter(n,Wn(2));
%         [A,B,C,D] = butter(n,Wn);
%         sos = zp2sos(z,p,k);
%         figure;freqz(sos,512,1/TR)
        
        if doFilterType == 1 % Band-pass filter
            [filtb,filta] = butter(n,Wn);
            [zfi,pfi,kfi] = butter(n,Wn);
        elseif doFilterType == 2 % High-pass
            [filtb,filta] = butter(n,Wn(1),'high');
            [zfi,pfi,kfi] = butter(n,Wn(1),'high');
        elseif doFilterType == 3 % Low-pass
            [filtb,filta] = butter(n,Wn(2));
            [zfi,pfi,kfi] = butter(n,Wn(2));
        end
        sosfi = zp2sos(zfi,pfi,kfi);
    end
end

%% Detect outliers in movement and/or data
if doRegressMotion
    M       = spm_en(M);
    dM      = [zeros(1,6); diff(M(:,1:6),1,1)];    % First-order derivatives
    [r1, ~] = y_FD_Jenkinson(SS.f_rp,SS.f_in);% give the scan-scan displacement and add to movement parameters
    r1      = r1(trim_initvols+1:end,:);
    M       = [M r1];
    rms     = M(:,7);

    if ~isempty(SpikeMovAbsThr)  % Absolute translation threshold (don't need one for rotation?)
        aspk = find(rms>SpikeMovAbsThr);
        fprintf('%d spikes in absolute translation differences (based on threshold of %4.2f)\n',length(aspk),SpikeMovAbsThr)
    else
        aspk = [];
    end

    if ~isempty(SpikeMovRelThr)  % Relative (SD) threshold for translation and rotation
        rspk = find(rms > (mean(rms) + SpikeMovRelThr*std(rms)));
        fprintf('%d spikes in relative translation or rotation differences (based on threshold of %4.2f SDs)\n',length(rspk),SpikeMovRelThr)
    else
        rspk = [];
    end

    if ~isempty(SpikeDatRelThr)  % Relative (SD) threshold across all ROIs (dangerous? Arbitrary?)
        dY   = [zeros(1,Nr); diff(Y,1,1)];
        rms  = sqrt(mean(dY.^2,2));
        dspk = find(rms > (mean(rms) + SpikeDatRelThr*std(rms)));
        fprintf('%d spikes in mean data across ROIs (based on threshold of %4.2f SDs)\n',length(dspk),SpikeDatRelThr)
    else
        dspk = [];
    end

    % -----------------------------------------------
    % Create delta-function regressors for each spike
    % -----------------------------------------------
    spk = unique([aspk; rspk; dspk]); lspk = spk;
    for q = 1:(SpikeLag-1)
        lspk = [lspk; spk + q];
    end
    spk = unique(lspk);

    if ~isempty(spk)
        RSP = zeros(Ns,length(spk));
        n = 0;
        for p = 1:length(spk)
            if spk(p) <= Ns
                n=n+1;
                RSP(spk(p),n) = 1;
            end
        end
        fprintf('%d unique spikes in total\n',length(spk))
        %RSP = spm_en(RSP,0);
        strSPK= string(strcat('spk',num2str([1:size(RSP,2)]','%02d')));
    else
        RSP = [];
        strSPK = [];
    end
    
    % ----------------------------------------
    % Create expansions of movement parameters
    % ----------------------------------------
    if doMotionOpt==2 %Satterthwaite, 2012
        aM  = [M(:,1:6) dM M(:,1:6).^2 dM.^2];

    elseif doMotionOpt==1
        bf = eye(5);  % artifacts can last up to 5 TRs = 10s, according to Power et al (2013)
        bf = [bf; diff(eye(5))];
        U=[];
        for c=1:6;
            U(c).u = M(:,c);
            U(c).name{1}='c';
        end;
        aM = spm_Volterra(U,bf',2);
        aM = spm_en(aM,0);
    elseif doMotionOpt == 3 % 6rp
        aM = M;
    end    
end
 
%% Add Global? (Note: may be passed by User in S.C anyway)
% (recommended by Rik and Power et al, 2013, though will entail negative
% correlations, which can be problem for some graph-theoretic measures)
if doRegressGlobal
    C = [C SS.global];
    nameReg = [nameReg; 'GlobalSignal']; 
end

% Add WM signal
if doRegressWM
    C = [C SS.wm];
    nameReg = [nameReg; 'WMSignal']; 
end

% Add CSF signal
if doRegressCSF
    C = [C SS.csf];
    nameReg = [nameReg; 'CSFSignal'];
end
 
% Remove dummy entry in nameReg, if regressors have been addded
if size(nameReg,1)>1
    nameReg(nameReg=='') = [];
end


%% Include the derivative and the squared terms of the confounds

C = spm_en(C,0); % Eucledian normalisation

% Include derivative terms of confounds
if doDerivatives
    dC = [zeros(1,size(C,2)); diff(C,1,1)]; % First-oder derivative
    C = [C dC];
    nameReg = [nameReg; strcat(nameReg,'deriv')];
end


% Include squared terms of confounds
if doSquaredTerms
    C = [C C.^2];
    nameReg = [nameReg; strcat(nameReg,'squared')];
end



%% Combine all confounds (assumes more data than confounds, ie Ns > Nc) and perform dimension reduction (cf PCA of correlation matrix)

% Add CompCorr regressors
if doRegressCompCor
    C = [C CompCor];
    nameReg = [nameReg; string(strcat('Compcor',num2str([1:size(CompCor,2)]',('%02d'))))];
end

% Possible of course that some dimensions tiny part of SVD of X (so excluded) by happen to correlate highly with y...
% ...could explore some L1 (eg LASSO) or L2 regularisation of over-parameterised model instead,
% but LASSO takes ages (still working on possible L2 approach with spm_reml)
if doRegressMotion
    if doMotionOpt==1
        [U,scores] = spm_svd(aM,0);
        scores     = diag(scores).^2; scores = full(cumsum(scores)/sum(scores));
        Np         = find(scores > svd_thr); Np = Np(1);
        nmov       = Np;
        X0r        = full(U(:,1:Np));
        fprintf('%d modes left (from %d movement terms) - %4.2f%% variance of correlation explained\n',Np,size(aM,2),100*scores(Np))
    else
        nmov=24;
        X0r=aM;
    end
    strRP = string(strcat('rp',num2str([1:size(X0r,2)]','%02d')));
else
    X0r = [];
    RSP = [];
    spk = [];
    nmov = [];
    strRP = [];
end

X0r =[C X0r RSP]; 
nameReg = [nameReg; ...
            strRP;...
            strSPK];

nrem=size(RSP,2);

if doFilter & doFilterMethod == 1
    X0r = [X0r K];
    nameReg = [nameReg; string(strcat('dct',num2str([1:size(K,2)]',('%03d'))))];
end

% Add Detrend terms (linear and quadratic)
dlin  = [1:Ns]';
if doDetrend
    X0r = [ zscore([dlin dlin.^2]) X0r];
    nameReg = ['Linear'; 'Quadratic';nameReg];
end

% Z-scale covariates/regressors
if doZsoreCovariates 
    X0r  = zscore(X0r);
end

% Add constant term. If X0r is empty add constant term
if isempty(X0r)
    X0r = ones(Ns,1);
    nameReg = 'Constant';
else
    X0r = [ones(Ns,1) X0r];
    nameReg = ['Constant'; nameReg];
end


Nc  = size(X0r,2)+nrem;
if Nc >= Ns; error('Not enough dfs (scans) to estimate');
else fprintf('%d confounds for %d scans (%d left)\n',Nc,Ns,Ns-Nc); end



%% Create adjusted data and compute correlations, in case user wants for other metrics, eg, MI, and standardised data, if requested
% Calculate residual

if doRegressCovariates 
    R  = eye(Ns) - X0r*pinv(X0r);
    aY = R*Y;
elseif doFilterMethod == 1 % This is in case no covariate regression and Hallquist filtering was requested
    R  = eye(Ns) - K*pinv(K);
    aY = R*Y;
elseif doSeedFC
else % No post-processing
    aY = Y;
end

% Do Butterworth filtering, sequential
if doFilterMethod == 2
    aY = filtfilt(filtb,filta,aY);
%     aY = sosfilt(sosfi,aY,1);
end

nospike=setdiff(1:Ns, spk);
aY=aY(nospike,:);

if doZscoreYr
    aY = zscore(aY);
end

SS.Y    = Y;
SS.Yr   = aY;
SS.X0r  = X0r;
SS.Nr   = Nr; % Number of ROIs
SS.Ns   = Ns; % Number of timepoints (after triming dummy events)
SS.Nc   = Nc;
SS.nmov = nmov;
SS.nrem = nrem;
SS.nameRegressors = nameReg;


%% Compute and save SD/RSFA maps
if doSaveSD
   kat_fmri_postproc_glm_rsfa(SS, aY,V);% 
end

%% Compute and save relative CVR maps as in Liu et al 2017 Neuroimage
if doSaveCVR
    kat_fmri_postproc_glm_relCVR(SS,aY,V); % INCOMPLETE
end

%% Do GLMs for all target ROIs for a given source ROI
if doSaveFCnode2node
    SS = kat_fmri_postproc_glm_fc_node2node(SS);
end

%% Do GLMs across all voxels for a given source ROIs
if doSaveFCnode2voxel
    SS = kat_fmri_postproc_glm_fc_node2voxel(SS); 
end

%% Write out result for 4d images
if doSave4D 
    mkdir(fileparts(SS.f_out));
%     Info.V = V;
    if doAddConstant % May be needed by some spm functions for saving nii images
        kat_WriteOutNii(aY,SS.f_out,Info,doAddConstant);% Needs updating with spm_write_vol
    else
%         WriteOutNii(aY, SS.f_out, Info); % Old version using BWT
        dim  = V(1).dim;
        nvol = size(V,1)-doTrimVolumes;
        
        %-Check that number of volumes is correct (after removing dummy volumes)
        nvolCh2 = size(aY,1);
        if nvol ~= nvolCh2
            error('Number of volumes is conflicting.');
        end
        
        %-Move from 2d to 4d array
        aY4d = reshape(aY,[nvol dim(1) dim(2) dim(3)]);
        aY4d = permute(aY4d, [2 3 4 1]);
 
        %-mask image (if mask is provided)
        if isfield(SS,'f_mask')
            Ym = logical(spm_read_vols(spm_vol(SS.f_mask)));
            aY4d = aY4d.*repmat(Ym,1, 1, 1, nvol);
        end
        
        %-write the 4D nifty file
        for ivol=1:nvol
            V(ivol).fname   = SS.f_out;%regexprep(SS.f_out,'.nii','_1.nii');
            spm_write_vol(V(ivol),aY4d(:,:,:,ivol));
        end
    end
    return;
end
