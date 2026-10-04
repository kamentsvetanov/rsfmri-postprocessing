function [SS] = kat_fmri_postproc_glm_fc_node2node(SS)
% Auxiliary function to fmri_processing_GLM to generate functional
% connectivity based on Rik's GLM function

% try pflag               = SS.do.PartialFC;             catch pflag                 = 0; end
% try svd_thr = SS.svd_thr;  catch svd_thr = .99; end
    
% --------------------
% Read ROIs
% --------------------
Vroi = spm_vol(SS.f_rois);
Yroi = spm_read_vols(Vroi);
dim  = Vroi.dim;


% Check if ROI image is 3D or 4D, e.g. having more than on ROIs

if size(Yroi,4)==1
    % Move from 3D image with unique int intesity value associated with
    % each roi to 4D image, where each volume represents a mask of the RIO,
    % e.g. ROI with intensity value 2 becomes volume 2 in the new 4D image
    Yroi    = int32(Yroi);
    nvols   = numel(unique(Yroi))-1; % exclude 0 valued voxels
    Yroi4d  = zeros([size(Yroi) nvols]);
    for ivol = 1:nvols
        Yroi4d(:,:,:,ivol) = Yroi==ivol;
    end
    Yroi = Yroi4d;
end
Yroi = permute(Yroi,[4 1 2 3]);
Yroi = logical(Yroi(:,:));


% -------------------
% Unpack SS variables
% -------------------
Ns          = SS.Ns;
Nr          = size(Yroi,1);
Nc          = SS.Nc;
Nv          = size(Yroi,2);               
TR          = SS.TR;
Y           = SS.Y;

if SS.do.zscoreY
    Y = zscore(Y);
end

X0r         = SS.X0r;
doPreWhiten = SS.do.PreWhiten;
doWantZval  = SS.do.wantZval;
doPartialFC = SS.do.PartialFC;
svd_thr     = SS.svd_thr;

Ar    = [1:Nr];
Zmat  = zeros(Nr,Nv);  % Matrix of Z statistics for each pairwise regression
Bmat  = zeros(Nr,Nv);
Pmat  = zeros(Nr,Nv);


pZmat = [];
pBmat = [];


%% Create comprehensive model of residual autocorrelation (to counter Eklund et al, 2012, Neuroimage)
T     = (0:(Ns - 1))*TR;                    % time
d     = 2.^(floor(log2(TR/4)):log2(64));    % time constants (seconds)
Q    = {};                                  % dictionary of components
for i = 1:length(d)
    for j = 0:1
        Q{end + 1} = toeplitz((T.^j).*exp(-T/d(i)));
    end
end



%[V h]   = rik_reml(CY,X0r,Q,1,0,4);   % if could assume that error did not depend on seed timeseries
%W       = spm_inv(spm_sqrtm(V));

lNpc = 0;  % Just for printing pflag output below
for n = 1:Nr
%     tic
%     fprintf('%d \n',n);
    Yr = Y;        % Yr contains all other timeseries, so ANCOVA below run on Nr-2 timeseries in one go
    Xr = mean(Yr(:,Yroi(n,:)),2);
    Xtmp(:,n) = Xr;

%     Xr = spm_en(Xr);     % Normalise regressors so Betas can be compared directly across (seed) ROIs (not necessary if Standardised Y already)?
    X  = [Xr X0r ones(Ns,1)]; % add a constant
%     X  = spm_en(X);     % Normalise regressors so Betas can be compared directly across (seed) ROIs (not necessary if Standardised Y already)?

    %        if lassoflag == 1  % takes too long (particularly for cross-validation to determine lambda)
    %            for pn = 1:length(Ir)
    %                Ir = setdiff(Ar,pn);
    %
    %                Yr = Y(:,Ir);        % Yr contains all other timeseries, so ANCOVA below run on Nr-2 timeseries in one go
    %                Xr = Y(:,pn);
    %                Xr = spm_en(Xr,0);     % Normalise regressors so Betas can be compared directly across (seed) ROIs (not necessary is Standardised Y already)?
    %
    %                fprintf('l');
    %                [lB,lfit] = lasso(X0r,Yr(:,pn),'CV',10,'Options',opts);
    %                keepX     = find(lB(:,lfit.Index1SE));
    %                X         = [Xr X0r(:,keepX)];
    %                Nc        = length(keepX);
    %
    %                [V h]   = rik_reml(CY,X,Q,1,0,4); % rik_reml is just version of spm_reml with fprintf commented out to speed up
    %                W       = spm_inv(spm_sqrtm(V));
    %
    %                %% Estimate T-value for regression of first column (seed timeseries)
    %                [T(pn),dfall,Ball] = spm_ancova(W*X,speye(Ns,Ns),W*Yr(:,pn),[1 zeros(1,Nc)]');
    %                df(pn) = dfall(2);
    %                B(pn)  = Ball(1);
    %            end
    %            fprintf('\n');
    %        else

    %% Estimate autocorrelation of error (pooling across ROIs) and prewhitening matrix
    W = [];
    if doPreWhiten
%         Y = detrend(Y,'constant');
        Y = zscore(Y);
        CY = cov(Y');  % Correct to only do this after any standardisation?
        [V h]   = rik_reml(CY,X,Q,1,0,4);   % rik_reml is just version of spm_reml with fprintf commented out to speed up
        W       = spm_inv(spm_sqrtm(V));
    else
        W       = speye(Ns);
    end

    %% Estimate T-value for regression of first column (seed timeseries)

    if doWantZval
        try 
            [T,df,B]   = spm_ancova(W*X,speye(Ns,Ns),W*Yr,[1 zeros(1,Nc+1)]'); % Adding 1 to NC for constant term

            % check for nan ,e.g. ROI timeseries are empty; edit kat
            if sum(isnan(T))>0
                disp(SS.SubID)
                idx = ~isnan(T); 
                T = T(idx); 
                Pmat(n,:) =nan;
                Pmat(n,idx) = 1-spm_Tcdf(abs(T),df(2));
                P          = spm_Tcdf(abs(T),df(2));
                P(P>0.9999999999999999)=0.9999999999999999;
                Zmat(n,:) =nan; % Edit kat
                Zmat(n,idx) = norminv(P);  % Needs Matlab Stats toolbox, but has larger range of Z (so not so many "Inf"s)
                Zmat(n,idx)=Zmat(n,idx).*(T./(abs(T)));
            else
                Pmat(n,idx) = 1-spm_Tcdf(abs(T),df(2));
                P           = spm_Tcdf(abs(T),df(2));
                P(P>0.9999999999999999)=0.9999999999999999;
                Zmat(n,idx) = norminv(P);  % Needs Matlab Stats toolbox, but has larger range of Z (so not so many "Inf"s)
                Zmat(n,idx)=Zmat(n,idx).*(T./(abs(T)));
            end
%             if sum(isinf(Zmat(n,Ir)))>0
%                 disp(SS.SubID)
%             end
        catch 
            Zmat(n,Ir) =nan; 
            Pmat(n,Ir) =nan; 
        end

    else
        try
            B = pinv(X)*Yr;
        catch
            %keyboard
        end
    end
    Bmat(n,:) = B(1,:);
%     toc
    fprintf('.')
end

fprintf('\n Done. \n');
SS.Zmat     = Zmat;
SS.Bmat     = Bmat;
SS.Pmat     = Pmat;

%% Write results to images
mkdir(fileparts(SS.f_out));

%-Move from 2d to 3d array
aY4d = reshape(Bmat,[Nr dim(1) dim(2) dim(3)]);
aY4d = permute(aY4d, [2 3 4 1]);

%-mask image (if mask is provided)
if isfield(SS,'f_mask')
    Ym = logical(spm_read_vols(spm_vol(SS.f_mask)));
    aY4d = aY4d.*repmat(Ym,1, 1, 1, Nr);
end

%-write the 4D nifty file
for iroi=1:Nr
    f_out           = SS.f_out;%regexprep(SS.f_out,'.nii','_1.nii');
    f_out           = regexprep(f_out,'.nii',[SS.namerois{iroi} '.nii']); % Append ROI name
    Vroi(1).fname   = f_out;
    spm_write_vol(Vroi(1),aY4d(:,:,:,iroi));

    % Apply smoothing
    try smoothFilt = SS.smooth; catch smoothFilt = []; end
    for iSmooth = 1:numel(smoothFilt)
        sfilt = smoothFilt(iSmooth);
        f_new = regexprep(Vroi(1).fname,'.nii',sprintf('_s%d.nii',sfilt));
        spm_smooth(Vroi(1).fname,f_new,[sfilt sfilt sfilt]);
        % Mask again
        Vtemp = spm_vol(f_new);
        Y = spm_read_vols(Vtemp);
        Ymasked = Y.*Ym;
        spm_write_vol(Vtemp,Ymasked);
    end
    delete(Vroi(1).fname);
end

% fprintf('\n')

return


