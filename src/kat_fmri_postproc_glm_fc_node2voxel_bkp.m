function [SS] = kat_fmri_postproc_glm_fc_node2voxel(SS)
% Auxiliary function to fmri_processing_GLM to generate seed-based
% activation maps based on Rik's GLM function

% -------------------
% Unpack SS variables
% -------------------
Ns          = SS.Ns;
Nr          = SS.Nr;
Nc          = SS.Nc;
TR          = SS.TR;
Y           = SS.Y;
X0r         = SS.X0r;
doWantZval  = SS.do.wantZval;


Ar    = [1:Nr];
Zmat  = zeros(Nr);  % Matrix of Z statistics for each pairwise regression
Bmat  = zeros(Nr);
Pmat  = zeros(Nr);

if SS.do.zscoreY
    Y = zscore(Y);
end


%% Create comprehensive model of residual autocorrelation (to counter Eklund et al, 2012, Neuroimage)
T     = (0:(Ns - 1))*TR;                    % time
d     = 2.^(floor(log2(TR/4)):log2(64));    % time constants (seconds)
Q    = {};                                  % dictionary of components
for i = 1:length(d)
    for j = 0:1
        Q{end + 1} = toeplitz((T.^j).*exp(-T/d(i)));
    end
end


% Loop through all ROIs/Seed Regions
lNpc = 0;  % Just for printing pflag output below

for n = 1:Nr
    

    Xr = mean(Y(:,Yroi(n,:)),2);% Seed region 
    X  = [Xr X0r ones(Ns,1)]; % add a constant

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


        if doWantZval
        try 
            [T,df,B]   = spm_ancova(W*X,speye(Ns,Ns),W*Yr,[1 zeros(1,Nc+1)]'); % Adding 1 to NC for constant term

            % check for nan ,e.g. ROI timeseries are empty; edit kat
            if sum(isnan(T))>0
                disp(SS.SubID)
                idx = ~isnan(T); 
                T = T(idx); 
                Pmat(n,Ir) =nan;
                Pmat(n,Ir(idx)) = 1-spm_Tcdf(abs(T),df(2));
                P          = spm_Tcdf(abs(T),df(2));
                P(P>0.9999999999999999)=0.9999999999999999;
                Zmat(n,Ir) =nan; % Edit kat
                Zmat(n,Ir(idx)) = norminv(P);  % Needs Matlab Stats toolbox, but has larger range of Z (so not so many "Inf"s)
                Zmat(n,Ir(idx))=Zmat(n,Ir(idx)).*(T./(abs(T)));
            else
                Pmat(n,Ir) = 1-spm_Tcdf(abs(T),df(2));
                P          = spm_Tcdf(abs(T),df(2));
                P(P>0.9999999999999999)=0.9999999999999999;
                Zmat(n,Ir) = norminv(P);  % Needs Matlab Stats toolbox, but has larger range of Z (so not so many "Inf"s)
                Zmat(n,Ir)=Zmat(n,Ir).*(T./(abs(T)));
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
    Bmat(n,Ir) = B(1,:);
end


fprintf('\n Done. \n');
SS.Zmat     = Zmat;
SS.Bmat     = Bmat;
SS.Pmat     = Pmat;
% SS.pflag    = doPartialFC;
SS.partial_Zmat = pZmat;
SS.partial_Bmat = pBmat;
% fprintf('\n')
