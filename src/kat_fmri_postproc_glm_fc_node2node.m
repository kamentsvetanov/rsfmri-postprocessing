function [SS] = kat_fmri_postprocessing_GLM_fc(SS)
% Auxiliary function to fmri_processing_GLM to generate functional
% connectivity based on Rik's GLM function

% try pflag               = SS.do.PartialFC;             catch pflag                 = 0; end
% try svd_thr = SS.svd_thr;  catch svd_thr = .99; end
    

% -------------------
% Unpack SS variables
% -------------------
Ns          = SS.Ns;
Nr          = SS.Nr;
Nc          = SS.Nc;
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
Zmat  = zeros(Nr);  % Matrix of Z statistics for each pairwise regression
Bmat  = zeros(Nr);
Pmat  = zeros(Nr);

if doPartialFC
    pZmat = zeros(Nr);  % Matrix of Z statistics for each pairwise partial regression
    pBmat = zeros(Nr);  % Matrix of Z statistics for each pairwise partial regression
else
    pZmat = [];
    pBmat = [];
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



%[V h]   = rik_reml(CY,X0r,Q,1,0,4);   % if could assume that error did not depend on seed timeseries
%W       = spm_inv(spm_sqrtm(V));

lNpc = 0;  % Just for printing pflag output below
for n = 1:Nr
%     tic
%     fprintf('%d \n',n);
    Ir = setdiff(Ar,n);

    Yr = Y(:,Ir);        % Yr contains all other timeseries, so ANCOVA below run on Nr-2 timeseries in one go
    Xr = Y(:,n);
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


    %% Estimate T-value for PARTIAL regression of first column (seed timeseries) - takes ages!
    if doPartialFC
       
        % SVD only once with all timeseries, as a signle node has
        % negligible impact on it
        % Ok to assume error and hence SVD unaffected by addition of Xr in Yr?
%         YrEN    = spm_en(Yr,0);       % Is necessary for SVD below
%         [U,ss]  = spm_svd(zscore([Yr X0r]),0);
%         ss      = diag(ss).^2; ss = full(cumsum(ss)/sum(ss));
%         Np      = find(ss > svd_thr); Np = Np(1);
        
        
        for pn = 1:length(Ir)

            pIr = setdiff(Ar,[n Ir(pn)]);
            pY  = spm_en(Y(:,pIr),0);       % Is necessary for SVD below
%             pY  = Y(:,pIr);
            % SVD
            if svd_thr < 1
                [U,ss] = spm_svd([pY X0r],0);
                ss     = diag(ss).^2; ss = full(cumsum(ss)/sum(ss));
                Np    = find(ss > svd_thr); Np = Np(1);
                XY0   = full(U(:,1:Np));
%                 Npc   = Np;
            else
                pIr = setdiff(Ar,[n Ir(pn)]);
                pY  = Y(:,pIr);
                XY0   = [pY X0r];
            end
           

            XY0 = [XY0 ones(Ns,1)];  % Reinsert mean
            Npc = size(XY0,2);

            %                if lassoflag == 1
            %                    fprintf('l')
            %                    [lB,lfit] = lasso(XY0,Yr(:,pn),'CV',10,'Options',opts);
            %                    keepX     = find(lB(:,lfit.Index1SE));
            %                    X         = [Xr XY0(:,keepX)];
            %                    Nc        = length(keepX);
            %
            %                    [V h]   = rik_reml(CY,X,Q,1,0,4); % rik_reml is just version of spm_reml with fprintf commented out to speed up
            %                    W       = spm_inv(spm_sqrtm(V));
            %
            %                    %% Estimate T-value for regression of first column (seed timeseries)
            %                    [T,df,B]   = spm_ancova(W*X,speye(Ns,Ns),W*Yr(:,pn),[1 zeros(1,Nc)]');
            %                else

            if Npc >= (Ns-1)  % -1 because going to add Xr below
                warning('Not enough dfs (scans) to estimate - just adjusting data and ignoring loss of dfs');
                R = eye(Ns) - pY*pinv(pY);    % Residual-forming matrix
                [T,df,B] = spm_ancova(W*X,speye(Ns,Ns),R*W*Yr(:,pn),[1 zeros(1,Npc)]');  % Ok to assume error and hence W unaffected by addition of pY in X?
            else

                if lNpc ~= Npc  % Just to reduce time taken to print to screen
                    fprintf('   partial for seed region %d and target region %d: %d confounds for %d scans (%d left)\n',n,pn,Npc,Ns,Ns-Npc);
                end
                lNpc = Npc;

                X  = [Xr XY0];

                %% Commented out below because ok to assume error and hence W unaffected by addition of pY in X?
                %                   [V h]   = rik_reml(CY,X,Q,1,0,4);   % rik_reml is just version of spm_reml with fprintf commented out to speed up
                %                   W       = spm_inv(spm_sqrtm(V));

                %% Estimate T-value for regression of first column (seed timeseries)
                [T,df,B]   = spm_ancova(W*X,speye(Ns,Ns),W*Yr(:,pn),[1 zeros(1,Npc)]');
            end

            try
                pZmat(n,Ir(pn)) = norminv(spm_Tcdf(T,df(2)));
            catch
                pZmat(n,Ir(pn)) = spm_invNcdf(spm_Tcdf(T,df(2)));
            end

            pBmat(n,Ir(pn)) = B(1);
        end
    end
%     toc
    fprintf('.')
end

fprintf('\n Done. \n');
SS.Zmat     = Zmat;
SS.Bmat     = Bmat;
SS.Pmat     = Pmat;
% SS.pflag    = doPartialFC;
SS.partial_Zmat = pZmat;
SS.partial_Bmat = pBmat;
% fprintf('\n')

return





% 
% function [C] = spm_trace(A,B)
% % fast trace for large matrices: C = spm_trace(A,B) = trace(A*B)
% % FORMAT [C] = spm_trace(A,B)
% %
% % C = spm_trace(A,B) = trace(A*B) = sum(sum(A'.*B));
% %__________________________________________________________________________
% % Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% 
% % Karl Friston
% % $Id: spm_trace.m 4805 2012-07-26 13:16:18Z karl $
% 
% % fast trace for large matrices: C = spm_trace(A,B) = trace(A*B)
% %--------------------------------------------------------------------------
% C = sum(sum(A'.*B));


%% Similar to above ROI correlation estimation

% 
% %add the spikes back in for the GLM
% X0r=[X0r RSP];
% 
% 
% 
% %% Code for estimating corration betwen ROIs; not relevant to my needs
% %% Pool data covariance over ROIs (assuming enough of them!), unless specified
% if isempty(CY)
%     CY = cov(Y');  % Correct to only do this after any standardisation?
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% Main Loop
% %% Do GLMs for all target ROIs for a given source ROI
% Ar    = [1:Nr];
% Zmat  = zeros(Nr);  % Matrix of Z statistics for each pairwise regression
% Bmat  = zeros(Nr);
% 
% if pflag
%     pZmat = zeros(Nr);  % Matrix of Z statistics for each pairwise partial regression
%     pBmat = zeros(Nr);  % Matrix of Z statistics for each pairwise partial regression
% else
%     pZmat = [];
%     pBmat = [];
% end
% 
% %[V h]   = rik_reml(CY,X0r,Q,1,0,4);   % if could assume that error did not depend on seed timeseries
% %W       = spm_inv(spm_sqrtm(V));
% 
% lNpc = 0;  % Just for printing pflag output below
% for n = 1:Nr
% %     Ir = setdiff(Ar,n);
%     
% %     Yr = Y(:,Ir);        % Yr contains all other timeseries, so ANCOVA below run on Nr-2 timeseries in one go
% %     Xr = Y(:,n);
%     
%     %       Xr = spm_en(Xr);     % Normalise regressors so Betas can be compared directly across (seed) ROIs (not necessary if Standardised Y already)?
%     X  = [Xr X0r];
%     
%     %        if lassoflag == 1  % takes too long (particularly for cross-validation to determine lambda)
%     %            for pn = 1:length(Ir)
%     %                Ir = setdiff(Ar,pn);
%     %
%     %                Yr = Y(:,Ir);        % Yr contains all other timeseries, so ANCOVA below run on Nr-2 timeseries in one go
%     %                Xr = Y(:,pn);
%     %                Xr = spm_en(Xr,0);     % Normalise regressors so Betas can be compared directly across (seed) ROIs (not necessary is Standardised Y already)?
%     %
%     %                fprintf('l');
%     %                [lB,lfit] = lasso(X0r,Yr(:,pn),'CV',10,'Options',opts);
%     %                keepX     = find(lB(:,lfit.Index1SE));
%     %                X         = [Xr X0r(:,keepX)];
%     %                Nc        = length(keepX);
%     %
%     %                [V h]   = rik_reml(CY,X,Q,1,0,4); % rik_reml is just version of spm_reml with fprintf commented out to speed up
%     %                W       = spm_inv(spm_sqrtm(V));
%     %
%     %                %% Estimate T-value for regression of first column (seed timeseries)
%     %                [T(pn),dfall,Ball] = spm_ancova(W*X,speye(Ns,Ns),W*Yr(:,pn),[1 zeros(1,Nc)]');
%     %                df(pn) = dfall(2);
%     %                B(pn)  = Ball(1);
%     %            end
%     %            fprintf('\n');
%     %        else
%     
%     %% Estimate autocorrelation of error (pooling across ROIs) and prewhitening matrix
%     if PreWhiten
%         [V h]   = rik_reml(CY,X,Q,1,0,4);   % rik_reml is just version of spm_reml with fprintf commented out to speed up
%         W       = spm_inv(spm_sqrtm(V));
%     else
%         W       = speye(Ns);
%     end
%     
%     %% Estimate T-value for regression of first column (seed timeseries)
%     
%     if WantZval
%         [T,df,B]   = spm_ancova(W*X,speye(Ns,Ns),W*Yr,[1 zeros(1,Nc)]');
%         try
%             Zmat(n,Ir) = norminv(spm_Tcdf(T,df(2)));  % Needs Matlab Stats toolbox, but has larger range of Z (so not so many "Inf"s)
%         catch
%             Zmat(n,Ir) = spm_invNcdf(spm_Tcdf(T,df(2)));
%         end
%     else
%         try
%             B = pinv(X)*Yr;
%         catch
%             %keyboard
%         end
%     end
%     Bmat(n,Ir) = B(1,:);
% 
%     
%     %% Estimate T-value for PARTIAL regression of first column (seed timeseries) - takes ages!
%     if pflag
%         for pn = 1:length(Ir)
%             
%             pIr = setdiff(Ar,[n Ir(pn)]);
%             pY  = spm_en(Y(:,pIr),0);       % Is necessary for SVD below
%             
%             % SVD
%             if svd_thr < 1
%                 [U,S] = spm_svd([pY X0],0);
%                 S     = diag(S).^2; S = full(cumsum(S)/sum(S));
%                 Np    = find(S > svd_thr); Np = Np(1);
%                 XY0   = full(U(:,1:Np));
%                 Npc   = Np;
%             else
%                 XY0   = [pY X0];
%             end
%             
%             XY0 = [XY0 ones(Ns,1)];  % Reinsert mean
%             Npc = size(XY0,2);
%             
%             %                if lassoflag == 1
%             %                    fprintf('l')
%             %                    [lB,lfit] = lasso(XY0,Yr(:,pn),'CV',10,'Options',opts);
%             %                    keepX     = find(lB(:,lfit.Index1SE));
%             %                    X         = [Xr XY0(:,keepX)];
%             %                    Nc        = length(keepX);
%             %
%             %                    [V h]   = rik_reml(CY,X,Q,1,0,4); % rik_reml is just version of spm_reml with fprintf commented out to speed up
%             %                    W       = spm_inv(spm_sqrtm(V));
%             %
%             %                    %% Estimate T-value for regression of first column (seed timeseries)
%             %                    [T,df,B]   = spm_ancova(W*X,speye(Ns,Ns),W*Yr(:,pn),[1 zeros(1,Nc)]');
%             %                else
%             
%             if Npc >= (Ns-1)  % -1 because going to add Xr below
%                 warning('Not enough dfs (scans) to estimate - just adjusting data and ignoring loss of dfs');
%                 R = eye(Ns) - pY*pinv(pY);    % Residual-forming matrix
%                 [T,df,B] = spm_ancova(W*X,speye(Ns,Ns),R*W*Yr(:,pn),[1 zeros(1,Npc)]');  % Ok to assume error and hence W unaffected by addition of pY in X?
%             else
%                 
%                 if lNpc ~= Npc  % Just to reduce time taken to print to screen
%                     fprintf('   partial for seed region %d and target region %d: %d confounds for %d scans (%d left)\n',n,pn,Npc,Ns,Ns-Npc);
%                 end
%                 lNpc = Npc;
%                 
%                 X  = [Xr XY0];
%                 
%                 %% Commented out below because ok to assume error and hence W unaffected by addition of pY in X?
%                 %                   [V h]   = rik_reml(CY,X,Q,1,0,4);   % rik_reml is just version of spm_reml with fprintf commented out to speed up
%                 %                   W       = spm_inv(spm_sqrtm(V));
%                 
%                 %% Estimate T-value for regression of first column (seed timeseries)
%                 [T,df,B]   = spm_ancova(W*X,speye(Ns,Ns),W*Yr(:,pn),[1 zeros(1,Npc)]');
%             end
%             
%             try
%                 pZmat(n,Ir(pn)) = norminv(spm_Tcdf(T,df(2)));
%             catch
%                 pZmat(n,Ir(pn)) = spm_invNcdf(spm_Tcdf(T,df(2)));
%             end
%             
%             pBmat(n,Ir(pn)) = B(1);
%         end
%     end
%     fprintf('.')
% end
% fprintf('\n')
% 
% return
% 
% %figure,imagesc(Zmat); colorbar
% %figure,imagesc(pZmat); colorbar
% 
% 