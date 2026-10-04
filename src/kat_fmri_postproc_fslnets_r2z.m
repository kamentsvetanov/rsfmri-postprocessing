function netmats = kat_fmri_postprocessing_fslnets_r2z(ts);
% convert from r to z using FSLNets method
% Fisher transform to Z statistics (of correlation coeffieicnt) to account
% for autocorrelation in the signal using the FSLnets package (Smith et al., 2011).

N = ts.Nnodes;
ts.Nsubjects = 1;
netmats = ts.netmats;
MethodType = 1;


 % quick crappy estimate of median AR(1) coefficient
    arone=[];
    for s=1:ts.Nsubjects
      grot=ts.ts((s-1)*ts.NtimepointsPerSubject+1:s*ts.NtimepointsPerSubject,:);
      for i=1:N
        g=grot(:,i);  arone=[arone sum(g(1:end-1).*g(2:end))/sum(g.*g)];
      end
    end
    arone=median(arone);

    % create null data using the estimated AR(1) coefficient
    clear grot*; grotR=[];
    for s=1:ts.Nsubjects
      for i=1:N
        grot(1)=randn(1);
        for t=2:ts.NtimepointsPerSubject
          grot(t)=grot(t-1)*arone+randn(1);
        end
        grotts((s-1)*ts.NtimepointsPerSubject+1:s*ts.NtimepointsPerSubject,i)=grot;
      end
      if MethodType==1
        grotr=corr(grotts((s-1)*ts.NtimepointsPerSubject+1:s*ts.NtimepointsPerSubject,:));
      else
        grotr=-inv(cov(grotts((s-1)*ts.NtimepointsPerSubject+1:s*ts.NtimepointsPerSubject,:)));
        grotr=(grotr ./ repmat(sqrt(abs(diag(grotr))),1,N)) ./ repmat(sqrt(abs(diag(grotr)))',N,1);
      end
      grotR=[grotR; grotr(eye(N)<1)];
    end
    grotZ=0.5*log((1+grotR)./(1-grotR));
    RtoZcorrection=1/std(grotZ);

    netmats=real(0.5*log((1+netmats)./(1-netmats))*RtoZcorrection);
  end