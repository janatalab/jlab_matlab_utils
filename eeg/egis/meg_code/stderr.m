function errors = stderr(coh,ntrial);

rtcoh = sqrt(coh);
errors = sqrt(2)*(ones(size(coh,1),size(coh,2))-coh)./rtcoh;
errors = errors/sqrt(ntrial);
