function ref_trialdata = rereference(trialdata)

ref_trialdata = zeros(size(trialdata,1),size(trialdata,2) - 1);

ref_trialdata = trialdata(:,1:size(trialdata,2) - 1);

ref_signal = trialdata(:,size(trialdata,2))*ones(1,size(ref_trialdata,2));

ref_trialdata = ref_trialdata - ref_signal;

