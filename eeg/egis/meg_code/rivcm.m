
cm1 = hot(32);
cm1 = cm1(1:28,:);			% ramesh uses 1:56 of 64

cm2 = jet(64);
deltab = cm2(2,3) - cm2(1,3);
newblue = 0:deltab:(cm2(1,3)-deltab);
cm2 = [ [zeros(length(newblue'),2) newblue']; cm2];
cm2 = cm2(32:-1:1,:);

cm = [cm2; cm1];