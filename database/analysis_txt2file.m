function analysis_txt2file(as,fname)
% Writes textual responses from the analysis structure out to a file.
%
% analysis_txt2file(as,fname)
%

% 08/20/05 PJ
% 10/25/05 PJ  -- Fixed handling of empty stimidx field

td = as.txt;

nsubs = length(td.subid);
nquest = length(td.qid);

max_stims = size(td.data,3);

fid = fopen(fname,'wt');
fprintf(fid,'Question\tSubject\tStimID\tResponse\n');
for iq = 1:nquest
  fprintf(fid,'\n%s\n',td.qtxt{iq});
  for isub = 1:nsubs
    fprintf(fid,'\n\t%s\n', td.subid{isub});
    for istim = 1:max_stims
      if td.niter(isub,iq) && ~isempty(td.data{isub,iq,istim})
	if ~isempty(td.stimidx)
	  stim_id = as.stims.ids(td.stimidx(isub,istim));
	  stim_str = sprintf('%d', stim_id);
	else
	  stim_str = '';
	end
	fprintf(fid,'\t\t%s\t%s\n', stim_str, td.data{isub,iq,istim});
      end
    end % for istim
  end % for isub=
end % for iquest
fclose(fid);

return
