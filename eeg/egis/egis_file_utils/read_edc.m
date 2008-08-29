function edit_codes = read_edc(fname, fhdr, chdr, nchan)
% read_edc.m  -- reads an edit codes file
%
% fhdr: EGIS format specific vector
% chdr: EGIS format specific vector
%

% Modification history:
% 7/7/95		Written by Petr Janata
%
% 7/15/01 PJ - Modified to read .edc format file without relying on any EGIS
% format information

mode = 'GENERIC';

if nargin > 1 & nargin < 4
  mode = 'EGIS';
end

if nargin == 4
  mode = 'GENERIC';
end

if nargin < 1
  nchan = input('How many channels: ');
end

% Initialize some things
switch mode
  case 'EGIS'
    ses_hdr_offsets_v;
    edit_codes = ones(sum(chdr(:,NTrials)), fhdr(NChan));
  case 'GENERIC'
    edit_codes = ones(1,nchan);
end

fid = fopen(fname);
if fid == -1
  error(sprintf('Error opening edit codes file: %s', fname))
end

disp(['Reading edit codes file <' fname '>']);
done = 0; 
lines_read = 0;

while ~done
  theline = fgetl(fid);
  if theline == -1
    done = 1;
	disp(['Finished reading edit codes file (' int2str(lines_read) ' lines)']);
  else
    lines_read = lines_read + 1;
    theline = sscanf(theline, '%d');
    nelem = length(theline);
    thecell = theline(1); 
    thetrial = theline(2);
    max_trials(thecell) = thetrial;
	
    switch mode
      case 'EGIS'
	if (thecell <= 0) | (thecell > fhdr(NCells))
	  warning(['Invalid cell (' int2str(thecell) ') specified']);
	end
	if (thetrial > chdr(thecell, NTrials))
	  warning(['Invalid trial number (' int2str(thetrial) ') for cell ' ...
		int2str(thecell) ]);
	end
	
	trial_offset = sum(chdr(1:current_cell, NTrials)) - chdr(current_cell, ...
	    NTrials) + thetrial;
	
	if (theline(3) >= 2)
	  edit_codes(trial_offset, :) = zeros(1, size(edit_codes,2));
	end
	if (nelem > 3) & (theline(3) == 1)
	  for i = 4:nelem
	    if (theline(i) < 0) | (theline(i) > fhdr(NChan))
	      disp(['Abnormality in line ' int2str(lines_read + 1) ' of edit codes file']);
	    elseif theline(i) ~= 0
	      edit_codes(trial_offset, theline(i)) = 0;
	    end
	  end
	end % case EGIS
      case 'GENERIC'
	trial_offset = sum(max_trials(1:thecell));
	if theline(3) >=2
	  edit_codes(trial_offset, 1:nchan) = 0;
	else
	  edit_codes(trial_offset, 1:nchan) = 1;
	end
	if (nelem > 3) & (theline(3) == 1)
	  edit_codes(trial_offset, theline(4:end)) = 0;
	end
    end % switch mode
  end % if theline == -1
end % while ~done

fclose(fid);
