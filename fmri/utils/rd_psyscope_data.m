function data = rd_psyscope_data(fname)
% data = rd_psyscope_data(fname);
%
% Reads event timing information from psyscope data file
%

% 03/24/00 PJ

fid = fopen(fname,'rt');

if fid == -1
  error(sprintf('Could not open: <%s>'))
end

data = struct('trial',[],'condition',[],'time',[],'state',[]);

%
% Look for start of data
%

found_start = 0;

while ~found_start
  line = fgetl(fid);
  if findstr(line,'Trial') & findstr(line,'Condition')
    disp('Found start of data ...')
    found_start = 1;
  end
end

done = 0;
nevents = 0;

disp('Scanning for events ...')

while ~done
  line = fgets(fid);
  if line == -1
    done = 1;
  else
    vals = sscanf(line,'%d',3);

    if vals
      nevents = nevents + 1;
      data.trial(nevents) = vals(1);
      data.time(nevents) = vals(2);
      data.state(nevents) = vals(3);
    end
  end
end

disp(sprintf('Found %d events',nevents))

fclose(fid);