% convert_biosemi_locs.m
%
% Creates .ced files from Biosemi files

fpath = '/afs/cmb.ucdavis.edu/share/matlab/janata/eeg/electrode_files';

files = {'biosemi_128ch.txt', ...
      'biosemi_256ch.txt'};

nfiles = length(files);

for ifile = 1:nfiles
  fname = fullfile(fpath, files{ifile});
  fprintf('Converting file: %s\n', fname);
  
  eloc = readlocs(fname,'filetype','custom',...
	'format',{'labels','custom1','custom2','custom3'});

  nchan = length(eloc);
  for ichan = 1:nchan
    % Convert phi coordinates
    eloc(ichan).sph_phi = -abs(eloc(ichan).custom1)+90;
    
    % Convert the theta coordinates
    eloc(ichan).sph_theta = eloc(ichan).custom2-90*sign(eloc(ichan).custom1);
  end
  
  % Remove the custom fields
  rmfield(eloc,{'custom1','custom2','custom3'});
  
  % Save as .ced file
  [fpath,fstub] = fileparts(fname);
  outfname = fullfile(fpath, sprintf('%s.ced',fstub));
  fprintf('Writing new electrode location file to: %s\n', outfname);
  writelocs(eloc,outfname,'filetype','chanedit');
end % for ifile=
