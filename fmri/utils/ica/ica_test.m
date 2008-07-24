% ica_test.m
%
% My first attempt at running ICA on fMRI data
%

% 02/01/01

GET_FILE_INFO = 0;
LOAD_RAW_DATA = 0;
COMPUTE_ICA = 0;
SAVE_ICA_DATA = 0;
LOAD_ICA_DATA = 0;
COMPUTE_COMPONENTS = 0;
WRITE_COMPONENTS = 0;
WRITE_Z = 0;

rootdir = '/data1/prime_long/';
outdir = './testimages/';

sub_id = '25jan01JP';
run_id = 1;

if GET_FILE_INFO
  % Specify the files to load
  P = spm_get('Files',fullfile(rootdir,sub_id,'epi', sprintf('run%d',run_id)),'sn*.img');
  
  nvol = size(P,1);

  % Get the dimension info for the first file
  disp(sprintf('Getting info for %d files', nvol))
  V = spm_vol(P);
end

if LOAD_RAW_DATA
  % Load the data
  clear indata
  disp(sprintf('Loading data'))
  indata = zeros(nvol,prod(V(1).dim(1:3)));

  for ivol = 1:nvol
    fid = fopen(deblank(P(ivol,:)),'rb','l');
    indata(ivol,:) = fread(fid,inf,'int16')';
    fclose(fid);
  end % for ivol
end

if COMPUTE_ICA
  % Run through the ICA algorithm
  disp('Running ICA algorithm')
  [weights sphere] = runica(indata);
end

if SAVE_ICA_DATA
  disp('Saving data')
  save(sprintf('%s_ica.mat',sub_id), 'weights','sphere')
end

if LOAD_ICA_DATA
  disp('Loading ICA data')
  load(sprintf('%s_ica.mat',sub_id))
end

if COMPUTE_COMPONENTS
  disp('Computing components')
  components = weights * indata;
end

if WRITE_COMPONENTS
  disp('Writing components')
  for ivol = 1:nvol
    outname = fullfile(outdir,sprintf('%s_comp%04d',sub_id,ivol));
    Vo(ivol) = V(ivol);
    Vo(ivol).fname = outname;
    Vo(ivol).descrip = sprintf('ICA component: %d', ivol);

    tmp = components(ivol,:);
    if WRITE_Z
      tmp = (tmp-mean(tmp))/std(tmp);
      max_z(ivol) = max(tmp);
      min_z(ivol) = min(tmp);
    end
    Vo(ivol) = spm_write_vol(Vo(ivol), reshape(tmp,Vo(ivol).dim(1:3)));
  end
end