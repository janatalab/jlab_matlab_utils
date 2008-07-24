function [xCon] = reset_xCon(stats_dir)
%  [xCon] = reset_xCon(stats_dir);
%
%  Clobbers and initializes an xCon.mat file that is created during an SPM99
%  analysis.
%
%  stats_dir is the name of the analysis directory that contains the SPM.mat
%  and xCon.mat files.
%
%  This script is particularly useful when using the SPM99 batch mode in a
%  situation where you change the contrasts you are computing and you want to
%  flush out all of the previous contrasts without having to run the whole
%  analysis over again, i.e. you just want to evaluate the updated set of
%  contrasts. 
%
%  Lines that initialize the xCon structure were pulled from the spm_spm.m
%  script.

% 1/1/04 Petr Janata

xCon_fname = fullfile(stats_dir,'xCon.mat');
if exist(xCon_fname)
  fprintf('Deleting ... %s\n', xCon_fname);
      delete(xCon_fname);
end

load(fullfile(stats_dir,'SPM.mat'),'xX');

F_iX0 = struct(	'iX0',		[],...
    'name',		'all effects');
xCon  = spm_FcUtil('Set',F_iX0(1).name,'F','iX0',F_iX0(1).iX0,xX.xKXs);
for i = 2:length(F_iX0)
  xcon = spm_FcUtil('Set',F_iX0(i).name,'F','iX0',F_iX0(i).iX0,xX.xKXs);
  xCon = [xCon xcon];
end

save(xCon_fname,'xCon');
