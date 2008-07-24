% design_batch.m

global rootpath

type = {'defaults_edit','model','contrasts','headers',...
'means','realign','coreg','normalize','smooth'};


analyses = struct('type',[],'index',[],'work_dir',[],'mfile',[]);

work_dir = {rootpath, ...
      };

mfile = ...
    { './design_model.m' ...
      };

na = 1;

analyses.type(na) = 2;
analyses.index(na) = 1;
analyses.work_dir(na) = 1;
analyses.mfile(na) = 1;

