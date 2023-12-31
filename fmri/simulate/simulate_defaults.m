%
% Defaults _ copied from ts1, initially created for Schubert batch analyses
% 

%
% Defaults for analysis 'defaults_edit'
%

%
% NOTE: It is better to set normalization defaults in the script that handles
% all the normalization batching.  If they are set here, any values set in the
% normalization batch script are not used!
%

% Set general defaults

defaults_edit = struct( ...
   'type_area', [1 2], ...
   'index', 	[1 1] ...
);

defaults_edit(2) = struct( ...
    'type_area', 4, ...
    'index', 1 ...
    );

type_area = {'Misc','Printing','Hdr','Statistics', ...
      'Normalisation','RealignCoreg','Reset'};
	   
Misc = struct( ...
    'log_to_file', 1, ...
    'log_file_name', sprintf('analysis_log_%s', datestr(now,1)), ...
    'cmdline', 1, ...
    'grid', 0.1 ...
    );

Printing(1) = struct( ...
   'printing_mode', 1,  ...
   'postscript_filename', sprintf('realign_%s_%s', datestr(now,1), datestr(now,15)),  ...
   'postscript_type', 1,  ...
   'default_printer', 0,  ...
   'printer_name', '',  ...
   'post_type', '-dps',  ...
   'graphics_filename', '',  ...
   'graph_type', 3,  ...
   'print_string', '' ...
);

Printing(2) = Printing(1);
Printing(2).postscript_filename = sprintf('coregister_%s', datestr(now,1));
Printing(2).postscript_type = 2;

Hdr = struct( ...
   'defaults_area',  3, ...
   'image_size_voxels',  '10 20 30', ...
   'voxel_size_mm',  '4 4 4', ...
   'scale',  1, ...
   'data_type',  16, ...
   'offset',  0, ...
   'origin_voxels',  '1 1 1', ...
   'description',  'oufff' ...
   );

RealignCoreg = struct( ...
    'defaults_area', 4, ...
    'separate_combine' , -11, ...
    'create' , 1, ...
    'adjust', 0, ...
    'mask', 1, ...
    'reg_quality', 1.0 ...   
    );

Normalisation(1) = struct( ...
   'defaults',  1, ...
   'estimates',  1, ...
   'custom_estimates',  ones(1,12), ...
   'custom_norm',  -1, ...
   'nonlin_func_nb', 13, ...		%  7=4x5x4
   'func_nb',  0, ...
   'nonlin_ite_nb',  12, ...
   'nonlin_regular',  0.01, ...
   'mask_brain',  0, ...
   'mask_object_brain',  0 ...
);

Normalisation(2).defaults = 0;
Normalisation(2).bounding_box = [1];
Normalisation(2).voxel_sizes = 0;
Normalisation(2).voxel_sizes_custom = [3.5 3.5 5];

Statistics = struct( ...
   'fMRI_T', 	   16, ...
   'fMRI_T0', 	   1, ...
   'F_threshold',  1 ...
);

