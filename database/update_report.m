function r = update_report(r, report_str)
%

if ~isfield(r,'msglog')
  r.msglog = {};
end

r.msglog = [r.msglog;{report_str}];

if isfield(r,'report_on_fly') && ~isempty(r.report_on_fly)
  fprintf('%s\n',report_str)
end

if isfield(r,'logfile_fid')
  fprintf(r.logfile_fid,'%s\n',report_str);
end
