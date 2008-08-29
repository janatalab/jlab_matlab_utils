function [fid,fname]=put_fid(perm,tempfilename, title)
%function [fid,fname]=put_fid(perm,tempfilename, title)
%
%This function prompts for a filename and pathname using the ui.
%It check first to see if a non-null string is returned from the ui.
%If a valid filename is entered, the function attempts to open the
%file with the specified permission, and returns a fid and 
%pathname+filename string.
%
%Argument perm is a file permission string (c.f. fopen). Arguments
%mask and title are passed to uigetfile (c.f.).
%

if nargin < 1
	error('Function requires a file permision string as the first argument.');
elseif nargin == 1
	tempfilename='*.*';
	title='Save As:';
elseif nargin == 2
	title='Save As:';
end

[fname, pathname] = uiputfile(tempfilename,title);
if (fname == '')|(fname == 0)
   fid= -1;
   disp('No filename selected. You have to click on a name')
   return;
end

fname = [pathname fname];

[fid, message]=fopen(fname,perm);
if fid == -1
  disp(message)
end
