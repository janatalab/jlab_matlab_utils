function ha = add_fighdr(pp)
% ha = add_fighdr(pp);
%
% Adds an axes to a figure for the purpose of adding header information.
% Returns a handle to the axes.
%

% Need to add some field checking and default handling
try tfs = pp.titlefontsize; catch tfs = 14; end

vert_offset = 0.95;
ha = axes('position', [0 vert_offset 1 1-vert_offset], 'units','norm');
axis('off')

% Add the title
text(0.5,0.5, pp.title, 'fontsize', tfs, 'horizontalalign', 'center');

% Add a timestamp
text(0,1, datestr(now), 'horizontalalign','left','verticalalign','top')
