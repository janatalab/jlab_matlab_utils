function scale_fig_axes(fignum,scale_factor)
% scale_fig_axes(fignum,scale_factor)
%
% Multiplies all axes positions in the figure by the scale factor

% 02/01/07 Petr Janata

kids = get(fignum,'children');
nkids = length(kids);

for ikid = 1:nkids
  if strcmp(get(kids(ikid),'type'),'axes')
    old_pos = get(kids(ikid),'position');
    new_pos = old_pos*scale_factor;  % scale
    % Adjust offsets within page
    new_pos([1 2]) = new_pos([1 2])+(old_pos([1 2])-new_pos([1 2]))/2;
    set(kids(ikid),'position',new_pos)
  end
end % for ikid=
