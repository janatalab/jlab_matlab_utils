function ss = get_stable_sets(cdata, scs)
% ss = get_stable_sets(cdata, scs)
%
% Returns sets of components that belong together.
%
% Figures out all of the conjunction chains working through the cdata
% matrices.  
%
% scs - stable components in each session
%

% 10/19/05 Petr Janata - extracted from earlier script

nsess = length(scs);

stable_ics = find(scs{1});
nstable = length(stable_ics);

ss = {};
for istable = 1:nstable
  start_ic_idx = stable_ics(istable);
  fprintf('Conjoining ICs starting with session 1 IC %d\n', start_ic_idx);

  curr_ic_idx = start_ic_idx;

  curr_chain_idx = 1;
  ss{istable} = [];
  ss{istable}(1,:) = [1 curr_ic_idx];
  for isess = 1:nsess
    jsess = mod(isess,nsess)+1;
    tmp = cdata{isess,jsess}(curr_ic_idx,:)' .* repmat(scs{jsess},1,length(curr_ic_idx));
    curr_ic_idx = find(any(tmp,2));
    nic = length(curr_ic_idx);
    if nic
      start_idx = size(ss{istable},1)+1;
      stop_idx = start_idx+nic-1;
      ss{istable}(start_idx:stop_idx,:) = [ones(nic,1)*jsess,curr_ic_idx(:)];
    else
      break
    end
  end % for isess
  
  % Now get the unique ICs
  ss{istable} = unique(ss{istable},'rows');
  
end % for istable