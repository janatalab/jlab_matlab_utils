function fsf = setup_fsl_con(fsf, ev_names, conlist, Fconlist)
% fsf = setup_fsl_con(fsf, conlist, Fconlist);
%
% Populates fields in the fsf contrast structures given the contrast
% information in conlist and Fconlist.
%
% ev_names - names of explanatory variables
% conlist - cell array in which each element is specifies the name of the ev to
%           use for a t-test
% Fconlist - cell array in which each element is a cell array containing a list
%            of EVs to use for the F-test

% 11/03/05 Petr Janata

fsf.con = create_con('real');  % specify structure

nev = length(ev_names);

ncon = length(conlist);
conmtx = zeros(ncon,nev);
for icon = 1:ncon
  con = create_con('real');
  con.conname = conlist{icon};
  
  targ_col = strmatch(conlist{icon},ev_names,'exact');
  conmtx(icon,targ_col) = 1;
  
  con.con_vect = conmtx(icon,:);
  
  fsf.con(icon) = con;
end

fsf.ncon_real = sum(strcmp('real',{fsf.con.type}));
fsf.ncon_orig = sum(strcmp('orig',{fsf.con.type}));

fsf.Ftest = [];

nF = length(Fconlist);

fsf.nftests_real = nF;

% I don't currently handle orig vs real contrasts very well
Fconmtx = zeros(ncon,nF);
for ifcon = 1:nF
  [conmask, con_idxs] = ismember(Fconlist{ifcon}, {fsf.con.conname});
  if all(conmask)
    for icon = 1:length(con_idxs)
      if isempty(fsf.con(con_idxs(icon)).ftest)
	fsf.con(con_idxs(icon)).ftest = zeros(1,nF);
      end
      fsf.con(con_idxs(icon)).ftest(ifcon) = 1;
    end
  end
end
