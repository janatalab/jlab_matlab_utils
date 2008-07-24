function fh = ensemble_enum_hist(subfunc_name)
% fh = quest_resp_dist(subfunc_name)
%
% Calculates and plots the distributions of the response tallies in the
% different categories of an enum. 
%

% 01/24/07 Petr Janata - adapted from quest_resp_dist.m

if nargin < 1
  subfunc_name = 'calc';
end

% Return a function handle to the nested function that we are going to use
switch subfunc_name
  case 'summary'
    fh = {@calc,@total};
  case 'calc'
    fh = {@calc};
  case 'total'
    fh = {@total};
  case 'per_category'
    fh = {@per_category};
  case 'all'
    fh = {@calc,@total,@per_category};
end

function out = calc(ds,params)
  % Figure out how many question IDs we are dealing with
  qids = params.qids;
  nqid = length(qids);
  
  % Get the column constants
  FD = set_form_col_const(ds.vars);
  
  for iqid = 1:nqid
    % Make the mask for the qid
    
    
  end % for iqid
end  % calc()


if isfield(p,'queue')
 proc_queue = p.queue;
else
  proc_queue = {'calc','report'};
end

if isfield(p, 'is_bitmask')
  is_bitmask = p.is_bitmask;
else
  is_bitmask = 0;
end

pp = p.pp;  % Printing parameters

nq = length(proc_queue);
for iq = 1:nq
  proc_str = proc_queue{iq};
  
  switch proc_str
    case 'calc'

      nsub = size(resp_data,1);

      % Find the index into the data matrix which corresponds to this question ID
      idx = find(as.num.qid == p.qid);
      
      % Extract information about the response categories
      enum_cats = mysql_resolve_enum(as.num.dfid(idx),p.CONN_ID);
      enum_cats = enum_cats{1};
      ncat = length(enum_cats);

      % Extract the relevant data and arrange such that subjects are in columns
      tmpmtx = squeeze(resp_data(:,idx,:))';
      
      % If we are dealing with bitmask data, we have to do something a bit
      % different because the categories are represented directly by category
      % bits being turned off and on, rather than by the number that
      % corresponds to the category.  Thus, we arrive at catcnt in different
      % ways.
      if is_bitmask
	catcnt = zeros(ncat,nsub);
	
	% Replace NaNs with zeros
	
	nstim = zeros(1,nsub);
	% Loop over subjects
	for isub = 1:nsub
	  bitmask = data2bitmask(tmpmtx(:,isub),ncat);
	  nstim(isub) = sum(any(bitmask,2));
	  
	  catcnt(:,isub) = sum(bitmask)';
	end
      else
	% For each subject, determine the number of stimuli that fall into each enum
	% category.
	catcnt = hist(tmpmtx,1:ncat);
      end

      % For each subject, determine the proportion of stims that the subject heard
      % that fell into each category. If we are dealing with a bitmask question
      % in which multiple responses are possible, divide by the number of stims
      % rather than the number of responses so that the proportions reflect the
      % number of stims rather than the number of responses.

      if is_bitmask
	catprop = catcnt./repmat(nstim,ncat,1);
      else
	catprop = catcnt./repmat(sum(catcnt),ncat,1); 
      end				      

      % Summary statistics
      catprop = catprop';  % subjects in rows

    case 'report'
      p.nfig = p.nfig+1;
      figure(p.nfig), clf
      subplot(2,1,1)
      h = bar('v6',1:ncat,mean(catprop));
      add_errorbars(h,(std(catprop)/sqrt(size(catprop,1)))');
      set(gca,'xticklabel',enum_cats)
      title(as.num.qtxt{idx})
      
      subplot(2,1,2)
      hist(tmpmtx(:),1:ncat)
      set(gca,'xticklabel',enum_cats)
      title('Overall counts within each category')

      pp.pagehdr.title = 'Distributions of responses in each category';
      autobio_add_fighdr(pp);

      % Print the figure to a file
      if pp.write2file
	fprintf('Printing figure to file: %s\n', pp.figfname);
	print(pp.figfname, pp.printargs{:})
      end

      % Plot subject distributions of observed proportions
      p.nfig = p.nfig+1;
      figure(p.nfig), clf
      hist_scale = 0:0.1:1;
      for icat = 1:ncat
	subplot(ncat,1,icat)
	hist_vals = hist(catprop(:,icat),hist_scale);
	bar(hist_scale,hist_vals)
	set(gca,'xlim',[-0.05 1.05])
	title(enum_cats{icat})
  
	if isfield(p,'add_nsub') && p.add_nsub
	  text(0.9,0.8,sprintf('N=%d', nsub),'units','norm');
	end  
      end % for icat=
  
      pp.pagehdr.title = sprintf('Distributions of responses to the question:\n%s', as.num.qtxt{idx});
      autobio_add_fighdr(pp);

      % Print the figure to a file
      if pp.write2file
	fprintf('Printing figure to file: %s\n', pp.figfname);
	print(pp.figfname, pp.printargs{:})
      end
  end % switch proc_str
end % for iq

out.catprop = catprop;
out.catcnt = catcnt;
out.enum_cats = enum_cats;
