function an_st = ensemble_qid_x_qid_chi2(data,params)

% Calculates Pearson's chi-squared goodness of fit test for data in an ensemble_qid_x_qid_table
% 
% this script could probably be a bit fancier, but at the moment provides
% simple functionality to calculate goodness of fit between two qids tabulated
% using ensemble_qid_x_qid_table.m
%
% see also http://en.wikipedia.org/wiki/Pearson%27s_chi-square_test 
%
% ARGUMENTS
%  data - the results of ensemble_qid_x_qid_table
%  params - at this point, nothing from params is used in this script
% RETURNS
%  an_st - 
%
% FB 11/20/2007

an_st = {};
na = 0;

dcols = set_var_col_const(data{1}.vars);
eqxq_tbls = data{1}.data{dcols.count};

il = length(eqxq_tbls(:,1,1));
jl = length(eqxq_tbls(1,:,1));

df = (il-1)*(jl-1);
chi2 = 0;

ctbl = zeros(il,jl);
ctbli = zeros(il);
ctblj = zeros(jl);
cttl = sum(sum(sum(eqxq_tbls)));

for ii=1:il
    ctbli(ii) = sum(sum(eqxq_tbls(ii,:,:)));
    for ij=1:jl
        ctblj(ij) = sum(sum(eqxq_tbls(:,ij,:)));
        ctbl(ii,ij) = sum(eqxq_tbls(ii,ij,:));
        chi2e = (ctbli(ii)*ctblj(ij))/cttl;
        chi2 = chi2 + (ctbl(ii,ij)-chi2e)^2/chi2e;
    end
end

na=na+1;
an_st{na} = init_analysis_struct;
an_st{na}.type = 'ensemble_qid_x_qid_chi2';
an_st{na}.vars = {'chi2','df'};
an_cols = set_var_col_const(an_st{na}.vars);
an_st{na}.data = {};
an_st{na}.data{an_cols.chi2} = chi2;
an_st{na}.data{an_cols.df} = df;

function an_st = init_analysis_struct
  an_st.type = '';
  an_st.vars = {};
  an_st.data = {};
  an_st.meta = [];
  
