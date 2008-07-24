function con=create_con(mode)
% con=create_con(mode);
%
% Creates a contrast structure for specification of FEAT model
%
% mode can either be 'orig' or 'real'

% 04/22/05 PJ

try mode(1);
catch mode='orig';
end

con.type = mode;

% Display images for contrast_real 1
con.conpic=1;

% Title for contrast_real 1
con.conname='';

% Real contrast_real vector 1
% Ultimately has EV elements
con.con_vect=[];

% F-test 1 element 1
con.ftest=[];
