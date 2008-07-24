% MYSQL - Interact with a MySQL database server
%
%   If no output arguments are given on the left, then display results.
%   If arguments are given, then return requested data silently.
%
%   mysql( 'open', host, user, password )
%      Open a connection with specified parameters, or defaults if none
%            host:  default is local host. Use colon for port number
%            user:  default is Unix login name.
%        password:  default says connect without password.
%
%       Examples: mysql('open','arkiv')       %  connect on port 0
%                 mysql('open','arkiv:2215')
%
%   mysql('close')
%      Close the current connection
%
%   mysql('use',db)  or   mysql('use db')
%      Set the current database to db   Example:  mysql('use cme')
%
%   mysql('status')  or  mysql()  or mysql   (no arguments)
%      Display information about the connection and the server.
%      Return    0     if connection is open and functioning
%             nonzero  if something is not correct (see code for details)
%
%   mysql( query )
%      Send the given query or command to the MySQL server
%
%      With no output arguments on the left side, display the result
%      If arguments are given on the left, then each argument
%          is set to the column of the returned query.
%      Dates and times are converted to Matlab format: dates are
%          serial day number, and times are fraction of day.
%      String variables are returned as cell arrays.
%
%      Example:
%      [ t, p ] = mysql('select time,price,askbid from cme.sp
%                  where date="1997-04-30" and expir like "1997-06-%"');
%         (but be sure to put quoted text all on one input line)
%      Returns time and price for trades on the June 1997 contract
%      that occured on April 30, 1997.
%
%   Multiple connections:  The program can maintain up to 10
%   independent connections. Any command may be preceded by a
%   connection id -- an integer from 0 to 9 -- to apply the
%   command to that connection. Default id is zero.
%   Example:
%      mysql('open','host1')    %  connection 0 to host 1
%      mysql(5,'open','host2')  %  connection 5 to host 2
%      mysql                    %  status of all connections
