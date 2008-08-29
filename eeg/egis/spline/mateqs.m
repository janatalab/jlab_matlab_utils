function [p,q,error_check]= mateqs2(w,x,y,z,v,k,kinv,a,ainv,e)

n = size(x,2);
mm = 10;
et = e';
kv = kinv*v';

ev = et*kv;

q = ainv*ev;

eq = e*q;

keq = kinv*eq;

p = kv - keq;

kp = k*p;

verror = v' - (kp +eq);

etp = et*p;

error_check(1) = sum(verror);
error_check(2) = sum(etp);
