function F=Fcritical(p,ndf,ddf)

     r=0.001;
     if(ddf>250)
        ddf=250;
     end

     fmax=sqrt(ndf)*p^-0.5;
     f=[0:r:fmax];
     g=betainc((ddf./(ddf+ndf.*f)),ddf/2,ndf/2);
     h=abs(g-p);
     t=find(h==min(h));
     F=f(t);

