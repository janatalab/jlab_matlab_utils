% Function to Return Image Volume Local Maxima
% That are, on average, Separated by One FWHM of the Image
% Autocorrelation Function.  Coordinates Given
% in Voxel Units.

% LocalMax = [Voxel Value, X location, Y location, Z location];

function LocalMax=FindLocalMax(image, XDIM, YDIM, ZDIM, Sx, Sy, Sz, Resels, Thresh)
     %fid=fopen(image,'r');
     %X=fread(fid,'float32');
     %fclose(fid);
     
     %X=reshape(image,XDIM,YDIM,ZDIM);

     LocalMax=[];

     LMax=0.;
     tmp=zeros(1,4);
     for x=ceil(Sx/2):ceil(Sx):XDIM-ceil(Sx/2),
	 for y=ceil(Sy/2):ceil(Sy):YDIM-ceil(Sy/2),
	     for z=ceil(Sz/2):ceil(Sz):ZDIM-ceil(Sz/2),
			
		 for i=(-ceil(Sx/2)+1):(ceil(Sx/2)-1),
		     for j=(-ceil(Sy/2)+1):(ceil(Sy/2)-1),
			 for k=(-ceil(Sz/2)+1):(ceil(Sz/2)-1),
							      
			     Index=((z+k-1)*YDIM+(y+j))*XDIM+(x+i);
			     if( (image(Index)>=Thresh) & (image(Index)>LMax) )
                               LMax=image(Index);
                               tmp=[LMax, (x+i), (y+j), (z+k), Index];
                             end

                         end
                     end
                 end

	         if(sum(tmp(2:4))>0)
                    LocalMax=[LocalMax;tmp];
                 end

                 LMax=0.;
                 tmp=[0 0 0 0];

              end
          end
      end

% Sort in Descending Order of the First Column
[a, I]=sort(-LocalMax(:,1));
LocalMax=LocalMax(I,:);


return;


