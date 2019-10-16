function out = spacing_fn(xx,yy)

%spacing_fn specifies the average spacing of a ganglion cell from its neighbors when the coordinates of the ganglion 
%cell is at (xx,yy) degrees relative to the fovea. The data used are from Drasdo et al (2007).

%Foveal spacing
s0 = 1/120;

%Eccentricities at which the spacing is half that of the fovea. en=nasal, et=temporal, es=superior, and ei=inferior
%visual field.
en = 1.633922;
et = 1.6666;
es = 1.126081;
ei = 1.488036;

%Spacings for the 4 quadrants
Q1 = s0.*(sqrt(xx.^2./et^2 + yy.^2./es^2) + 1);
Q2 = s0.*(sqrt(xx.^2./en^2 + yy.^2./es^2) + 1);
Q3 = s0.*(sqrt(xx.^2./en^2 + yy.^2./ei^2) + 1);
Q4 = s0.*(sqrt(xx.^2./et^2 + yy.^2./ei^2) + 1);

%Piece together the 4 quadrants
Q=zeros(size(xx));
Q((xx>=0)&(yy>=0)) = Q1((xx>=0)&(yy>=0));
Q((xx<0)&(yy>=0)) = Q2((xx<0)&(yy>=0));
Q((xx<0)&(yy<0)) = Q3((xx<0)&(yy<0));
Q((xx>=0)&(yy<0)) = Q4((xx>=0)&(yy<0));

%Output
out = Q;

