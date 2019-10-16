function out = create_GC_mosaic(sizeXY)

%create_GC_mosaic creates a ganglion cell mosaic that fully covers a square with edge length sizeXY, in degrees. The
%fovea is at the center of this square. The ganglion cell mosaic represents the locations of ganglion cells in the right
%eye. The distributions are governed by the half eccentricities defined in spacing_fn. The output is a structure that
%includes the field GC_mosaic, which is the precomputed ganglion cell mosaic used in retina_V1_model. Save GC_mosaic as
%GC_mosaic to use in retina_V1_model. The mosaic is automatically displayed.

%OUTPUTS (fields of output structure):
%all_x = floating point x coordinates of all created ganglion cells, including those outside the square region.
%all_y = floating point y coordinates of all created ganglion cells, including those outside the square region.
%x_float = floating point x coordinates of ganglion cells within square region.
%y_float = floating point y coordinates of ganglion cells within square region.
%x_int = integer x coordinates of ganglion cells within square region.
%y_int = integer y coordinates of ganglion cells within square region.
%GC_mosaic = integer ganglion cell mosaic. Save GC_mosaic as GC_mosaic to use in retina_V1_model.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CREATE RINGS OF GANGLION CELLS

%Maximum radius of a circle that inscribes the square region
max_radius = sqrt((sizeXY/2)^2+(sizeXY/2)^2);

%We start with two initial points.
x1 = 0; y1 = 0; x2 = 0; y2 = spacing_fn(x1,y1) + y1; x0 = x2; y0 = y2;

%These vectors store the data.
x_vector = [x1 x2]; y_vector = [y1 y2];

%These variables store data about the current contour.
x_contour = x1; y_contour = y1;
contourIndex = 1; sameContour = true;

%radius keeps track of the size of the map.
radius = y2;

%Loop for creating map begins.
while radius < max_radius      
        
    if sameContour  
        IC = intersect_circle(x1,y1,x2,y2);
        %We find which of a larger number of randomly generated points
        %centered at (x2,y2) has maximum distance from all points within a
        %given radius of (x2,y2).  This prevents artifacts.
        tr = spacing_fn(x2,y2); 
        tx = x1; ty = y1;
        fx = find(x_contour==x1);
        if length(x_contour)>fx
            for tt = fx+1:min(fx+2,length(x_contour))
                if sqrt((x_contour(tt)-x2)^2+(y_contour(tt)-y2)^2) < spacing_fn(x_contour(tt),y_contour(tt))+tr
                    tx=cat(2,tx,x_contour(tt)); ty=cat(2,ty,y_contour(tt));
                end
            end
        elseif length(x_contour) == fx
            tx=cat(2,tx,x_contour(1)); ty=cat(2,ty,y_contour(1));
        end
        
        %We find the intersections between the circle centered at (x2,y2) and all other circles of the points in tx,ty.
        qx=0; qy=0;
        for ii = 1:length(tx)
            ICQ = intersect_circle(x2,y2,tx(ii),ty(ii));
            if sqrt(ICQ.xx1^2+ICQ.yy1^2)>sqrt(qx^2+qy^2)
                if ((y2<=0) && (atan2(ICQ.yy1,ICQ.xx1)>atan2(y2,x2)) && (atan2(ICQ.yy1,ICQ.xx1)<atan2(y2,x2)+pi)) ||...
                        ((y2>0) && (mod(atan2(ICQ.yy1,ICQ.xx1),2*pi)>mod(atan2(y2,x2),2*pi)) &&...
                        (mod(atan2(ICQ.yy1,ICQ.xx1),2*pi)<mod(atan2(y2,x2),2*pi)+pi))
                    qx=ICQ.xx1; qy=ICQ.yy1;
                end
            end
            if sqrt(ICQ.xx2^2+ICQ.yy2^2)>sqrt(qx^2+qy^2)
                if ((y2<=0) && (atan2(ICQ.yy2,ICQ.xx2)>atan2(y2,x2)) && (atan2(ICQ.yy2,ICQ.xx2)<atan2(y2,x2)+pi)) ||...
                        ((y2>0) && (mod(atan2(ICQ.yy2,ICQ.xx2),2*pi)>mod(atan2(y2,x2),2*pi)) &&...
                        (mod(atan2(ICQ.yy2,ICQ.xx2),2*pi)<mod(atan2(y2,x2),2*pi)+pi))
                    qx=ICQ.xx2; qy=ICQ.yy2;
                end
            end
        end
        %We add some noise.
        x2 = qx+randn(1)/15*spacing_fn(qx,qy); y2 = qy+randn(1)/15*spacing_fn(qx,qy);
        
        %Now, we choose (x1,y1) to be the nearest point on the previous contour.
        DC = sqrt((x_contour-x2).^2+(y_contour-y2).^2);
        x1 = x_contour(DC==min(DC)); y1 = y_contour(DC==min(DC));
        
        %We check whether this should be the last point on the contour.
        if y0<=0
            if (sqrt((y0-y2)^2+(x0-x2)^2)<...
                    1.5*spacing_fn(x0,y0)) && (mod(atan2(y0,x0),2*pi)-mod(atan2(y2,x2),2*pi)>0) &&...
                    (mod(atan2(y0,x0),2*pi)-mod(atan2(y2,x2),2*pi)<pi)
                sameContour = false;                
            end
        elseif y0>0
            if (sqrt((y0-y2)^2+(x0-x2)^2)<...
                    1.5*spacing_fn(x0,y0)) && (atan2(y0,x0)-atan2(y2,x2)>0) && (atan2(y0,x0)-atan2(y2,x2)<pi)
                sameContour = false;                
            end
        end
                
    elseif ~sameContour
        %We store coordinates of all points on the previous contour.
        x_contour = x_vector(contourIndex+1:end); y_contour = y_vector(contourIndex+1:end);
                
        %We space the points on the previous contour more evenly.
        EC = circle_spacer(x_contour,y_contour);
        x_contour = EC.CX; y_contour = EC.CY;
        x_vector(contourIndex+1:end) = x_contour; y_vector(contourIndex+1:end) = y_contour;
        
        contourIndex = length(x_vector);    
        
        newRand = floor(length(x_contour)/3);
        if newRand ~= 1
            x1 = x_contour(newRand); y1 = y_contour(newRand);
            x2 = x_contour(newRand-1); y2 = y_contour(newRand-1);
        else
            x1 = x0; y1 = y0; 
        end
        
        IC = intersect_circle(x1,y1,x2,y2);
        
        %We choose the point that is further away from the origin.
        D1 = sqrt(IC.xx1^2+IC.yy1^2); D2 = sqrt(IC.xx2^2+IC.yy2^2);
        if D1>D2             
            x2 = IC.xx1; y2 = IC.yy1; 
        elseif D2>D1 
            x2 = IC.xx2; y2 = IC.yy2; 
        else
            error('D1=D2')
        end          
        x0 = x2; y0 = y2; %This is the first point on the next contour.        
        
        sameContour = true;
    end   
    
    %We store the data.
    radius = min(sqrt(x_contour.^2 + y_contour.^2));
    if radius > max_radius
        x_vector = x_vector(1:end-1); y_vector = y_vector(1:end-1);
    else
        x_vector = cat(2,x_vector,x2); y_vector = cat(2,y_vector,y2);        
    end         
            
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INTEGER COORDINATE MOSAIC

%Create a ganglion cell mosaic with integer coordinates at 120 pixels per degree
ppd = 120; %fixed pixels per degree
xv = round(x_vector*ppd); yv = round(y_vector*ppd); %pixel coordinates

%Find the data that lies in the desired square.
xv(x_vector>sizeXY/2)=NaN; yv(x_vector>sizeXY/2)=NaN;
xv(x_vector<-sizeXY/2)=NaN; yv(x_vector<-sizeXY/2)=NaN;
xv(y_vector>sizeXY/2)=NaN; yv(y_vector>sizeXY/2)=NaN;
xv(y_vector<-sizeXY/2)=NaN; yv(y_vector<-sizeXY/2)=NaN;
xv=xv(~isnan(xv)); yv=yv(~isnan(yv));

%Create the ganglion cell mosaic
GCmapX = xv+floor(sizeXY/2*ppd)+1; GCmapY = yv+floor(sizeXY/2*ppd)+1;
SFbitMap = zeros(sizeXY*ppd+1);
idx = sub2ind(size(SFbitMap),GCmapX,GCmapY);
SFbitMap(idx)=1;
SFbitMap = rot90(SFbitMap); %Satisfies e2 = half eccentricities for right eye, defined in spacing_fn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT

out.all_x = x_vector; out.all_y = y_vector; %floating point coordinates of all created ganglion cells
out.x_float = xv; out.y_float = yv; %floating point coordinates of all ganglion cells within square region
out.x_int = cast(xv,'int32'); out.y_int = cast(yv,'int32'); %integer coordinates within square region
out.GC_mosaic = cast(SFbitMap,'int32'); %ganglion cell mosaic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIGURES

figure(1)
plot(out.all_x,out.all_y,'k.') 
axis equal
title('All created ganglion cells'); xlabel('degrees'); ylabel('degrees')

figure(2)
plot(out.x_int,out.y_int,'k.')
axis([min(out.x_int),max(out.x_int),min(out.y_int),max(out.y_int)])
title('Ganglion cell mosaic'); xlabel('pixels'); ylabel('pixels')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SUBROUTINE

function out_IC = intersect_circle(a,b,c,d)

%intersect_circle solves for the intersection of two circles whose centers are (a,b) and (c,d)

r = spacing_fn(a,b); 
s = spacing_fn(c,d);

e = c-a; %difference in x coordinates
f = d-b; %difference in y coordinates
p = sqrt(e^2+f^2); %distance between centers
k = (p^2+r^2-s^2)/(2*p); %distance from center 1 to line joining points of intersection

%Solutions:
xx1 = a+e*k/p+(f/p)*sqrt(r^2-k^2);
yy1 = b+f*k/p-(e/p)*sqrt(r^2-k^2);

xx2 = a+e*k/p-(f/p)*sqrt(r^2-k^2);
yy2 = b+f*k/p+(e/p)*sqrt(r^2-k^2);

%Output
out_IC.xx1 = xx1;
out_IC.xx2 = xx2;
out_IC.yy1 = yy1;
out_IC.yy2 = yy2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SUBROUTINE

function out_C = circle_spacer(cx,cy)

%circle_spacer takes in a vector of coordinates (cx,cy) and outputs a modified circle in which the points are more
%equally spaced.  It is important that any interval that is too small (if there is any) be the interval between the last
%element (cx(end),cy(end)) and the first element (cx(1),cy(1)).

[Th,R] = cart2pol(cx,cy);

%This is the angle between the 1st and last coordinates.
if (Th(1)<-pi/2) && (Th(end)>pi/2)
    Th_1e = mod(Th(1),2*pi)-Th(end);
else
    Th_1e = Th(1)-Th(end);
end

%This is the angle between the 1st and 2nd coordinates.
if (Th(2)<-pi/2) && (Th(1)>pi/2)
    Th_21 = mod(Th(2),2*pi)-Th(1);
else
    Th_21 = Th(2)-Th(1);
end

%This is the increment we subtract from each interval except between Th(1)
%and Th(end).
inc = Th_21-Th_1e;

%We subtract the proper increments.
I = 1:length(cx);
Th = Th-(I-1).*inc./length(cx);
Th(Th<-pi) = Th(Th<-pi) + 2*pi;

[CX,CY] = pol2cart(Th,R);

%Output
out_C.CX = CX;
out_C.CY = CY;