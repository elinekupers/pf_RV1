function out = envelope(T)

%envelope fits the best fitting 2D Gaussian to a target image T.

[X,Y] = meshgrid(-floor(size(T,2))/2:floor(size(T,2))/2-1,-floor(size(T,1))/2:floor(size(T,1))/2-1);
xdata = zeros(size(X,1),size(Y,2),2);
xdata(:,:,1) = X; xdata(:,:,2) = Y;
x0 = [1,0,50,0,50,0];

TT = abs(T);
Q = lsqcurvefit(@D2GaussFunctionRot,x0,xdata,TT);
envelope =  D2GaussFunctionRot(Q,xdata);

%Output
out = envelope./max(envelope(:)); %normalize the envelope

