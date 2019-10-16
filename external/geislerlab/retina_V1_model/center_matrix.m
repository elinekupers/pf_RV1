function out = center_matrix(A,B)

%center_matrix centers matrix A in matrix B such that the centers of both matrices align. A must fit within B.

%Matrix centers
szA = size(A); szB = size(B); cxB = center(szB(1)); cyB = center(szB(2));

%Center A in B, creating matrix C
C = B;
C(cxB-floor(szA(1)/2):cxB-floor(szA(1)/2)+szA(1)-1,cyB-floor(szA(2)/2):cyB-floor(szA(2)/2)+szA(2)-1) = A;

%Output
out = C;


