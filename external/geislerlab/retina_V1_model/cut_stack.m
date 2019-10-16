function out = cut_stack(stack, cx, cy, sx, sy)

%cut_stack takes stack with dimensions [n,m,k] and outputs a new stack with dimensions [sx,sy,k]. Each image in the 
%output stack is centered at (cx,cy) in the corresponding layer of the input stack.

%Select the desired portion of the stack
fl_x = cx-center(sx)+1; fl_y = cy-center(sy)+1;
cl_x = cx+ceil(sx/2)-1; cl_y = cy+ceil(sy/2)-1;

%Output
out = stack(fl_x:cl_x,fl_y:cl_y,:);