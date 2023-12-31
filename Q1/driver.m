clear all; clc;

kappa = 1.0; % isotropic homogeneous heat conductivity

% manufactured solution and source term
exact   = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);

f = @(x,y) -2*x*(x-1)-2*y*(y-1);

% quadrature rule
n_int_xi  = 10;              % number of quadrature points in xi-direction
n_int_eta = 10;              % number of quadrature points in eta-direction
n_int_h   = 10;              % number of quadrature points on the Neumann boundary
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);
[xih,weighth]     = Gauss(n_int_h,-1,1); % Used for the Neumann boundary integration

% FEM mesh settings
n_en = 4;                   % 4-node quadrilateral element

ele = [10,20,30,100];        % Different number of elements
for i = 1:length(ele)
n_el_x = ele(i);                 % number of element in x-direction
n_el_y = ele(i);                 % number of element in y-direction
n_el   = n_el_x * n_el_y;   % total number of element in 2D domain

n_np_x = n_el_x + 1;        % number of node points in x-direction
n_np_y = n_el_y + 1;        % number of node points in y-direction
n_np   = n_np_x * n_np_y;   % total number of node points in 2D domain

% generate the coordinates of the nodal points
x_coor = zeros(n_np, 1);
y_coor = zeros(n_np, 1);

hh_x = 1 / n_el_x;          % mesh size in the x-direction
hh_y = 1 / n_el_y;          % mesh size in the y-direction

for ny = 1 : n_np_y
  for nx = 1 : n_np_x
    index = (ny-1)*n_np_x + nx; % nodal index for the (nx, ny) node
    x_coor(index) = (nx-1) * hh_x;
    y_coor(index) = (ny-1) * hh_y;
  end
end

% setup the IEN array for element with local node numbering as
% a=4 ------- a=3
% |           |
% |           |
% |           |
% |           |
% a=1 ------- a=2
IEN = zeros(n_en, n_el);
for ey = 1 : n_el_y
  for ex = 1 : n_el_x
    ee = (ey-1)*n_el_x + ex;
    IEN(1,ee) = (ey-1)* n_np_x + ex;
    IEN(2,ee) = (ey-1)* n_np_x + ex + 1;
    IEN(3,ee) =  ey   * n_np_x + ex + 1;
    IEN(4,ee) =  ey   * n_np_x + ex;
  end
end
% 
% % Q(a) the whole D boundary
% % ID array
% ID = zeros(n_np, 1);
% counter = 1;
% for ny = 2 : n_np_y - 1
%   for nx = 2 : n_np_x - 1
%     ID( (ny-1)*n_np_x + nx ) = counter;
%     counter = counter + 1;
%   end
% end

% Q(b) the D boundary and the N boundary
% ID = 0 if nodes on D boundary, so the L and R sides are zeros
% ID array
ID = zeros(n_np, 1);
counter = 1;
for ny = 1 : n_np_y
  for nx = 2 : n_np_x - 1
    ID( (ny-1)*n_np_x + nx ) = counter;
    counter = counter + 1;
  end
end

% Q(a) D boundary on the whole sides
% n_eq = n_np - n_np_x * 2 - n_np_y * 2 + 4;

% Q(b) D boundary on right and left sides, N boundary on top and bottom
%sides, the 4 corner points are in D boundary only
% the number of equations = number of nodes - number of D nodes
n_eq = n_np - n_np_y * 2;

% Construct two matrixes to store the D and N boundary data
% 1 ----0-----1                 1 ----1-----1
% 1 ----0-----1                 0 ----0-----0
% 1    g_b    1                 0    h_b    0
% 1 ----0-----1                 0 ----0-----0
% 1 ----0-----1                 0 ----0-----0
% 1 ----0-----1                 1 ----1-----1
g_b = ones(n_np,1); % set g = 1 in this case
for ny = 2 : n_np_y - 1  
    for nx = 2 : n_np_x - 1 
        g_b((ny-1) * n_np_x + nx) = 0; % Make the L and R sides to be 1, the remaind to be 0
    end
end

h_b = ones(n_np,1); % set h = 1 in this case
for ny = 2 : n_np_y - 1 
    for nx = 1 : n_np_x
        h_b((ny-1) * n_np_x + nx) = 0; % Make the T and B sides to be 1, the remaind to be 0
    end
end

LM = ID(IEN);

% Start the assembly procedure
K = spalloc(n_eq, n_eq, 9*n_eq);
F = zeros(n_eq, 1);

for ee = 1 : n_el
   k_ele = zeros(n_en, n_en);
   f_ele = zeros(n_en, 1);
%    g_ele = zeros(n_en, 1); %the g boundary on each element
%    h_ele = zeros(n_en, 1); %the h boundary on each element
   Nah   = zeros(n_en, 1); %Na * h on the Neumann boundary of each element
   kg    = zeros(n_en, 1); %k_ab * g_b on the Dirichlet boundary of each element

   x_ele = x_coor( IEN(1:n_en, ee) );
   y_ele = y_coor( IEN(1:n_en, ee) );
   
   g_ele = g_b(IEN(1:n_en,ee));
   h_ele = h_b(IEN(1:n_en,ee));
   
   %Q(a) whole Dirichlet boundary
   %Dirichlet boundary on each element
%    for bb = 1 : n_en
%        g_ele(bb) = g_function(x_ele(bb),y_ele(bb));
%    end

%    % Q(b) Dirichlet boundary on L and R sides
%    % Dirichlet boundary on each element
%    for bb = 1 : n_en
%        g_ele(bb) = g_function_b(x_ele(bb),y_ele(bb));
%    end
   
%    % Q(b) Neumann boundary on T and B sides
%    % Neumann boundary on each element
%    for bb = 1 : n_en
%       h_ele(bb) = h_function_b(x_ele(bb),y_ele(bb));
%    end

   % loop over quadrature points   
   for ll = 1 : n_int
     % we need the geometric mapping at each quadrature points
     x_l = 0.0; y_l = 0.0;
     dx_dxi = 0.0; dx_deta = 0.0;
     dy_dxi = 0.0; dy_deta = 0.0;
     for aa = 1 : n_en
        x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
        y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
        [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));  
        dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
        dx_deta = dx_deta + x_ele(aa) * Na_eta;
        dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
        dy_deta = dy_deta + y_ele(aa) * Na_eta;
     end

     detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

     % we loop over a and b to assemble the element stiffness matrix and load
     % vector
     for aa = 1 : n_en
       [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));

       % See page 147 of Sec. 3.9 of the textbook for the shape function
       % routine details
       Na_x = (Na_xi * dy_deta    - Na_eta * dy_dxi) / detJ;
       Na_y = (Na_xi * (-dx_deta) + Na_eta * dx_dxi)  / detJ;

       f_ele(aa) = f_ele(aa) + weight(ll) * detJ * f(x_l, y_l) * Quad(aa, xi(ll), eta(ll));
       
       for bb = 1 : n_en
         [Nb_xi, Nb_eta] = Quad_grad(bb, xi(ll), eta(ll));
         Nb_x = (Nb_xi * dy_deta    - Nb_eta * dy_dxi) / detJ;
         Nb_y = (Nb_xi * (-dx_deta) + Nb_eta * dx_dxi)  / detJ;

         k_ele(aa,bb) = k_ele(aa,bb) + weight(ll) * detJ * kappa * ( Na_x * Nb_x + Na_y * Nb_y);

       end % end of bb-loop
     end % end of aa-loop
   end % end of quadrature loop
   
   % Calculate the Neumann boundary integration (refer to L7 notes, local/element perspective)
   % Na * h on the Neumann boundary of each element
   % Refer to the 1-D Problem to help construct the following codes
   for ll_h = 1 : n_int_h
       % Bottom boundary of one element
       dx_dxi_bottom = 0.0;
       for aa = 1 : 2
           dx_dxi_bottom = dx_dxi_bottom + x_ele(aa) * Quad_grad(aa,xih(ll_h),-1); % -1 makes the zero in Quad
       end
       for aa = 1 : 2
           Nah(aa) = Nah(aa) + weighth(ll_h) * Quad(aa,xih(ll_h),-1) * h_ele(aa) * dx_dxi_bottom;
       end
       
       % Top boundary of one element
       dx_dxi_top = 0.0;
       for aa = 3 : 4
           dx_dxi_top = dx_dxi_top + x_ele(aa) * Quad_grad(aa,xih(ll_h),1); % +1 makes the zero in Quad
       end
       for aa = 3 : 4
           Nah(aa) = Nah(aa) + weighth(ll_h) * Quad(aa,xih(ll_h),1) * h_ele(aa) * dx_dxi_top;
       end
       
   end
   
   % Refer to the f_ele formulation in L7 notes
   % Plus the Neumann boundary terms as follows
   for aa = 1 : n_en 
       f_ele(aa) = f_ele(aa) + Nah(aa);
   end
   
   % global assembly
   for aa = 1 : n_en
     PP = LM(aa, ee);
     if PP > 0
       F(PP) = F(PP) + f_ele(aa);
       for bb = 1 : n_en
           QQ = LM(bb, ee);
           if QQ > 0
              K(PP, QQ) = K(PP, QQ) + k_ele(aa, bb);
           else
           % do something for non-zero g boundary condition
           F(PP) = F(PP) - k_ele(aa,bb) * g_ele(bb); % Refer to L7 notes, local/element perspective
           end
       end
     end
   end
end % end of element loop

% solve the linear system
d_temp = K \ F;

disp = zeros(n_np, 1);

% insert the solution vector back with the g-data
for ii = 1 : n_np
  index = ID(ii);
  if index > 0
    disp(ii) = d_temp(index);
  else
    % g boundary data where ID = 0; 
    % the value should be the same with g_function
    disp(ii) = g_b(ii); 
  end
end

% save the solution to file
save("FEM_solution", "disp", "n_el_x", "n_el_y", "exact", "exact_x", "exact_y");
save(sprintf('FEM_solution_%d.mat', ele(i)), "disp", "n_el_x", "n_el_y", "exact", "exact_x", "exact_y");

end

% EOF