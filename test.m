%% EIT project scratch work

% PDEtool function circleg
% returns pdegeom-compatible unit circle boundary points
% reference: 
% https://www.coe.ufrj.br/~borre/pub/Matlab_6.5/CD2/help/toolbox/pde/2exampl2.html
%[p,e,t] = initmesh('circleg','Hmax',0.1);

% define geometry
R1 = [3,4,-0.5,0.5,0.5,-0.5,0.5,0.5,-0.5,-0.5]';
C1 = [1,0,0,1]';
C1 = [C1;zeros(length(R1) - length(C1),1)];
gm = [R1,C1];
sf = 'R1+C1';

ns = char('R1','C1');
ns = ns';
[g, bt] = decsg(gm,sf,ns);

[p,e,t] = initmesh(g, 'Hmax',0.8);

% select only boundary elements from EXTERNAL boundary, 
% see: https://se.mathworks.com/matlabcentral/answers/293821-pdetool-identify-boundary-nodes
e = e(:,find(e(7,:) == 0));

% uncomment once prototype script complete
%[p,e,t]=refinemesh(g,p,e,t);


%% EIT quantities
% sigma : conductivity
% u : voltage
% g : (known/specified/applied) current on boundary
% note, voltage u can be measured on boundary but not inside 
% Neumann-to-Dirichlet map Gamma : maps boundary g to boundary u
% ND map depends on conductivity sigma
% goal of EIT is to infer/approximate sigma by observing g and ND map ?

% specified g current (Neumann bc)
% g = @(x, y) 0; % (DUMMY)
g = @(x, y) cos(4*atan2(y,x)); % g = cos(theta)

% 1-parameter g current (will I use this?)
g1 = @(t) cos(t);

% model conductivity (underlying)
% dummy sigma below (box shape with higher conductivity)
sigma = @(x,y) 1 + 2*((-0.5 < x) & (x < 0.5) & (-0.5 < y) & (y < 0.5));

% how to loop over edges?
% nr of edges = nr of cols in e matrix
Nedge = size(e,2);

% edge b vector needs to be folded into true b vector
% because this problem does not have zero dirichlet bc
bedge = zeros(1,Nedge+1); % dummy
gvec = zeros(1,Nedge+1);

% true b vector
Nnod = size(p,2); % number of nodes
%b = sparse(Ndof,1); % sparse may be efficient
b = zeros(Nnod,1); % placeholder zero vector

% number of triangles
Ntri = size(t,2);

% POSSIBLY INCORRECTLY SPECIFIED BOUNDARY CONDITION
% TRY PLOTTING OR CHANGING

boundaryX = zeros(1,Nedge);
boundaryY = zeros(1,Nedge);
% loop over edges
for i=1:Nedge
    % midpoint rule evaluation of integral
    % simplicity of basis functions -> intersect in 1/2
    
    % coordinates of current edge
    x0 = p(1, e(1,i));
    y0 = p(2, e(1,i));
    x1 = p(1, e(2,i));
    y1 = p(2, e(2,i));

    % midpoints
    xm = (x0+x1)/2;
    ym = (y0+y1)/2;

    boundaryX(i) = xm; % debug
    boundaryY(i) = ym;

    % length of edge
    elen = sqrt((x1-x0)^2 + (y1-y0)^2);

    % contributions of current edge to bedge vector
    % add to correct b vector
    bedge(i) = bedge(i) + g(xm,ym)*elen/2;
    bedge(i+1) = bedge(i+1) + g(xm,ym)*elen/2;

    % add to gvec
    gvec(i) = g(x0,y0);
end
bedge(1) = bedge(1) + bedge(Nedge+1);
bedge = bedge(1:Nedge);
gvec(Nedge+1) = g(x1,y1);

% fold bedge into global b vector
%b(e(1,1:Nedge)) = bedge(:);
for i=1:Nedge
    b(e(1,i)) = b(e(1,i)) + bedge(i);
end

% initial plotting
% display mesh
figure(1)
pdemesh(p,e,t)
hold on
scatter(boundaryX,boundaryY,'black')
hold off

% look at sigma
X = -1:0.05:1;
Y = -1:0.05:1;
for i=1:size(X,2)
    for j=1:size(Y,2)
        Z(i,j) = sigma(X(i),Y(j));
    end
end

figure(2)
surf(Z)

%t = t(1:3, :); % row 4 only 1s

% basis functions
%phi = 
%dphi = 
% copypasted shamelessly from mycourses
% Quadrature weights and positions for reference triangle
wi = [1/6 1/6 1/6]';
% Each column corresponds to sample point coordinates
qi = [1/2 0; 1/2 1/2; 0 1/2]';

% Values of basis functions at each quadrature point
% Each column corresponds to values of a single basis function
phis = [1-qi(1, :)-qi(2, :); qi(1, :); qi(2, :)]';

% Values of derivatives. As the elements are linear, the derivatives are
% the same at all quadrature points.
% Each column corresponds to x- and y-derivatives of a single basis
% function (regardless of quadrature point)
dphis = [-1 1 0; -1 0 1];

% matrix A
A = zeros(Nnod, Nnod);

% loop over triangles
for i=1:Ntri
    
    % nodes of triangle i
    n = t(1:3, i);

    % midpoint of triangle i
    mpx = (p(1,n(1)) + p(1,n(2)) + p(1,n(2)))/3;
    mpy = (p(2,n(1)) + p(2,n(2)) + p(2,n(2)))/3;

    % lookup sigma
    s = sigma(mpx, mpy);

    % copypasted, see P44.m
    % Affine map
    p1 = p(:, n(1));
    p2 = p(:, n(2));
    p3 = p(:, n(3));
    
    Mk = [p2-p1, p3-p1];
    bk = p1;

    % We need the inverse of Mk to integrate the bilinear form.
    IMk = inv(Mk);

    % local A
    Aloc = zeros(3,3);
    for ii=1:3
        for jj=1:3
            % sigma implemented as simple lookup
            % POSSIBLY INCORRECT ASSEMBLY? SINGULAR MATRIX WHEN A\b
            Aloc(ii,jj) = Aloc(ii,jj) + s*abs(det(Mk))*(IMk'*dphis(:,ii)).'*(IMk'*dphis(:,jj))/2;
        end
    end
    
    % Distribute to global A matrix
    for ii = 1:3 % Loop over basis functions of one element
        for jj = 1:3 % Loop over basis functions of one element
            A(n(ii), n(jj)) = A(n(ii), n(jj)) + Aloc(ii, jj);
        end
    end

end

% parameters for voltage approximation
uh = zeros(Nnod, 1);
uh = A \ b; % why is this so singular?

% when evaluated in node points (p matrix) uh gives the function vals
X = p(1,:);
Y = p(2,:);
t = t(1:3,:);
figure(3)
patch(X(t), Y(t), uh(t), uh(t))
colormap(turbo)
colorbar

figure(4)
trisurf(t', X', Y', uh)
colormap(turbo)
colorbar