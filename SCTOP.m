%%%% TOPOLOGY OPTIMIZATION OF LATTICE MATERIAL %%%%
function SCTOP(nelx,nely,nelz,edgeN,volfrac,penal,rmin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nelx    : Elements number along x axis
% nely    : Elements number along y axis
% nelz    : Elements number along z axis
% edgeN   : Sub-cells number along each direction of an element
% volfrac : Maximum allowed volume fraction
% penal   : Penalization factor for SIMP-based method
% rmin    : Filter radius
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lx = 0.1; ly = 0.1; lz = 0.1; % Unit cell length in each direction
%% BASE MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9*E0;
nu = 0.3;
%% MESH FOR SUB-CELLS
nx = nelx*edgeN;
ny = nely*edgeN;
nz = nelz*edgeN;
%% PREPARE FILTER
iH = ones(nx*ny*nz*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nz
    for i1 = 1:nx
        for j1 = 1:ny
            e1 = (k1-1)*nx*ny + (i1-1)*ny+j1;
            for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nz)
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nx)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),ny)
                        e2 = (k2-1)*nx*ny + (i2-1)*ny+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIAL DESIGN VARIABLES
x = initialDesign(ny,nx,nz);
x_f = zeros(ny,nx,nz);
%% CHECK THE SYMMETRY OF THE INITIAL DESIGN
x_sym = (x + flip(x, 2)) / 2;
x_sym = (x_sym + flip(x_sym, 1)) / 2;
x_sym = (x_sym + flip(x_sym, 3)) / 2;
if(max(abs(x_sym(:)-x(:)))~=0) 
    error('Initial design asymmetry.');
end
%% OPTIMIZATION PARAMETERS
eta = 0.5;
beta = 4;
loop = 0;
change = 1;
maxloop = 200;
%% START ITERATION
iniStart_t = tic;
while change > 0.01 && loop < maxloop
    iniStart = tic;
    loop = loop+1;
    x_f(:) = (H*x(:))./Hs; % filtering
    rho = (tanh(beta*eta)+tanh(beta*(x_f-eta)))/(tanh(beta*eta)+tanh(beta*(1-eta))); % Heaviside projection
    drho = beta*(1-tanh(beta*(x_f-eta)).^2)/(tanh(beta*eta)+tanh(beta*(1-eta))); % derivative of Heaviside projection
    [DH,dDH] = homo3d(rho,lx,ly,lz,E0,Emin,nu,penal,edgeN); % homogenization
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    weight = [0,0,0,0,0,0,0,0,0,-1,-1,-1]; % weight coefficient w_ij for maximum shear
    f = weight*[DH(1,2); DH(2,1); DH(2,3); DH(3,2); DH(1,3); DH(3,1);...
        DH(1,1); DH(2,2); DH(3,3); DH(4,4); DH(5,5); DH(6,6)];
    df = zeros(ny,nx,nz);
    for i = 1:nx*ny*nz
        df(i) = weight*[dDH{i}(1,2); dDH{i}(2,1); dDH{i}(2,3);...
            dDH{i}(3,2);dDH{i}(1,3);dDH{i}(3,1);
            dDH{i}(1,1);dDH{i}(2,2);dDH{i}(3,3);
            dDH{i}(4,4);dDH{i}(5,5);dDH{i}(6,6)] * drho(i);
    end
    df(:) = H*(df(:)./Hs);
    %% CONSTRAINT FUNCTION AND SENSITIVITY ANALYSIS
    dg = zeros(ny,nx,nz);
    dg(:) = H*(drho(:)./Hs);
    %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    l1 = 0; l2 = 1e9; move = 0.1;
    while (l2-l1)/(l2+l1) > 1e-4 && l2 > 1e-40
        lmid = 0.5*(l2+l1);
        xnew = max(0.001,max(x-move,min(1,min(x+move,x.*(max(1e-10,-df./dg/lmid)).^0.3))));
        x_f(:) = (H*xnew(:))./Hs;
        rho = (tanh(beta*eta)+tanh(beta*(x_f-eta)))/(tanh(beta*eta)+tanh(beta*(1-eta)));
        if sum(rho(:)) > volfrac*nx*ny*nz, l1 = lmid; else l2 = lmid; end
    end
    change = max(abs(x(:)-xnew(:)));
    x = xnew;
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4g Vol.:%7.3f maxrho.:%7.3f ch.:%7.3f Time.:%10.3g \n',loop,f, ...
        mean(rho(:)),max(rho(:)),change,toc(iniStart));
end
toc(iniStart_t)
%% PLOT ISOSURFACE
plot_3d(rho,ny,nx,nz);
end

%% FUNCTIONS
function plot_3d(rho,ny,nx,nz)
isovals = shiftdim(reshape(rho,ny,nx,nz),2);
isovals = smooth3(isovals,'box',1);
patch(isosurface(isovals,0.5),'FaceColor',[0 1 0],'EdgeColor','none');
patch(isocaps(isovals,0.5),'FaceColor',[0 1 0],'EdgeColor','none');
view(3); axis equal tight off; camlight;
camlight headlight;
lighting gouraud;
material shiny;
end

%% INITIAL DESIGN 1
function x = initialDesign(ny,nx,nz)
x = 0.4*ones(ny,nx,nz);
for k = 1:nz
    for i = 1:nx
        for j = 1:ny
            if sqrt((i-nx/2-0.5)^2+(j-ny/2-0.5)^2+(k-nz/2-0.5)^2) < min([nx,ny,nz])/3
                x(j,i,k) = 0.001;
            end
        end
    end
end
end

% ======================================================================= %
% A compact and efficient MATLAB code for topology optimization of material
% structures
% Developed by: Chao Li et al.
% Email: chao.li-6@student.uts.edu.au
%
% Main reference:
% Liu, K., Tovar, A.: An efficient 3d topology optimization code 
% written in matlab. Structural and multidisciplinary optimization 50(6), 
% 1175–1196 (2014)
%
% **************************   Disclaimer   ***************************** %
% The authors reserve all rights for the programs. The programs may be
% distributed and used for academic and educational purposes. The authors
% do not guarantee that the code is free from errors, and they shall not be
% liable in any event caused by the use of the program.
% ======================================================================= %