%-------------------------------------------------------------------------
%
% Solve the time-independent Schroedinger equation to get eigenstates and
% -energies. The Hamiltonian is first constructed in pseudospectral
% representation using the DVR approximation for the potential energy
% operator. Diagonalizing the Hamiltonian yields the eigenfunctions
% directly in coordinate space.
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2011 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function eigen()

global hamilt psi space time

%% Initial setup, output etc.

util.disp(' ')
util.disp('*******************************************************')
util.disp('Solve TISE by direct diagonalization of the Hamiltonian')
util.disp('*******************************************************')

%%  Set default values for hamilt.eigen and psi.eigen
if ~isfield (psi,'eigen')
    psi.eigen = [];
end

if ~isfield (psi.eigen,'start')
    psi.eigen.start = 0; % ground state
end

if ~isfield (psi.eigen,'stop')
    psi.eigen.stop = 10; % 10th excited state
end

if ~isfield(hamilt,'eigen')
    hamilt.eigen = []; 
end

if ~isfield(hamilt.eigen, 'cutoff')
    hamilt.eigen.cutoff = 0; % no cutoff
end

if ~isfield(hamilt.eigen, 'storage')
    hamilt.eigen.storage = 'f'; % full
end

hamilt.eigen.size = space.size.n_tot * hamilt.coupling.n_eqs;
   
%% Output parameters
util.disp(['Size of matrix                : ' int2str(hamilt.eigen.size) ' x ' int2str(hamilt.eigen.size)])
if hamilt.eigen.cutoff == 0
    util.disp('No cutoff value');
else
    util.disp(['Cut-off absolute values below : ' num2str(hamilt.eigen.cutoff)])
end
if hamilt.eigen.storage == 'f'
    util.disp('Storing full matrix');
elseif hamilt.eigen.storage == 's'
    util.disp('Storing sparse matrix');
else
    util.error('Wrong storage scheme. hamilt.eigen.storage has to be "f" (full) or "s" (sparse)');
end
util.disp('   ')
util.disp(['First eigenstate investigated : ' num2str(psi.eigen.start)])
util.disp(['Last  eigenstate investigated : ' num2str(psi.eigen.stop)])

% Do some dummy setup to avoid crashes. Over the course of the rewrite, parts of
% the code have been adapted to the TDSE only. When reimplementing the TISE, we
% either have to settle this completely, or as a quick fix, create some dummy
% variables.
time.sub.delta = 1e-10;
time.main.n = psi.eigen.stop - psi.eigen.start + 1;
time.main.grid   = (psi.eigen.start : psi.eigen.stop)';
time.main.total  = psi.eigen.stop;


%% Set up the Hamiltonian matrix in DVR space
util.disp (' ')
util.disp('Start setting up matrix ...')
tim = cputime;

if hamilt.eigen.storage == 'f'
    hamilt.eigen.matrix = zeros(hamilt.eigen.size);
else
    hamilt.eigen.matrix = sparse(hamilt.eigen.size, hamilt.eigen.size);
end

% Create kinetic energy matrix from both the grid's internal kinetic energy
% operators and external kinetic energy operators. If we solve the TISE for
% coupled equations, the resulting matrix should have a blockwise-diagonal form;
% we fill only the left upper square and copy the result to the other blocks.
totsize = space.size.n_tot;
for k = 1:space.size.n_dim
    util.disp([num2str(k) '. kinetic energy'])
    hamilt.eigen.matrix(1:totsize, 1:totsize) = ...
        hamilt.eigen.matrix(1:totsize,1:totsize) + kinetic2dvr(space.dof{k});
end

if isfield(hamilt, 'kin')
    for ii = 1:numel(hamilt.kin)
        util.disp([num2str(ii) '. custom kinetic energy operator']);
        hamilt.eigen.matrix(1:totsize, 1:totsize) = ...
            hamilt.eigen.matrix(1:totsize,1:totsize) + kinetic2dvr(hamilt.kin{ii});
    end
end

for m = 2:hamilt.coupling.n_eqs
    hamilt.eigen.matrix((m-1)*totsize+1:m*totsize, (m-1)*totsize+1:m*totsize) = ...
        hamilt.eigen.matrix(1:totsize, 1:totsize);
end

% Add the potential energy and the diabatic coupling elements
% Since the potential energy is not given as a diagonal matrix, but as a grid, we have
% to make a matrix out of it.
util.disp('potential energy')

for m = 1:hamilt.coupling.n_eqs
    for n = m:hamilt.coupling.n_eqs
        if isempty(hamilt.pot.grid_ND{m,n})
            continue;
        end

        pot = reshape(hamilt.pot.grid_ND{m,n}, space.size.n_tot, 1);
        if hamilt.eigen.storage == 's'
            pot = sparse(pot);
        end

        hamilt.eigen.matrix((m-1)*totsize+1:m*totsize, (n-1)*totsize+1:n*totsize) = ...
            hamilt.eigen.matrix((m-1)*totsize+1:m*totsize, (n-1)*totsize+1:n*totsize) ...
            + diag(pot);

        if m ~= n
            % we know that there are only diabatic couplings
            hamilt.eigen.matrix((n-1)*totsize+1:n*totsize, (m-1)*totsize+1:m*totsize) = ...
                diag(pot);
        end
    end
end

% Cut-off
hamilt.eigen.matrix(abs(hamilt.eigen.matrix) < hamilt.eigen.cutoff) = 0;

util.disp (['Finished after [CPU seconds] : ' num2str(cputime-tim)])
util.disp (' ')

%% Enable symmetry restrictions. So far, this is restricted to a 1-D grid.
if isfield(hamilt.eigen, 'symmetry')

    % Total number (N) of grid points should be even!
    N = space.size.n_tot;
    if mod(N,2) 
        util.error ('Symmetry adaption only for even number of grid points')
    end
    
    % Mapping the full problem to a symmetry-adapted basis
    switch class(space.dof{1})
        
        case {'grid.legendre','grid.hermite'}
            main = eye ( N/2);
            anti = fliplr(main);
            switch lower(hamilt.eigen.symmetry)
                case 'g' % N/2 rows in transformation matrix
                    util.disp ('Even (g) parity (A_g irrep of C_i group) - GAU-LEG-HER')
                    hamilt.eigen.transform = [main anti];
                case 'u' % N/2 rows in transformation matrix
                    util.disp ('Even (g) parity (A_g irrep of C_i group) - GAU-LEG-HER')
                    hamilt.eigen.transform = [main -anti];
                otherwise
                    util.error ('Symmetry should be "g" (=even) or "u" (=odd) for GAU-LEG-HER')
            end
            
        case 'grid.fft' % Note absence of grid point at space.dof{1}.r_max
            switch hamilt.eigen.symmetry
                case 'g' % N/2+1 rows in transformation matrix
                    util.disp ('Even (g) parity (A_g irrep of C_i group) - FFT')
                    main = eye ( N/2-1 );
                    anti = fliplr(main);
                    init = zeros(1,N); init(1)=1;
                    last = zeros(1,N); last(N/2+1)=1;
                    fill = zeros(N/2-1,1);
                    hamilt.eigen.transform = [init; fill main fill anti; last];
                case 'u' % N/2-1 rows in transformation matrix
                    util.disp ('Odd (u) parity (A_u irrep of C_i group) - FFT')
                    main = eye ( N/2-1 );
                    anti = fliplr(main);
                    fill = zeros(N/2-1,1);
                    hamilt.eigen.transform = [fill main fill -anti];
                case 'A1'
                    util.disp ('1st (A_1) irrep of C_2v group - FFT')
                    main = eye ( N/4-1 );
                    anti = fliplr(main);
                    init = zeros(1,N); init(1)=1; init(N/2+1)=1;
                    last = zeros(1,N); last(N/4+1)=1;last(3*N/4+1)=1;
                    fill = zeros(N/4-1,1);
                    hamilt.eigen.transform = [init; fill main fill anti fill main fill anti; last];
                case 'B1'
                    util.disp ('2nd irrep of C_2v group - FFT')
                    main = eye ( N/4-1 );
                    anti = fliplr(main);
                    init = zeros(1,N); init(1)=1; init(N/2+1)=-1;
                    fill = zeros(N/4-1,1);
                    hamilt.eigen.transform = [init; fill main fill -anti fill -main fill anti];
                case 'B2'
                    util.disp ('3rd irrep of C_2v group - FFT')
                    main = eye ( N/4-1 );
                    anti = fliplr(main);
                    last = zeros(1,N); last(N/4+1)=1;last(3*N/4+1)=-1;
                    fill = zeros(N/4-1,1);
                    hamilt.eigen.transform = [fill main fill anti fill -main fill -anti; last];
                case 'A2'
                    util.disp ('4th irrep of C_2v group - FFT')
                    main = eye ( N/4-1 );
                    anti = fliplr(main);
                    fill = zeros(N/4-1,1);
                    hamilt.eigen.transform = [fill main fill -anti fill main fill -anti];
                otherwise
                    util.error ('Symmetry should be "g" (=even) or "u" (=odd) or A1,A2,B1,B2 for FFT')
            end
            
    end
    
    % normalize rows of transformation matrices
    for j=1:size(hamilt.eigen.transform,1)
        hamilt.eigen.transform(j,:) = hamilt.eigen.transform(j,:) ...
            / norm(hamilt.eigen.transform(j,:));
    end
    
    % for more than one surface, we just create a block-diagonal matrix, where each
    % block looks the same, and transforms the wave function on the corresponding
    % surface to the reduces representation.
    if hamilt.coupling.n_eqs > 1
        newtrafo = zeros(space.size.n_tot/2 + hamilt.coupling.n_eqs, space.size.n_tot);
        for m = 1:hamilt.coupling.n_eqs
            nrows = size(hamilt.eigen.transform,1);
            newtrafo((m-1)*nrows+1:m*nrows, (m-1)*N+1:m*N) = hamilt.eigen.transform;
        end
        hamilt.eigen.transform = newtrafo;
    end
    
    % transform Hamiltonian matrix
    hamilt.eigen.matrix = hamilt.eigen.transform*hamilt.eigen.matrix*hamilt.eigen.transform';
end


%% Find eigenvalues and eigenvectors of Hamiltonian matrix; Sort them!
tim = cputime;

if hamilt.eigen.storage == 'f'
    util.disp('Start diagonalizing full matrix ...')
    [ hamilt.eigen.eig_vecs, hamilt.eigen.eig_vals ] = eig ( hamilt.eigen.matrix );
else
    util.disp(['Density of matrix :' num2str(nnz(hamilt.eigen.matrix) / space.size.n_tot^2)])
    util.disp('Start diagonalizing sparse matrix ...')
    opts.disp = 0;
    [ hamilt.eigen.eig_vecs, hamilt.eigen.eig_vals ] = ...
        eigs ( hamilt.eigen.matrix, psi.eigen.stop+1, 'sm', opts );
end

util.disp('Sorting eigenvalues...')

[sorted, order] = sort(diag(hamilt.eigen.eig_vals));
hamilt.eigen.eig_vals = sorted;
sortvecs = zeros(size(hamilt.eigen.eig_vecs));
for k= 1:length(order)
    sortvecs(:, k) = hamilt.eigen.eig_vecs(:, order(k));
end
hamilt.eigen.eig_vecs = sortvecs;


util.disp (['Finished after [CPU seconds] : ' num2str(cputime-tim)])


%% List eigenvalues
util.disp(' ')
util.disp('Table of eigenvalues')
for ii=psi.eigen.start : psi.eigen.stop
    util.disp([ int2str(ii) ' : ' num2str(util.real(hamilt.eigen.eig_vals(ii+1)))])
end
util.disp(' ')
