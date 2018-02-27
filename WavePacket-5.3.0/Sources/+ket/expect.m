%--------------------------------------------------------------------------
%
% For given wavefunction (psi) and their represention (space) and for
% given Hamiltonian operator (hamilt), this function returns structures
% (expect and uncert) containing all relevant expectation values and 
% their uncertainties
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function expect ( step )
global expect hamilt psi space uncert

% Loop over (coupled) wavefunctions
for m = 1:hamilt.coupling.n_eqs

    %-----------------------------------------------------------------
    % Expectation values and uncertainties of INDIVIDUAL wavefunctions
    %-----------------------------------------------------------------
        
    % Population 
    expect.ind.pop{m}(step) = sum ( abs(psi.dvr.grid_ND{m}(:)).^2 .* space.dvr.weight_ND(:));

    % If population exceeds pre-specified threshold
    if expect.ind.pop{m}(step) > eps

        % additional multiplicative operators
        if isfield(space, 'amo')
            for p = 1:length(space.amo)
                expect.ind.amo {m}(step,p) = sum (  space.amo{p}.grid_ND{m,m}(:)   .*abs(psi.dvr.grid_ND{m}(:)).^2 .* space.dvr.weight_ND(:) );
            end
        end

        % Transform our wave function in the pure spectral basis. That way,
        % we do never have to care about weights (they are always 0) when
        % getting the "momentum" expectation values.
        
        % Components of position and momentum vectors
        for k=1:space.size.n_dim
            expect.ind.dvr {m}(step,k) = sum ( space.dvr.grid_ND{k}(:)   .*abs(psi.dvr.grid_ND{m}(:)).^2  .* space.dvr.weight_ND(:) ) / expect.ind.pop{m}(step);
            expect.ind.dvr2{m}(step,k) = sum ( space.dvr.grid_ND{k}(:).^2.*abs(psi.dvr.grid_ND{m}(:)).^2  .* space.dvr.weight_ND(:) ) / expect.ind.pop{m}(step);

            % Now integrate in the spectral basis to get the "momentum" expectation values.
            % momentum applies the momentum operator to the input.
            psi.mom = momentum(space.dof{k}, psi.dvr.grid_ND{m});
            expect.ind.fbr {m}(step,k) = sum ( conj(psi.dvr.grid_ND{m}(:)) .* psi.mom(:) .* space.dvr.weight_ND(:) ) / expect.ind.pop{m}(step);
            expect.ind.fbr {m}(step,k) = util.real(expect.ind.fbr {m}(step,k));

            psi.mom = momentum(space.dof{k}, psi.mom);
            expect.ind.fbr2 {m}(step,k) = sum ( conj(psi.dvr.grid_ND{m}(:)) .* psi.mom(:) .* space.dvr.weight_ND(:) ) / expect.ind.pop{m}(step);
            expect.ind.fbr2 {m}(step,k) = util.real(expect.ind.fbr2 {m}(step,k));
        end

        %  Potential/kinetic/all energy (in position/momentum representation)
        expect.ind.pot {m}(step) = sum ( hamilt.pot.grid_ND{m,m}(:)   .*abs(psi.dvr.grid_ND{m}(:)).^2 .* space.dvr.weight_ND(:)) / expect.ind.pop{m}(step);
        expect.ind.pot2{m}(step) = sum ( hamilt.pot.grid_ND{m,m}(:).^2.*abs(psi.dvr.grid_ND{m}(:)).^2 .* space.dvr.weight_ND(:)) / expect.ind.pop{m}(step);
    end

    % See right below... Here we just zero them.
    expect.ind.kin {m}(step) = 0;
    expect.ind.kin2{m}(step) = 0;
end

% Calculate the kinetic energy independently for each dimension and sum up
% (our operator is a sum of operators). Due to the workings of kinetic
% (does all coupled SE's at once), we have to permute the oder of the loops
% If we want to be at least a bit efficient.

for k = 1:space.size.n_dim
    % Calculate T|Psi> and put it in psi.dvr.new_ND; If the kinetic energy was
    % not yet initialised, it is so after this call.
    space.dof{k} = kinetic(space.dof{k}, false);
    for m = 1:hamilt.coupling.n_eqs
        if expect.ind.pop{m}(step) > eps
            expect.ind.kin{m}(step) = expect.ind.kin{m}(step) + sum( conj(psi.dvr.new_ND{m}(:)) .* psi.dvr.grid_ND{m}(:) .* space.dvr.weight_ND(:)) / expect.ind.pop{m}(step);
        end
    end

    % Apply T another time on psi to get <E_kin^2>
    kinetic(space.dof{k}, true);
    for m = 1:hamilt.coupling.n_eqs
        if expect.ind.pop{m}(step) > eps
            expect.ind.kin2{m}(step) = expect.ind.kin2{m}(step) + sum( conj(psi.dvr.new_ND{m}(:)) .* psi.dvr.grid_ND{m}(:) .* space.dvr.weight_ND(:)) / expect.ind.pop{m}(step);
        end
    end
end

% Now if T = T1 + T2 + ..., kin2 contains T1^2, T2^2, ...
% what is missing are the cross terms Ti Tj (i ~= j), which we have to add here
for k = 1:space.size.n_dim
    for l = 1:space.size.n_dim
        if k == l
            continue
        end
        kinetic(space.dof{k}, false);
        kinetic(space.dof{l}, true);
        for m = 1:hamilt.coupling.n_eqs
            if expect.ind.pop{m}(step) > eps
                expect.ind.kin2{m}(step) = expect.ind.kin2{m}(step) + ...
                    sum( conj(psi.dvr.new_ND{m}(:)) .* psi.dvr.grid_ND{m}(:) ...
                    .* space.dvr.weight_ND(:)) / expect.ind.pop{m}(step);
            end
        end
    end
end

% For external kinetic energy operators, the same thing has to be done as well.
if isfield(hamilt, 'kin')
    % First, get T_i|Psi> and T_i^2|Psi>
    for n = 1:length(hamilt.kin)
        hamilt.kin{n} = apply(hamilt.kin{n}, false);
        for m = 1:hamilt.coupling.n_eqs
            if expect.ind.pop{m}(step) > eps
                expect.ind.kin{m}(step) = expect.ind.kin{m}(step) + sum( conj(psi.dvr.new_ND{m}(:)) .* psi.dvr.grid_ND{m}(:) .* space.dvr.weight_ND(:)) / expect.ind.pop{m}(step);
            end
        end

        apply(hamilt.kin{n}, true);
        for m = 1:hamilt.coupling.n_eqs
            if expect.ind.pop{m}(step) > eps
                expect.ind.kin2{m}(step) = expect.ind.kin{m}(step) + sum( conj(psi.dvr.new_ND{m}(:)) .* psi.dvr.grid_ND{m}(:) .* space.dvr.weight_ND(:)) / expect.ind.pop{m}(step);
            end
        end
    end

    % Then get the cross-terms. There are three terms here: grid+external, external+grid,
    % external+external
    for n = 1:length(hamilt.kin)
        for k = 1:space.size.n_dim
            apply(hamilt.kin{n}, false);
            kinetic(space.dof{k}, true);
            for m = 1:hamilt.coupling.n_eqs
                if expect.ind.pop{m}(step) > eps
                    expect.ind.kin2{m}(step) = expect.ind.kin2{m}(step) + ...
                        sum( conj(psi.dvr.new_ND{m}(:)) .* psi.dvr.grid_ND{m}(:) ...
                        .* space.dvr.weight_ND(:)) / expect.ind.pop{m}(step);
                end
            end

            kinetic(space.dof{k}, false);
            apply(hamilt.kin{n}, true);
            for m = 1:hamilt.coupling.n_eqs
                if expect.ind.pop{m}(step) > eps
                    expect.ind.kin2{m}(step) = expect.ind.kin2{m}(step) + ...
                        sum( conj(psi.dvr.new_ND{m}(:)) .* psi.dvr.grid_ND{m}(:) ...
                        .* space.dvr.weight_ND(:)) / expect.ind.pop{m}(step);
                end
            end
        end
    end

    for n = 1:length(hamilt.kin)
        for o = 1:length(hamilt.kin)
            if n == o
                continue;
            end
            apply(hamilt.kin{n}, false);
            apply(hamilt.kin{n}, true);
            for m = 1:hamilt.coupling.n_eqs
                if expect.ind.pop{m}(step) > eps
                    expect.ind.kin2{m}(step) = expect.ind.kin2{m}(step) + ...
                        sum( conj(psi.dvr.new_ND{m}(:)) .* psi.dvr.grid_ND{m}(:) ...
                        .* space.dvr.weight_ND(:)) / expect.ind.pop{m}(step);
                end
            end
        end
    end
end



% Sometimes the expectation value becomes a tiny bit complex (10^-18 or so).
% Remove the imaginary part here.
for m = 1:hamilt.coupling.n_eqs
    expect.ind.kin{m}(step) = util.real(expect.ind.kin{m}(step));
    expect.ind.kin2{m}(step) = util.real(expect.ind.kin2{m}(step));
end


% Now we can "continue" with the original loop
for m = 1:hamilt.coupling.n_eqs
    if expect.ind.pop{m}(step) > eps
        expect.ind.all{m}(step) = expect.ind.pot{m}(step) + expect.ind.kin{m}(step);
    
        % Uncertainties
        uncert.ind.dvr{m}(step,:) = sqrt ( expect.ind.dvr2{m}(step,:) - expect.ind.dvr{m}(step,:).^2  );
        uncert.ind.fbr{m}(step,:) = sqrt ( expect.ind.fbr2{m}(step,:) - expect.ind.fbr{m}(step,:).^2  );

        uncert.ind.pot{m}(step) = sqrt ( expect.ind.pot2{m}(step) - expect.ind.pot{m}(step).^2  );
        uncert.ind.kin{m}(step) = sqrt ( abs(expect.ind.kin2{m}(step) - expect.ind.kin{m}(step).^2));
        % the abs() shall remove spurious tiny imaginary contributions in the result.

        %---------------------------------------------------------------------------------------
        % Expectation values of TOTAL wavefunction as weighted sum over individual wavefunctions
        %---------------------------------------------------------------------------------------
        
        % Spatial projection function
        if isfield(space, 'amo')
            expect.tot.amo(step,:) = expect.tot.amo(step,:) + expect.ind.amo{m}(step,:) * expect.ind.pop{m}(step);
        end
            
        % Position and momentum
        expect.tot.dvr(step,:) = expect.tot.dvr(step,:) + expect.ind.dvr{m}(step,:) * expect.ind.pop{m}(step);
        expect.tot.fbr(step,:) = expect.tot.fbr(step,:) + expect.ind.fbr{m}(step,:) * expect.ind.pop{m}(step);

        % Potential/kinetic/all energy
        expect.tot.pot(step  ) = expect.tot.pot(step  ) + expect.ind.pot{m}(step  ) * expect.ind.pop{m}(step);
        expect.tot.kin(step  ) = expect.tot.kin(step  ) + expect.ind.kin{m}(step  ) * expect.ind.pop{m}(step);
        expect.tot.all(step  ) = expect.tot.all(step  ) + expect.ind.all{m}(step  ) * expect.ind.pop{m}(step);
        
    end
    
    % Population
    expect.tot.pop(step)   = expect.tot.pop(step) + expect.ind.pop{m}(step);
        
end
    
% Normalize spatial projection
if isfield(space, 'amo')
    expect.tot.amo(step,:) = expect.tot.amo(step,:) / expect.tot.pop(step);
end

% Normalize position/momentum
expect.tot.dvr(step,:) = expect.tot.dvr(step,:) / expect.tot.pop(step);
expect.tot.fbr(step,:) = expect.tot.fbr(step,:) / expect.tot.pop(step);

% Normalize potential/kinetic/all energy 
expect.tot.pot(step  ) = expect.tot.pot(step  ) / expect.tot.pop(step);
expect.tot.kin(step  ) = expect.tot.kin(step  ) / expect.tot.pop(step);
expect.tot.all(step  ) = expect.tot.all(step  ) / expect.tot.pop(step);

% Total energy including couplings: Should be a real number
ket.hamilt(0,0,0);
energy = 0;
for m = 1:hamilt.coupling.n_eqs
    energy = energy + sum ( conj(psi.dvr.grid_ND{m}(:)) .* psi.dvr.new_ND{m}(:) .* space.dvr.weight_ND(:) );
end
energy = energy / expect.tot.pop(step);
expect.tot.tot(step) = util.real(energy);
