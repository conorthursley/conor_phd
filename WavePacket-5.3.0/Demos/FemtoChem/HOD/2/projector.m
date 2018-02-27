% Copyright (C) 2010 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

% Projects the wave function "wavefunc" onto the requested exit channel.
%
% We project on the part of the space where the O-D bond is longer than the O-H
% bond. Furthermore, we require that the O-D bond has a minimum length to avoid
% contributions from the excited bound states.
function res = projector(wavefunc)

global space

prj = ones(size(wavefunc{1}));
prj(space.dvr.grid_ND{1} < space.dvr.grid_ND{2}) = 0;
prj(space.dvr.grid_ND{1} < 4) = 0;

res = cell(2,1);

res{2} = prj .* wavefunc{1};
res{1} = zeros(size(res{1}));	% the electronic ground state is rather uninteresting
