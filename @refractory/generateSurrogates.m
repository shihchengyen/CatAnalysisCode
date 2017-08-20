function surrogates = generateSurrogates(obj,varargin)
%@refractory/generateSurrogates Generate surrogate spike trains using q(t)
%   SURROGATES = generateSurrogates(OBJ,VARARGIN) generates surrogate
%   spike trains for all cells using the free firing rate, q(t),
%   and the recovery function, w(t). The SURROGATES structure
%   contains the following fields:
%      SURROGATES{CELLN,SETN}{REPN}(SPIKEN) - spike time
%         for set SETN, cell CELLN, repetition REPN and spike number
%         SPIKEN.
%
%   The optional input arguments are:
%      cells - number or array of numbers indicating which cells to
%              use (default obj.data.nCells).
%      sets - number of sets of surrogates to generate (default is
%             1000).
%
%   RASTERS = generateSurrogates(OBJ,'cells',1:obj.data.nCells,'set',1000)

Args = struct('cells',1:obj.data.nCells,'sets',1000);

Args = getOptArgs(varargin,Args);

% pre-allocate memory
surrogates = cell(length(Args.cells),Args.sets);

for i = 1:Args.sets
	for j = Args.cells
		% arrange cells in rows with each set being one column
		fprintf('Computing set %i cell %i\n',i,j);
		surrogates{j,i} = getSurrogateSpikes(obj.data,j,varargin{:});
	end
end
