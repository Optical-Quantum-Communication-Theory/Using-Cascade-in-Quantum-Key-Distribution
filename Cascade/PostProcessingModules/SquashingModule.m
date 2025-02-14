function [squashedExpectations, modParser] = SquashingModule(expectations,params,options,debugInfo)
%SQUASHINGMODULE Postprocessing module to perform a simple squashing map.
% The squashing map squashes a (4 or 6) x (16 or 64) table of conditional
% expectation values from detector patterns into a (4 or 6) x (5 or 7) 
% table of conditional expectations corresponding to an equivalent lossy 
% qubit channel.
% Inputs:
% * expectations: The original table of expectation values based on
%   detector patterns
% Outputs:
% * squashedExpectations: The new table of expectations squashed down to a
%   qubit with loss
% See also BB84DecoyChannelFunc, BB84DecoyKeyRateFunc
arguments
    expectations (:,:,:) {mustBeDetectorDataSize}
    params(1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;

%% module parser
modParser = moduleParser(mfilename);
modParser.parse(params);

params = modParser.Results;

% debugInfo = struct();

try
    %% extract dimensions from size of expectation table 
    % (note that the size has already been validated
    [dimA, dimB, ~] = size(expectations);
    
    %% Generate squashing map based on the dimensions

    squashingMapB = generateSquashingMap(dimB);

    if dimA > 6
        % this is entanglement-based, so Alice needs to be squashed too
        squashingMapA = generateSquashingMap(dimA);
    else
        squashingMapA = eye(dimA);
    end
    
    %% apply squashing map(s)
    squashedExpectations = pagemtimes(expectations, squashingMapB);
    squashedExpectations = pagemtimes(squashingMapA', squashedExpectations);

catch err
    ErrorHandling.handle(options.errorHandling, err, debugInfo);
    return
end

end

function squashMap = generateSquashingMap(dimension)
    if dimension == 16
        squashMap = squashingMap16();
    elseif dimension == 64
        squashMap = squashingMap64();
    else
        exception = MException("SquashingModule:dimension", ...
                        sprintf("Error in squashing map generation: unknown expectation table size given: %d", dimension));
                    throwAsCaller(exception);
    end
end

% transforms nx16 detector data into an nx5 p(A,B) table
function mapping = squashingMap16()
    mapping = zeros(16, 5);

    % detector order: 1 | 2 | 3 | 4 | vac
    % where the numbers refer to the signal in the corresponding row of Alice

    % most are just mapped to discarded cross clicks
    mapping(:, 5) = 1;

    % vacuum
    mapping = quickMap(mapping, [0,0,0,0], [0,0,0,0,1]);

    % single clicks
    mapping = quickMap(mapping, [1,0,0,0], [1,0,0,0,0]);
    mapping = quickMap(mapping, [0,1,0,0], [0,1,0,0,0]);
    mapping = quickMap(mapping, [0,0,1,0], [0,0,1,0,0]);
    mapping = quickMap(mapping, [0,0,0,1], [0,0,0,1,0]);

    % double clicks: map to 50/50 0 or 1 inside basis
    mapping = quickMap(mapping, [1,1,0,0], [.5,.5,0,0,0]);
    mapping = quickMap(mapping, [0,0,1,1], [0,0,.5,.5,0]);
end

% transforms nx64 detector data into an nx7 p(A,B) table
function mapping = squashingMap64()
    px = 1/3; % this squasing map only works with symmetric 6-state basis choice
    %probability of discarding
    pas = @(prob,N) (1-prob)^N/(2^(N-1)) + prob^N;
    pdis = 1 - 0.5*(pas(px, 3))/(1-pas(px, 3));

    % create mapping
    mapping = zeros(64, 7);

    % detector order: 1 | 2 | 3 | 4 | 5 | 6 | vac
    % where the numbers refer to the signal in the corresponding row of Alice

    % most are cross clicks
    ccdiscard = ones(1,7)*(1-pdis)/6;
    ccdiscard(7) = pdis;
    assert(sum(ccdiscard)==1, "ccdiscard does not sum to 1!");
    for i = 1 : 64
        mapping(i,:) = ccdiscard;
    end

    % vacuum
    mapping = quickMap(mapping, [0,0,0,0,0,0], [0,0,0,0,0,0,1]);

    % single clicks
    mapping = quickMap(mapping, [1,0,0,0,0,0], [1,0,0,0,0,0,0]);
    mapping = quickMap(mapping, [0,1,0,0,0,0], [0,1,0,0,0,0,0]);
    mapping = quickMap(mapping, [0,0,1,0,0,0], [0,0,1,0,0,0,0]);
    mapping = quickMap(mapping, [0,0,0,1,0,0], [0,0,0,1,0,0,0]);
    mapping = quickMap(mapping, [0,0,0,0,1,0], [0,0,0,0,1,0,0]);
    mapping = quickMap(mapping, [0,0,0,0,0,1], [0,0,0,0,0,1,0]);

    % double clicks: map to 50/50 0 or 1 inside basis
    mapping = quickMap(mapping, [1,1,0,0,0,0], [.5,.5,0,0,0,0,0]);
    mapping = quickMap(mapping, [0,0,1,1,0,0], [0,0,.5,.5,0,0,0]);
    mapping = quickMap(mapping, [0,0,0,0,1,1], [0,0,0,0,.5,.5,0]);
end


% helper function for squashExpectations
function mapping = quickMap(mapping,pattern,remapping)
    mapping(sub2indPlus([2,2,2,2],pattern+1),:) = remapping;
end

% validation function for size of joint expectations
function value = mustBeDetectorDataSize(expectations)
    [m, n, ~] = size(expectations);
    if (m == 4 || m == 6 || m == 16 || m == 64) && (n == 16 || n == 64)
        value = true;
    else
        value = false;
    end
end