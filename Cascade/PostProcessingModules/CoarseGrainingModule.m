function [newObs, newExp, modParser] = CoarseGrainingModule(expectations, params,options,debugInfo)
%COARSEGRAININGMODULE Postprocessing module to perform coarse graining on statistics.
% This module takes in a full table of expectation values and either keeps
% it the same (fullstat=1), trims it down to only events where Alice and
% Bob measure in the same basis (fullstat=0), or reduces it all the way to
% QBER statistics (fullstat=-1).
%
% Inputs:
% * expectations: The original table of expectation values
% * observablesJoint: The corresponding observables
% * fullstat: Denotes which statistics to keep:
%   1: keep all statistics (this module does nothing)
%   0: keep only statistics where Alice and Bob measure in the same basis
%   -1: use only QBER statistics
% * senderBasisCount (optional): the amount of bases Alice sends in
%   (defualt 2)
% * receiverBasisCount (optional): the amount of bases Bob measures in
%   (default 2)
% Outputs:
% * newExp: the resulting table of expectation values after
%   coarse graining
% * newObs: A cell array of the corresponding observables to the
%   coarse grained expectation values

arguments
    expectations (:,:,:) {mustBeNumeric}
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;

%% module parser
modParser = moduleParser(mfilename);
% modParser.addRequiredParam("expectationsJoint", @(x) mustBeNumeric(x));
modParser.addRequiredParam("observablesJoint",@(x) allCells(x,@ishermitian));
modParser.addRequiredParam("fullstat",@(x) mustBeMember(x, [-1,0,1]))
modParser.addOptionalParam("senderBasisCount", 2, @(x) mustBeMember(x, [2,3]));
modParser.addOptionalParam("receiverBasisCount", 2, @(x) mustBeMember(x, [2,3]));
modParser.parse(params);

params = modParser.Results;

extractParameters(params);

% debugInfo = struct();

expectationsJoint = expectations;

try
    % helper function to compute the sum of observables stored in a cell
    cellsum = @(X) sum(cat(3, X{:}), 3);
    [~, ~, nDecoy] = size(expectationsJoint);
    switch fullstat
        case -1 % QBER and gain statistics only
            jointObservables = cell(4,1);
            jointExpectations = zeros(4,1,nDecoy);
            % Z gain
            jointObservables{1} = cellsum({observablesJoint{1:2,1:2}});
            jointExpectations(1,1:nDecoy) = sum(expectationsJoint(1:2,1:2,:), 'all');
            % Z QBER
            jointObservables{2} = observablesJoint{1,2} + observablesJoint{2,1};
            jointExpectations(2,1:nDecoy) = expectationsJoint(1,2,:) + expectationsJoint(2,1,:);
            % X gain
            jointObservables{3} = cellsum({observablesJoint{3:4,3:4}});
            jointExpectations(3,1:nDecoy) = sum(expectationsJoint(3:4,3:4,:), 'all');
            % X QBER
            jointObservables{4} = observablesJoint{3,4} + observablesJoint{4,3};
            jointExpectations(4,1:nDecoy) = expectationsJoint(3,4,:) + expectationsJoint(4,3,:);
            if senderBasisCount == 3 && receiverBasisCount == 3
                % Y gain
                jointObservables{5} = cellsum({observablesJoint{3:4,3:4}});
                jointExpectations(5,1:nDecoy) = sum(expectationsJoint(3:4,3:4,:), 'all');
                % X QBER
                jointObservables{6} = observablesJoint{3,4} + observablesJoint{4,3};
                jointExpectations(6,1:nDecoy) = expectationsJoint(3,4,:) + expectationsJoint(4,3,:);
            end
        case 0 % same basis statistics
            % indices of observables to add from old table
            observableYIdx = [1, 1, 2, 2, 3, 3, 4, 4];
            observableXIdx = [1, 2, 1, 2, 3, 4, 3, 4];
            % locations in new table to put them
            newYIdx = [1, 1, 2, 2, 3, 3, 4, 4];
            newXIdx = [1, 2, 1, 2, 1, 2, 1, 2];
            if senderBasisCount == 3 && receiverBasisCount == 3
                observableYIdx = [observableYIdx, 5, 5, 6, 6];
                observableXIdx = [observableXIdx, 5, 6, 5, 6];
                newYIdx = [newYIdx, 5, 5, 6, 6];
                newXIdx = [newXIdx, 1, 2, 1, 2];
                jointObservables = cell(6,2);
                jointExpectations = zeros(6,2);
            else
                jointObservables = cell(4,2);
                jointExpectations = zeros(4,2);
            end
            % loop through the indices and add each one
            for iElement = 1 : numel(observableYIdx)
                y = observableYIdx(iElement);
                x = observableXIdx(iElement);
                newY = newYIdx(iElement);
                newX = newXIdx(iElement);
                jointObservables{newY,newX} = observablesJoint{y, x};
                jointExpectations(newY,newX,1:nDecoy) = expectationsJoint(y, x, :);
            end
        case 1 % full statistics
            jointObservables = observablesJoint;
            jointExpectations = expectationsJoint;
    end

    newObs = jointObservables;
    newExp = jointExpectations;

catch err
    ErrorHandling.handle(options.errorHandling, err, debugInfo);
    return
end

end