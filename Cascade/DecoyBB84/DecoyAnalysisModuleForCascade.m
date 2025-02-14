function [bipartiteExpectationsU, bipartiteExpectationsL, modParser] = DecoyAnalysisModuleForCascade(expectations, params, options,debugInfo)
% DecoyAnalysisModule A postprocessing module to perform decoy analysis 
% based on statistics given as bipartite expectations between Alice and 
% Bob. We always use the linear program method for doing deocy. 
% Inputs:
% * decoys: A cell array of the decoy intensities used
% * expectations: An n x m x d table of the expectations
%   produced by Alice and Bob's measurements, where d = the number of
%   decoys and n and m may be any value.
% * fullstat: The fullstat parameter chosen in the protocol. 
%       -1: use QBER and gain statistics only
%       0: use only statistics where Alice and Bob choose identical bases
%       1: use all statistics
% Outputs:
% * bipartiteExpectationsL: a table of lower bounds based on decoy 
%   analysis. Has the same dimensions as expectations.
% * bipartiteExpectationsU: a table of upper bounds  based on decoy 
%   analysis. Has the same dimensions as expectations.
% Options:
% * decoyTolerance: the looseness on decoy bounds that we allow for
%   numerical stability
% * decoySolver: the choice of solver to use for decoy analysis linear 
%   programs
% See also BB84DecoyKeyRateFunc, BB84DecoyPostprocessing
arguments
    expectations (:,:,:) {mustBeNumeric}
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end 

    optionsParser = makeGlobalOptionsParser("decoyAnalysis");
    optionsParser.addOptionalParam("decoyTolerance", 1e-12, @(x) x>=0);
    optionsParser.addOptionalParam("decoySolver", "SDPT3", @mustBeText);
    optionsParser.parse(options);
    options = optionsParser.Results;
    
    modParser = moduleParser(mfilename);
    
    modParser.addRequiredParam("decoys", @(x) mustBeCellOf(x, 'numeric'));
    modParser.addRequiredParam("fullstat", @(x) mustBeMember(x, [-1, 0, 1]));
    
    modParser.parse(params);
    params = modParser.Results;
    
    % debugInfo = struct();
    
    decoys = params.decoys;
    bipartiteExpectationsWCP = expectations;
    fullstat = params.fullstat;
    
    dA = size(bipartiteExpectationsWCP, 1);
    dB = size(bipartiteExpectationsWCP, 2);

    bipartiteExpectationsL = zeros(dA, dB);
    bipartiteExpectationsU = zeros(dA, dB);
    % loop through each squashed expectation value and run the decoy linear
    % programs on it
    parfor iDimA = 1 : dA
        for jDimB = 1 : dB
            decoyExpectations = bipartiteExpectationsWCP(iDimA,jDimB,:);
            [bipartiteExpectationsL(iDimA,jDimB), bipartiteExpectationsU(iDimA,jDimB)] = decoyLP(decoys, decoyExpectations, options);
        end
    end
    bipartiteExpectationsU = max(bipartiteExpectationsU,eps);
    bipartiteExpectationsL = max(bipartiteExpectationsL,0);
end



% function to run the decoy linear programs
function [Y1L, Y1U] = decoyLP(decoys, decoyExpectations, options)
    % parameters for decoy
    nPhoton = 10; % upper bound on photon statistics
    nDecoy = numel(decoys); % # of decoys
    Poisson = @(mu, n) exp(-mu)*mu^n/factorial(n); % anonymous function to query Poisson distributions
    decoyTolerance = options.decoyTolerance; % tolerance for error in decoy LP

    % objective function in linear program -- use this to make sure that
    % only Y1 (element 2) is maximized
    Obj = zeros(1, nPhoton+1);
    Obj(2) = 1; 

    % solve for upper bound
    try 
        cvx_begin quiet
            % initialize cvx solver
            cvx_solver(convertStringsToChars(options.decoySolver));

            variable Y(nPhoton + 1)
            maximize Obj*Y
            for k = 1 : nPhoton+1
                Y(k) <= 1;
                Y(k) >= 0;
            end
            for i = 1 : nDecoy
                C = zeros(1, nPhoton+1);
                Ptotal = 0;
                for k = 1 : nPhoton+1
                    P = Poisson(decoys{i}, k-1);
                    C(k) = P;
                    Ptotal = Ptotal+P;
                end
                C*Y <= decoyExpectations(i) + decoyTolerance;
                C*Y >= decoyExpectations(i) - decoyTolerance - (1-Ptotal);
            end
        cvx_end
        if strcmp(cvx_status, 'Infeasible') | strcmp(cvx_status, 'Failed')
            exception = MException("BB84DecoyChannelFunc:decoyError", ...
                sprintf("Error in decoy analysis linear program! Status: %s", cvx_status));
            throwAsCaller(exception);
        end
    catch ME
        throwAsCaller(ME);    
    end

    Y1U = Y(2);

    % solve for lower bound
    try
        cvx_begin quiet
            % initialize cvx solver
            cvx_solver(convertStringsToChars(options.decoySolver));

            variable Y(nPhoton+1)
            minimize Obj*Y
            for k = 1 : nPhoton+1
                Y(k) <= 1;
                Y(k) >= 0;
            end
            for i = 1 : nDecoy
                C = zeros(1, nPhoton+1);
                Ptotal = 0;
                for k = 1 : nPhoton+1
                    P = Poisson(decoys{i}, k-1);
                    C(k) = P;
                    Ptotal = Ptotal+P;
                end
                C*Y <= decoyExpectations(i) + decoyTolerance;
                C*Y >= decoyExpectations(i) - decoyTolerance - (1-Ptotal);
            end
        cvx_end
        if strcmp(cvx_status, 'Infeasible') | strcmp(cvx_status, 'Failed')
            exception = MException("BB84DecoyChannelFunc:decoyError", ...
                "Error in decoy analysis linear program! Status: %s", cvx_status);
            throwAsCaller(exception);
        end
    catch ME
        throwAsCaller(ME);    
    end

    Y1L = Y(2);
end