function [keyRate, modParser] = BB84DecoyKeyRateFuncForCascade(params,options,mathSolverFunc,debugInfo)
% BB84DECOYKEYRATEFUNC A key rate function for use with expectations
% intended for decoy analysis. Given 4x16 detector data conditioned on
% state sent by Alice, this function will squash the statistics, apply
% decoy analysis, construct Kraus operators and key map for the G and Z
% maps, and package observables and expectations appropriately for use with
% an asymptotic inequality solver. Has custom behavious required for the
% Cascade Project.
%
% Input parameters:
% * pz: The probability that Alice measures in the Z-basis (for BB84 protocols,
%   this is also the probability that Bob measures in the Z-basis). It
%   must be between 0 and 1.
% * decoys: a cell of the intensities used in decoy analysis. These are the
%   mean photon numbers that Alice can choose from when performing the
%   decoy protocol
% * deltaLeak: The amount of leakage from error correction.
% * observables: The joint observables from Alice and Bob's
%   measurments. These are organized in a 4x16 table where rows are the
%   state Alice sent and columns are Bob's detector patterns. The
%   observables must be hermitian and each must be the size dimA*dimB by
%   dimA*dimB. The observables assume the spaces are ordered A \otimes B.
%   They also should be positive semi-definite and should sum to identity,
%   but this is hard to check because of machine precision issues.
% * expectationsU, expectationsL: The upper and lower bounds for 
%   conditional expectations from Alice and Bob's measurements. These 
%   should line up with the corresponding observables at each entry.
% * fullstat: Can be -1, 0, or 1. -1 means we use QBER and gain statistics
%   only; 0 means we only use statistics where Alice and Bob are in the
%   same basis; and 1 means we use all expectations and observables.
% Outputs:
% * keyrate: Key rate of the QKD protocol.
% Options:
% * verboseLevel: (global option) See makeGlobalOptionsParser for details.
% * errorHandling: (global option) See makeGlobalOptionsParser for details.
% DebugInfo:
% * krausSum: sum_i K^\dagger_i*K_i. For a completely positive trace
%   non-increasing map, this sum should be <=I. 
%
% See also QKDKeyRateModule, BB84DecoyChannelFunc,
% BB84LossyDescriptionFunc, makeGlobalOptionsParser
arguments
    params (1,1) struct
    options (1,1) struct
    mathSolverFunc (1,1) function_handle
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;

%% module parser
modParser = moduleParser(mfilename);

modParser.addRequiredParam("observablesJoint",@(x) allCells(x,@ishermitian));
modParser.addRequiredParam("expectationsCond",@(x) all(x>=0 & x<=1,"all") );

modParser.addRequiredParam("krausOps", @(x) isCPTNIKrausOps(x));
% modParser.addRequiredParam("keyMap", @(x) mustBeAKeyMap(x));
modParser.addOptionalParam("keyMap",@(x)mustBeA(x,"KeyMapElement"))
% disp(modParser.Parameters)

modParser.addRequiredParam("decoys", @(x) mustBeCellOf(x, 'numeric'));

%modParser.addRequiredParam("f", @(x) x>=1);
%modParser.addRequiredParam("announcementsA", @(x) mustBeNumeric(x));
%modParser.addRequiredParam("announcementsB", @(x) mustBeNumeric(x));

modParser.addRequiredParam("fullstat", @(x) mustBeMember(x, [-1, 0, 1]));

modParser.addRequiredParam("pz",@(x) mustBeInRange(x,0,1));
modParser.addOptionalParam("darkCountRate", @(x) x>=0 && x<1);

modParser.addRequiredParam("rhoA", @(x) abs(trace(x)-1)<=eps);

modParser.parse(params);

params = modParser.Results;

% debugInfo = struct();
% debugInfo.mathSolver = struct();
debugMathSolver = debugInfo.addLeaves("mathSolver");

dimA = 2;
dimB = 3;

try

    mathSolverInput = struct();
    
    %% postprocessing step
    %% Dark count postprocessing
    % simulate dark counts in detectors using a postprocessing map
    %[modifiedExpectations, ~, debugInfo.darkCountPostprocessing] = DarkCountPostprocessingModule(params.expectationsJoint, params, options);
    


    %% Squashing map
    [squashedExpectations, ~] = SquashingModule(params.expectationsCond, params, options, debugInfo);
    %these are actually conditioned expectations.
    

    %% Decoy analysis
    debugDecoyAnalysis = debugInfo.addLeaves("decoyAnalysis");
    [totalCondExpectationsU, totalCondExpectationsL, ~] = DecoyAnalysisModuleForCascade(squashedExpectations, params, options, debugDecoyAnalysis);
    debugDecoyAnalysis.storeInfo("totalCondExpectationsU",totalCondExpectationsU);
    debugDecoyAnalysis.storeInfo("totalCondExpectationsL",totalCondExpectationsL);
    % note that these bounds are computed separately for each entry on the
    % table. Therefore, they commute with coarse-graining.
    
    
    deltaLeak = 0 ;

    %we do decoy on the conditional squashed Expectations..

    %% turn bipartite expectations from conditional to joint probabilities
    pz = params.pz;
    probList = [pz/2, pz/2, (1-pz)/2, (1-pz)/2];
    expectationsU = diag(probList)*totalCondExpectationsU;
    expectationsL = diag(probList)*totalCondExpectationsL;
        
    
    %% Coarse graining
    [observables, expectationsL, ~] = CoarseGrainingModule(expectationsL, params, options, debugInfo);
    [observables, expectationsU, ~] = CoarseGrainingModule(expectationsU, params, options, debugInfo);

    
   
    %% translate for the math solver
    %now we have all the parts we need to get a key rate from the a math
    %solver, but we need to put it into a form it can understand.
    %first we give it the kraus operators for the G map and the projection
    %operators for the key map (Z).
    mathSolverInput.krausOps = params.krausOps;
    mathSolverInput.keyProj = params.keyMap;
    % also include rhoA from the description if it was given
    if ~isequaln(params.rhoA,nan)
        mathSolverInput.rhoA = params.rhoA;
    end

    numObs = numel(observables);

    % inequality constraints = all of Bob's observations
    mathSolverInput.inequalityConstraints = arrayfun(@(x)...
        InequalityConstraint(observables{x},expectationsL(x),expectationsU(x)), 1:numObs);

    [relEnt,~] = mathSolverFunc(mathSolverInput,debugMathSolver);

    %check if the math solver caught an error and we need to return
    if isfield(debugMathSolver.info,"error")
        debugInfo.error = debugMathSolver.error;
        % debugInfo.mathSolver.debugInfo = rmfield(debugInfo.mathSolver.debugInfo,"error");
        keyRate = 0;
        return
    end
    
    %decoy pSignal to modify relative entropy
    pSignal = params.decoys{1}*exp(-params.decoys{1});
    
    keyRate = max(pSignal*relEnt-deltaLeak,0);
    
    if options.verboseLevel>=1
        fprintf("Key rate: %e\n",keyRate);
    end
    
    %set the upper bound as well for debug
    if isfield(debugMathSolver.info,"relEntUpperBound")
        debugInfo.keyRateUpperBound = max(pSignal*debugMathSolver.relEntUpperBound - deltaLeak,0);
    
        if options.verboseLevel>=2
            fprintf("Key rate upper bound: %e\n",debugInfo.keyRateUpperBound)
        end
    end
catch err
    ErrorHandling.handle(options.errorHandling, err, debugInfo);
    keyRate = 0;
    return
end

end



%% functions to perform decoy analysis