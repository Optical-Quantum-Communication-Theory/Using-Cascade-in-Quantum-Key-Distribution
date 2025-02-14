function [keyRate, modParser] = BB84AsymptoticQubitKeyRateFuncForCascade(params,options,mathSolverFunc,debugInfo)
% BB84ASYMPTOTICQUBITKEYRATEFUNC A simple key rate function for a qubit BB84 protocol with
% no loss that uses the Schmidt decomposition to turn Alice's 4d space of
% signals sent to a 2d space. Constructs the Kraus operatros, key map,
% constraints on Alice's system, and selects which joint expectations to
% use to constrain the math solver.
%
% Input parameters:
% * pz: The probability that Alice measures in the Z-basis (for this protocol,
%   it's also the probability that Bob measures in the Z-basis aswell). It
%   must be between 0 and 1.
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * f: error correction effiency. It 1 means we are correcting at the
%   Shannon limit. A more practicle value would be around 1.16.
% * observablesJoint: The joint observables from Alice and Bob's
%   measurments which they perform on the (idealy) max entangled state. The
%   observables must be hermitian and each must be the size dimA*dimB by
%   dimA*dimB. The observables assume the spaces are ordered A \otimes B.
%   They also should be positive semi-definite and should sum to identity,
%   but this is hard to check because of machine precision issues.
% * expectationsJoint: The joint expectations (as an array) from Alice and
%   Bob's measurements that line up with it's corresponding observable in
%   observablesJoint. These values should be betwen 0 and 1.
% * fullstat: Controls which parts of Alice and Bob's joint statistics are
%   used to bound the attack. There are 3 options:
%   * -1: Only the QBER and gain statistics are used.
%   * 0: Statistics where Alice and Bob chose the same basis choice are
%     used.
%   * 1: All joint statistics are used.
% Outputs:
% * keyrate: Key rate of the QKD protocol.
% Options:
% * verboseLevel: (global option) See makeGlobalOptionsParser for details.
% * errorHandling: (global option) See makeGlobalOptionsParser for details.
% DebugInfo:
% * krausSum: sum_i K^\dagger_i*K_i. For a completely positive trace
%   non-increasing map, this sum should be <=I. 
%
% See also QKDKeyRateModule, BB84DescriptionFunc, makeGlobalOptionsParser
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


%% modParser
modParser = moduleParser(mfilename);

modParser.addRequiredParam("observablesJoint",@(x) allCells(x,@ishermitian));
modParser.addRequiredParam("expectationsJoint",@(x) all(or(x>=0,x<=1),"all"));
modParser.addAdditionalConstraint(@isEqualSize,["observablesJoint","expectationsJoint"]);

modParser.addRequiredParam("krausOps", @(x) isCPTNIKrausOps(x));
% modParser.addRequiredParam("keyMap", @(x) mustBeAKeyMap(x));
modParser.addOptionalParam("keyMap",@(x)mustBeA(x,"KeyMapElement"))
% disp(modParser.Parameters)

modParser.addRequiredParam("dimA",@(x) mustBeInteger(x));
modParser.addRequiredParam("dimB", @(x) mustBeInteger(x));
modParser.addAdditionalConstraint(@(observables,dimA,dimB) allCells(observables,@(x) size(x,1) == dimA*dimB),["observablesJoint","dimA","dimB"])

modParser.addRequiredParam("pz",@(x) mustBeInRange(x,0,1));

%modParser.addRequiredParam("f", @(x) x>=1);
modParser.addOptionalParam("fullstat", 1, @(x) mustBeMember(x, [-1,0,1]));
%modParser.addRequiredParam("announcementsA", @(x) mustBeNumeric(x));
%modParser.addRequiredParam("announcementsB", @(x) mustBeNumeric(x));

modParser.addRequiredParam("rhoA", @(x) abs(trace(x)-1)<=eps);
modParser.parse(params);

params = modParser.Results;
% disp(params)

%% simple setup
% deprecatedDebugInfo = struct();
% debugInfo.mathSolver = struct();
debugMathSolver = debugInfo.addLeaves("mathSolver");

try

    mathSolverInput = struct();

    %% Postprocessing step
    % For postprocessing, we need to perform error correction and
    % (optionally) coarse grain the statistics.
    %% Error correction
    %[deltaLeak, ~, ~, debugInfo.errorCorrection] = ErrorCorrectionModule(params.expectationsJoint, params, options);
    deltaLeak = 0;
    
    %% Coarse graining
    [observables, expectations, ~] = CoarseGrainingModule(params.expectationsJoint, params, options, debugInfo);

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

    mathSolverInput.equalityConstraints = arrayfun(@(x)...
        EqualityConstraint(observables{x},expectations(x)),1:numObs);


    % now we call the math solver function on the formulated inputs, with
    % it's options.
    % [relEnt,debugInfo.mathSolver.modParser,debugInfo.mathSolver.debugInfo] = mathSolverFunc(mathSolverInput,mathSolverOptions);
    [relEnt,~] = mathSolverFunc(mathSolverInput,debugMathSolver);
    % relEnt = rand;

    %check if the math solver caught an error and we need to return a key
    %rate of 0.
    if isfield(debugMathSolver.info,"error")
        debugInfo.error = debugMathSolver.error;
        % debugInfo.mathSolver.debugInfo = rmfield(debugInfo.mathSolver.debugInfo,"error");
        keyRate = 0;
        return
    end

    %ensure we don't spit out a negative key rate.
    keyRate = max(relEnt-deltaLeak,0);

    if options.verboseLevel>=1
        fprintf("Key rate: %e\n",keyRate);
    end

    %set the upper bound as well for debuging
    if isfield(debugMathSolver.info,"relEntUpperBound")
        debugInfo.keyRateUpperBound = max(debugMathSolver.relEntUpperBound - deltaLeak,0);

        if options.verboseLevel>=2
            fprintf("Key rate upper bound: %e",debugInfo.keyRateUpperBound)
        end
    end

catch err
    ErrorHandling.handle(options.errorHandling, err, debugInfo);
    keyRate = 0;
    return
end
end