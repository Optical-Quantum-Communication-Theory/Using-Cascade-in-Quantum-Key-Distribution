function [newParams, modParser]= BB84ChannelFuncForCascade(params, options, debugInfo)
% BB84ChannelFunc a simple channel function for a qubit based BB84 protocol
% with no loss. This channel model allows for depolariztion and
% misalignment between Alice and Bob's detectors.Here, Shmidt decomposition
% was used to shrink Alice from a 4d space to a 2d space. With this, we can
% model Alice and Bob through a maximally entangled state, and apply the
% channel modely directly to them. 
%
% Used for the Cascade project. We also include a "replacement" channel
% that replaces the state at the output of Alice's lab with the H state
% with some fixed probability. This is interesting because it breaks
% certain symmetries in the statistics, and we need to do this for the
% Cascade project. 
%
% Input parameters:
% * dimA: The dimension of Alice's system. In this channel it is assumed to
%   be 2.
% * dimB: The dimension of Bob's system. In this channel it is assumed to
%   be 2.
% * observablesJoint: The joint observables from Alice and Bob's
%   measurments which they perform on the (idealy) max entangled state. The
%   observables must be hermitian and each must be the size dimA*dimB by
%   dimA*dimB. The observables assume the spaces are ordered A \otimes B.
%   They also should be positive semi-definite and should sum to identity,
%   but this is hard to check because of machine precision issues.
% * depolarization: The amount of depolarization applied to the signal
%   Alice sends to Bob. At maximum depolarization (depolariztion =1) a pure
%   qubit state is converted to a maximally mixed state. Depolarization
%   should be between 0 and 1.
% * misalignmentAngle: Angle Alice and Bob's bases are misaligned by around
%   the Y-axis. For example, Bob's detectors could be slightly rotated away
%   from the incoming signals. Although calculations are done on the Bloch
%   sphere, angles should not be given in that form. Give them as half, the
%   angle of rotation on the bloch sphere. Angles must be real numbers.
% Output parameters:
% * expectationsJoint: The joint epxectations for Alice and Bob's
%   measurement of the signals. Simply formed by taking the
%   observablesJoint and applying them to a simulated rhoAB.
% Options:
% * errorHandling: (global option) See makeGlobalOptionsParser for details.
% DebugInfo:
% * rhoAB: Alice and Bob's shared density matrix after the channel has
%   acted on it. Usefull for checking the channel has been applied
%   correctly.
%
% See also QKDChannelModule, makeGlobalOptionsParser
arguments
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser("BB84ChannelFunc");
optionsParser.parse(options);
options = optionsParser.Results;

%% module parser

modParser = moduleParser("BB84ChannelFunc");
modParser.addRequiredParam("observablesJoint",@(x) allCells(x,@ishermitian));

modParser.addRequiredParam("dimA",@(x) x==2);
modParser.addRequiredParam("dimB", @(x) x ==2);
modParser.addRequiredParam('replaceProb', @(x) mustBeInRange(x,0,1));


modParser.addOptionalParam("depolarization",0,@(x) mustBeInRange(x,0,1));
modParser.addOptionalParam("misalignmentAngle",0,@mustBeReal);

modParser.addAdditionalConstraint(@(observables,dimA,dimB) allCells(observables,@(x) size(x,1) == dimA*dimB),["observablesJoint","dimA","dimB"])

modParser.parse(params);

params = modParser.Results;

%% simple setup
newParams = struct();

try

    dimA = params.dimA;
    dimB = params.dimB;

    %% generate the density matrix shared between Alice and Bob.

    % generate the maximally entangled density matrix for these as Alice
    % and are the same dimension.
    rhoAB = 0;
    for index = 1:dimA
        rhoAB = rhoAB + kron(zket(dimA,index),zket(dimB,index));
    end

    rhoAB = (rhoAB*rhoAB')/dimA;
    
    rhoAB=PerturbationChannel(rhoAB, params.replaceProb);

    %depolarize
    rhoAB = depolarizationChannel(rhoAB,params.depolarization,dimA);

    %apply misalignment around the Y axis on the bloch sphere.
    %(Normally this formula has theta/2 to convert an angle on the bloch sphere
    %to a regular angle, but we can skip that as we are already feeding a
    %regular angle into the formula).
    rotMat = cos(params.misalignmentAngle)*eye(2)...
        - 1i*sin(params.misalignmentAngle)*[0,1i;-1i,0];
    rotMat = kron(eye(dimA),rotMat);
    rhoAB = rotMat*rhoAB*rotMat';


    %save transfered state to the debugInfo
    debugInfo.storeInfo("rhoAB",rhoAB);



    %% generate the expectation values
    expectationsJoint = zeros(size(params.observablesJoint));

    for index = 1:numel(params.observablesJoint)
        expectationsJoint(index) = trace(params.observablesJoint{index}'*rhoAB);
    end

    newParams.expectationsJoint = expectationsJoint;
catch err
    ErrorHandling.handle(options.errorHandling, err, debugInfo);
    return
end
end


% Applies a depolarization channel to
function rhoPrime = depolarizationChannel(rho, depol,dimA)
% Nielsen & Chuang p. 379
% depolarization for a 2d system.

% thanks Scott!
krausOps = cell(1,4);
pauliMatrices = {eye(2), [0,1;1,0], [0,1i;-1i,0], [1,0;0,-1]};
krausOps{1} = kron(eye(dimA), sqrt(1-3/4*depol)*eye(2));
for iOp = 2 : length(pauliMatrices)
    krausOps{iOp} = kron(eye(dimA), sqrt(depol)*pauliMatrices{iOp}/2);
end
rhoPrime = ApplyMap(rho, krausOps);
end

%
function rhoOut=PerturbationChannel(rhoIn, Prob)
    %implements rho -> (1-p) rho + p [0]
    
    op1=kron(eye(2),sqrt(1-Prob)*eye(2));
    op2=kron(eye(2),sqrt(Prob)*zket(2,1)*zket(2,1)');
    op3=kron(eye(2),sqrt(Prob)*zket(2,1)*zket(2,2)');

    rhoOut=ApplyMap(rhoIn,{op1,op2,op3});
end