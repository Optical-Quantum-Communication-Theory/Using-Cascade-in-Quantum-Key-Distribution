function [newParams,modParser] = BB84DescriptionFuncForCascade(params, options, debugInfo)
% BB84DescriptionFunc A simple description function for a qubit BB84
% protocol with no loss that uses the Schmidt decomposition to turn Alice's
% 4d space of signals sent to a 2d space. This is used for the Cascade
% project and includes additional announcements of the location of errors.
%
% Input parameters:
% * pz: The probability that ALice measures in the Z-basis (for this protocol,
%   it's also the probability that Bob measures in the Z-basis aswell). It
%   must be between 0 and 1.
% * feedback : 0 means W=X + Y is not announced. 1 means it is announced.
% Output parameters:
% * observablesJoint: The joint observables for Alice and Bob's measurement
%   of the signals.
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% Options:
% * errorHandling: (global option) See makeGlobalOptionsParser for details.
% DebugInfo:
% * krausSum: sum_i K^\dagger_i*K_i which should be <= identity for
%   a CPTNI map.
%
% See also QKDDescriptionModule, BB84KeyRateFunc, makeGlobalOptionsParser
arguments
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
modParser.addRequiredParam("pz",@(x) mustBeInRange(x,0,1));
modParser.addRequiredParam("feedback", @(x) mustBeMember(x,[0,1]));
modParser.parse(params)
params = modParser.Results;


%% simple setup
newParams = struct();
try

    pz = params.pz;
    px = 1 - pz;

    ketP = [1;1]/sqrt(2);
    ketM = [1;-1]/sqrt(2);
    ket0 = zket(2,1);
    ket1 = zket(2,2);

    %% dimension sizes of Alice and Bob
    dimA = 2;
    dimB = 2;
    newParams.dimA = dimA;
    newParams.dimB = dimB;

    %% rhoA for source replacement scheme
    newParams.rhoA = eye(dimA)/dimA;

    %% joint obserables
    POVMsA = {pz*diag([1,0]),pz*diag([0,1]),(1-pz)*(ketP*ketP'),(1-pz)*(ketM*ketM')};
    POVMsB = POVMsA;
    newParams.POVMA = POVMsA;
    newParams.POVMB = POVMsB;

    observablesJoint = cell(numel(POVMsA),numel(POVMsB));

    for indexA = 1:numel(POVMsA)
        for indexB = 1:numel(POVMsB)
            observablesJoint{indexA,indexB} = kron(POVMsA{indexA},POVMsB{indexB});
        end
    end

    newParams.observablesJoint = observablesJoint;

    %% Kraus Ops (for G map)
    % A: Alice's system, B: Bob's System, C: announcement register, R:
    % key register.
    % The Kraus operators are matrices from ABC \rightarrow RBC. Here we
    % used an isometry to shrink the Kraus operators from outputing on RABC
    % to just RBC. This lets us save time on computing eigen values later.
    % The factor of pz comes from a \sqrt(pz) from Alice's measurements(from shmidt
    % reduction) and \sqrt(pz) from Bob's measurements.

    if params.feedback==0

        % kraus operator for post-processing G map. The ordering of registers
        % is R, A, B, the two-dimensional announcement register (Alice's & Bob's announcement registers combined after sifting)
        krausOpZ = kron(kron(kron(zket(2,1),sqrt(pz)*(ket0*ket0'))+ kron(zket(2,2),sqrt(pz)*(ket1*ket1')), sqrt(pz) * eye(dimB)), [1;0]); % for Z basis
        krausOpX = kron(kron(kron(zket(2,1),sqrt(px)*(ketP*ketP'))+ kron(zket(2,2),sqrt(px)*(ketM*ketM')),sqrt(1-pz) * eye(dimB)),[0;1]); % for X basis
        krausOps = {krausOpZ, krausOpX};
     
        % components for the pinching Z map
        keyProj1 =kron(diag([1,0]), eye(dimA*dimB*2)); 
        keyProj2 = kron(diag([0,1]), eye(dimA*dimB*2));
        keyMap = {keyProj1, keyProj2};

        %Here we compute sum_i K^\dagger_i*K_i. Which should satisfy sum_i
        %K^\dagger_i*K_i <= I. A.K.A. the Kraus operators represent a
        %completely positive, trace non-increasing linear map.
        krausSum = 0;
        for index = 1:numel(krausOps)
            krausSum = krausSum+krausOps{index}'*krausOps{index};
        end
        debugInfo.storeInfo("krausSum",krausSum);

    elseif params.feedback == 1
        proj0=ket0*ket0';
        proj1=ket1*ket1'; %makers it easier to write stuff. 
        projplus=ketP*ketP';
        projminus=ketM*ketM';
         %kraus operator for post-processing G map. Ordering of registers is
        % R A B basis announcement, Feedback
        krausOpZ_F0=kron(kron(kron(kron(ket0,sqrt(pz)*proj0),sqrt(pz)*proj0)+kron(kron(ket1,sqrt(pz)*proj1),sqrt(pz)*proj1),ket0),ket0); %F=0
        krausOpZ_F1=kron(kron(kron(kron(ket0,sqrt(pz)*proj0),sqrt(pz)*proj1)+kron(kron(ket1,sqrt(pz)*proj1),sqrt(pz)*proj0),ket0),ket1); %F=1

        krausOpX_F0=kron(kron(kron(kron(ket0,sqrt(px)*projplus),sqrt(px)*projplus)+kron(kron(ket1,sqrt(px)*projminus),sqrt(px)*projminus),ket1),ket0); %F=0
        krausOpX_F1=kron(kron(kron(kron(ket0,sqrt(px)*projplus),sqrt(px)*projminus)+kron(kron(ket1,sqrt(px)*projminus),sqrt(px)*projplus),ket1),ket1); %F=1

        krausOps={krausOpZ_F0,krausOpZ_F1,krausOpX_F0,krausOpX_F1};
            
          % components for the pinching Z map
        keyProj1 =kron(diag([1,0]), eye(dimA*dimB*2*2)); 
        keyProj2 = kron(diag([0,1]), eye(dimA*dimB*2*2));
        keyMap = {keyProj1, keyProj2};
    end

    
    %% set key map, kraus ops, and announcements in new parameters
    newParams.krausOps = krausOps;
    newParams.keyMap = keyMap;
    %newParams.announcementsA = annA;
    %newParams.announcementsB = annB;
    %newParams.keyMapFunction = keyMapFunc;

catch err
    ErrorHandling.handle(options.errorHandling, err, debugInfo);
    return
end

end