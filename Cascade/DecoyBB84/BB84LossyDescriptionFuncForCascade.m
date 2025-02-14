function [newParams,modParser] = BB84LossyDescriptionFuncForCascade(params, options,debugInfo)
% BB84LossyDescriptionFuncForCascade A simple description function for a qubit BB84
% protocol with loss, using the Schmidt decomposition to turn Alice's
% 4d space of signals sent to a 2d space. This is used with
% BB84LossyKeyRateFunc.
%
% Input parameters:
% * pz: The probability that ALice measures in the Z-basis (for this protocol,
%   this is also the probability that Bob measures in the Z-basis). It
%   must be between 0 and 1.
% * fullstat: If true, protocol will use the full squashed expectations,
%   even when Alice and Bob's basis choices do not line up or when Bob
%   measures nothing. If false, protocol will only use those expectations
%   and observables wherein Alice and Bob chose the same basis.
% Output parameters:
% * observablesJoint: The joint observables for Alice and Bob's measurement
%   of the signals.
% Options:
% * errorHandling: (global option) See makeGlobalOptionsParser for details.
% DebugInfo:
% * krausSum: sum_i K^\dagger_i*K_i which should be <= identity for
%   a CPTNI map.
%
% See also QKDDescriptionModule, BB84DescriptionFunc, makeGlobalOptionsParser
arguments
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

modParser = moduleParser("BB84LossyDescriptionFunc");
modParser.addRequiredParam("pz",@(x) mustBeInRange(x,0,1));
modParser.addRequiredParam("feedback", @(x) mustBeMember(x, [0,1]));
modParser.parse(params)
params = modParser.Results;

pz = params.pz;

ketP = [1;1]/sqrt(2);
ketM = [1;-1]/sqrt(2);

newParams = struct();
try
    %% simple setup
    
    dimA = 2;
    dimB = 3;

    newParams.dimA = dimA;
    newParams.dimB = dimB;
    
    %% generate rhoA
    % in source replacement scheme, this is just the maximally mixed state
    newParams.rhoA = eye(dimA)/dimA;
    
    %% joint obserables
    POVMsA = {pz*diag([1,0]),pz*diag([0,1]),(1-pz)*(ketP*ketP'),(1-pz)*(ketM*ketM')};
    % include vacuum
    % block diagonal structure and order: 2x2 qubit, 1x1 vac
    POVMsB = {pz*diag([1,0,0]),pz*diag([0,1,0]),(1-pz)*([ketP;0]*[ketP;0]'),(1-pz)*([ketM;0]*[ketM;0]'), diag([0,0,1])};
    
    newParams.POVMA = POVMsA;
    newParams.POVMB = POVMsB;

    % add all joint observables; for coarse graining, look in key rate function
    observablesJoint = cell(numel(POVMsA),numel(POVMsB));
    for indexA = 1:numel(POVMsA)
        for indexB = 1:numel(POVMsB)
            observablesJoint{indexA,indexB} = kron(POVMsA{indexA},POVMsB{indexB});
        end
    end
    
    newParams.observablesJoint = observablesJoint;
    
    %% Kraus Ops (for G map)
    pz = params.pz;
    px = 1 - pz;

    if params.feedback==0
        % kraus operator for post-processing G map. The ordering of registers
        % is R, A, B, the two-dimensional announcement register (Alice's & Bob's announcement registers combined after sifting)
        krausOpZ = kron(kron(kron(zket(2,1),sqrt(pz)*diag([1,0]))+ kron(zket(2,2),sqrt(pz)*diag([0,1])), sqrt(pz) * diag([1,1,0]) ), [1;0]); % for Z basis
        krausOpX = kron(kron(kron(zket(2,1),sqrt(px)*ketP*ketP')+ kron(zket(2,2),sqrt(px)*ketM*ketM'),sqrt(px) * diag([1,1,0])  ),[0;1]); % for X basis
        krausOps = {krausOpZ, krausOpX};
     
        % components for the pinching Z map
        keyProj1 =kron(diag([1,0]), eye(dimA*dimB*2)); 
        keyProj2 = kron(diag([0,1]), eye(dimA*dimB*2));
        keyMap = {keyProj1, keyProj2};
    
    elseif params.feedback==1

        ket0=zket(2,1);
        ket1=zket(2,2);
        
        proj0Z=blkdiag(ket0*ket0',[0]);
        proj1Z=blkdiag(ket1*ket1',[0]);
        proj0X=blkdiag(ketP*ketP',[0]);
        proj1X=blkdiag(ketM*ketM',[0]);
    
    
   
    %  
         % kraus operator for post-processing G map. The ordering of registers
        % is R, A, B, the two-dimensional announcement register (Alice's &
        % Bob's announcement registers combined after sifting), W
        %We have additional F register that tells error or not!
        krausOpZW0 =kron( kron(kron(kron(zket(2,1),sqrt(pz)*ket0*ket0'),sqrt(pz)*proj0Z)+ kron(kron(zket(2,2),sqrt(pz)*ket1*ket1'),sqrt(pz)*proj1Z ),ket0),ket0); % for Z basis,F=0
        krausOpZW1 =kron( kron(kron(kron(zket(2,1),sqrt(pz)*ket0*ket0'),sqrt(pz)*proj1Z)+ kron(kron(zket(2,2),sqrt(pz)*ket1*ket1'),sqrt(pz)*proj0Z ),ket0),ket1); % for Z basis,F=1
        
        krausOpXW0 =kron( kron(kron(kron(zket(2,1),sqrt(px)*ketP*ketP'),sqrt(px)*proj0X)+ kron(kron(zket(2,2),sqrt(px)*ketM*ketM'),sqrt(px)*proj1X ),ket1),ket0); % for X basis,F=0
        krausOpXW1 =kron( kron(kron(kron(zket(2,1),sqrt(px)*ketP*ketP'),sqrt(px)*proj1X)+ kron(kron(zket(2,2),sqrt(px)*ketM*ketM'),sqrt(px)*proj0X ),ket1),ket1); % for X basis,F=1
    
        krausOps = {krausOpZW0, krausOpZW1, krausOpXW0, krausOpXW1};
     
        % components for the pinching Z map
        keyProj1 =kron(diag([1,0]), eye(dimA*dimB*2*2)); 
        keyProj2 = kron(diag([0,1]), eye(dimA*dimB*2*2));
        keyMap = {keyProj1, keyProj2};

    end

    
   
    
    krausSum = 0;
    for index = 1:numel(krausOps)
        krausSum = krausOps{index}'*krausOps{index};
    end
    debugInfo.storeInfo("krausSum",krausSum);
    
    

    %% automatic key projection / kraus operator generation 
    %keyMapFunc = @(x,a,b) (a==b)*x;
    %annA = [1,1,2,2];
    %annB = [1,1,2,2,0];

    %% set key map, kraus ops, and announcements in new parameters
    newParams.krausOps = krausOps;
    newParams.keyMap = keyMap;
    %newParams.announcementsA = annA;
    %newParams.announcementsB = annB;
    %newParams.keyMapFunction = keyMapFunc;
    %newParams.blockDimsA = [2];
    %newParams.blockDimsB = [2,1];

catch err
    ErrorHandling.handle(options.errorHandling, err, debugInfo);
    return
end
end