function expectations = coherentSourceChannelForCascade(signalIntensity, eta, detectorEfficiency, misalignmentAngle, pd, basesA, probDistB)
% coherentSourceChannel a function that simulates detection data for Bob
% given the signal Alice sent, based on parameters of the channel and 
% Alice's state. This function is what we use to generate data for decoy 
% analysis; with experimental data, this channel is replaced.
% The only change in Cascade is that misalignmentAngle is defined
% physically instead of on the block sphere 
% Input parameters:
% * signalIntensity: The mean photon number of the signal Alice sent.
% * eta: the transmissivity of the quantum channel; equivalent to 1 - loss.
%   Must be between 0 and 1 inclusive.
% * detectorEfficiency: the efficiency of Bob's detectors. Must be between
%   0 and 1 inclusive
% * misalignmentAngle: Angle Alice and Bob's bases are misaligned by around
%   the Y-axis. For example, Bob's detectors could be slightly rotated away
%   from the incoming signals. 
% * pd: The dark count rate of Bob's detectors
% * basesA: the number of bases Alice sends in (2 for BB84)
% * probDistB: A Dist of probabilities in which Bob measures in. For
%   symmetric BB84, this would be [1/2, 1/2].
% Outputs:
% * expectations: a table of conditional probabilities which has size (2*basesA) by
%   2^(length(probDistB)).
% See also BB84DecoyChannelFunc, BB84DecoyKeyRateFunc
arguments
    signalIntensity (1,1) double
    eta (1,1) double 
    detectorEfficiency (1,1) double
    misalignmentAngle (1,1) double 
    pd (1,1) double
    basesA (1,1) {mustBeInteger} 
    probDistB (1,:) {mustBeProbDist}
end
    % basis choice probabilities
    bProb = probDistB;
    basesB = numel(probDistB); % number of bases Bob measures in

    % total transmittance of the channel
    transmittance = eta*detectorEfficiency;
    intensity = signalIntensity*transmittance;

    % determine the table of Bob's clicks given what Alice sent
    bobClicksGivenAlice = zeros(2*basesA, 2*basesB);
    for basisA = 1 : basesA
        for keyA = 1 : 2
            % construct the equivalent coherent state that Alice sent
            % (loss is built into intensity)
            cohVec = coherentSourceVec(basisA, keyA, intensity);

            % misalignment rotation
            %this does block sphere rotations
            cohVec = coherentRotation(cohVec, 2*misalignmentAngle, [0,0,1]);

            % convert to linear index
            indexA = sub2indPlus(2*ones(1, basesA), [keyA, basisA]); % groups together same basis and alternates key bits
            for basisB = 1 : basesB
                % rotate from basisB to Z
                cohVecTemp = pauliBasis(basisB, false)'*cohVec;

                % account for how much intensity is sent to this detector's
                % basis
                cohVecTemp = sqrt(bProb(basisB))*cohVecTemp;

                for keyB = 1 : 2
                    % convert to linear index
                    indexB = sub2indPlus(2*ones(1, basesB), [keyB, basisB]);

                    % probability the detector doesn't click
                    probNoClick = exp(-cohVecTemp(keyB)*cohVecTemp(keyB)');

                    % probability of clicking/not clicking when we include
                    % dark counts:
                    probClick = 1-probNoClick*(1-pd);

                    bobClicksGivenAlice(indexA, indexB) = probClick;
                end
            end
        end
    end

    % now construct click pattern probabilities

    dimPatterns = 2*ones(1,2*basesB);
    numPatterns = prod(dimPatterns);

    expectations = zeros(2*basesA, numPatterns); 

    % toggle between clickProb and 1-clickProb
    clickProbSwitch = @(click, clickProb) (click==0).*(1-clickProb) + (click~=0).*clickProb;

    for basisA = 1 : basesA
        for keyA = 1 : 2
            indexA = sub2indPlus(2*ones(1,basesA), [keyA, basisA]);
            for patternIndex = 1 : numPatterns
                patternVec = ind2subPlus(dimPatterns, patternIndex)-1;
                expectations(indexA, patternIndex) = prod(clickProbSwitch(patternVec, bobClicksGivenAlice(indexA, :)));
            end
        end
    end
end

function coherentVec = coherentSourceVec(basis,key,intensity)
    %Prepare a coherent state polarization in the given basis and signal state
    coherentVec = sqrt(intensity)*pauliBasis(basis,false)*zket(2,key);
end

function vec = coherentRotation(vec, theta, axisZXY)
    vec = rotMatrix(theta,axisZXY)*vec;
end

function rotMat = rotMatrix(theta,axisZXY)
    %used some of Scott's code as a base for this
    %normalize the axis
    axisZXY = axisZXY/norm(axisZXY);
    
    PauliZ = [1,0;0,-1];
    PauliX = [0,1;1,0];
    PauliY = [0,-1i;1i,0];
    
    rotMat = cos(theta/2)*eye(2)-1i*sin(theta/2)...
        *(axisZXY(1)*PauliZ+axisZXY(2)*PauliX+axisZXY(3)*PauliY);
end