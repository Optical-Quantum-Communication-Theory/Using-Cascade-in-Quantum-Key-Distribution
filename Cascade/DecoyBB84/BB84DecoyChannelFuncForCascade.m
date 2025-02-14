function [newParameters, modParser]= BB84DecoyChannelFuncForCascade(params,options,debugInfo)
% BB84DECOYCHANNELFUNC A channel function for BB84 using decoy state
% analysis. Given a collection of decoy intensities, this channel produces
% a group of 4x16 tables of expectations, one for each decoy intensity,
% which are the conditional probability for each of Bob's 16 detector
% patterns given Alice's signal sent. We also optionally include the effect
% of the replacement channel that replaces the state at Alice's output with
% some probability. This is done to break symmetries, for the Cascade
% project.
%
% Input parameters:
% * decoys: a cell of the intensities used in decoy analysis. These are the
%   mean photon numbers that Alice can choose from when performing the
%   decoy protocol
% * eta: the transmissivity of the quantum channel; equivalent to 1 - loss.
%   Must be between 0 and 1 inclusive.
% * detectorEfficiency: the efficiency of Bob's detectors. Must be between
%   0 and 1 inclusive
% * misalignmentAngle: Angle Alice and Bob's bases are misaligned by around
%   the Y-axis. For example, Bob's detectors could be slightly rotated away
%   from the incoming signals. Angles must be real numbers and refer to
%   Bloch sphere rotations.
% Output parameters:
% * expectationsJoint: The conditional expectations (as a 3D array) from 
%   Alice and Bob's measurements. This should be organized as a 4 x 16 x n 
%   array, where n = the number of intensities used in the decoy protocol.
%   In each table, these should line up with the corresponding observables 
%   at each entry.
% Options:
% * errorHandling: (global option) See makeGlobalOptionsParser for details.
% DebugInfo:
% * None.
%
% See also QKDChannelModule,  makeGlobalOptionsParser
arguments
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser("BB84DecoyChannelFunc");
optionsParser.parse(options);
options = optionsParser.Results;

%% module parser
modParser = moduleParser("BB84DecoyChannelFunc");
modParser.addRequiredParam("decoys", @(x) mustBeCellOf(x, 'numeric'));
modParser.addRequiredParam("pz", @(x) mustBeInRange(x, 0, 1));
modParser.addOptionalParam("eta", 1, @(x) mustBeInRange(x, 0, 1));
modParser.addOptionalParam("detectorEfficiency", 1, @(x) mustBeInRange(x, 0, 1));
modParser.addOptionalParam("misalignmentAngle",0,@mustBeReal); %angle on the bloch sphere NOT the satellite.
modParser.addOptionalParam("darkCountRate", 0, @(x) mustBeInRange(x, 0, 1));
modParser.addRequiredParam("PerturbProb",@(x) mustBeInRange(x, 0 ,1));

modParser.parse(params);

params = modParser.Results;

extractParameters(params);

%% simulate the signal intensities passing through the quantum channel
decoys = params.decoys;
nDecoys = length(decoys);

numBasesA = 2; % number of states Alice sends
numBasesB = 2; % number of states Bob can detect

% generate the 4x16 table of detector data for Bob given what Alice sent.
% first row is Alice sending H, 2nd row is Alice sending V, 3rd is Alice
% sending D, and 4th is Alice sending A
rawExpectations = zeros(2*numBasesA, 2^(numBasesB*2), nDecoys);
for iDecoy = 1 : nDecoys
    % coherentSourceChannel function simulates the generation of this data
    rawExpectations(:,:,iDecoy) = coherentSourceChannelForCascade(decoys{iDecoy}, eta, detectorEfficiency, misalignmentAngle, darkCountRate, numBasesA, [pz, 1-pz]);
end

% don't need any debug info
% debugInfo = struct();

%perturb data 
rawExpectations = PerturbData(rawExpectations, PerturbProb);

newParameters.expectationsCond = rawExpectations;

end

% this function just extracts the struct of parameters and drops them in
% the workspace of the caller, so we don't need to do params.{variable}
% each time
function extractParameters(paramStruct)
    fn = fieldnames(paramStruct);
    for iField = 1 : numel(fn)
        fname = string(fn(iField));
        field = paramStruct.(fname);
        if isstruct(field)
            extractParameters(field);
            continue;
        end
        assignin('caller', fname, field)
    end
end




%corresponds to a channel where Eve randomly erases output state of Alice's
%lab and replaces it with the the first signal state.
function expectations= PerturbData(expectations, Prob)
    nDecoy=size(expectations,3);
    nRows=size(expectations,1);

    
    for i=1:nDecoy
        for rowIndex=2:nRows

            expectations(rowIndex,:,i)=(1-Prob)*expectations(rowIndex,:,i)+Prob*expectations(1,:,i);
            
        end
    end
end