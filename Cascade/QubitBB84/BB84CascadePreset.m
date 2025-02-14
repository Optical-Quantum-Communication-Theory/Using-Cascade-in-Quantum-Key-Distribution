function qkdInput = BB84CascadePreset()
% BB84CascadePreset A preset that describes a qubit based BB84 protocol
% using source replacement. In this implementation, there is no tranmission
% loss modeled. Alice and Bob have the same probability of choosing the
% Z-basis (and same for the X-basis). Using the singular value
% decomposition, Alice was shrunk from 4d to 2d. This way we can now model
% Alice and Bob as sharing a 2d max entangled state. We include additional
% announcements of the location of errors, which is required by the Cascade
% project. We also have the option to use the replacement channel in the
% channel model. 

qkdInput = QKDSolverInput();

%% Parameters

%Here we start by setting what are the initial parameters we will use with
%the protocol. These parameters help define what needs to be fixed, what we
%want to scan over (usualy for graphing), and what needs to be optimized to
%extract more key.
%
%How do we know what parameters can be set? Each module (the section after
%this) has a series of inputs they request. If they aren't given by a
%previously executed module, then it us up to the user to specify them in
%the preset.

% depolarization: Amount of depolarization Bob's state experiences during
% transmission. For qubits, depolarization corresponds to a shrinking of the
% state on the Bloch sphere. With 0 depolarization a pure state remains on
% the surface of the Bloch sphere. A completely depolarized state
% (depolarization =1) is maximally mixed at the center of the Bloch sphere.

%qkdInput.addScanParameter("depolarization",num2cell(linspace(0,0.2,11)));
qkdInput.addFixedParameter("depolarization",0.1);

% misalignmentAngle: The angle Bob's detectors are misaligned from Alice's.
% This causes a unitary rotation around the Y-axis on the Bloch sphere.
% Note, this parameter is optional, and does not need to be specified.
%qkdInput.addFixedParameter("misalignmentAngle",0);
qkdInput.addScanParameter("misalignmentAngle",num2cell(linspace(0,pi/12,10)));

% pz: The probability of Alice/Bob choosing the Z-basis for transmission / 
% measurement. For this BB84 protocol, if they don't choose the Z-basis
% then they transmit / measure in the X-basis.
% qkdInput.addOptimizeParameter("pz",struct("lowerBound",0.01,"upperBound",0.99,"initVal",0.8));
qkdInput.addFixedParameter("pz",1/2);
% qkdInput.addScanParameter("pz", {1/2, 1/3, 1/4, 1/5, 1/6, 1/8, 1/10, 1/25});

% fullstat: Controls which parts of Alice and Bob's joint statistics are
% used to bound the attack. There are 3 options:
% * -1: Only the QBER and gain statistics are used.
% * 0: Statistics where Alice and Bob chose the same basis choice are
%   used.
% * 1: All joint statistics are used.
% qkdInput.addFixedParameter("fullstat",1);
% qkdInput.addFixedParameter("feedback",1); %whether W (location of errors) is announced or not.
% qkdInput.addFixedParameter("replaceProb",0.2);
%qkdInput.addFixedParameter("f",1); %we don't include cost of EC in this . So this value is useless. It is included so we can
%use the canonical asymptotic BB84 solver. 

% files = {'simpleBB84_1.mat'   'simpleBB84_2.mat'    'simpleBB84_3.mat'    'simpleBB84_4.mat' 'simpleBB84_5.mat'    'simpleBB84_6.mat'    'simpleBB84_7.mat'    'simpleBB84_8.mat'    'simpleBB84_9.mat'  'simpleBB84_10.mat'   'simpleBB84_11.mat' };
% files = convertCharsToStrings(files);
% qkdInput.addScanParameter("filepath", num2cell(files'));%{"exportTest.mat"});

%% modules

% modules are at the core of the software's design philosophy. There are 5
% core types of modules used. KeyRate, mathSolver, channel, description,
% and optimizer modules. We give a brief description of them here but
% please open the function and the module for a more detailed description.
% (right click and either view help or open the function and modules directly)
%

% This module works with the keyRate module to describe the class of
% protocols that the keyRate module can solve. Usually code goes here
% because it is likely to change to work on flavors of other protocols or
% when parts of the protocol are usefull for a channel model to have access
% to.
% This description only provides the joint observables so that the channel
% module can use them. It's useful when we want to swap out what our
% measuremnts are.
descriptionModule = QKDDescriptionModule(@BB84DescriptionFuncForCascade);
qkdInput.setDescriptionModule(descriptionModule);

% The channel model provides the expectation values used to bound what the
% state Alice, Bob (and the Eve) share. Here, this channel module requires
% the joint observables, and produces the joint expectation values.
channelModule = QKDChannelModule(@BB84ChannelFuncForCascade);
% channelModule = QKDChannelModule(@DataImportChannel);
qkdInput.setChannelModule(channelModule);

% The key rate module contains the proof techniques to produce a safe lower
% bound on the key rate. A large function of the key rate module is to
% reformulate all the inputs into a form that can be sent to the mathSolver
% module to determine the minimum relative entropy.
keyRateModule = QKDKeyRateModule(@BB84AsymptoticQubitKeyRateFuncForCascade);
qkdInput.setKeyRateModule(keyRateModule);

% The optimizer module is designed to tweak perameters to increase your
% keyrate. It's not used in this protocol, though use cases can include
% tweaking the intensity of coherent pulses used by Alice.
optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",0),struct("verboseLevel",0));
qkdInput.setOptimizerModule(optimizerMod);


% The mathSolver module takes in a description of linear equality, linear
% ineqaulity, trace norm constraints, Kraus operators (for the G map) and
% projectors for the key map (also known as the Z map). It then determine
% the worst case senario and produces the minimum relaive entropy between
% the key and Eve's information.
mathSolverOptions = struct();
mathSolverOptions.initmethod = 1; %closesest to maximally mixed.
mathSolverOptions.maxiter = 300; %number of iterations that should be performed
mathSolverOptions.maxgap = 1e-6;
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

%% global options
% options used everywhere. Many modules may have options given that
% overwrite these, but these are the basic options used everywhere.
% From the documentation on makeGlobalOptionsParser:
% Currently the global options are:
% * cvxSolver (default "SDPT3"): String naming the solver CVX should use.
% * cvxPrecision (default "high"): String that CVX uses to set the solver
%   precision.
% * verboseLevel (default 1): Non-negative integer telling the program how
%   much information it should display in the command window. 0, minimum; 1
%   basic information; 2, full details, including CVX output.
% * errorHandling (default 1): Integer 1, 2, or 3, detailing how the
%   program should handle run time errors. 1: catch and warn the user. The
%   key rate for the point is set to 0. 2: catch but don't warn the user.
%   The key rate for the point is set to 0. 3: don't catch the error and
%   let it up the stack.
qkdInput.setGlobalOptions(struct("errorHandling",ErrorHandling.DontCatch,"verboseLevel",1,"cvxSolver","mosek"));
end