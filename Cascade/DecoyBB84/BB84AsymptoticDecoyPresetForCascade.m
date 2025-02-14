function qkdInput = BB84AsymptoticDecoyPresetForCascade()
% BB84ASYMPTOTICDECOYPRESET A preset for lossy BB84 using decoy states.
% Instead of computing Alice and Bob's joint state directly (which would
% require photon number truncation), we directly simulate detector data
% based on channel parameters, then squash and perform decoy analysis on
% the detector data in order to generate expectations values to definte our
% acceptance set. This has custom parameters that are required for the
% Cascade project.

qkdInput = QKDSolverInput();

%% Parameters
qkdInput.addScanParameter("eta", num2cell((1:-0.1:0.5)));
%qkdInput.addFixedParameter("eta", 1);
qkdInput.addFixedParameter("misalignmentAngle", 0.2475); %(0.2255 corresponds to QBER of 0.5 for single photons)
%0.1419 is qber 0.02
%qkdInput.addScanParameter("misalignmentAngle", num2cell(linspace(0,0.225,10)));
%qkdInput.addFixedParameter("eta", 1);
% decoys should be in one group, which can be created with these lines:
qkdInput.addFixedParameter("GROUP_decoys_1", 0.5);
qkdInput.addFixedParameter("GROUP_decoys_2", 0.1);
qkdInput.addFixedParameter("GROUP_decoys_3", 0.001);
% end decoys

qkdInput.addFixedParameter("pz",1/2);
qkdInput.addFixedParameter("f",1);

qkdInput.addFixedParameter("darkCountRate", 0);



% description is the same as the lossy qubit description sinRce we squash
% Bob's detector data down to a lossy qubit equivalent
descriptionModule = QKDDescriptionModule(@BB84LossyDescriptionFuncForCascade);
qkdInput.setDescriptionModule(descriptionModule);

% channel model is very different from normal qubit
channelModule = QKDChannelModule(@BB84DecoyChannelFuncForCascade);
qkdInput.setChannelModule(channelModule);

% Key rate module performs squashing and decoy analysis
keyRateOptions = struct();
keyRateOptions.decoyTolerance = 1e-14;
keyRateOptions.decoySolver = "mosek";
keyMod = QKDKeyRateModule(@BB84DecoyKeyRateFuncForCascade, keyRateOptions);
qkdInput.setKeyRateModule(keyMod);

optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",0),struct("verboseLevel",0));
qkdInput.setOptimizerModule(optimizerMod);

% math solver options
mathSolverOptions = struct();
mathSolverOptions.initMethod = 1;
mathSolverOptions.maxIter = 250;
mathSolverOptions.maxGap = 1e-6;
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

% global options
qkdInput.setGlobalOptions(struct("errorHandling",ErrorHandling.DontCatch,"verboseLevel",1,"cvxSolver","mosek","cvxPrecision", "medium"));