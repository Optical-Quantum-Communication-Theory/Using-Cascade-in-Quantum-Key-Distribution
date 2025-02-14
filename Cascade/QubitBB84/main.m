%pick the preset file to use In this case we start with a basic qubit BB84
%protocol with no loss. Feel free to open up the module and look at what it
%sets!
%qkdInput = BB84FiniteQubitPreset();
qkdInput = BB84CascadePreset();
% qkdInput = BB84SimpleTestPreset();


filenames = {};

%
% replaceProb - defines perturbation factor, e.g. 0.0 or 0.2
%
% fullstat - defines grain type
% ** -1: Coarse-grained
% ** 0: Sifted Fine-grained
% ** 1: Fine-grained
%
% feedback
% ** 0: F
% ** 1: F'
%

for replaceProb = [0.2]
    for fullstat = -1:1
        if fullstat == -1
            grain = "CoarseGrained";
        elseif fullstat == 0
            grain = "SiftedFineGrained";
        else
            grain = "FineGrained";
        end

        for feedback = 0:1
            if feedback == 0
                fb = "WithoutFeedback";
            else
                fb = "WithFeedback";
            end

            qkdInput.addFixedParameter("replaceProb",replaceProb);
            qkdInput.addFixedParameter("fullstat",fullstat);
            qkdInput.addFixedParameter("feedback",feedback);

            if mod(replaceProb, 1) == 0
                % Integer format if whole number
                replaceProbStr = sprintf("%d", round(replaceProb));
            else
                % Decimal format otherwise
                replaceProbStr = strrep(sprintf("%.1f", replaceProb), '.', '+');
            end


            filename = sprintf(grain+fb+"AndPerturbOf%s.mat", replaceProbStr);
            filenames{end+1} = filename;
            disp(filename);

            %run the QKDSolver with this input
            results = MainIteration(qkdInput);
            save(filename,"results","qkdInput");
        end 
    end
end


%% plot the result
% scaleStyle = "dB-log";
scaleStyle = "linear";
% QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle",scaleStyle,"yScaleStyle",scaleStyle)
QKDPlot.plotParametersFromFiles(filenames,"misalignmentAngle","_keyRate_","xScaleStyle",scaleStyle,"yScaleStyle",scaleStyle)