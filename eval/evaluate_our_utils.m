
function evaluate_our_utils(label_classes, binary_labels, binary_outputs, scalar_outputs)
% clear classes
oldFolder = cd(fileparts(which(mfilename)));
mod = py.importlib.import_module('evaluate_our_utils');
% py.importlib.reload(mod);
res = py.evaluate_our_utils.evaluate_12ECG_score_mod(label_classes, binary_labels, binary_outputs, scalar_outputs);

output_string = sprintf('AUROC|AUPRC|Accuracy|F-measure|Fbeta-measure|Gbeta-measure|Challenge metric\n%.3f|%.3f|%.3f|%.3f|%.3f|%.3f|%.3f',...
                         res{1}, res{2}, res{3}, res{4}, res{5}, res{6}, res{7});
disp(output_string)
cd(oldFolder)
end
