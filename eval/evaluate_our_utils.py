#import matlab
import numpy as np
from importlib import reload  

import sys
sys.path.insert(0,'../../evaluation-2020-master')

import copy

import evaluate_12ECG_score
reload(evaluate_12ECG_score)


def evaluate_12ECG_score_mod(label_classes, labels, binary_outputs, scalar_outputs):
	label_classes = list(np.asarray(np.asarray(label_classes, dtype=np.int), dtype=np.str))
	labels = np.asarray(labels, dtype=np.bool)
	binary_outputs = np.asarray(binary_outputs, dtype=np.bool)
	scalar_outputs = np.asarray(scalar_outputs, dtype=np.float64)
	output_classes = copy.deepcopy(label_classes)

	# Define the weights, the SNOMED CT code for the normal class, and equivalent SNOMED CT codes.
	weights_file = '../../evaluation-2020-master/weights.csv'
	normal_class = '426783006'
	equivalent_classes = [['713427006', '59118001'], ['284470004', '63593006'], ['427172004', '17338001']]

	# # Find the label and output files.
	# print('Finding label and output files...')
	# label_files, output_files = evaluate_12ECG_score.find_challenge_files(label_directory, output_directory)
	#
	# # Load the labels and outputs.
	# print('Loading labels and outputs...')
	# label_classes, labels = evaluate_12ECG_score.load_labels(label_files, normal_class, equivalent_classes)
	# output_classes, binary_outputs, scalar_outputs = evaluate_12ECG_score.load_outputs(output_files, normal_class, equivalent_classes)

	label_classes, labels = prepare_labels(label_classes, labels, normal_class, equivalent_classes)
	output_classes, binary_outputs, scalar_outputs = prepare_outputs(output_classes, binary_outputs, scalar_outputs, normal_class, equivalent_classes)

	# Organize/sort the labels and outputs.
	classes, labels, binary_outputs, scalar_outputs = evaluate_12ECG_score.organize_labels_outputs(label_classes, output_classes, labels, binary_outputs, scalar_outputs)

	# Load the weights for the Challenge metric.
	weights = evaluate_12ECG_score.load_weights(weights_file, classes)

	# Only consider classes that are scored with the Challenge metric.
	indices = np.any(weights, axis=0) # Find indices of classes in weight matrix.
	classes = [x for i, x in enumerate(classes) if indices[i]]
	labels = labels[:, indices]
	scalar_outputs = scalar_outputs[:, indices]
	binary_outputs = binary_outputs[:, indices]
	weights = weights[np.ix_(indices, indices)]

	# Evaluate the model by comparing the labels and outputs.
	auroc, auprc = evaluate_12ECG_score.compute_auc(labels, scalar_outputs)
	accuracy = evaluate_12ECG_score.compute_accuracy(labels, binary_outputs)
	f_measure = evaluate_12ECG_score.compute_f_measure(labels, binary_outputs)
	f_beta_measure, g_beta_measure = evaluate_12ECG_score.compute_beta_measures(labels, binary_outputs, beta=2)
	challenge_metric = evaluate_12ECG_score.compute_challenge_metric(weights, labels, binary_outputs, classes, normal_class)

	# Return the results.
	return auroc, auprc, accuracy, f_measure, f_beta_measure, g_beta_measure, challenge_metric


def prepare_outputs(classes, binary_outputs, scalar_outputs, normal_class, equivalent_classes_collection):
	# For each set of equivalent class, use only one class as the representative class for the set and discard the other classes in the set.
	# The binary output for the representative class is positive if any of the classes in the set is positive.
	# The scalar output is the mean of the scalar outputs for the classes in the set.
	num_recordings = len(binary_outputs)
	remove_classes = list()
	remove_indices = list()
	for equivalent_classes in equivalent_classes_collection:
		equivalent_classes = [x for x in equivalent_classes if x in classes]
		if len(equivalent_classes) > 1:
			representative_class = equivalent_classes[0]
			other_classes = equivalent_classes[1:]
			equivalent_indices = [classes.index(x) for x in equivalent_classes]
			representative_index = equivalent_indices[0]
			other_indices = equivalent_indices[1:]

			binary_outputs[:, representative_index] = np.any(binary_outputs[:, equivalent_indices], axis=1)
			scalar_outputs[:, representative_index] = np.nanmean(scalar_outputs[:, equivalent_indices], axis=1)
			remove_classes += other_classes
			remove_indices += other_indices

	for x in remove_classes:
		classes.remove(x)
	binary_outputs = np.delete(binary_outputs, remove_indices, axis=1)
	scalar_outputs = np.delete(scalar_outputs, remove_indices, axis=1)

	# If any of the outputs is a NaN, then replace it with a zero.
	binary_outputs[np.isnan(binary_outputs)] = 0
	scalar_outputs[np.isnan(scalar_outputs)] = 0

	# If the binary outputs for the normal class and one or more other classes are positive, then make the binary output for the normal class negative.
	# If the binary outputs for all classes are negative, then make the binary output for the normal class positive.
	normal_index = classes.index(normal_class)
	for i in range(num_recordings):
		num_positive_classes = np.sum(binary_outputs[i, :])
		if binary_outputs[i, normal_index] == 1 and num_positive_classes > 1:
			binary_outputs[i, normal_index] = 0
		elif num_positive_classes == 0:
			binary_outputs[i, normal_index] = 1
	return classes, binary_outputs, scalar_outputs


def prepare_labels(classes, labels, normal_class, equivalent_classes_collection):
	num_recordings = len(labels)
	# For each set of equivalent class, use only one class as the representative class for the set and discard the other classes in the set.
	# The label for the representative class is positive if any of the labels in the set is positive.
	remove_classes = list()
	remove_indices = list()
	for equivalent_classes in equivalent_classes_collection:
		equivalent_classes = [x for x in equivalent_classes if x in classes]
		if len(equivalent_classes) > 1:
			representative_class = equivalent_classes[0]
			other_classes = equivalent_classes[1:]
			equivalent_indices = [classes.index(x) for x in equivalent_classes]
			representative_index = equivalent_indices[0]
			other_indices = equivalent_indices[1:]

			labels[:, representative_index] = np.any(labels[:, equivalent_indices], axis=1)
			remove_classes += other_classes
			remove_indices += other_indices

	for x in remove_classes:
		classes.remove(x)
	labels = np.delete(labels, remove_indices, axis=1)

	# If the labels for the normal class and one or more other classes are positive, then make the label for the normal class negative.
	# If the labels for all classes are negative, then make the label for the normal class positive.
	normal_index = classes.index(normal_class)
	for i in range(num_recordings):
		num_positive_classes = np.sum(labels[i, :])
		if labels[i, normal_index] == 1 and num_positive_classes > 1:
			labels[i, normal_index] = 0
		elif num_positive_classes == 0:
			labels[i, normal_index] = 1

	return classes, labels

# if __name__ == '__main__':
# 	import scipy.io
# 	mat = scipy.io.loadmat('test_data.mat')
# 	YTest_Performance = mat['YTest_Performance']
# 	classes = mat['classes'].ravel()
# 	testPredReal = mat['testPredReal']
# 	testPred_Performance = mat['testPred_Performance']
# 	auroc, auprc, accuracy, f_measure, f_beta_measure, g_beta_measure, challenge_metric = evaluate_12ECG_score_mod(classes, YTest_Performance, testPred_Performance, testPredReal)
# 	output_string = 'AUROC,AUPRC,Accuracy,F-measure,Fbeta-measure,Gbeta-measure,Challenge metric\n{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f}'.format(auroc, auprc, accuracy, f_measure, f_beta_measure, g_beta_measure, challenge_metric)
# 	print(output_string)