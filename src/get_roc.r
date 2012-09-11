results_folder = "~/active_site/active_site/src/test_data/"
score_file = paste(results_folder, "scores.csv", sep = '')
class_file = paste(results_folder, "true_classes.csv", sep = '')
scores = read.table(score_file, sep=',', header=F)
classes = read.table(class_file, sep=',', header=F)
scores = as.vector(as.matrix(scores[1,]))
classes = as.vector(as.matrix(classes[1,]))


library(ROCR)

pred = prediction(scores, classes)
perf = performance(pred,"tpr","fpr")
plot(perf)
