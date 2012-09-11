sink("/dev/null")
args<-commandArgs(TRUE)

suppressMessages(library(ROCR))

input_file = args[1]
input = read.table(input_file, header = FALSE, sep = ',', quote = '')

output_file = args[2]

iteration = args[3]
iter_string = as.character(iteration)


obj_val = as.numeric(args[4])
obj_str = format(obj_val, width = 3, digits = 3)

print(input_file)
print(output_file)

classes = input[,1]
scores = input[,2]

pred = prediction(scores, classes)
perf = performance(pred,"tpr","fpr")


jpeg(file = output_file, width = 480, height = 960)

attach(mtcars)
layout(matrix(c(1,2), 2, 1, byrow = T), widths = c(1), heights = c(1,1))


roc = performance(pred,"auc")
roc = roc@y.values
roc = roc[[1]]

roc_title_string = paste("iter:", iter_string, "obj:", obj_str, "roc", format(roc, digits=4), Sys.time(), sep = " ")

plot(perf, main = roc_title_string, ylim=c(0,1))

perf = performance(pred,"prec","rec")
xx = perf@x.values
yy = perf@y.values
xx = unlist(xx)
yy = unlist(yy)
yy[1] = 1
#f = approxfun(xx,yy)
#area_prec_rec = integrate(f,0,1)$value
area_prec_rec = 0

precrec_title_string = paste("iter", iter_string, "obj:", obj_str, "prec_rec", format(area_prec_rec, digits=4), Sys.time(), sep = " ")

plot(perf, main = precrec_title_string, ylim=c(0,1))


alt_input_file = args[4]
data_alt = 



dev.off()
sink()