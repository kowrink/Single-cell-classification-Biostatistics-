library(caret)
library(tictoc)
library(dplyr)
library(tidyverse)
seq_data=read.csv("sc_celseq2.count.csv")
meta_data=read.csv("sc_celseq2.metadata.csv",stringsAsFactors = FALSE)
#singlecell=colnames(seq_data)
seq_data=data.frame(t(seq_data))
#seq_data=cbind(singlecell,seq_data)
singlecell2=(meta_data$cell_line)
#meta_data=cbind(singlecell2,rownames(meta_data))
#df=merge(meta_data,seq_data,by.x="V2",by.y="singlecell")
#set.seed(2019)
#traindata=seq_data[sample(nrow(seq_data), 180), ]
#testdata=seq_data[-traindata,]
tic()
rf_10 <- train(x = as.matrix(seq_data),
                        y = factor(singlecell2),
                        method = "ranger",
                        num.trees = 200,
                        importance = "impurity",
                        trControl = trainControl(method = "oob"))
toc()
rf_10$finalModel %>%
  # extract variable importance metrics
  ranger::importance() %>%
  # convert to a data frame
  enframe(name = "variable", value = "varimp") %>%
  top_n(n = 20, wt = varimp) %>%
  # plot the metrics
  ggplot(aes(x = fct_reorder(variable, varimp), y = varimp)) +
  geom_col() +
  coord_flip() +
  labs(x = "RNA",
       y = "Variable importance (higher is more important)")

library(Rtsne)
library(ggplot2)
library(plotly)

tsne <- Rtsne(as.matrix(seq_data), perplexity = 50, pca = FALSE)

tsne_plot <- tsne$Y %>%
  as.data.frame() %>%
  mutate(word = row.names(seq_data)) %>%
  
  ggplot(aes(x = V1, y = V2, label = word,color=singlecell2)) + 
  geom_text(size = 3)
tsne_plot

library(e1071)
svm_vector = svm(factor(singlecell2) ~ ., data = tsne$Y, scale = FALSE, kernel = "radial", cost = 5)
yprec = predict(svm_vector, tsne$Y)
datapre=data.frame(cbind(yprec,tsne$Y,singlecell2),stringsAsFactors=FALSE)
datapre$singlecell2[data$singlecell2 == "H1975"] <- "1"
datapre$singlecell2[data$singlecell2 == "H2228"] <- "2"
datapre$singlecell2[data$singlecell2 == "HCC827"] <- "3"
ggplot(datapre, aes(x=as.numeric(V2), y=as.numeric(V3), color=singlecell2) )+ geom_point()
ggplot(datapre, aes(x=as.numeric(V2), y=as.numeric(V3), color=yprec)) + geom_point()
table(datapre$yprec,datapre$singlecell2)