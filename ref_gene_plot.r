library(ggplot2)

load("human_refuse_20190711.RData")

dim(ref_exp)
ref_exp[1:5,1:2]
data <- ref_exp[rownames(ref_exp)=='ACE2',]
test <- sort(data, decreasing=T)[1:50]
test <- as.data.frame(test)
colnames(test)<-'Data'
test$Sample <- rownames(test)
test$Type <- ifelse(stringi::stri_detect(rownames(test), fixed='Fetal'), 'Adult', 'Fetal')
test$Tissue <- ifelse(stringi::stri_detect(rownames(test), fixed='Kidney'), 'Kidney', 'Others')
test_adult <- test[test$Type=='Adult',]
test_fetal <- test[test$Type=='Fetal',]


ggplot(test_fetal, aes(x=reorder(Sample, Data), y=Data, fill=Tissue)) + 
  geom_bar(stat='identity', width=0.5) +
  labs(title = 'Fetal', y='Expression Level', x='Sample') +
  theme(axis.text.x=element_text(angle=90, size=7, vjust=0.5),
        axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=18))

ggplot(test_adult, aes(x=reorder(Sample, Data), y=Data, fill=Tissue)) + 
  geom_bar(stat='identity', width=0.5) +
  labs(title = 'Adult', y='Expression Level', x='Sample') +
  theme(axis.text.x=element_text(angle=90, size=7, vjust=0.5),
        axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=18))
