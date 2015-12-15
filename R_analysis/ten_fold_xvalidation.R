require(caret)
setwd("/home/rtraborn/Daphnia/CAGE/TCO/motif_analysis/Dp_homer_tag_directory")
tss.dataset <- read.table(file="tss.txt",header=FALSE)
colnames(tss.dataset) <- c("PeakID","chr","start","end","strand", "NormalizedTagCount","focusRatio","findPeaksScore","FoldChangevsLocal","p-valuevsLocal","DispersionRatio","PeriodicRatio")
flds <- createFolds(tss.dataset$PeakID, k = 10, list = TRUE, returnTrain = FALSE)

names(flds)[1] <- "fold1"
names(flds)[2] <- "fold2"
names(flds)[3] <- "fold3"
names(flds)[4] <- "fold4"
names(flds)[5] <- "fold5"
names(flds)[6] <- "fold6"
names(flds)[7] <- "fold7"
names(flds)[8] <- "fold8"
names(flds)[9] <- "fold9"
names(flds)[10] <- "fold10"

#creating the individual fold dataset

fold_1 <- tss.dataset[flds$fold1,]
fold_2 <- tss.dataset[flds$fold2,]
fold_3 <- tss.dataset[flds$fold3,]
fold_4 <- tss.dataset[flds$fold4,]
fold_5 <- tss.dataset[flds$fold5,]
fold_6 <- tss.dataset[flds$fold6,]
fold_7 <- tss.dataset[flds$fold7,]
fold_8 <- tss.dataset[flds$fold8,]
fold_9 <- tss.dataset[flds$fold9,]
fold_10 <- tss.dataset[flds$fold10,]

#performing each individual validation

#validation 1
train1 <- rbind(fold_2,fold_3,fold_4,fold_5,fold_6,fold_7,fold_8,fold_9,fold_10)
test1 <- fold_1

#validation 2
train2 <- rbind(fold_1,fold_3,fold_4,fold_5,fold_6,fold_7,fold_8,fold_9,fold_10)
test2 <- fold_2

#validation 3
train3 <- rbind(fold_1,fold_2,fold_4,fold_5,fold_6,fold_7,fold_8,fold_9,fold_10)
test3 <- fold_3

#validation 4
train4 <- rbind(fold_1,fold_2,fold_3,fold_5,fold_6,fold_7,fold_8,fold_9,fold_10)
test4 <- fold_4

#validation 5
train5 <- rbind(fold_1,fold_2,fold_3,fold_4,fold_6,fold_7,fold_8,fold_9,fold_10)
test5 <- fold_5

#validation 6
train6 <- rbind(fold_1,fold_2,fold_3,fold_4,fold_5,fold_7,fold_8,fold_9,fold_10)
test6 <- fold_6

#validation 7
train7 <- rbind(fold_1,fold_2,fold_3,fold_4,fold_5,fold_6,fold_8,fold_9,fold_10)
test7 <- fold_7

#validation 8
train8 <- rbind(fold_1,fold_2,fold_3,fold_4,fold_5,fold_6,fold_7,fold_9,fold_10)
test8 <- fold_8

#validation 9
train9 <- rbind(fold_1,fold_2,fold_3,fold_4,fold_5,fold_6,fold_7,fold_8,fold_10)
test9 <- fold_9

#validation 10
train10 <- rbind(fold_1,fold_2,fold_3,fold_4,fold_5,fold_6,fold_7,fold_8,fold_9)
test10 <- fold_10

train_list <- list(train1,train2,train3,train4,train5,train6,train7,train8,train9, train10)
test_list <- list(test1,test2,test3,test4,test5,test6,test7,test8,test9,test10)

str(train_list)
str(test_list)

for(i in seq_along(train_list)) {
    write.table(train_list[[i]], paste("train_list_",i, ".txt", sep = ""),
                col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
}


for(i in seq_along(test_list)) {
    write.table(test_list[[i]], paste("test_list_",i, ".txt", sep = ""),
                col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
}

q()
