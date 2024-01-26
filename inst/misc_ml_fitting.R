library(CINfits)
library(caret)
library(randomForest)
library(dplyr)
library(tidyr)
library(Biobase)

# convert depth to abs cn
depthtocn<-function(x,purity,seqdepth) #converts readdepth to copy number given purity and single copy depth
{
    (x/seqdepth-2*(1-purity))/purity
}

meta.old <- read.table("C:/Users/Phil/Downloads/OneDrive_1_26-01-2024/fit_QC_predownsample_withsmooth.tsv",header = T,sep = "\t")

qc.data <- meta.old %>%
            select(SAMPLE_ID,ploidy,purity,use) %>%
            mutate(const = SAMPLE_ID != lag(SAMPLE_ID)) %>%
            mutate(const = ifelse(is.na(const),TRUE,const)) %>%
            group_by(SAMPLE_ID) %>%
            filter(const == TRUE)

rds.rel <- readRDS("C:/Users/Phil/Downloads/OneDrive_1_26-01-2024/britroc_smoothed_copyNumbersSegmented.rds")

samples <- qc.data[which(qc.data$SAMPLE_ID %in% colnames(rds.rel)),]
rds.rel <- rds.rel[,which(colnames(rds.rel) %in% qc.data$SAMPLE_ID)]

pData(rds.rel)$purity <- samples$purity[match(pData(rds.rel)$name,samples$SAMPLE_ID)]
pData(rds.rel)$ploidy <- samples$ploidy[match(pData(rds.rel)$name,samples$SAMPLE_ID)]

abs_profiles <- rds.rel[fData(rds.rel)$use,]
# For each
for(sample in pData(rds.rel)$name){
    # Index and subselect sample
    ind <- which(colnames(rds.rel)==sample)
    relcn <- rds.rel[,ind]
    to_use <- fData(relcn)$use #
    relcn <- relcn[to_use,]
    smooth.bool <- FALSE
    # Extract cn and ploidy
    copynumber <- assayDataElement(relcn,"copynumber")
    rel_ploidy <- mean(copynumber,na.rm=T)
    ploidy <- pData(relcn)$ploidy
    purity <- pData(relcn)$purity
    cellploidy <- ploidy*purity+2*(1-purity)
    seqdepth <- rel_ploidy/cellploidy

    # Extract CN and Segs
    cn <- assayDataElement(relcn,"copynumber")
    seg <- assayDataElement(relcn,"segmented")

    # Convert to abs
    abs_cn <- depthtocn(cn,purity,seqdepth)
    abs_seg <- depthtocn(seg,purity,seqdepth)
    assayDataElement(relcn,"copynumber") <- abs_cn
    assayDataElement(relcn,"segmented") <- abs_seg
    # Add to abs RDS
    assayDataElement(abs_profiles,"copynumber")[,ind] <- abs_cn
    assayDataElement(abs_profiles,"segmented")[,ind] <- abs_seg
}

abs_segs <- CINSignatureQuantification:::getSegTable(abs_profiles)
abs_segs[,2:4] <- apply(abs_segs[,2:4],MARGIN = 2,FUN = as.numeric)

stats <- as.data.frame(calculateCINStats(abs_segs))

m.data <- stats %>%
    mutate(SAMPLE_ID = rownames(.)) %>%
    left_join(.,qc.data,by = "SAMPLE_ID") %>%
    select(-c(ploidy.y,const)) %>%
    rename("ploidy" = "ploidy.x") %>%
    tibble::column_to_rownames("SAMPLE_ID")

data <- m.data
data$use <- as.factor(data$use)

set.seed(0990)
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3))
train <- data[ind==1,]
test <- data[ind==2,]

rf <- randomForest(use ~ ., data=train, importance = T, proximity = T)
print(rf)

p1 <- predict(rf, train)
confusionMatrix(p1, train$use)

p2 <- predict(rf, test)
confusionMatrix(p2, test$use)

plot(rf)

t <- tuneRF(train[,-6], train[,6],
            stepFactor = 0.5,
            plot = F,
            ntreeTry = 500,
            trace = TRUE,
            improve = 0.05)

hist(treesize(rf),
     main = "No. of Nodes for the Trees",
     col = "green")

varImpPlot(rf,
           sort = T,
           n.var = 4,
           main = "Variable Importance")
importance(rf)

MDSplot(rf, train$use)
