library(e1071)
library(caTools)
library(tictoc)

main_path <- "/data/amohammed/nci"

# Import drug activity data set & separate meta data from data
raw_drug_activity <- read.csv(file="raw_drug_activity.csv",row.names=1)
drug_meta <- raw_drug_activity[,c(1:9)]
drug_activity <- raw_drug_activity[,-c(1:9)]

# Import module eigengene data
raw_MEs <- read.csv(file="MEs.csv",row.names=1)


# ML.func(drug_activity= Drug activity dataset,
#         dim= Dimension of ME combinations,
#         act= Threshold for active classification ,
#         percent_training= Percentage of data used for training,
#         ncycles= # of times model is trained and tested before taking average)

ML.func <- function (drug_activity,dim,act,percent_training,ncycles) {
  
  # Store drug names in a vector called drugz
  drugz <- row.names(drug_activity)
  
  # Create a list of data frames that contain combinations of module eigengenes
  # based on the number of dimensions specified. Must be at least 2.
  EG_combinations <- combn(raw_MEs,dim,simplify = FALSE)
  
  # Store the names of each module eigengene combination in a character vector 
  # called ME_comb_names
  ME_comb_names <- {}
  for (i in 1:length(EG_combinations)) {
    cnames <- colnames(EG_combinations[[i]])
    ME_comb_names[i] <- capture.output(cat(cnames, sep= "_"))
  }
  
  # Create an empty data frame that will be used to store the results of the 
  # drug being analyzed
  results <- data.frame(matrix(ncol = 5, nrow = length(EG_combinations)))
  colnames(results) <- c("Accuracy", "Sensitivity", "Specificity", "PPV", "NPV")
  
  # Create an empty list that will be used to store each list of the results data 
  # frames for each drug. This will be a nested list.
  results_list <- {}
  
  for (drug in drugz) {
    
    # Get the activity values for the drug being analyzed & store it in an 
    # array called raw_act
    raw_act <- t(drug_activity[drug,])
    
    # Remove rows that contain NA values 
    id <- which(!is.na(raw_act[,1]))
    activity <<- raw_act[id,1]
    MEs <- raw_MEs[id,]
    EG_combz <- combn(MEs,dim,simplify = FALSE)
    
    # Classify each activity value as active or inactive based on the value of 
    # activity threshold "act" & store the result in a vector called indicator. 
    # 1 for active, 0 for inactive. 
    indicator <- rep(0,length(activity))
    id <- which(activity>act)
    indicator[id] <- 1
    indicator <- as.factor(indicator)
    
    # Train and test SVMs for each module eigengene combination
    for (i in 1:length(EG_combz)) {
      dat <- data.frame(EG_combz[[i]], y = indicator)
      results[i,] <- SVM_(dat,dim,percent_training,ncycles)
    }
    
    # Use ME_comb_names to name the rows of the results data frame 
    row.names(results) <- ME_comb_names
    
    # Sort the results from greatest to least based on accuracy
    results_sorted <- results[order(results$Accuracy, decreasing = TRUE),,drop=FALSE]
    
    # Create a copy of the sorted results data frame that includes the names of 
    # each module eigengene combination as a column. That copy will be exported 
    # as a sheet in an excel spreadsheet. This is done to prevent the loss of 
    # row names.
    ME_comb_names_sorted <- data.frame(MEs = row.names(results_sorted))
    results_sorted_xlsx <- cbind(ME_comb_names_sorted,results_sorted)
    
    # Extract the row of the ME combination that had the best accuracy, along 
    # with the meta data for the drug it is associated with
    best_combination_xlsx <- results_sorted_xlsx[1,]
    drug_name <- drug_meta[drug,"Drug.NAme"]
    Mechanism <- drug_meta[drug,"Mechanism.of.action"]
    best_combination_xlsx <- cbind(data.frame(NSC= drug),
                                   data.frame(Drug_Name= drug_name),
                                   as.data.frame(Mechanism), 
                                   best_combination_xlsx)
    
    # Make a list of the data frames for each drug & and store it in the nested
    # list results_list
    df_list <- list(results,results_sorted,results_sorted_xlsx,best_combination_xlsx)
    names(df_list) <- c("results","results_sorted","results_sorted_xlsx","best_combination_xlsx")
    results_list[[drug]] <- df_list
  }
  
  return(results_list)
}






# Function for training and testing SVMs 
SVM_ <- function (dat,dim,percent_training,ncycles) {
  
  # Empty objects for collecting average stats
  accuracy_all <- {}
  sensitivity_all <- {}
  specificity_all <- {}
  PPV_all <- {}
  NPV_all <- {}
  
  # In order to avoid the "Model is empty!" error, a while loop is used to continue
  # shuffling the training set until there are at least 2 unique values in the 
  # indicator column. If the values of that column are all 0s or all 1s, it will
  # trigger the "Model is empty!" error.
  training <- data.frame(y=c(0,0))
  
  # SVMs are trained and tested the number of times specified in ncycles for 
  # the inputted module eigengene combination. Once the value of ncycles is 
  # reached, the average of each stat is calculated.
  for(k in 1:ncycles) {
    
    while (length(unique(training$y)) < 2) {
      # The sample.int() function from the 'caTools' package is used to
      # randomly subset a percentage of data to be used as the training set.
      # This shuffle will continue until there are at least 2 unique values 
      # in the indicator column y of the training set.
      id <- sample.int(nrow(dat), floor(length(activity) * percent_training/100))
      training <- dat[id,]
    }
    
    # The data that was not included in the training set is placed in the 
    # test set
    test <- dat[-id,]
    
    # Here is where the SVM model is trained
    classifier <- svm(formula = y ~ ., data = training, 
                      type = 'C-classification', kernel = 'linear')
    
    # Here is where the model is tested
    ncols <- dim + 1
    y_pred <- predict(classifier, newdata = test[,-ncols])
    
    # Creation of a confusion matrix based on the test results
    cm <- table(test[, ncols], y_pred)
    
    # The stats from each iteration is stored
    accuracy_all[k] <- (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
    sensitivity_all[k] <- cm[2,2]/(cm[2,2] + cm[2,1])
    specificity_all[k] <- cm[1,1]/(cm[1,1] + cm[1,2])
    PPV_all[k] <- cm[2,2]/(cm[2,2] + cm[1,2])
    NPV_all[k] <- cm[1,1]/(cm[1,1] + cm[2,1])
  }
  
  # NA values that resulted from dividing by zero are converted to 0
  sensitivity_all[is.na(sensitivity_all)] <- 0
  specificity_all[is.na(specificity_all)] <- 0
  PPV_all[is.na(PPV_all)] <- 0
  NPV_all[is.na(NPV_all)] <- 0
  
  # The average value of each stat is calculated and placed in a data frame 
  # called results
  acc <- mean(accuracy_all)
  sen <- mean(sensitivity_all)
  spe <- mean(specificity_all)
  PPV <- mean(PPV_all)
  NPV <- mean(NPV_all)
  results <- data.frame(Accuracy= acc, Sensitivity= sen, Specificity= spe, PPV= PPV, NPV= NPV)
  
  return(results)
}




# tic()
# test <- ML.func(drug_activity[1,],2,6,80,100)
# toc()

best.comb.xlsx <- function (results_df,drugnames) {
  out <- data.frame(matrix(ncol = ncol(results_df[[drugnames[1]]][["best_combination_xlsx"]]), 
                           nrow = nrow(drug_activity)))
  colnames(out) <- colnames(results_df[[drugnames[1]]][["best_combination_xlsx"]])
  for (i in 1:length(drugnames)) {
    out[i,] <- results_df[[i]][["best_combination_xlsx"]]
  }
  return(out)
}



drugnames <- row.names(raw_drug_activity)

################################################################################

SVM_results_2D_6_80_100 <- ML.func(drug_activity,2,6,80,100)
saveRDS(SVM_results_2D_6_80_100, file.path(main_path, "SVM_results_2D_6_80_100.rds"))

SVM_results_2D_6.5_80_100 <- ML.func(drug_activity,2,6.5,80,100)
saveRDS(SVM_results_2D_6.5_80_100, file.path(main_path, "SVM_results_2D_6.5_80_100.rds"))

SVM_results_2D_7_80_100 <- ML.func(drug_activity,2,7,80,100)
saveRDS(SVM_results_2D_7_80_100, file.path(main_path, "SVM_results_2D_7_80_100.rds"))

SVM_results_2D_7.5_80_100 <- ML.func(drug_activity,2,7.5,80,100)
saveRDS(SVM_results_2D_7.5_80_100, file.path(main_path, "SVM_results_2D_7.5_80_100.rds"))

SVM_top_2D_6_80_100 <- best.comb.xlsx(SVM_results_2D_6_80_100,drugnames)
SVM_top_2D_6.5_80_100 <- best.comb.xlsx(SVM_results_2D_6.5_80_100,drugnames)
SVM_top_2D_7_80_100 <- best.comb.xlsx(SVM_results_2D_7_80_100,drugnames)
SVM_top_2D_7.5_80_100 <- best.comb.xlsx(SVM_results_2D_7.5_80_100,drugnames)

dfs <- list(SVM_top_2D_6_80_100,SVM_top_2D_6.5_80_100,SVM_top_2D_7_80_100,SVM_top_2D_7.5_80_100)
names(dfs) <- c("act= 6", "act= 6.5", "act= 7", "act= 7.5")
writexl::write_xlsx(dfs,path=file.path(main_path, "SVM_2D_100_cycles.xlsx"))

################################################################################

SVM_results_2D_6_80_500 <- ML.func(drug_activity,2,6,80,500)
saveRDS(SVM_results_2D_6_80_500, file.path(main_path, "SVM_results_2D_6_80_500.rds"))

SVM_results_2D_6.5_80_500 <- ML.func(drug_activity,2,6.5,80,500)
saveRDS(SVM_results_2D_6.5_80_500, file.path(main_path, "SVM_results_2D_6.5_80_500.rds"))

SVM_results_2D_7_80_500 <- ML.func(drug_activity,2,7,80,500)
saveRDS(SVM_results_2D_7_80_500, file.path(main_path, "SVM_results_2D_7_80_500.rds"))

SVM_results_2D_7.5_80_500 <- ML.func(drug_activity,2,7.5,80,500)
saveRDS(SVM_results_2D_7.5_80_500, file.path(main_path, "SVM_results_2D_7.5_80_500.rds"))

SVM_top_2D_6_80_500 <- best.comb.xlsx(SVM_results_2D_6_80_500,drugnames)
SVM_top_2D_6.5_80_500 <- best.comb.xlsx(SVM_results_2D_6.5_80_500,drugnames)
SVM_top_2D_7_80_500 <- best.comb.xlsx(SVM_results_2D_7_80_500,drugnames)
SVM_top_2D_7.5_80_500 <- best.comb.xlsx(SVM_results_2D_7.5_80_500,drugnames)

dfs <- list(SVM_top_2D_6_80_500,SVM_top_2D_6.5_80_500,SVM_top_2D_7_80_500,SVM_top_2D_7.5_80_500)
names(dfs) <- c("act= 6", "act= 6.5", "act= 7", "act= 7.5")
writexl::write_xlsx(dfs,path=file.path(main_path, "SVM_2D_500_cycles.xlsx"))

################################################################################

SVM_results_2D_6_80_1000 <- ML.func(drug_activity,2,6,80,1000)
saveRDS(SVM_results_2D_6_80_1000, file.path(main_path, "SVM_results_2D_6_80_1000.rds"))

SVM_results_2D_6.5_80_1000 <- ML.func(drug_activity,2,6.5,80,1000)
saveRDS(SVM_results_2D_6.5_80_1000, file.path(main_path, "SVM_results_2D_6.5_80_1000.rds"))

SVM_results_2D_7_80_1000 <- ML.func(drug_activity,2,7,80,1000)
saveRDS(SVM_results_2D_7_80_1000, file.path(main_path, "SVM_results_2D_7_80_1000.rds"))

SVM_results_2D_7.5_80_1000 <- ML.func(drug_activity,2,7.5,80,1000)
saveRDS(SVM_results_2D_7.5_80_1000, file.path(main_path, "SVM_results_2D_7.5_80_1000.rds"))


SVM_top_2D_6_80_1000 <- best.comb.xlsx(SVM_results_2D_6_80_1000,drugnames)
SVM_top_2D_6.5_80_1000 <- best.comb.xlsx(SVM_results_2D_6.5_80_1000,drugnames)
SVM_top_2D_7_80_1000 <- best.comb.xlsx(SVM_results_2D_7_80_1000,drugnames)
SVM_top_2D_7.5_80_1000 <- best.comb.xlsx(SVM_results_2D_7.5_80_1000,drugnames)

dfs <- list(SVM_top_2D_6_80_1000,SVM_top_2D_6.5_80_1000,SVM_top_2D_7_80_1000,SVM_top_2D_7.5_80_1000)
names(dfs) <- c("act= 6", "act= 6.5", "act= 7", "act= 7.5")
writexl::write_xlsx(dfs,path=file.path(main_path, "SVM_2D_1000_cycles.xlsx"))

################################################################################

SVM_results_2D_6_80_10000 <- ML.func(drug_activity,2,6,80,10000)
saveRDS(SVM_results_2D_6_80_10000, file.path(main_path, "SVM_results_2D_6_80_10000.rds"))

SVM_results_2D_6.5_80_10000 <- ML.func(drug_activity,2,6.5,80,10000)
saveRDS(SVM_results_2D_6.5_80_10000, file.path(main_path, "SVM_results_2D_6.5_80_10000.rds"))

SVM_results_2D_7_80_10000 <- ML.func(drug_activity,2,7,80,10000)
saveRDS(SVM_results_2D_7_80_10000, file.path(main_path, "SVM_results_2D_7_80_10000.rds"))

SVM_results_2D_7.5_80_10000 <- ML.func(drug_activity,2,7.5,80,10000)
saveRDS(SVM_results_2D_7.5_80_10000, file.path(main_path, "SVM_results_2D_7.5_80_10000.rds"))

SVM_top_2D_6_80_10000 <- best.comb.xlsx(SVM_results_2D_6_80_10000,drugnames)
SVM_top_2D_6.5_80_10000 <- best.comb.xlsx(SVM_results_2D_6.5_80_10000,drugnames)
SVM_top_2D_7_80_10000 <- best.comb.xlsx(SVM_results_2D_7_80_10000,drugnames)
SVM_top_2D_7.5_80_10000 <- best.comb.xlsx(SVM_results_2D_7.5_80_10000,drugnames)

dfs <- list(SVM_top_2D_6_80_10000,SVM_top_2D_6.5_80_10000,SVM_top_2D_7_80_10000,SVM_top_2D_7.5_80_10000)
names(dfs) <- c("act= 6", "act= 6.5", "act= 7", "act= 7.5")
writexl::write_xlsx(dfs,path=file.path(main_path, "SVM_2D_10000_cycles.xlsx"))

