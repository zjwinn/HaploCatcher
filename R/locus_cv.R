#' Haplotype Prediction: Cross Validation of KNN and RF Models
#'
#' This function performs the analysis featured in Winn et al 2022 where genome wide markers are used to train machine learning models to identify if genotypes have or do not have specific alleles of a QTL/gene. This function is used to perform cross validation where a random partition of the total available data is used to train a model and a reserved testing partition is used to validate. This function is used for a single round of cross validation.
#' @param geno_mat An imputed, number-coded, genotypic matrix which has n rows of individuals and m columns of markers. Row names of the matrix should be representative of genotypic IDs and column names should be representative of marker IDs. Missing data is not allowed. Numeric coding of genotypes can vary as long as it remains consistant among markers.
#' @param gene_file A dataframe containing at least three columns labeled as follows: 'Gene', 'FullSampleName', and 'Call'. The 'Gene' column contains the name of the gene for which the observation belongs to. The 'FullSampleName' column contains the genotypic ID which corresponds exactly to the column name in the genotypic matrix. The 'Call' column contains the marker call which corresponds to the gene for that genotype. Other information may be present in this dataframe beyond these columns, but the three listed columns above are obligatory.
#' @param gene_name A character string which matches the name of the gene which you are trying to perform cross validation for. This character string must be present in your gene_file 'Gene' column.
#' @param marker_info A dataframe containing the following three columns: 'Marker', 'Chromosome', and 'BP_Position'. The 'Marker' column contains the names of the marker which are present in the genotypic matrix. The 'Chromosome' column contains the corresponding chromosome/linkage group to which the marker belongs. The 'Position' column contains the physical or centimorgan position of the marker. All markers present in the genotypic matrix must be listed in this dataframe. If physical or centimorgan positions are unavailable for the listed markers, a numeric dummy variable ranging from one to n number of markers may be provided instead.
#' @param chromosome A character string which matches the name of the chromosome upon which the gene resides. This chromosome name must be present in the marker_info file.
#' @param ncor_markers A numeric variable which represents the number of markers the user want to use in model training. Correlation among markers to the gene call is calculated and the top n markers specified are retained for training. The default setting is 50 markers.
#' @param n_neighbors A numeric variable which represents the number of neighbors to use in KNN. Default is 50.
#' @param percent_testing A numeric variable which ranges such that x|0<x<1. This means that this number can be neither zero nor one. This number represents the percent of the total data available the user wants to retain to validate the model. The default setting is 0.20.
#' @param percent_training A numeric variable which ranges such that x|0<x<1. This means that the number can be neither zero nor one. This number represents the percent of the total data available the user wants to retain for training of the model.The default setting is 0.80.
#' @param include_hets A logical variable which determines if the user wishes to include heterozygous calls or not. Default is FALSE.
#' @param include_models A logical variable which determines if the user wishes to include the trained models in the results object for further testing. Warning: the models are quite large and running this will result in a very large results object. Default is FALSE.
#' @param verbose A logical variable which determines if the user wants text feedback. Default is TRUE.
#' @param graph A logical variable which determines if the user wants plots displayed. Default is FALSE.
#'
#' @return This function returns a list of list which contains the following list objects: 'confu', 'preds', 'models', and 'data'. The 'confu' list contains the confusion matrix objects for both the random forest and k-nearest neighbors models. The 'preds' list contains the predictions made by the separate models. The 'models' contains the two caret model objects for both the random forest and k-nearest neighbors models. The 'data' list contains the training and test data frames made by the function.
#'
#' @export
#'
#' @examples
#'
#' #read in the genotypic data matrix
#' data("geno_mat")
#'
#' #read in the marker information
#' data("marker_info")
#'
#' #read in the gene compendium file
#' data("gene_comp")
#'
#' #run the function without hets for a very limited number of markers and neighbors
#' #due to requirements by cran, this must be commented out
#' #to run, place this code in the console and remove comments
#' #fit<-locus_cv(geno_mat=geno_mat, #the genotypic matrix
#' #             gene_file=gene_comp, #the gene compendium file
#' #             gene_name="sst1_solid_stem", #the name of the gene
#' #             marker_info=marker_info, #the marker information file
#' #             chromosome="3B", #name of the chromosome
#' #             ncor_markers=2, #number of markers to retain
#' #             n_neighbors=1, #number of neighbors
#' #             percent_testing=0.2, #percentage of genotypes in the validation set
#' #             percent_training=0.8, #percentage of genotypes in the training set
#' #             include_hets=FALSE, #include hets in the model
#' #             include_models=TRUE, #include models in the final results
#' #             verbose=TRUE, #allows text output
#' #             graph=TRUE) #allows graph output
#'
#' @importFrom randomForest randomForest
#' @importFrom lattice qq
#' @importFrom ggplot2 ggplot
#' @importFrom caret train
#' @importFrom foreach %dopar%

locus_cv<-function(geno_mat, #genotypic matrix
                   gene_file, #locus info file
                   gene_name, #name of gene in gene_file
                   marker_info, #marker information file
                   chromosome, #chromosome of gene
                   ncor_markers=50, #number of markers to keep
                   n_neighbors=50, #number of nearest neighbors
                   percent_testing=0.2, #percent partition testing
                   percent_training=0.8, #percent partition training
                   include_hets=FALSE, #include heterozygous calls
                   include_models=FALSE, #include the model object form caret
                   verbose=TRUE, #produce text
                   graph=FALSE #produce graphs
                   ){
  #check if
  if(!gene_name %in% gene_file$Gene){

    stop("'gene_name' not found in 'geno_file'!")

  }

  #make dataset
  classification<-gene_file[gene_file$Gene==gene_name,]

  #check if
  if(length(classification$FullSampleName[base::duplicated(classification$FullSampleName)])>0){

    base::stop("There are duplicated individuals in the 'gene_file' or there are missing FullSampleNames in the 'gene_file'!")

  }else if(base::class(base::try(base::ncol(gene_file[,c("Gene", "FullSampleName", "Call")])))=="try-error"){

    base::stop("The 'gene_file' does not have the columns 'Gene', 'FullSampleName', and 'Call'!")

  }

  #check if
  if(include_hets==TRUE){

    if(length(unique(classification$Call))<=2){
      stop("You specified 'include_hets=TRUE' yet there appears to only be two catigorizations for your call... is your data formatted properly?")
    }else if(length(unique(classification$Call))>3){
      stop("There appears to be more than 3 catigories in your gene_file... only biallelic loci are supported at this time! Please properly format your file!")
    }

  }else if(include_hets==FALSE){

    if(length(unique(classification$Call))>3){
      stop("There appears to be more than 3 catigories in your gene_file... only biallelic loci are supported at this time! Please properly format your file!")
    }

  }

  #check if
  if(include_hets==FALSE){

    if(verbose==TRUE){base::print("Note: Removing heterozygous calls from the dataframe")}
    classification<-classification[!classification$Call %in% classification$Call[grep("het_", classification$Call)], ]

  }else if(!is.logical(include_hets)){

    base::stop("Argument 'include_hets' must be a logical argument! (TRUE or FALSE)")

  }else{

    if(verbose==TRUE){base::print("Note: User has requested heterozygous calls remain in the dataset")}

  }

  #check if
  if(percent_training>1 | percent_training<0 | !base::is.numeric(percent_training)){

    base::stop("'percent_training' is a non-integer, numeric variable bound between 0 and 1!")

  }

  #define the training
  training<-classification[base::sample(1:nrow(classification), size = percent_training*nrow(classification)),]

  #check if
  if(percent_testing>1 | percent_testing<0 | !base::is.numeric(percent_testing)){

    base::stop("'percent_testing' is a non-integer, numeric variable bound between 0 and 1!")

  }else if(percent_training+percent_testing>1){

    base::stop("'percent_testing'+'percent_training' cannot be greater than 1!")

  }

  #define the test
  test<-classification[!classification$FullSampleName %in% training$FullSampleName,]

  #check if
  if(base::unique(base::duplicated(rownames(geno_mat)))>1){

    base::stop("There are duplicated individuals in the 'geno_mat'!")

  }else if(!base::is.numeric(geno_mat)){

    base::stop("The 'geno_mat' contains non-numeric markers!")

  }


  #pull genetic information for individuals
  geno_matrix<-geno_mat[base::rownames(geno_mat) %in% classification$FullSampleName,]

  #check if
  if(base::class(base::try(base::ncol(marker_info[,c("Marker", "Chromosome", "BP_Position")])))=="try-error"){

    base::stop("columns are misnamed or missing in the 'marker_info'!")

  }else if(!chromosome %in% marker_info$Chromosome){

    base::stop("'chromosome' is not defined in the 'marker_info' chromosome column!")

  }

  #subset markers
  selected_markers<-marker_info[marker_info$Chromosome==chromosome,]

  #pull markers
  geno_matrix<-geno_matrix[,base::colnames(geno_matrix) %in% selected_markers$Marker]

  #check if
  if(base::length(base::unique(training$FullSampleName %in%  base::rownames(geno_matrix)))>1){

    base::stop("Individuals in the training parition are not found in the 'geno_mat'!")

  }

  #pull training individuals
  geno_sub_training<-geno_matrix[base::rownames(geno_matrix) %in% training$FullSampleName,]

  #order the markers
  geno_sub_training<-geno_sub_training[base::order(base::rownames(geno_sub_training)), ]

  #order the fullsamplenames in training
  training<-training[base::order(training$FullSampleName),]

  #cbind
  training<-base::cbind(training, geno_sub_training)

  #make a correlation set
  corr<-training[,c("Call", base::colnames(geno_sub_training))]
  corr$Call<-base::match(corr$Call, unique(corr$Call))
  corr[,1:base::ncol(corr)]<-base::lapply(corr[,1:base::ncol(corr)], base::as.numeric)

  #perform correlation
  corr<-base::suppressWarnings(stats::cor(corr))
  corr<-base::round(corr, 10)
  corr<-base::data.frame("|r|"=corr[-1,1],
                         check.names = FALSE)
  corr$`|r|`<-base::abs(corr$`|r|`)
  corr$marker<-base::rownames(corr)
  corr$BP_Position<-selected_markers$BP_Position
  corr$MBP_Position<-base::round(base::as.numeric(corr$BP_Position)/1000000, 2)

  #plot results
  if(graph==TRUE){

    base::plot(x=corr$MBP_Position,
               y=corr$`|r|`,
               main = base::paste("Correlational Study of markers on", chromosome, "for", gene_name),
               xlab = "Megabasepair (Mbp) Position",
               ylab = "Absolute Value of Correlation (|r|)")

  }

  #plot threshold
  corr<-corr[base::order(corr$`|r|`, decreasing = TRUE),]
  corr<-corr[1:ncor_markers,]

  #plot threshold line
  if(graph==TRUE){

    graphics::abline(h=min(corr$`|r|`), col = "red", lty=2)
    graphics::legend("bottomleft",legend=base::paste("Top", ncor_markers, "correlated markers thresold"),  col = "red", lty = 2 )

  }

  #pull markers in training
  training<-training[,
                     c("FullSampleName",
                       "Call",
                       corr$marker)]

  #subset markers and bind them to testing set
  geno_sub_testing<-geno_matrix[base::rownames(geno_matrix) %in% test$FullSampleName,
                                base::colnames(geno_matrix) %in% corr$marker]
  geno_sub_testing<-geno_sub_testing[base::order(rownames(geno_sub_testing)),]
  test<-test[base::order(test$FullSampleName),]

  #cbind
  test<-base::cbind(test,geno_sub_testing)

  #select
  test<-test[,
             c("FullSampleName",
               "Call",
               corr$marker)]

  #format datasets
  training[,1:2]<-base::lapply(training[,1:2], base::as.factor)
  test[,1:2]<-base::lapply(test[,1:2], as.factor)
  training[,3:base::ncol(training)]<-base::lapply(training[,3:base::ncol(training)], base::as.numeric)
  test[,3:base::ncol(test)]<-base::lapply(test[,3:base::ncol(test)], base::as.numeric)

  #show frequency in the testing population
  if(verbose==TRUE){

    a<-base::data.frame(base::table(training$Call)/ base::nrow(training))
    base::colnames(a)=c("Call", "Frequency")
    print(knitr::kable(a, caption = "Frequency of Calls in Training", digits = 2))
    a<-base::data.frame(base::table(test$Call)/ base::nrow(test))
    base::colnames(a)=c("Call", "Frequency")
    print(knitr::kable(a, caption = "Frequency of Calls in Test", digits = 2))

  }

  #make sure they are dataframes
  training<-base::as.data.frame(training)
  test<-base::as.data.frame(test)

  #make folds
  fold<-caret::createFolds(training[,2], k = 5)

  #make train control
  cont_train<-caret::trainControl(classProb = TRUE,
                                  verboseIter = FALSE,
                                  savePredictions = TRUE,
                                  index = fold,
                                  number = 5,
                                  repeats = 1000,
                                  method = "repeatedcv")

  #send message
  if(verbose==TRUE){print("Note: Running Models...")}

  if(n_neighbors<=2){

    #make alternate grid
    grid_tune<-base::expand.grid(k = 1)

  }else{

    #make grid
    grid_tune<-base::expand.grid(k = base::seq(from = 1,
                                               to = n_neighbors,
                                               by = 2))
  }

  #fit model
  fit_1<-caret::train(Call ~ .,
                      data = training[,-1],
                      method = "knn",
                      tuneGrid = grid_tune,
                      trControl = cont_train)

  #check number of markers
  if(ncor_markers<5){

    #make grid
    grid_tune<-base::expand.grid(mtry=1)

    }else{

    #make grid
    grid_tune<-base::expand.grid(mtry = c(1,
                                          base::seq(from = 0,
                                                    to = (base::ncol(training)-2),
                                                    by = 5)[2:length(base::seq(from = 0,
                                                                               to = (base::ncol(training)-2),
                                                                               by = 5))]))
    }

  #fit model
  fit_2<-caret::train(Call ~ .,
                      data = training[,-1],
                      method = "rf",
                      tuneGrid = grid_tune,
                      trControl = cont_train)

  #send message
  if(verbose==TRUE){base::print("Note: Done!")}
  if(verbose==TRUE){base::print("Note: Displaying cross validation tuning results...")}

  #summary of models
  if(verbose==TRUE){

    base::print(fit_1)
    base::print(fit_2)

  }

  #send message
  if(verbose==TRUE){base::print("Note: Producing predictions and calculating confusion matricies on the test population")}

  #predict
  pred_1<-stats::predict(fit_1, test[,-1])
  pred_1<-base::data.frame(FullSampleName=test[,1],
                           Model = "K-Nearest Neighbors",
                           Gene=gene_name,
                           Observed_Call=test[,2],
                           Predicted_Call=pred_1)
  pred_2<-stats::predict(fit_2, test[,-1])
  pred_2<-base::data.frame(FullSampleName=test[,1],
                           Model = "Random Forest",
                           Gene=gene_name,
                           Observed_Call=test[,2],
                           Predicted_Call=pred_2)

  #calculate confusion matraix
  confu_1<-caret::confusionMatrix(pred_1$Observed_Call,
                                  pred_1$Predicted_Call)
  confu_2<-caret::confusionMatrix(pred_2$Observed_Call,
                                  pred_2$Predicted_Call)

  #show tables
  if(verbose==TRUE){

    base::print(knitr::kable(confu_1$table, caption = "Confusion matrix of K-Nearest Neighbors predictions"))
    base::print(knitr::kable(confu_2$table, caption = "Confusion matrix of Random Forest predictions"))

  }

  #show results of models
  a<-base::data.frame(Parameters=names(confu_1$overall),
                      `K-Nearest Neighbors`=confu_1$overall,
                      `Random Forest`=confu_2$overall,
                      check.names = FALSE,
                      row.names = NULL)
  a$Parameters=c("Accuracy",
                 "Kappa",
                 "Accuracy_Lower_CI",
                 "Accuracy_Upper_CI",
                 "Accuracy_Null",
                 "Accuracy_P_Value",
                 "Mcnemar_P_Value")

  if(verbose==TRUE){

    base::print(knitr::kable(a, caption = "Overall Accuracy Parameters", digits = 3))

  }

  #show results of models
  if(include_hets==FALSE){

    a<-base::data.frame(Parameters=names(confu_1$byClass),
                        `K-Nearest Neighbors`=confu_1$byClass,
                        `Random Forest`=confu_2$byClass,
                        check.names = FALSE,
                        row.names = NULL)
    a<-a[c(1,2,5,6,11),]
    a$Parameters[5]="Balanced_Accuracy"
    base::rownames(a)=NULL

    if(verbose==TRUE){

      base::print(knitr::kable(a, caption = "By-Class Accuracy Parameters", digits = 3))

    }

  }

  if(include_hets==TRUE){

    a<-base::as.data.frame(base::rbind(confu_1$byClass, confu_2$byClass))
    a$Model=c(base::rep("K-Nearest Neighbors", 3), base::rep("Random Forest",3))
    a$Class=base::rownames(a)
    a$Class=base::gsub("Class..", "", a$Class)
    a$Class=base::gsub("[.1]", "", a$Class)
    a<-a[,c("Model",
            "Class",
            "Sensitivity",
            "Specificity",
            "Precision",
            "Recall",
            "Balanced Accuracy")]
    base::colnames(a)[7]="Balanced_Accuracy"
    rownames(a)=NULL

    if(verbose==TRUE){

      base::print(knitr::kable(a, caption = "By-Class Accuracy Parameters", digits = 3))

    }

  }

  #make results object
  confu<-base::list(knn=confu_1,
                    rf=confu_2)
  preds<-base::list(knn=pred_1,
                    rf=pred_2)
  models<-base::list(knn=fit_1,
                     rf=fit_2)
  data<-base::list(training=training,
                   test=test)

  if(include_models==TRUE){

    if(verbose==TRUE){base::print("Note: User has request that models remain in the results object")}
    results<-base::list(data_frames=data,
                        trained_models=models,
                        test_predictions=preds,
                        confusion_matrices=confu)

  }else if(include_models==FALSE){

    if(verbose==TRUE){base::print("Note: User has request that models are omitted from the results object")}
    results<-base::list(data_frames=data,
                        test_predictions=preds,
                        confusion_matrices=confu)

  }

  #retrun the results
  return(results)
}
