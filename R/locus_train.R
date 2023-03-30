#' Haplotype Prediction: Training Models for Use in Forward Prediction
#'
#' This function is used to train a model for use in forward prediction of lines which have no record
#'
#' @param geno_mat An imputed, number-coded, genotypic matrix which has n rows of individuals and m columns of markers. Row names of the matrix should be representative of genotypic IDs and column names should be representative of marker IDs. Missing data is not allowed. Numeric coding of genotypes can vary as long as it remains consistant among markers.
#' @param gene_file A dataframe containing at least three columns labeled as follows: 'Gene', 'FullSampleName', and 'Call'. The 'Gene' column contains the name of the gene for which the observation belongs to. The 'FullSampleName' column contains the genotypic ID which corresponds exactly to the column name in the genotypic matrix. The 'Call' column contains the marker call which corresponds to the gene for that genotype. Other information may be present in this dataframe beyond these columns, but the three listed columns above are obligatory.
#' @param gene_name A character string which matches the name of the gene which you are trying to perform cross validation for. This character string must be present in your gene_file 'Gene' column.
#' @param marker_info A dataframe containing the following three columns: 'Marker', 'Chromosome', and 'BP_Position'. The 'Marker' column contains the names of the marker which are present in the genotypic matrix. The 'Chromosome' column contains the corresponding chromosome/linkage group to which the marker belongs. The 'Position' column contains the physical or centimorgan position of the marker. All markers present in the genotypic matrix must be listed in this dataframe. If physical or centimorgan positions are unavailable for the listed markers, a numeric dummy variable ranging from one to n number of markers may be provided instead.
#' @param chromosome A character string which matches the name of the chromosome upon which the gene resides. This chromosome name must be present in the marker_info file.
#' @param ncor_markers A numeric variable which represents the number of markers the user want to use in model training. Correlation among markers to the gene call is calculated and the top n markers specified are retained for training. The default setting is 50 markers.
#' @param n_neighbors A numeric variable which represents the number of neighbors to use in KNN. Default is 50.
#' @param include_hets A logical variable which determines if the user wishes to include heterozygous calls or not. The default setting is FALSE.
#' @param verbose A logical variable which determines if the user wants text feedback. Default setting is TRUE.
#' @param graph A logical variable which determines if the user wants graphs displayed. default setting is FALSE.
#' @param set_seed A numeric variable that is used to set a seed for reproducible results if the user is running the function once for use in the "locus_pred" function. If the user wishes to run the function many times with a random seed and decide the outcome by voting, use the function "locus_voting" instead. The default setting is NULL.
#' @param models_request A character string which defines what models are to be ran. K-nearest neighbors is abbreviated as "knn" and random forest is "rf". If both models are desired, use the text string "all". Default setting is "all".
#'
#' @return This function returns a list of list which contains: "seed", "models_request" ,"models", and "data". The "seed" object is the seed set by the user. If no seed was provided this will appear as a character stating "no_seed_set". The "models_request" item hold the models requested. The "models" object contains the trained models. The "data" object contains the data used to train the models.
#'
#' @export
#'
#' @examples
#'
#' #set seed for reproducible sampling
#' set.seed(022294)
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
#' #Note: in practice you would have something like a gene file
#' #that does not contain any lines you are trying to predict.
#' #However, this is for illustrative purposes on how to run the function
#'
#' #sample data in the gene_comp file to make a traning population
#' train<-gene_comp[gene_comp$FullSampleName %in%
#'                    sample(gene_comp$FullSampleName,
#'                           round(length(gene_comp$FullSampleName)*0.8),0),]
#'
#' #pull vector of names, not in the train, for forward prediction
#' test<-gene_comp[!gene_comp$FullSampleName
#'                 %in% train$FullSampleName,
#'                 "FullSampleName"]
#'
#' #run the function with hets
#' fit<-locus_train(geno_mat=geno_mat, #the genotypic matrix
#'                  gene_file=train, #the gene compendium file
#'                  gene_name="sst1_solid_stem", #the name of the gene
#'                  marker_info=marker_info, #the marker information file
#'                  chromosome="3B", #name of the chromosome
#'                  ncor_markers=2, #number of markers to retain
#'                  n_neighbors=3, #number of neighbors
#'                  include_hets=FALSE, #include hets in the model
#'                  verbose = FALSE, #allows for text and graph output
#'                  set_seed = 022294, #sets a seed for reproduction of results
#'                  models = "knn") #sets what models are requested
#'
#' #predict the lines in the test population
#' pred<-locus_pred(locus_train_results=fit,
#'                  geno_mat=geno_mat,
#'                  genotypes_to_predict=test)
#'
#' #see predictions
#' head(pred)
#'
#' @importFrom randomForest randomForest
#' @importFrom lattice qq
#' @importFrom ggplot2 ggplot
#' @importFrom caret train
#' @importFrom foreach %dopar%

locus_train<-function(geno_mat, #numeric genotypic matrix
                      gene_file, #locus info file
                      gene_name, #name of the gene in the gene_file
                      marker_info, #marker information file
                      chromosome, #chromosome of gene
                      ncor_markers=50, #number of correlated markers to use
                      n_neighbors=50, #number of neighbors
                      include_hets=FALSE, #include hets in the model
                      verbose=FALSE, #allows for text outputs
                      set_seed=NULL, #sets a seed for reproduction of results
                      models_request="all", #sets what models are requested
                      graph=FALSE #allows for graphical outputs
                      ){

  #set seed
  base::set.seed(set_seed)

  #check if
  if(!gene_name %in% gene_file$Gene){

    stop("Error:'gene_name' not found in 'geno_file'!")

  }

  #make dataset
  classification<-gene_file[gene_file$Gene==gene_name,]

  #check if
  if(length(classification$FullSampleName[base::duplicated(classification$FullSampleName)])>0){

    base::stop("Error: There are duplicated individuals in the 'gene_file'!")

  }else if(base::class(base::try(base::ncol(gene_file[,c("Gene", "FullSampleName", "Call")])))=="try-error"){

    base::stop("Error: The 'gene_file' does not have the columns 'Gene', 'FullSampleName', and 'Call'!")

  }

  #check if
  if(include_hets==FALSE){

    if(verbose==TRUE){base::print("Note: Removing heterozygous calls from the dataframe")}
    classification<-classification[!classification$Call %in% classification$Call[grep("het_", classification$Call)], ]

  }else if(!is.logical(include_hets)){

    base::stop("Error: Argument 'include_hets' must be a logical argument! (TRUE or FALSE)")

  }else{

    if(verbose==TRUE){base::print("Note: User has requested heterozygous calls remain in the dataset")}

  }

  #pull genetic information for individuals
  geno_matrix<-geno_mat[base::rownames(geno_mat) %in% classification$FullSampleName,]

  #check if
  if(base::class(base::try(base::ncol(marker_info[,c("Marker", "Chromosome", "BP_Position")])))=="try-error"){

    base::stop("Error: columns are misnamed or missing in the 'marker_info'!")

  }else if(!chromosome %in% marker_info$Chromosome){

    base::stop("Error: 'chromosome' is not defined in the 'marker_info' chromosome column!")

  }

  #check if
  if(!models_request %in% c("knn","rf","all")){

    base::stop("Error: 'models_request' is incorrect, use either 'knn', 'rf', or 'all'!")

  }

  #subset markers
  selected_markers<-marker_info[marker_info$Chromosome==chromosome,]

  #pull markers
  geno_matrix<-geno_matrix[,base::colnames(geno_matrix) %in% selected_markers$Marker]

  #order the markers
  geno_matrix<-geno_matrix[base::order(base::rownames(geno_matrix)), ]

  #order the fullsamplenames in training
  training<-classification[base::order(classification$FullSampleName),]

  #cbind
  training<-base::cbind(training, geno_matrix)

  #check if
  if(base::unique(training$FullSampleName %in%  base::rownames(geno_matrix))>1){

    base::stop("Error: Individuals in the training parition are not found in the 'geno_mat'!")

  }

  #make a correlation set
  corr<-training[,c("Call", base::colnames(geno_matrix))]
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

  #format datasets
  training[,1:2]<-base::lapply(training[,1:2], base::as.factor)
  training[,3:base::ncol(training)]<-base::lapply(training[,3:base::ncol(training)], base::as.numeric)

  #show frequency in the testing population
  if(verbose==TRUE){

    a<-base::data.frame(base::table(training$Call)/ base::nrow(training))
    base::colnames(a)=c("Call", "Frequency")
    print(knitr::kable(a, caption = "Frequency of Calls in Training", digits = 2))

  }

  #make sure they are dataframes
  training<-base::as.data.frame(training)

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

  #check number of markers
  if(n_neighbors<=2){

    #make alternate grid
    grid_tune<-base::expand.grid(k =1)

  }else{

    #make grid
    grid_tune<-base::expand.grid(k = base::seq(from = 1,
                                               to = n_neighbors,
                                               by = 2))
  }


  #fit model
  if(models_request=="all" | models_request=="knn"){

    fit_1<-caret::train(Call ~ .,
                        data = training[,-1],
                        method = "knn",
                        tuneGrid = grid_tune,
                        trControl = cont_train)

  }

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
  if(models_request=="all" | models_request=="rf"){

    fit_2<-caret::train(Call ~ .,
                        data = training[,-1],
                        method = "rf",
                        tuneGrid = grid_tune,
                        trControl = cont_train)

  }

  #send message
  if(verbose==TRUE){base::print("Note: Done!")}
  if(verbose==TRUE){base::print("Note: Displaying cross validation tuning results...")}

  #summary of models
  if(verbose==TRUE){

    if(models_request=="all" | models_request=="knn"){base::print(fit_1)}
    if(models_request=="all" | models_request=="rf"){base::print(fit_2)}

  }

  #make results object
  seed<-base::ifelse(base::is.null(set_seed), "no_seed_set", set_seed)

  if(models_request=="all"){

    models<-base::list(knn=fit_1,
                       rf=fit_2)

  }else if(models_request=="knn"){

    models<-fit_1

  }else if(models_request=="rf"){

    models<-fit_2

  }

  data<-training
  results<-base::list(seed=seed,
                      models_request=models_request,
                      trained_models=models,
                      data=data)

  #retrun the results
  return(results)
}
