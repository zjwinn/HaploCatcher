#' Haplotype Prediction: Permutation Cross Validation of KNN and RF Models
#'
#'This function performs the analysis featured in Winn et al 2022 where genome wide markers are used to train machine learning models to identify if genotypes have or do not have specific alleles of a QTL/gene. This function is used to perform cross validation where a random partition of the total available data is used to train a model and a reserved testing partition is used to validate. This function is used for permutation based cross validation.
#'
#' @param n_perms A numeric variable defining the number of permutations to perform. This value may range from one to infinity. Default is 30.
#' @param geno_mat An imputed, number-coded, genotypic matrix which has n rows of individuals and m columns of markers. Row names of the matrix should be representative of genotypic IDs and column names should be representative of marker IDs. Missing data is not allowed. Numeric coding of genotypes can vary as long as it remains consistent among markers.
#' @param gene_file A dataframe containing at least three columns labeled as follows: 'Gene', 'FullSampleName', and 'Call'. The 'Gene' column contains the name of the gene for which the observation belongs to. The 'FullSampleName' column contains the genotypic ID which corresponds exactly to the column name in the genotypic matrix. The 'Call' column contains the marker call which corresponds to the gene for that genotype. Other information may be present in this dataframe beyond these columns, but the three listed columns above are obligatory.
#' @param gene_name A character string which matches the name of the gene which you are trying to perform cross validation for. This character string must be present in your gene_file 'Gene' column.
#' @param marker_info A dataframe containing the following three columns: 'Marker', 'Chromosome', and 'BP_Position'. The 'Marker' column contains the names of the marker which are present in the genotypic matrix. The 'Chromosome' column contains the corresponding chromosome/linkage group to which the marker belongs. The 'Position' column contains the physical or centimorgan position of the marker. All markers present in the genotypic matrix must be listed in this dataframe. If physical or centimorgan positions are unavailable for the listed markers, a numeric dummy variable ranging from one to n number of markers may be provided instead.
#' @param n_neighbors A numeric variable which represents the number of neighbors to use in KNN. Default is 50.
#' @param chromosome A character string which matches the name of the chromosome upon which the gene resides. This chromosome name must be present in the marker_info file.
#' @param ncor_markers A numeric variable which represents the number of markers the user want to use in model training. Correlation among markers to the gene call is calculated and the top n markers specified are retained for training. The default setting is 50 markers.
#' @param percent_testing A numeric variable which ranges such that x|0<x<1. This means that this number can be neither zero nor one. This number represents the percent of the total data available the user wants to retain to validate the model. The default setting is 0.20.
#' @param percent_training A numeric variable which ranges such that x|0<x<1. This means that the number can be neither zero nor one. This number represents the percent of the total data available the user wants to retain for training of the model.The default setting is 0.80.
#' @param include_hets A logical variable which determines if the user wishes to include heterozygous calls or not. The default setting is FALSE.
#' @param include_models A logical variable which determines if the user wishes to include the trained models in the results object for further testing. Warning: the models are quite large and running this will result in a very large results object. The default setting is FALSE.
#' @param verbose A logical variable which determines if the user wants plots displayed and text feedback from each permutation. Regardless of this parameter, the function will display the name of the gene which is being cross validated and the current progress of the permutations. Default setting is FALSE.
#' @param parallel A logical variable which determines if the user wants the cross validation performed in parallel. Default is FALSE. If the user defines that parallel is TRUE, all visual and textual feedback will not be rendered.
#' @param n_cores A numerical vector which denotes the number of cores used for parallel processor. If "parallel" option is TRUE and n_cores is not specified, then the number of available cores minus one will be assigned to processing.
#'
#' @return
#' This function returns a list of list with the following objects: "Overall_Parameters", "By_Class_Parameters", "Overall_Summary", "By_Class_Summary", and "Raw_Permutation_Info". The "Overall_Parameters" data frame contains all the relevant parameters for each permutation overall. The "By_Class_Parameters" data frame contains all the relevant parameters for each permutation by class.The "Overall_Summary" data frame contains all the relevant parameters overall summarized across permutations. The "By_Class_Summary" data frame contains all the relevant parameters by class summarized across permutations. The "Raw_Permutation_Info" is a list of list which contains each permutations model info as described in the "locus_cv" function.
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
#' #run permutational analysis - commented out for package specifications
#' #to run, copy and paste without '#' into the console
#'
#' #fit<-locus_perm_cv(n_perms = 10, #the number of permutations
#' #                   geno_mat=geno_mat, #the genotypic matrix
#' #                   gene_file=gene_comp, #the gene compendium file
#' #                   gene_name="sst1_solid_stem", #the name of the gene
#' #                   marker_info=marker_info, #the marker information file
#' #                   chromosome="3B", #name of the chromosome
#' #                   ncor_markers= 25, #number of markers to retain
#' #                   n_neighbors = 25, #number of nearest-neighbors
#' #                   percent_testing=0.2, #percentage of genotypes in the validation set
#' #                   percent_training=0.8, #percentage of genotypes in the training set
#' #                   include_hets=FALSE, #excludes hets in the model
#' #                   include_models=FALSE, #excludes models in results object
#' #                   verbose = FALSE) #excludes text
#'
#'
#' @importFrom randomForest randomForest
#' @importFrom lattice qq
#' @importFrom ggplot2 ggplot
#' @importFrom caret train
#' @importFrom foreach %dopar%

locus_perm_cv<-function(n_perms=30, #number of permutations
                        geno_mat, #the genotypic matrix
                        gene_file, #the gene compendium file
                        gene_name, #the name of the gene
                        marker_info, #the marker information file
                        chromosome, #name of the chromosome
                        ncor_markers=50, #number of markers to retain
                        n_neighbors=50, #number of nearest neighbors to use
                        percent_testing=0.2, #percentage of genotypes in the validation set
                        percent_training=0.8, #percentage of genotypes in the training set
                        include_hets=FALSE, #include hets in the model
                        include_models=FALSE, #include models,
                        verbose=FALSE, #if the function should put out text
                        parallel=FALSE, #if the function should be ran in parallel
                        n_cores=NULL #number of cores to use
                        ){

  #check inputs
  if(!is.logical(parallel)){

    stop("Option 'parallel' must be a logical argument of TRUE or FALSE!")

  }

  if(!is.logical(include_hets)){

    stop("Option 'include_hets' must be a logical argument of TRUE or FALSE!")

  }

  if(!is.logical(include_models)){

    stop("Option 'include_models' must be a logical argument of TRUE or FALSE!")

  }

  if(!is.logical(verbose)){

    stop("Option 'verbose' must be a logical argument of TRUE or FALSE!")

  }

  if(!is.numeric(n_cores)==TRUE & is.null(n_cores)==FALSE){

    stop("Option 'n_cores' must be either an numeric argument or set to NULL!")

  }

  if(parallel==TRUE){

    #set cores
    if(is.null(n_cores)){
      n_cores=parallel::detectCores()-1
    }

    #assign cluster

    cluster<-parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cluster)

    #announce
    if(verbose==TRUE){
      print("Running permutations in parallel! All text output suppressed!")
    }

    #set verbose false
    verbose=FALSE

    results<-foreach::foreach(i=1:n_perms) %dopar% {

      #run model
      a<-HaploCatcher::locus_cv(geno_mat = geno_mat,
                                gene_file = gene_file,
                                gene_name = gene_name,
                                marker_info = marker_info,
                                chromosome = chromosome,
                                ncor_markers = ncor_markers,
                                n_neighbors = n_neighbors,
                                percent_testing = percent_testing,
                                percent_training = percent_training,
                                include_hets = include_hets,
                                include_models = include_models,
                                verbose = verbose)

      #put in a list object
      return(a)

    }

    #stop the cluster
    parallel::stopCluster(cluster)

    #re-register sequential
    foreach::registerDoSEQ()

    #rename for the rest of the function
    for(i in 1:length(results)){
      names(results)[[i]]=base::paste("Permutation_", i, sep = "")
    }


  }else{

    #make an object for the results
    results<-base::list()

    #announce
    if(verbose==TRUE){
      base::print(base::paste("Conducting permutational cross validation for", gene_name))
      base::print("0% complete!")
    }

    for(i in 1:n_perms){

      #run function
      a<-HaploCatcher::locus_cv(geno_mat = geno_mat,
                                gene_file = gene_file,
                                gene_name = gene_name,
                                marker_info = marker_info,
                                chromosome = chromosome,
                                ncor_markers = ncor_markers,
                                n_neighbors = n_neighbors,
                                percent_testing = percent_testing,
                                percent_training = percent_training,
                                include_hets = include_hets,
                                include_models = include_models,
                                verbose = verbose)

      #put in a list object
      results[[base::paste("Permutation_", i, sep = "")]]<-a

      #calc percent and announce
      a<-i
      a<-base::paste(base::round((a/n_perms)*100,0), "%", sep = "")
      if(verbose==TRUE){base::print(base::paste(a, "complete!"))}

      #remove
      remove(a)

    }

    }

    #make object for results
    return_results<-base::list()

    #summarize results
    if(include_hets==TRUE){

      a<-base::c()

      #pull summaries
      for(i in base::names(results)){

        b<-base::as.data.frame(base::t(results[[i]]$confusion_matrices$knn$overall))
        b$Model="K-Nearest Neighbors"
        b$Permutation=base::gsub("Permutation_", "", i)
        c<-base::as.data.frame(base::t(results[[i]]$confusion_matrices$rf$overall))
        c$Model="Random Forest"
        c$Permutation=base::gsub("Permutation_", "", i)
        a<-rbind(a,b,c)
        remove(b,c)

      }

      base::colnames(a)<-c("Accuracy",
                           "Kappa",
                           "Accuracy_Lower_CI",
                           "Accuracy_Upper_CI",
                           "Accuracy_Null",
                           "Accuracy_P_Value",
                           "Mcnemar_P_Value",
                           "Model",
                           "Permutation")

      a<-a[,c("Permutation",
              "Model",
              "Accuracy",
              "Kappa",
              "Accuracy_Lower_CI",
              "Accuracy_Upper_CI",
              "Accuracy_Null",
              "Accuracy_P_Value",
              "Mcnemar_P_Value")]

      a$Permutation<-base::as.numeric(a$Permutation)

      return_results[["Overall_Parameters"]]<-a
      remove(a)

      a<-base::c()

      #pull summaries for by class
      for(i in base::names(results)){

        b<-base::as.data.frame(results[[i]]$confusion_matrices$knn$byClass)
        b$Model="K-Nearest Neighbors"
        b$Class=base::gsub("Class: ", "", base::rownames(b))
        b$Permutation=base::gsub("Permutation_", "", i)
        base::rownames(b)=NULL
        c<-base::as.data.frame(results[[i]]$confusion_matrices$rf$byClass)
        c$Model="Random Forest"
        c$Class=base::gsub("Class: ", "", base::rownames(c))
        c$Permutation=base::gsub("Permutation_", "", i)
        base::rownames(c)=NULL
        c$Permutation=base::gsub("Permutation_", "", i)
        a<-rbind(a,b,c)
        remove(b,c)

      }

      a<-a[,c(14,12,13,1,2,5,6,11)]
      base::colnames(a)[8]="Balanced_Accuracy"
      base::rownames(a)=NULL
      a$Permutation<-base::as.numeric(a$Permutation)

      return_results[["By_Class_Parameters"]]<-a
      remove(a)

      #summarize
      a<-stats::aggregate(. ~ Model, data=return_results$Overall_Parameters[,2:4], base::mean)
      base::colnames(a)[2:3]=base::paste("Mean_", base::colnames(a)[2:3], sep = "")
      b<-stats::aggregate(. ~ Model, data=return_results$Overall_Parameters[,2:4], base::min)
      base::colnames(b)[2:3]=base::paste("Min_", base::colnames(b)[2:3], sep = "")
      c<-stats::aggregate(. ~ Model, data=return_results$Overall_Parameters[,2:4], base::max)
      base::colnames(c)[2:3]=base::paste("Max_", base::colnames(c)[2:3], sep = "")
      d<-stats::aggregate(. ~ Model, data=return_results$Overall_Parameters[,2:4], stats::sd)
      base::colnames(d)[2:3]=base::paste("SD_", base::colnames(d)[2:3], sep = "")
      a<-base::cbind(a, b[,2:3],c[,2:3],d[,2:3])
      a<-a[,c("Model",
              "Mean_Accuracy",
              "Min_Accuracy",
              "Max_Accuracy",
              "SD_Accuracy",
              "Mean_Kappa",
              "Min_Kappa",
              "Max_Kappa",
              "SD_Kappa")]

      return_results[["Overall_Summary"]]<-a

      remove(a,b,c,d)

      #summarize
      a<-stats::aggregate(. ~ Model+Class, data=return_results$By_Class_Parameters[,2:base::ncol(return_results$By_Class_Parameters)], base::mean)
      base::colnames(a)[3:base::ncol(a)]=base::paste("Mean_", base::colnames(a)[3:base::ncol(a)], sep = "")
      b<-stats::aggregate(. ~ Model+Class, data=return_results$By_Class_Parameters[,2:base::ncol(return_results$By_Class_Parameters)], base::min)
      base::colnames(b)[3:base::ncol(b)]=base::paste("Min_", base::colnames(b)[3:base::ncol(b)], sep = "")
      c<-stats::aggregate(. ~ Model+Class, data=return_results$By_Class_Parameters[,2:base::ncol(return_results$By_Class_Parameters)], base::max)
      base::colnames(c)[3:base::ncol(c)]=base::paste("Max_", base::colnames(c)[3:base::ncol(c)], sep = "")
      d<-stats::aggregate(. ~ Model+Class, data=return_results$By_Class_Parameters[,2:base::ncol(return_results$By_Class_Parameters)], stats::sd)
      base::colnames(d)[3:base::ncol(d)]=base::paste("SD_", base::colnames(d)[3:base::ncol(d)], sep = "")
      a<-base::cbind(a, b[,3:base::ncol(b)],c[,3:base::ncol(c)],d[,3:base::ncol(d)])
      a<-a[,c("Model",
              "Class",
              "Mean_Sensitivity",
              "Min_Sensitivity",
              "Max_Sensitivity",
              "SD_Sensitivity",
              "Mean_Specificity",
              "Min_Specificity",
              "Max_Specificity",
              "SD_Specificity",
              "Mean_Precision",
              "Min_Precision",
              "Max_Precision",
              "SD_Precision",
              "Mean_Recall",
              "Min_Recall",
              "Max_Recall",
              "SD_Recall",
              "Mean_Balanced_Accuracy",
              "Min_Balanced_Accuracy",
              "Max_Balanced_Accuracy",
              "SD_Balanced_Accuracy")]

      return_results[["By_Class_Summary"]]<-a
      remove(a,b,c,d)

    }else if(include_hets==FALSE){

      a<-base::c()

      #pull summaries for overall
      for(i in base::names(results)){

        b<-base::as.data.frame(base::t(results[[i]]$confusion_matrices$knn$overall))
        b$Model="K-Nearest Neighbors"
        b$Permutation=base::gsub("Permutation_", "", i)
        c<-base::as.data.frame(base::t(results[[i]]$confusion_matrices$rf$overall))
        c$Model="Random Forest"
        c$Permutation=base::gsub("Permutation_", "", i)
        a<-rbind(a,b,c)
        remove(b,c)

      }

      base::colnames(a)<-c("Accuracy",
                           "Kappa",
                           "Accuracy_Lower_CI",
                           "Accuracy_Upper_CI",
                           "Accuracy_Null",
                           "Accuracy_P_Value",
                           "Mcnemar_P_Value",
                           "Model",
                           "Permutation")

      a<-a[,c("Permutation",
              "Model",
              "Accuracy",
              "Kappa",
              "Accuracy_Lower_CI",
              "Accuracy_Upper_CI",
              "Accuracy_Null",
              "Accuracy_P_Value",
              "Mcnemar_P_Value")]

      a$Permutation<-base::as.numeric(a$Permutation)

      return_results[["Overall_Parameters"]]<-a
      remove(a)

      a<-base::c()

      #pull summaries for by class
      for(i in base::names(results)){

        b<-base::as.data.frame(base::t(results[[i]]$confusion_matrices$knn$byClass))
        b$Model="K-Nearest Neighbors"
        b$Permutation=base::gsub("Permutation_", "", i)
        c<-base::as.data.frame(base::t(results[[i]]$confusion_matrices$rf$byClass))
        c$Model="Random Forest"
        c$Permutation=base::gsub("Permutation_", "", i)
        a<-rbind(a,b,c)
        remove(b,c)

      }

      a<-a[,c(13,12,1,2,5,6,11)]
      base::colnames(a)[7]="Balanced_Accuracy"
      base::rownames(a)=NULL
      a$Permutation<-base::as.numeric(a$Permutation)

      return_results[["By_Class_Parameters"]]<-a
      remove(a)

      #summarize
      a<-stats::aggregate(. ~ Model, data=return_results$Overall_Parameters[,2:4], base::mean)
      base::colnames(a)[2:3]=base::paste("Mean_", base::colnames(a)[2:3], sep = "")
      b<-stats::aggregate(. ~ Model, data=return_results$Overall_Parameters[,2:4], base::min)
      base::colnames(b)[2:3]=base::paste("Min_", base::colnames(b)[2:3], sep = "")
      c<-stats::aggregate(. ~ Model, data=return_results$Overall_Parameters[,2:4], base::max)
      base::colnames(c)[2:3]=base::paste("Max_", base::colnames(c)[2:3], sep = "")
      d<-stats::aggregate(. ~ Model, data=return_results$Overall_Parameters[,2:4], stats::sd)
      base::colnames(d)[2:3]=base::paste("SD_", base::colnames(d)[2:3], sep = "")
      a<-base::cbind(a, b[,2:3],c[,2:3],d[,2:3])
      a<-a[,c("Model",
              "Mean_Accuracy",
              "Min_Accuracy",
              "Max_Accuracy",
              "SD_Accuracy",
              "Mean_Kappa",
              "Min_Kappa",
              "Max_Kappa",
              "SD_Kappa")]

      return_results[["Overall_Summary"]]<-a

      remove(a,b,c,d)

      #summarize
      a<-stats::aggregate(. ~ Model, data=return_results$By_Class_Parameters[,2:base::ncol(return_results$By_Class_Parameters)], base::mean)
      base::colnames(a)[2:base::ncol(a)]=base::paste("Mean_", base::colnames(a)[2:base::ncol(a)], sep = "")
      b<-stats::aggregate(. ~ Model, data=return_results$By_Class_Parameters[,2:base::ncol(return_results$By_Class_Parameters)], base::min)
      base::colnames(b)[2:base::ncol(b)]=base::paste("Min_", base::colnames(b)[2:base::ncol(b)], sep = "")
      c<-stats::aggregate(. ~ Model, data=return_results$By_Class_Parameters[,2:base::ncol(return_results$By_Class_Parameters)], base::max)
      base::colnames(c)[2:base::ncol(c)]=base::paste("Max_", base::colnames(c)[2:base::ncol(c)], sep = "")
      d<-stats::aggregate(. ~ Model, data=return_results$By_Class_Parameters[,2:base::ncol(return_results$By_Class_Parameters)], stats::sd)
      base::colnames(d)[2:base::ncol(d)]=base::paste("SD_", base::colnames(d)[2:base::ncol(d)], sep = "")
      a<-base::cbind(a, b[,2:base::ncol(b)],c[,2:base::ncol(c)],d[,2:base::ncol(d)])
      a<-a[,c("Model",
              "Mean_Sensitivity",
              "Min_Sensitivity",
              "Max_Sensitivity",
              "SD_Sensitivity",
              "Mean_Specificity",
              "Min_Specificity",
              "Max_Specificity",
              "SD_Specificity",
              "Mean_Precision",
              "Min_Precision",
              "Max_Precision",
              "SD_Precision",
              "Mean_Recall",
              "Min_Recall",
              "Max_Recall",
              "SD_Recall",
              "Mean_Balanced_Accuracy",
              "Min_Balanced_Accuracy",
              "Max_Balanced_Accuracy",
              "SD_Balanced_Accuracy")]

      return_results[["By_Class_Summary"]]<-a
      remove(a)

    }

    #put the raw results in the object
    return_results[["Raw_Permutation_Info"]]<-results

    #return the results
    return(return_results)
}


