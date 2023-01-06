#' Haplotype Prediction: Using Trained Models to Make Predictions
#'
#' @param locus_train_results This object is a the results object from the "locus_train" function.
#' @param geno_mat This is a genotypic matrix with genotypes of individuals you have genotyped and characterized for the locus/QTL/gene of interest and individuals that have only been genotyped with a genome wide marker panel. It is important to note that the markers in the genome wide panel *must* be shared between training population and the population you wish to forward predict. In the case of genotyping-by-sequencing markers, the training and test populations should be discovered and produced together into a genotype file. All marker platforms, however, are compatable, as long as the training and forward prediction population share the same markers genome-wide.
#' @param genotypes_to_predict This is a character vector of genotypes in the geno_mat which the user wishes to predict. If this object contains names in the training population, they will be omitted in the results to avoid bias.
#'
#' @return a data frame with two columns: genotype names and predictions.
#'
#' @export
#'
#' @examples
#'
#'
#'
geno_pred<-function(geno_train_results,
                    geno_mat,
                    genotypes_to_predict){

  a<-geno_train_results
  b<-geno_mat
  c<-genotypes_to_predict

  if(a$models_request=="all"){

    #pull models
    knn<-a$trained_models$knn
    rf<-a$trained_models$rf

    #subset to only have the indiviudals who you want to predict
    b<-b[base::rownames(b) %in% c, base::colnames(b) %in% a$selected_markers$selected_markers]

    #predict
    pred_1=stats::predict(knn, b)
    pred_2=stats::predict(rf, b)

    #make a prediction object
    pred<-base::data.frame(FullSampleName=rownames(b),
                           Prediction_KNN=pred_1,
                           Prediction_RF=pred_2)

    #make results
    results<-pred

    #return
    return(results)

  }else if(a$models_request=="knn"){

    #pull models
    knn<-a$trained_models$knn

    #subset to only have the indiviudals who you want to predict
    b<-b[base::rownames(b) %in% c, base::colnames(b) %in% a$selected_markers$selected_markers]

    #predict
    pred_1=stats::predict(knn, b)

    #make a prediction object
    pred<-base::data.frame(FullSampleName=rownames(b),
                           Prediction_KNN=pred_1)

    #make results
    results<-pred

    #return
    return(results)

  }else if(a$models_request=="rf"){

    #pull models
    rf<-a$trained_models$rf

    #subset to only have the indiviudals who you want to predict
    b<-b[base::rownames(b) %in% c, base::colnames(b) %in% a$selected_markers$selected_markers]

    #predict
    pred_2=stats::predict(rf, b)

    #make a prediction object
    pred<-base::data.frame(FullSampleName=rownames(b),
                           Prediction_RF=pred_2)

    #make results
    results<-pred

    #return
    return(results)

  }else{

    base::stop("Error: Is the 'geno_train_results' correct? Check your inputs!")

  }

}
