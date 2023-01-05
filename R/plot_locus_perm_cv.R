#' Visualize Permutation CV Results
#'
#' This function takes a list of list object from the "locus_perm_cv" function and creates a summary graphic of the accuracy, kappa, sensitivity, and specificity of the models ran. If heterozygous individuals were left in the cross validation scheme, by-class sensitivity and specificity will be displayed; otherwise, displayed parameters will be of the overall model.
#'
#' @param results An object of class "list" which is derived from the "locus_perm_cv" function.
#' @param individual_images A logical argument that defines if the user wants both the composite image and the full image. Default setting is FALSE.
#'
#' @return
#' Prints a ggplot2 image
#'
#' @export
#'
#' @examples
#'#read in the genotypic data matrix
#'data("geno_mat")
#'
#'#read in the marker information
#'data("marker_info")
#'
#'#read in the gene compendium file
#'data("genecomp")
#'
#'#run permutational analysis
#'fit<-locus_perm_cv(n_perms = 3, #the number of permutations
#'                   geno_mat=geno_mat, #the genotypic matrix
#'                   gene_file=genecomp, #the gene compendium file
#'                   gene_name="sst1_solid_stem", #the name of the gene
#'                   marker_info=marker_info, #the marker information file
#'                   chromosome="3B", #name of the chromosome
#'                   ncor_markers=50, #number of markers to retain
#'                   percent_testing=0.2, #percentage of genotypes in the validation set
#'                   percent_training=0.8, #percentage of genotypes in the training set
#'                   include_hets=TRUE, #includes hets in the model
#'                   include_models=TRUE, #includes models in results object
#'                   verbose = TRUE) #includes text/plots
#'
#'#plot results
#'plot_locus_perm_cv(fit,
#'                   individual_images=TRUE)
#'
#'@importFrom patchwork plot_layout

plot_locus_perm_cv<-function(results,
                             individual_images=FALSE){

  if(base::is.list(results)==FALSE|
     base::length(base::names(results))<5|
     base::length(base::names(results))>5){

    base::stop("Object provided to function is not a list or is not a list derived from the function 'locus_perm_cv'!")

  }

  #this line is so that the devtools::check passes
  Accuracy<-Kappa<-Sensitivity<-Specificity<-Model<-Class<-NULL

  if(base::is.null(results$By_Class_Parameters$Class)){

    a<-results$Overall_Parameters
    b<-ggplot2::ggplot(data = a, ggplot2::aes(y=Accuracy, x=Model, fill=Model))+
      ggplot2::geom_boxplot()+
      ggplot2::theme_bw()+
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())+
      ggplot2::labs(title = "Overall Accuracy")
    c<-ggplot2::ggplot(data = a, ggplot2::aes(y=Kappa, x=Model, fill=Model))+
      ggplot2::geom_boxplot()+
      ggplot2::theme_bw()+
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     legend.position = "none")+
      ggplot2::labs(title = "Overall Kappa")
    a<-results$By_Class_Parameters
    d<-ggplot2::ggplot(data = a, ggplot2::aes(y=Sensitivity, x=Model, fill=Model))+
      ggplot2::geom_boxplot()+
      ggplot2::theme_bw()+
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     legend.position = "none")+
      ggplot2::labs(title = "Overall Sensitivity")
    e<-ggplot2::ggplot(data = a, ggplot2::aes(y=Specificity, x=Model, fill=Model))+
      ggplot2::geom_boxplot()+
      ggplot2::theme_bw()+
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     legend.position = "none")+
      ggplot2::labs(title = "Overall Specificity")

    if(individual_images==TRUE){
      base::print(b)
      base::print(c)
      base::print(d)
      base::print(e)
    }

    f<-((b+c)/(d+e))+
      patchwork::plot_layout(guides = "collect")
    base::print(f)

  }else{
    a<-results$Overall_Parameters
    b<-ggplot2::ggplot(data = a, ggplot2::aes(y=Accuracy, x=Model, fill=Model))+
      ggplot2::geom_boxplot()+
      ggplot2::theme_bw()+
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())+
      ggplot2::labs(title = "Overall Accuracy")
    c<-ggplot2::ggplot(data = a, ggplot2::aes(y=Kappa, x=Model, fill=Model))+
      ggplot2::geom_boxplot()+
      ggplot2::theme_bw()+
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     legend.position = "none")+
      ggplot2::labs(title = "Overall Kappa")
    a<-results$By_Class_Parameters
    d<-ggplot2::ggplot(data = a, ggplot2::aes(y=Sensitivity, x=Model, fill=Model))+
      ggplot2::facet_grid(~base::paste("Class =", Class))+
      ggplot2::geom_boxplot()+
      ggplot2::theme_bw()+
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     legend.position = "none")+
      ggplot2::labs(title = "By-Class Sensitivity")
    e<-ggplot2::ggplot(data = a, ggplot2::aes(y=Specificity, x=Model, fill=Model))+
      ggplot2::facet_grid(~base::paste("Class =", Class))+
      ggplot2::geom_boxplot()+
      ggplot2::theme_bw()+
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     legend.position = "none")+
      ggplot2::labs(title = "By-Class Specificity")

    if(individual_images==TRUE){
      base::print(b)
      base::print(c)
      base::print(d)
      base::print(e)
    }

    f<-((b+c)/(d+e))+
      patchwork::plot_layout(guides = "collect")
    base::print(f)
  }

   remove(a,b,c,d,e,f)

}
