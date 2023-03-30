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
#'
#' #refer to vignette for an in depth look at the plot_locus_perm_cv function
#' vignette("An_Intro_to_HaploCatcher", package = "HaploCatcher")
#'
#' @importFrom patchwork plot_layout
#' @importFrom ggplot2 ggplot

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
    a$Class<-ifelse(a$Class %in% a$Class[grep("non_",a$Class, ignore.case = TRUE)], "-/-",
                    ifelse(a$Class %in% a$Class[grep("het_", a$Class, ignore.case = TRUE)], "+/-", "+/+"))
    a$Class<-factor(a$Class, levels = c("+/+", "+/-", "-/-"))
    d<-ggplot2::ggplot(data = a, ggplot2::aes(y=Sensitivity, x=Model, fill=Model))+
      ggplot2::facet_grid(~Class)+
      ggplot2::geom_boxplot()+
      ggplot2::theme_bw()+
      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     legend.position = "none")+
      ggplot2::labs(title = "By-Class Sensitivity")
    e<-ggplot2::ggplot(data = a, ggplot2::aes(y=Specificity, x=Model, fill=Model))+
      ggplot2::facet_grid(~Class)+
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
