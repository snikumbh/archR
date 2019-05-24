# Date created: 22/03/2019
# Author: snikumbh 
# 

plot_separate_basis_coef_maps <- function(result_object, name_prefix, given_results_path){
  
  
  
  if(identical("NMFfit", class(result_object)[1])){
    #my_palette <- colorRampPalette(c("white", "yellow", "orange", "blue"))(n = 299)
    print("Plotting NMF result...")
    full_fname <- file.path(given_results_path, paste0(name_prefix, "_basismap.pdf"))
    pdf(full_fname)
    basismap(result_object)
    basismap(result_object)
    dev.off()
    print("basismap PDF: %s\n", full_fname)
    #
    #
    full_fname <- file.path(given_results_path, paste0(name_prefix, "_coefmap.pdf"))
    pdf(full_fname)
    coefmap(result_object, Colv = NA, col=my_palette)
    coefmap(result_object, col=my_palette)
    dev.off()
    print("coefmap PDF: %s\n", full_fname)
  }
  else if(identical("CogapsResult", class(result_object)[1])){
    my_palette <- colorRampPalette(c("white", "yellow", "orange", "blue"))(n = 299)
    print("Plotting CoGAPS result...")
    full_fname <- file.path(given_results_path, paste0(name_prefix, "_basismap.pdf"))
    pdf(full_fname)
    feature_mat <- getFeatureLoadings(result_object)
    # rownames(feature_mat) <- seq(1:dim(sample_mat)[1])
    heatmap.2(feature_mat, Colv=NA, Rowv = NA, density.info = "none", dendrogram = "none", trace = "none", col=my_palette, scale = "none")
    dev.off()
    #
    #
    full_fname <- file.path(given_results_path, paste0(name_prefix, "_coefmap.pdf"))
    pdf(full_fname)
    sample_mat <- getSampleFactors(result_object)
    rownames(sample_mat) <- seq(1:dim(sample_mat)[1])
    # colnames(sample_mat) <- rep("P")
    heatmap.2(sample_mat, Colv=NA, Rowv = NA, density.info = "none", dendrogram = "none", trace = "none", col=my_palette, scale = "none")
    dev.off()
  }else{
    
    print("do nothing")
    
  }
  # # basismap(given_data_mat.nmf, Rowv= NA, Colv = NA, scale='none', subsetRow = TRUE)
  # basismap(given_data_mat.nmf)
  # basis_mat <- basis(given_data_mat.nmf)
  # # basis_mat <- sparsify_mat(basis_mat)
  # my_palette <- colorRampPalette(c("white", "yellow", "orange", "blue"))(n = 299)
  # heatmap.2(basis_mat, Colv=NA, Rowv = NA, density.info = "none", dendrogram = "none", trace = "none", col=my_palette, scale = "none")
  # dev.off()
  # sprintf("basismap PDF: %s\n", full_fname)
  # #
  # full_fname <- file.path(given_results_path, paste0(name_prefix, "_coefmap.pdf"))
  # pdf(full_fname)
  # coefmap(given_data_mat.nmf, Colv = NA, col=my_palette)
  # coefmap(given_data_mat.nmf, col=my_palette)
  # dev.off()
  # sprintf("coefmap PDF: %s\n", full_fname)
}