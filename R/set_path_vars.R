set_path_vars <- function(project_home_path){
  
  src_path <- file.path(project_home_path, 'src')
  results_path <- file.path(project_home_path, 'results')
  figures_path <- results_path #file.path(results_path, 'figures')
  data_path <- file.path(project_home_path, 'data')
  
  synthetic_data_results_path <- file.path(results_path, 'synthetic-data-results')
  synthetic_data_path <- file.path(data_path, 'synthetic-data')
  
  zebrafish_data_results_path <- file.path(results_path, 'zebrafish-data-results')
  zebrafish_data_path <- file.path(data_path, 'zebrafish-911-shifting-promoters')
  
  drosophila_data_results_path <- file.path(results_path, 'drosophila-data-results')
  drosophila_data_path <- file.path(data_path, 'drosophila-nplb')
  
  
  paths_list <- list(project_home_path = project_home_path, 
                     src_path = src_path, 
                     results_path = results_path, 
                     figures_path = figures_path, 
                     data_path = data_path, 
                     synthetic_data_path = synthetic_data_path, 
                     synthetic_data_results_path = synthetic_data_results_path, 
                     zebrafish_data_path = zebrafish_data_path,
                     zebrafish_data_results_path = zebrafish_data_results_path,
            		     drosophila_data_path = drosophila_data_path,
            		     drosophila_data_results_path = drosophila_data_results_path
                     )
  return(paths_list)
}
