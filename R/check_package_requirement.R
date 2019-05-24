check_package_requirement <- function(list_of_packages){
  
  # Install the packages in list_of_packages, if not already installed
  packages_to_install <- setdiff(list_of_packages, rownames(installed.packages()))
  if (length(packages_to_install) > 0) {
    install.packages(setdiff(list_of_packages, rownames(installed.packages())))  
  }else{
    sprintf("All required packages available!\n")
  }  
  return (NULL)
}