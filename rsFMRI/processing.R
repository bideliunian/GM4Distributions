##################
#  Read ADHD Data
#   Site: NYU
##################
library(stringr)
library(RNifti)

file_path <- "D:/Research/GM4Distributions/processed_data/NYU_preproc_filtfix/NYU"
atlas_path <- "D:/Research/GM4Distributions/processed_data/ADHD200_CC200_TCs_filtfix/templates"
saving_path <- "D:/Research/GM4Distributions/processed_data/NYU_processed"

image_atlas <- readNifti(paste(atlas_path, "ADHD200_parcellate_200.nii.gz", sep="/"))
array_atlas <- image_atlas[,,]

# image <- readNifti(paste(file_path, "0010001/sfnwmrda0010001_session_1_rest_1.nii.gz", sep="/"))

# extract ScanDir ID
phenotypic <- read.csv(paste(file_path, "NYU_phenotypic.csv", sep="/"))

phenotypic_pass <- phenotypic[phenotypic$QC_Rest_1 == "1", ]
scandir_id <- phenotypic_pass$ScanDir.ID
dx_group <- phenotypic_pass$DX # ADHD has DX!=0, control has DX==0
scandir_id <- str_pad(scandir_id, 7, pad="0")

n <- length(scandir_id)
nodes_list <- as.integer(names(table(array_atlas))[-1])
p <- length(nodes_list)
tau <- 172

# processed data

for(i in 1:n){
  processed_list <- vector("list", length = p)
  names(processed_list) <- names(table(array_atlas))[-1]
  id <- scandir_id[i]
  image_full <- readNifti(paste0(file_path, "/", id, "/sfnwmrda", id, "_session_1_rest_1.nii.gz"))
  array_full <- image_full[,,,]
  
  for (j in 1:p) {
    node <- nodes_list[j]
    indices <- which(array_atlas == node)
    array_selected <- array(dim = c(length(indices), tau))
    # array_selected <- apply(X = array_full, FUN = function(x) x[which(array_atlas == j)], MARGIN = 4)
    for (k in 1:tau) {
      array_selected[, k] <- array_full[,,,k][indices]
    }
    processed_list[j] <- list(array_selected)
    print(paste(j, "-th node, number of voxels", dim(array_selected)[1]))
  }
  
  save(processed_list, file = paste(saving_path, "/", id, ".RData", sep = ""))
  print(paste(i, "-th subject finished"))
  
  gc()
}

