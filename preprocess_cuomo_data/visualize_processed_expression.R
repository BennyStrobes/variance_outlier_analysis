args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(cowplot)
library(umap)
library(ggplot2)
library(RColorBrewer)
options(bitmapType = 'cairo', device = 'pdf')


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


######################################
# Make correlation heatmap correlating covariates with loadings
#######################################
make_covariate_loading_correlation_heatmap <- function(covariates, loadings) {
	# Covariates columns to consider
	#print(summary(covariates[, valid_covariates]))
	# Remove unimportant columns
	#covariates$experiment_day = factor(paste0(covariates$experiment, "_", covariates$day))
	#covariates$frac_alt_reads = covariates$n_alt_reads/covariates$n_total_reads
    loadings <- as.matrix(loadings)
    #valid_covariates <- 2:85
    #covs <- covariates[,valid_covariates]
    valid_covariates <- c(6, 8, 10, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 31, 32, 33, 36, 37, 38, 40, 41, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,55, 56, 57, 58, 59, 60, 69, 71, 77, 78, 79, 80, 81,83, 87, 89)
 	covariate_type <- c("num", "cat", "cat", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "cat", "num", "num", "num", "num", "num", "num", "num", "cat", "num", "num")

 	num_cov = length(valid_covariates)
 	for (iter in 1:num_cov) {
 		print(paste0(colnames(covariates)[valid_covariates[iter]], " ", covariate_type[iter]))
 	}
 	print(length(valid_covariates))
 	print(length(covariate_type))
 	print(colnames(covariates))
 	cov_names <- colnames(covariates)[valid_covariates]
 	print(cov_names)
 	num <- length(cov_names)
 	#print(cov_names)


  	covs <- covariates[,valid_covariates]
	#covs <- covariates


    # Initialize PVE heatmap
    factor_colnames <- paste0("Factor", 1:(dim(loadings)[2]))
    factor_rownames <- colnames(covs)
    pve_map <- matrix(0, dim(covs)[2], dim(loadings)[2])
    colnames(pve_map) <- factor_colnames
    rownames(pve_map) <- colnames(covs)


    # Loop through each PC, COV Pair and take correlation
    num_pcs <- dim(loadings)[2]
    num_covs <- dim(covs)[2]
    for (num_pc in 1:num_pcs) {
        for (num_cov in 1:num_covs) {
            pc_vec <- loadings[,num_pc]
            cov_vec <- covs[,num_cov]
            if (covariate_type[num_cov] == "cat") {
            #print(cov_vec[1:10])
            	lin_model <- lm(pc_vec ~ factor(cov_vec))
        	} else {
        		lin_model <- lm(pc_vec ~ cov_vec)
        	}
            pve_map[num_cov, num_pc] <- summary(lin_model)$adj.r.squared
        }
    }
    
    print("hello4")
    ord <- hclust( dist(scale(pve_map), method = "euclidean"), method = "ward.D" )$order

    melted_mat <- melt(pve_map)
    colnames(melted_mat) <- c("Covariate", "Loading","PVE")

    melted_mat$Covariate = factor(melted_mat$Covariate, levels=rownames(pve_map)[ord])
    melted_mat$Loading = factor(melted_mat$Loading, levels=factor_colnames)
	 #  Use factors to represent covariate and pc name
   	# melted_mat$Covariate 
    # melted_mat$Covariate <- factor(melted_mat$Covariate, levels = rownames(pve_map)[ord])
    #melted_mat$PC <- substr(as.character(melted_mat$PC),3,5)
    #melted_mat$PC <- factor(melted_mat$PC, levels=paste0("", 1:(length(unique(melted_mat$PC)))))

    
    #levels(melted_mat$PC) = paste0("PC", 1:(length(levels(melted_mat$PC))))
    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=Covariate, y=Loading)) + geom_tile(aes(fill=PVE)) + scale_fill_gradient2(midpoint=-.05, guide="colorbar")
    heatmap <- heatmap + labs(y="",fill="VE")
    heatmap <- heatmap + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
    # Save File
    return(heatmap)

}


#############################
# Command line args
#############################
processed_expression_dir <- args[1]
visualization_dir <- args[2]



# Input files
covariate_file <- paste0(processed_expression_dir, "cell_mata_data_day_0.txt")
covariates <- read.table(covariate_file, header=TRUE, sep="\t", comment.char="",quote="")

pca_loading_file <- paste0(processed_expression_dir, "pca_loadings_day_0.txt")
pca_loadings <- read.table(pca_loading_file, header=FALSE)




######################################
# Make Heatmap correlating known covariates with PCs
#######################################
output_file <- paste0(visualization_dir, "pca_loading_covariate_heatmap.pdf")
heatmap <- make_covariate_loading_correlation_heatmap(covariates, pca_loadings) 
ggsave(heatmap, file=output_file, width=10.2, height=5.5, units="in")