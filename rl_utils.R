suppressMessages(library(nichenetr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(circlize))

#' Utility function for NicheNet 
#' Perform Recept-Ligand interaction predictions for clusters of sender genes & receiver genes
#'
#' @param count_matrix scRNA-seq Gene x Cell expression matrix (filtered, low-normaliztion transformed)
#' @param geneset_oi marker genes of interest from cluster of receiver genes
#' @param es expressed genes from cluster of sender genes
#' @param er expressed genes from cluster of receiver genes
#' @param from_count filter genes that exist in the whole count matrix instead of genes of intereest or expressed genes
#' @param N number of potential "ligand" from sender genes
#' @param cutoff
#' @param prior_models list of prior networks and matrices from NicheNet
#' @param fig_path figure path to save circos plots
#' @return result_list list of (1). ligand-target network (2). ligand-receptor network; (3). top upstraem ligands; (4). top receptors
#' 
#' @example 
#' results <- get_RL_pred(mat, rc_marker_gsoi, sc_expressed_genes, rc_expressed_genes, N=20)
#' res_lr_network <- results$p_ligand_receptor_network  # retrieve learnt lr-network
#' res_ligands <- results$best_upstream_ligands  # retrieve best predicted ligand
#' res_receptors <- results$top_receptors
#'
get_RL_pred <- function(count_matrix, geneset_oi, es=NULL, er=NULL, from_count=FALSE, N=NULL, cutoff=0.25, prior_models=NULL, name_l, name_r, fig_path){
	
	if (!is.null(prior_models)) {
		lr_network <- prior_models$lr_network
		weighted_networks <- prior_models$weighted_networks
		ligand_target_matrix <- prior_models$ligand_target_matrix
	} else {
		# Download prior models
		lr_network <- readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
		weighted_networks <- readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
		ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
	}
	
	if (from_count){
		#only keep genes that are found in the ligand_target_matrix
		background_expressed_genes = rownames(count_matrix) %>% .[. %in% rownames(ligand_target_matrix)]
		expressed_genes_sender = rownames(count_matrix) %>% .[. %in% colnames(ligand_target_matrix)]
		expressed_genes_receiver = rownames(count_matrix)
	}
	else{
		background_expressed_genes = er %>% .[. %in% rownames(ligand_target_matrix)]
		expressed_genes_sender = es %>% .[. %in% colnames(ligand_target_matrix)]
		expressed_genes_receiver = er
	}
	
	#choose from the expressed sender genes, the ones that can serve as ligands based on the lr network
	#choose from the expressed receiver genes, the ones tat can serve as receptors based on the lr network
	#combine to form the expressed (possible) lr network
	
	ligands = lr_network %>% pull(from) %>% unique()
	expressed_ligands = intersect(ligands,expressed_genes_sender)
	#print("-------------")
	#print(expressed_ligands)
	
	receptors = lr_network %>% pull(to) %>% unique()
	expressed_receptors = intersect(receptors,expressed_genes_receiver)
	#print("-------------")
	#print(expressed_receptors)
	
	lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)
	
	#get the potential ligands from the lr network
	potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
	
	#print("=================")
	#print(potential_ligands)
	
	#print("getting potential ligands")
	#score the potential ligands according to their power to predict the geneset_oi
	ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
	
	#print("top ligand activities")
	#print(ligand_activities %>% arrange(-pearson))
	
	#choose how many of the potential ligands to consider moving forward
	top_N = length(potential_ligands)
	if (!is.null(N)){
		top_N = min(top_N, N)
	}
	#print(top_N)
	
	best_upstream_ligands = ligand_activities %>% top_n(top_N, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
	#print(best_upstream_ligands)
	
	#show histogram of ligand activity scores
	#p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
	#	geom_histogram(color="black", fill="darkorange")  + 
	#	# geom_density(alpha=.1, fill="orange") +
	#	geom_vline(aes(xintercept=min(ligand_activities %>% top_n(top_N, pearson) %>% pull(pearson))), color="red",
	#						 linetype="dashed", size=1) + 
	#	labs(x="ligand activity (PCC)", y = "# ligands") +
	#	theme_classic()
	#p_hist_lig_activity
	
	
	# ----------------------------------------------------------
	# get the ligand-receptor network of the top-ranked ligands
	# ----------------------------------------------------------
	lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>%
		distinct(from,to)
	best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
	
	# get the weights of the ligand-receptor interactions as used in the NicheNet model
	lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
	
	# Visualize receptor-ligand networks via circos
	visualize_RL_pred(lr_network_top_df = lr_network_top_df, cell_type_r = name_r, cell_type_l = name_l, fig_path = fig_path )

	# convert to a matrix
	lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
	lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
	#print(lr_network_top_df$to)
	
	# perform hierarchical clustering to order the ligands and receptors
	dist_receptors = dist(lr_network_top_matrix, method = "binary")
	hclust_receptors = hclust(dist_receptors, method = "ward.D2")
	order_receptors = hclust_receptors$labels[hclust_receptors$order]
	
	
	dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
	hclust_ligands = hclust(dist_ligands, method = "ward.D2")
	order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
	
	#visualize ligand receptor predictions
	#print(order_receptors)
	vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
	p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>%
		make_heatmap_ggplot("Top-ranked ligands","Receptors expressed by receptor cells", color = "mediumvioletred",
                            x_axis_position = "top",legend_title = "Prior interaction potential")
	#p_ligand_receptor_network
	
	return_list = list("pred_lr_network" = vis_ligand_receptor_network,
	                   "p_ligand_receptor_network" = p_ligand_receptor_network,
                     "best_upstream_ligands" = best_upstream_ligands,
                     "top_receptors" = lr_network_top_df$to)
	
	return(return_list)
}


#' Utility function for NicheNet
#' Visualize and save circos plots for predicted RL pairs of given two clustersaa
#'
#'
visualize_RL_pred <- function(lr_network_top_df, 
                              cell_type_r, 
                              cell_type_l, 
                              fig_path, 
                              cutoff_percentile = 0.8,
                              width_same_cell_same_ligand_type = 0.8,
                              width_different_cell = 6,
                              width_ligand_receptor = 15,
                              width_same_cell_same_receptor_type = 0.8) {
  lr_network_top_df <- lr_network_top_df %>% rename(ligand = from, receptor = to)
  lr_network_top_df <- cbind(lr_network_top_df, receptor_type = cell_type_r)
  lr_network_top_df <- cbind(lr_network_top_df, ligand_type = cell_type_l)

  cutoff_weight <- lr_network_top_df$weight %>% quantile(cutoff_percentile) # cutoff weight for RL pairs to be included in the network visualization
  
  # Specify dataframe & color parameters for the circos plot
  grid_col_ligand <- setNames(object = "lawngreen", nm = cell_type_l)
  grid_col_receptor <- setNames(object = "darkred", nm = cell_type_r)
  grid_col_tbl_ligand <- tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
  grid_col_tbl_receptor <- tibble(receptor_type = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)

  circos_links <- lr_network_top_df %>% mutate(ligand = paste(ligand, " "))
  circos_links <- circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_receptor)
  links_circle <- circos_links %>% select(ligand, receptor, weight)
  
  ligand_color <- circos_links %>% distinct(ligand, color_ligand_type)
  grid_ligand_color <- ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
  receptor_color <- circos_links %>% distinct(receptor, color_receptor_type)
  grid_receptor_color <- receptor_color$color_receptor_type %>% set_names(receptor_color$receptor)
  
  grid_col <- c(grid_ligand_color, grid_receptor_color)
  transparency <- circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 
  
  # Preparation of gaps & formatting for circos plot
  ligand_order <- circos_links$ligand %>% unique()
  receptor_order <- circos_links$receptor %>% unique()
  order <- c(ligand_order, receptor_order)
  
  gaps <- c(
    rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == cell_type_l) %>% distinct(ligand) %>% nrow()-1)),
    width_ligand_receptor,
    rep(width_same_cell_same_receptor_type, times = (circos_links %>% filter(receptor_type == cell_type_r) %>% distinct(receptor) %>% nrow()-1)),
    width_ligand_receptor
  )
  
  # Circos plot & save to file
  # check fig directory exists, create directory if not
  fig_path = paste0(getwd(), fig_path)
  if (!dir.exists(fig_path)) {
    dir.create(fig_path)
  }
  
  tiff(filename = paste0(fig_path, "circos_rl_", cell_type_l, '_', cell_type_r, ".tiff"), width = 15, height = 15, units = "in", res = 300)
  circos.par(gap.degree = gaps)
  chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_weight,annotationTrack = "grid", 
               preAllocateTracks = list(track.height = 0.075))
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
  }, bg.border = NA) #
  circos.clear()
  dev.off()
}
