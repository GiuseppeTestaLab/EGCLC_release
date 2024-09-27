# some adapted from https://github.com/CompEpigen/methrix/blob/master/R/methrix_plot.R

violin_plot <- function(meth_levels, stat = "median"){
  
  df_plot <- as.data.frame(meth_levels)
  plot.data <- suppressWarnings(data.table::melt(df_plot))
  colnames(plot.data) <- c("variable", "Meth")
  
  p <- ggplot2::ggplot(plot.data, ggplot2::aes(x = variable, y = Meth,fill = variable)) + ggplot2::geom_violin(alpha = 0.8) + ggplot2::theme_classic(base_size = 14) + labs(title=stat,
    y = "CpG methylation levels (%)") + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 12,colour = "black", angle = 90, vjust = 0.5, hjust=1), 
          axis.text.y = element_text(size = 12, colour = "black"), legend.title = element_blank()) +stat_summary(fun = stat, geom = "crossbar", width = 0.5,colour = "red")
  p
}


do_PCA <- function(meth_levels){
  #remove NA
  meth <- meth_levels[complete.cases(meth_levels), , drop = FALSE]
  meth_pca <- prcomp(x = t(meth), retx = TRUE)
  
  pc_vars <- meth_pca$sdev^2/sum(meth_pca$sdev^2)
  names(pc_vars) <- colnames(meth_pca$x)
  pc_vars <- round(pc_vars, digits = 2)
  
  pca_res <- list(PC_matrix = meth_pca$x, var_explained = pc_vars)
  
  return(pca_res)
}


plot_PCA <- function(pca_res, anno = NULL, col_anno = NULL, shape_anno = NULL,
                     pc_x = "PC1", pc_y = "PC2", show_labels = FALSE, custom_colors = NULL, point_size = 4) {
  
  X <- Y <- color_me <- shape_me <- row_names <- NULL
  
  pc_vars <- pca_res$var_explained
  pca_res <- as.data.frame(pca_res$PC_matrix)
  pca_res$row_names <- rownames(pca_res)
  
  x_lab <- paste0(pc_x, " [", pc_vars[pc_x]*100, " %]")
  y_lab <- paste0(pc_y, " [", pc_vars[pc_y]*100, " %]")
  
  if (!is.null(col_anno) || !is.null(shape_anno)) {
    pd <- as.data.frame(anno)
    pd <- pd[rownames(pca_res), , drop = FALSE]
    pca_res <- cbind(pca_res, pd)
  }
  
  if (!is.null(col_anno)) {
    col_anno_idx <- which(colnames(pca_res) == col_anno)
    if (length(col_anno_idx) == 0) {
      stop(paste0(col_anno, " not found in provided object"))
    } else {
      colnames(pca_res)[col_anno_idx] <- "color_me"
    }
  }
  
  if (!is.null(shape_anno)) {
    shape_anno_idx <- which(colnames(pca_res) == shape_anno)
    if (length(shape_anno_idx) == 0) {
      stop(paste0(shape_anno, " not found in provided object"))
    } else {
      colnames(pca_res)[shape_anno_idx] <- "shape_me"
    }
  }
  
  pc_x_idx <- which(colnames(pca_res) == pc_x)
  pc_y_idx <- which(colnames(pca_res) == pc_y)
  colnames(pca_res)[c(pc_x_idx, pc_y_idx)] <- c("X", "Y")
  
  if (all(c("color_me", "shape_me") %in% colnames(pca_res))) {
    pca_gg <- ggplot(data = pca_res, aes(x = X, y = Y, color = color_me,
                                         shape = shape_me, label = row_names)) + 
      geom_point(size = point_size) +
      xlab(pc_x) + ylab(pc_y) + labs(color = col_anno, shape = shape_anno)
  } else if ("color_me" %in% colnames(pca_res)) {
    pca_gg <- ggplot(data = pca_res, aes(x = X, y = Y, color = color_me,
                                         label = row_names)) + 
      geom_point(size = point_size) + xlab(pc_x) + ylab(pc_y) +
      labs(color = col_anno)
  } else if ("shape_me" %in% colnames(pca_res)) {
    pca_gg <- ggplot(data = pca_res, aes(x = X, y = Y, shape = shape_me,
                                         label = row_names)) + 
      geom_point(size = point_size) + xlab(pc_x) + ylab(pc_y) +
      labs(shape = shape_anno)
  } else {
    pca_gg <- ggplot(data = as.data.frame(pca_res), aes(x = X, y = Y,
                                                        label = row_names)) + 
      geom_point(size = point_size, fill = "black", color = "gray70") + xlab(pc_x) + ylab(pc_y)
  }
  
  pca_gg <- pca_gg + xlab(label = x_lab) + ylab(label = y_lab) + theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(colour = "black", size = 12),
          axis.text.y = element_text(colour = "black", size = 12))
  
  if (!is.null(col_anno) && !is.null(custom_colors)) {
    pca_gg <- pca_gg + scale_color_manual(values = custom_colors)
  } else {
    pca_gg <- pca_gg + scale_color_brewer(palette = "Dark2")
  }
  
  if (show_labels) {
    pca_gg <- pca_gg + ggrepel::geom_label_repel(aes(label = row_names),
                                                 box.padding   = 0.35, 
                                                 point.padding = 0.5,
                                                 segment.color = 'grey50', show.legend = FALSE) 
  }
  
  pca_gg
}

#https://compepigen.github.io/methrix_docs/articles/04-DMRs.html

DMRbarplot <- function(DMR, comparison, y_lim = NULL, interactive = FALSE){
  count <- c(nrow(DMR[[comparison]][DMR[[comparison]]$diff.Methy<0,]),
             nrow(DMR[[comparison]][DMR[[comparison]]$diff.Methy>0,]))
  length <- c(sum(DMR[[comparison]][DMR[[comparison]]$diff.Methy<0,"length"]), 
              sum(DMR[[comparison]][DMR[[comparison]]$diff.Methy<0,"length"])) 
  names(count) <-c(paste0("Higher methylation in ", strsplit(names(DMR)[comparison], "vs")[[1]][2]), paste0("Higher methylation in ", strsplit(names(DMR)[comparison], "vs")[[1]][1])) 
  names(length) <-c(paste0("Higher methylation in ", strsplit(names(DMR)[comparison], "vs")[[1]][2]), paste0("Higher methylation in ", strsplit(names(DMR)[comparison], "vs")[[1]][1])) 
  data <- data.frame(Direction=c("Gain", "Loss"), Count=count, Length=length)
  
  g <- ggplot(data=data)+geom_col(aes(x=factor(1), y=Count, fill=Direction), position = "dodge")+theme_bw()+theme(axis.title.x = element_blank(), axis.text.x = element_blank())+scale_fill_brewer(palette = "Dark2")+ggtitle(paste0("Methylation of ", strsplit(names(DMR)[comparison], "vs")[[1]][2], " vs ", strsplit(names(DMR)[comparison], "vs")[[1]][1]))+scale_y_continuous(limits = y_lim, labels = comma)
  
  if (interactive == TRUE){
    ggplotly(g)
  }else {g}
}

DMRbarplot_length <- function(DMR, comparison, y_lim = NULL){
  count <- c(nrow(DMR[[comparison]][DMR[[comparison]]$diff.Methy<0,]),
             nrow(DMR[[comparison]][DMR[[comparison]]$diff.Methy>0,]))
  length <- c(sum(DMR[[comparison]][DMR[[comparison]]$diff.Methy<0,"length"]), 
              sum(DMR[[comparison]][DMR[[comparison]]$diff.Methy<0,"length"])) 
  names(count) <-c(paste0("Higher methylation in ", strsplit(names(DMR)[comparison], "vs")[[1]][2]), paste0("Higher methylation in ", strsplit(names(DMR)[comparison], "vs")[[1]][1])) 
  names(length) <-c(paste0("Higher methylation in ", strsplit(names(DMR)[comparison], "vs")[[1]][2]), paste0("Higher methylation in ", strsplit(names(DMR)[comparison], "vs")[[1]][1])) 
  data <- data.frame(Direction=c("Gain", "Loss"), Count=count, Length=length)
  
  p <- ggplot(data=data)+geom_col(aes(x=factor(1), y=Length, fill=Direction), position = "dodge")+theme_bw()+theme(axis.title.x = element_blank(), axis.text.x = element_blank())+scale_fill_brewer(palette = "Dark2")+ggtitle(paste0("Length DMRs ", strsplit(names(DMR)[comparison], "vs")[[1]][2], " vs ", strsplit(names(DMR)[comparison], "vs")[[1]][1]))+scale_y_continuous(limits = y_lim, labels = comma)
  
  ggplotly(p)
  
}

DMRbarplot_all <- function(DMR, comparisons_ordered, label_y, y_lim = NULL, interactive = FALSE){

  data <- data.frame("Comparison" = comparisons_ordered, "LOSS" = rep(NA, length(comparisons_ordered)), "GAIN" = rep(NA, length(comparisons_ordered)))
  
  for (comparison in comparisons_ordered) {
    data[data$Comparison==comparison ,]$LOSS <- nrow(DMR[[comparison]][DMR[[comparison]]$diff.Methy>0,]) 
    data[data$Comparison==comparison ,]$GAIN <- nrow(DMR[[comparison]][DMR[[comparison]]$diff.Methy<0,])
  } 
  
  data_long <- tidyr::pivot_longer(data, cols = c(GAIN, LOSS), names_to = "Methylation", values_to = "Count")
  data_long$Comparison <- factor(data_long$Comparison, level = rev(comparisons_ordered)) #rev to have them ordered in the plot as preferred
  data_long$Count_plot <- ifelse(data_long$Methylation == "GAIN", -1 * data_long$Count, data_long$Count)
  
  p <- ggplot(data_long, aes(x = Count_plot, y = Comparison, fill = Methylation)) +
    geom_bar(stat = "identity", position = "identity", width = 0.5) +
    #geom_text(aes(label = ifelse(Count_plot != 0, as.character(Count), "")), vjust = 1.5) +
    #scale_fill_manual(values = c("LOSS" = "#609CFF", "GAIN" = "#F8766D")) + 
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "Number of DMRs", y = "Comparison") +
    theme_minimal() +
    theme(legend.position = "top")  +
    scale_x_continuous(limits = c(-y_lim,y_lim), labels = abs) +
    scale_y_discrete(labels = label_y) 
  
  if (interactive == TRUE){
    ggplotly(p)
  }else {p}
}
