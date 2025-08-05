#' workhorse function for \code{brainconn()}
#'
#' returns a ggraph object of plotted brain connectivity matrix
#' @author Sidhant Chopra
#' @import ggraph
#' @import ggplot2
#' @import grid


brainconn2=function (atlas, background = "ICBM152", view = "top", conmat = NULL, 
                     node.size = 4, node.color = "network", all.nodes = FALSE, 
                     edge.color = "black", edge.alpha = 0.8, edge.width = 1, 
                     edge.color.weighted = FALSE, labels = FALSE, show.legend = TRUE, 
                     thr = NULL, uthr = NULL, scale.edge.width = NULL, label.size = 1.5, 
                     label.edge.weight = FALSE, background.alpha = 1, bg_xmax = 0, 
                     bg_xmin = 0, bg_ymax = 0, bg_ymin = 0) 
{
  ifelse(is.character(atlas), data <- get(atlas), data <- atlas)
  if (background != "ICBM152" && view == "ortho") {
    stop("Custom background image detected, view cannot be 'ortho', please select top,\n        bottom, left, right, front or back.")
  }
  if (!is.null(thr)) {
    conmat[conmat < thr] <- 0
  }
  if (!is.null(uthr)) {
    conmat[conmat > thr] <- 0
  }
  if (view == "ortho") {
    ortho_list <- list()
    ortho_views <- c("top", "left", "front")
    for (v in 1:3) {
      view <- ortho_views[v]
      bg <- paste0("ICBM152_", view)
      m <- get(bg)
      w <- matrix(rgb(m[, , 1], m[, , 2], m[, , 3], m[, 
                                                      , 4] * background.alpha), nrow = dim(m)[1])
      background <- rasterGrob(w)
      nparc <- dim(data)[1]
      if (!exists("conmat")) {
        conmat <- matrix(0L, nrow = nparc, ncol = nparc)
      }
      conmat <- as.matrix(conmat)
      rownames(conmat) <- colnames(conmat)
      ifelse(isSymmetric.matrix(conmat) == TRUE, directed <- FALSE, 
             directed <- TRUE)
      if (all.nodes == FALSE && directed == FALSE) {
        include.vec <- vector(length = dim(data)[1])
        for (i in 1:dim(conmat)[1]) {
          ifelse(any(conmat[i, ] != 0), include.vec[i] <- 1, 
                 include.vec[i] <- 0)
        }
        data <- data[as.logical(include.vec), , drop = F]
        conmat <- conmat[which(rowSums(conmat, na.rm = T) != 
                                 0), which(colSums(conmat, na.rm = T) != 0), 
                         drop = F]
      }
      if (all.nodes == FALSE && directed == TRUE) {
        include.vec <- vector(length = dim(data)[1])
        for (i in 1:dim(conmat)[1]) {
          ifelse(any(conmat[i, ] != 0) | any(conmat[, 
                                                    i] != 0), include.vec[i] <- 1, include.vec[i] <- 0)
        }
      }
      if (all.nodes == TRUE) {
        include.vec <- vector(length = dim(data)[1])
        include.vec <- rep(1, length = dim(data)[1])
      }
      ifelse(v == 1, show.legend <- T, show.legend <- F)
      ortho_list[[v]] <- build_plot2(conmat = conmat, data = data, 
                                     background = background, node.size = node.size, 
                                     view = view, node.color = node.color, thr = thr, 
                                     uthr = uthr, edge.color = edge.color, edge.alpha = edge.alpha, 
                                     edge.width = edge.width, scale.edge.width = scale.edge.width, 
                                     show.legend = show.legend, labels = labels, 
                                     label.size = label.size, include.vec = include.vec, 
                                     edge.color.weighted = edge.color.weighted, label.edge.weight = label.edge.weight)
      if (is.environment(edge.color) == T) {
        ortho_list[[v]] <- ortho_list[[v]] + edge.color
      }
    }
    right_col <- plot_grid(ortho_list[[2]], ortho_list[[3]], 
                           nrow = 2, rel_heights = c(1, 1.45))
    p <- plot_grid(ortho_list[[1]], right_col, rel_widths = c(1.8, 
                                                              1.2))
    return(p)
  }
  if (background == "ICBM152") {
    bg <- paste0("ICBM152_", view)
    m <- get(bg)
    w <- matrix(rgb(m[, , 1], m[, , 2], m[, , 3], m[, , 
                                                    4] * background.alpha), nrow = dim(m)[1])
  }
  if (background != "ICBM152") {
    m <- OpenImageR::readImage(background)
  }
  w <- matrix(rgb(m[, , 1], m[, , 2], m[, , 3], m[, , 4] * 
                    background.alpha), nrow = dim(m)[1])
  background <- rasterGrob(w)
  nparc <- dim(data)[1]
  if (!exists("conmat")) {
    conmat <- matrix(0L, nrow = nparc, ncol = nparc)
  }
  conmat <- as.matrix(conmat)
  rownames(conmat) <- colnames(conmat)
  ifelse(isSymmetric.matrix(conmat) == TRUE, directed <- FALSE, 
         directed <- TRUE)
  if (all.nodes == FALSE && directed == FALSE) {
    include.vec <- vector(length = dim(data)[1])
    for (i in 1:dim(conmat)[1]) {
      ifelse(any(conmat[i, ] != 0), include.vec[i] <- 1, 
             include.vec[i] <- 0)
    }
    data <- data[as.logical(include.vec), , drop = F]
    conmat <- conmat[which(rowSums(conmat, na.rm = T) != 
                             0), which(colSums(conmat, na.rm = T) != 0), drop = F]
  }
  if (all.nodes == FALSE && directed == TRUE) {
    include.vec <- vector(length = dim(data)[1])
    for (i in 1:dim(conmat)[1]) {
      ifelse(any(conmat[i, ] != 0) | any(conmat[, i] != 
                                           0), include.vec[i] <- 1, include.vec[i] <- 0)
    }
  }
  if (all.nodes == TRUE) {
    include.vec <- vector(length = dim(data)[1])
    include.vec <- rep(1, length = dim(data)[1])
  }
  p <- build_plot2(conmat = conmat, data = data, background = background, 
                   node.size = node.size, view = view, node.color = node.color, 
                   thr = thr, uthr = uthr, edge.color = edge.color, edge.alpha = edge.alpha, 
                   edge.width = edge.width, scale.edge.width = scale.edge.width, 
                   show.legend = show.legend, labels = labels, label.size = label.size, 
                   include.vec = include.vec, edge.color.weighted = edge.color.weighted, 
                   label.edge.weight = label.edge.weight, bg_xmax = bg_xmax, 
                   bg_xmin = bg_xmin, bg_ymax = bg_ymax, bg_ymin = bg_ymin)
  return(p)
}


build_plot2 <- function(conmat,
                        data,
                        data.row=NULL,
                        data.col=NULL,
                        background,
                        node.size,
                        node.color="network",
                        thr=NULL,
                        uthr=NULL,
                        view,
                        edge.color,
                        edge.alpha,
                        edge.width,
                        show.legend,
                        label.size,
                        labels,
                        include.vec=NULL,
                        scale.edge.width,
                        edge.color.weighted,
                        label.edge.weight,
                        bg_xmin=0,
                        bg_ymin=0,
                        bg_xmax=0,
                        bg_ymax=0,
                        ...) {
  
  
  if (view =="top"){
    x.mni<-data$x.mni
    y.mni<-data$y.mni
    depth <- data$z.mni
    xmax = 70     + bg_xmax
    xmin = -75    + bg_xmin
    ymax = 73     + bg_ymax
    ymin = -107   + bg_ymin
  }
  
  if (view =="bottom"){
    x.mni<-data$x.mni*-1
    y.mni<-data$y.mni
    depth <- data$z.mni*-1
    xmax = 70     + bg_xmax
    xmin = -70    + bg_xmin
    ymax = 73     + bg_ymax
    ymin = -107   + bg_ymin
  }
  
  if (view =="front"){
    x.mni<-data$x.mni
    y.mni<-data$z.mni
    depth <- data$y.mni
    xmax = 70     + bg_xmax
    xmin = -70    + bg_xmin
    ymax = 80     + bg_ymax
    ymin = -48    + bg_ymin
  }
  
  
  if (view =="back"){
    x.mni<-data$x.mni*-1
    y.mni<-data$z.mni
    depth <- data$y.mni*-1
    xmax = 70    + bg_xmax
    xmin = -70   + bg_xmin
    ymax = 80    + bg_ymax
    ymin = -48   + bg_ymin
  }
  
  
  if (view =="left"){
    x.mni<-data$y.mni*-1
    y.mni<-data$z.mni
    depth <- data$x.mni
    xmax = 103   + bg_xmax
    xmin = -72   + bg_xmin
    ymax = 77    + bg_ymax
    ymin = -50   + bg_ymin
  }
  
  ##fix below
  if (view =="right"){
    x.mni<-data$y.mni
    y.mni<-data$z.mni
    depth <- data$x.mni*-1
    xmax = 103   + bg_xmax
    xmin = -140  + bg_xmin
    ymax = 77    + bg_ymax
    ymin = -50   + bg_ymin
  }
  
  
  #is matrix directed (i.e. symetric)
  ifelse(isSymmetric.matrix(conmat)==TRUE,
         directed <- FALSE,
         directed <- TRUE)
  
  #is matrix weighed
  ifelse(all(conmat %in% c(0,1))==TRUE,
         weighted <- FALSE,
         weighted <- TRUE)
  
  #should edges be colored by weight
  #  ifelse(edge.color=="weight", edge.color.weighted <- T, edge.color.weighted <- F)
  
  if (!exists("conmat")) stop(print("Please enter a valid connectivity matrix"))
  
  
  if(directed == F) {
    conmat[upper.tri(conmat)] <- 0 #only take bottom tri of matrix to stop the edge labels being plotted twice
    layout <- create_layout(graph = conmat, layout ="stress", circular=TRUE)
    layout$x <- x.mni
    layout$y <- y.mni
  }
  
  
  
  if(directed == T) {
    layout <- create_layout(graph = conmat, layout ="stress", circular=TRUE)
    layout$x <- x.mni
    layout$y <- y.mni
    layout$facet <- include.vec
  }
  
  
  #make graph
  
  if(directed == T && weighted==F){
    p <- ggraph(layout) +
      annotation_custom(background, xmax = xmax ,xmin = xmin , ymax = ymax , ymin = ymin ) +
      geom_edge_parallel(color=edge.color,
                         edge_width = edge.width,
                         edge_alpha = edge.alpha,
                         arrow = arrow(length = unit(3, 'mm')),
                         end_cap = circle((node.size/2)+0.6, 'mm'))  +
      #  ggraph::geom_edge_loop0(aes(strength=node.size*3), color=edge.color, edge_width = edge.width, arrow = arrow(length = unit(1, 'mm'))) +
      coord_fixed(xlim = c(-70,70), ylim = c(-107,73))
  }
  
  if(directed == T && weighted==T && edge.color.weighted==F && label.edge.weight==F){p <- ggraph(layout) +
    annotation_custom(background, xmax = xmax ,xmin = xmin , ymax = ymax , ymin = ymin ) +
    geom_edge_parallel(aes(width=weight),
                       color=edge.color,
                       edge_alpha = edge.alpha,
                       arrow = arrow(length = unit(3, 'mm')),
                       end_cap = circle(node.size/2, 'mm'))  +
    geom_edge_loop0(aes(strength=node.size*3, width=weight),
                    color=edge.color,
                    edge_alpha = edge.alpha,
                    arrow = arrow(length = unit(3, 'mm'))) +
    coord_fixed(xlim = c(-70,70), ylim = c(-107,73))
  
  }
  
  if(directed == T && weighted==T && edge.color.weighted==F && label.edge.weight==T){p <- ggraph(layout) +
    annotation_custom(background, xmax = xmax ,xmin = xmin , ymax = ymax , ymin = ymin ) +
    geom_edge_parallel(aes(width=weight, label=round(weight,3)),
                       color=edge.color,
                       edge_alpha = edge.alpha,
                       arrow = arrow(length = unit(3, 'mm')),
                       end_cap = circle(node.size/2, 'mm'),
                       angle_calc = 'along',
                       alpha = 0,
                       label_dodge = unit(2.5, 'mm'),
                       label_size = 2,
                       fontface = "bold")  +
    geom_edge_loop0(aes(strength=node.size*3, width=weight, label=round(weight,3)),
                    color=edge.color,
                    edge_alpha = edge.alpha,
                    arrow = arrow(length = unit(3, 'mm')),
                    angle_calc = 'none',
                    alpha = 0,
                    label_dodge = unit(6, 'mm'),
                    label_size = 2,
                    vjust = -1,
                    fontface = "bold") +
    coord_fixed(xlim = c(-70,70), ylim = c(-107,73))
  
  }
  
  if(directed == T && weighted==T && edge.color.weighted==T && label.edge.weight==F){
    p <- ggraph(layout) +
      annotation_custom(background, xmax = xmax ,xmin = xmin , ymax = ymax , ymin = ymin ) +
      geom_edge_parallel(aes(color=log(weight)),
                         edge_alpha = edge.alpha,
                         edge_width = edge.width,
                         arrow = arrow(length = unit(3, 'mm')),
                         end_cap = circle(node.size/2, 'mm')) +
      geom_edge_loop(aes(strength=node.size*3, color=log(weight)),
                     edge_width = edge.width,
                     edge_alpha = edge.alpha,
                     arrow = arrow(length = unit(3, 'mm'))) +
      coord_fixed(xlim = c(-70,70), ylim = c(-107,73))
    
  }
  
  
  
  if(directed == T && weighted==T && edge.color.weighted==T && label.edge.weight==T){p <- ggraph(layout) +
    annotation_custom(background, xmax = xmax ,xmin = xmin , ymax = ymax , ymin = ymin ) +
    geom_edge_parallel(aes(color=log(weight), label=round(log(weight),3)),
                       #color=edge.color,
                       edge_alpha = edge.alpha,
                       edge_width = edge.width,
                       arrow = arrow(length = unit(3, 'mm')),
                       end_cap = circle(node.size/2, 'mm'),
                       angle_calc = 'along',
                       alpha = 0,
                       label_dodge = unit(2.5, 'mm'),
                       label_size = 2,
                       fontface = "bold") +
    geom_edge_loop(aes(strength=node.size*3, color=log(weight), label=round(log(weight),3)),
                   edge_width = edge.width,
                   edge_alpha = edge.alpha,
                   arrow = arrow(length = unit(3, 'mm')),
                   angle_calc = 'none',
                   alpha = 0,
                   label_dodge = unit(6, 'mm'),
                   label_size = 2,
                   vjust = -1,
                   fontface = "bold") +
    coord_fixed(xlim = c(-70,70), ylim = c(-107,73))
  
  }
  
  
  
  if(directed == F && weighted==F){p <- ggraph(layout, circular = FALSE) +
    annotation_custom(background, xmax = xmax ,xmin = xmin , ymax = ymax , ymin = ymin ) +
    geom_edge_link(color=edge.color,
                   edge_width = edge.width,
                   edge_alpha = edge.alpha) +
    coord_fixed(xlim = c(-70,70), ylim = c(-107,73))
  }
  
  
  
  
  if(directed == F && weighted==T && edge.color.weighted==F && label.edge.weight==F){
    p <- ggraph(layout, circular = FALSE) +
      annotation_custom(background, xmax = xmax ,xmin = xmin , ymax = ymax , ymin = ymin ) +
      geom_edge_link(aes(width=weight),
                     color=edge.color,
                     edge_alpha = edge.alpha) +
      coord_fixed(xlim = c(-70,70), ylim = c(-107,73))
  }
  
  if(directed == F && weighted==T && edge.color.weighted==F && label.edge.weight==T){
    p <- ggraph(layout, circular = FALSE) +
      annotation_custom(background, xmax = xmax ,xmin = xmin , ymax = ymax , ymin = ymin ) +
      geom_edge_link(aes(width=weight, label=round(weight,3)),
                     color=edge.color,
                     edge_alpha = edge.alpha,
                     angle_calc = 'along',
                     alpha = 0,
                     label_dodge = unit(2.5, 'mm'),
                     label_size = 2,
                     fontface = "bold") +
      coord_fixed(xlim = c(-70,70), ylim = c(-107,73))
  }
  if(directed == F && weighted==T && edge.color.weighted==T && label.edge.weight==F){
    flat_matrix <- as.vector(conmat) # 过滤掉零值
    non_zero_elements <- flat_matrix[flat_matrix != 0]
    # a1=max(non_zero_elements)-min(non_zero_elements)
    # range=c(max(non_zero_elements),
    #         0.8*(max(non_zero_elements)-1),
    #         1,
    #         0.8*(min(non_zero_elements)-1),
    #         min(non_zero_elements))
    rounded_value=1
    range=c(1*rounded_value,
            0.8*rounded_value,
            0.6*rounded_value,
            0.4*rounded_value,
            0,
            -0.4*rounded_value,
            -0.6*rounded_value,
            -0.8*rounded_value,
            -1*rounded_value)
    color_scale = c("#FF0000", "#F94D4B", "#F08E88", "#F1C2B9", "#fdffb3", "#D3EDDD", "#94BBEC", "#4579BF", "#0000FF")
    color_scale = c("darkred", "#F94D4B", "#FF0000", "#F08E88", "#fdffb3", "#D3EDDD", "#4579BF", "#0000FF","darkblue")
    
    # 定义颜色端点1
    color_ends <- c("darkred", "#FF0000","#F1C2B9", "white", "#94BBEC", "#0000FF", "darkblue")
    
    # 创建平滑渐变
    color_scale <- colorRampPalette(color_ends)(10)
    
    
    # 定义颜色端点2
    color_ends <- c("darkred", "#FF0000","#F08E88", "#fdffb3", "white","#94BBEC", "#4579BF", "#0000FF", "darkblue")
    
    # 创建平滑渐变
    color_scale <- colorRampPalette(color_ends)(10)
    
    # 定义颜色端点3
    color_ends <- c("#A50F15", "#DE2D26", "#FB6A4A", "#FCAE91", 
                    "#FEE5D9",  # 中性浅色
                    "#EFF3FF", "#BDD7E7", "#6BAED6", "#2171B5", "#08306B")
    color_scale <- colorRampPalette(color_ends)(11)
    # if(min(non_zero_elements)>1){
    #   # range=c(max(non_zero_elements),
    #   #         0.6*(max(non_zero_elements)-1),
    #   #         0.4*(max(non_zero_elements)-1),
    #   #         1)
    #   range=c(1*rounded_value,
    #           0.8*rounded_value,
    #           0.6*rounded_value,
    #           0.4*rounded_value)
    #   color_scale=c("#FF0000", "#F94D4B", "#F08E88", "#F1C2B9")
    # }
    # if(max(non_zero_elements)<1){
    #   # range=c(1,
    #   #         0.4*(min(non_zero_elements)-1),
    #   #         0.6*(min(non_zero_elements)-1),
    #   #         min(non_zero_elements))
    #   range=c(1*rounded_value,
    #           0.5*rounded_value,
    #           0)
    #   color_scale=c("#94BBEC",
    #                 "#5C9BEB",
    #                 "blue")
    # }
    values_to_plot = scales::rescale(range)
    
    
    #color_scale = c("darkred", "#F94D4B", "#FF0000", "#F08E88", "#fdffb3", "#D3EDDD", "#4579BF", "#0000FF","darkblue")
    #color_scale <- rev(c("#000080", "#0000FF", "#4169E1", "#87CEEB", "#F0F8FF",
                        # "#FFE4E1", "#F08E88", "#FF6347", "#FF0000"))
    
    p <- ggraph(layout, circular = FALSE) +
      annotation_custom(background, xmax = xmax ,xmin = xmin , ymax = ymax , ymin = ymin ) +
      geom_edge_link(aes(colour=log(weight)),
                     edge_width = 0.4,
                     edge_alpha = edge.alpha) +
      # scale_edge_color_gradient2(low = "blue", mid = "#fdffb3", high = "red", midpoint = 0,limits=c(-0.5,0.5)) +
      scale_edge_color_gradientn(colors = color_scale,
                                 values = values_to_plot,
                                 na.value = "white",
                                 limits=c(-0.13,0.13)) +
      # scale_edge_color_stepsn(
      #   # colors = c("#FF0000", "#F1C2B9", "#fdffb3", "#D3EDDD", "#94BBEC", "#0000FF"),
      #   colors = c("#0000FF", "#4579BF","#94BBEC", "#D3EDDD", "#fdffb3", "#F1C2B9", "#F08E88", "#FF0000"),
      #   breaks  = c(-0.1, -0.005,-0.001, 0,0.001, 0.005, 0.1),
      #   limits = c(-1, 1)
      # ) +
      # scale_edge_color_stepsn(
      #   colors = c("darkblue","#0000FF", "#4579BF","#94BBEC", "#D3EDDD", "#fdffb3", "#F1C2B9", "#F94D4B", "#FF0000","darkred"),
      #   breaks = c(-1,-0.1, -0.005,-0.001,-0.0005, 0,0.0005,0.001, 0.005, 0.1,1),  # 扩展断点以包含无限
      #   values = scales::rescale(
      #     c(-1, -0.1, -0.005, 0, 0.005, 0.1, 1),
    #     to = c(0, 1)
    #   ),  # 将断点标准化到 [0,1] 范围
    #   limits = c(-1, 1),  # 根据数据实际范围调整，示例设为 [-1,1]
    #   na.value = "white",
    #   guide = guide_colorsteps(even.steps = T)  # 禁用均匀步长
    # ) +
    coord_fixed(xlim = c(-70,70), ylim = c(-107,73))
  }
  if(directed == F && weighted==T && edge.color.weighted==T && label.edge.weight==T){
    p <- ggraph(layout, circular = FALSE) +
      annotation_custom(background, xmax = xmax ,xmin = xmin , ymax = ymax , ymin = ymin ) +
      geom_edge_link(aes(colour=log(weight), label=round(log(weight),3)),
                     edge_width = edge.width,
                     edge_alpha = edge.alpha,
                     angle_calc = 'along',
                     alpha = 0,
                     label_dodge = unit(2.5, 'mm'),
                     label_size = 2,
                     fontface = "bold") +
      coord_fixed(xlim = c(-70,70), ylim = c(-107,73))
  }
  if (weighted==T && !is.null(scale.edge.width)){
    p <- p + scale_edge_width(range = scale.edge.width)
  }
  if(view=="left") {
    p <- p + coord_fixed(xlim = c(-64,98), ylim = c(-44,76)) }
  if(view=="right") {
    p <- p + coord_fixed(xlim = c(-98,64), ylim = c(-44,76)) }
  if(directed == F){
    # ifelse(node.color=="network",
    #        p <- p + geom_node_point(size=node.size, aes(colour=as.factor(data$network))),
    #        p <- p + geom_node_point(size=node.size, colour=node.color))
    nature_colors <- c(
      Vis         = "#4E79A7",  # 沉稳蓝 (视觉网络)
      SomMot      = "#59A14F",  # 森林绿 (感觉运动)
      DorsAttn    = "#E15759",  # 暗红色 (背侧注意)
      SalVentAttn = "#B07AA1",  # 灰紫色 (腹侧注意)
      Limbic      = "#F28E2B",  # 橙黄色 (边缘系统)
      Cont        = "#76B7B2",  # 青绿色 (控制网络)
      Default     = "#EDC948",  # 浅金色 (默认模式网络)
      Subcortical = "#A0A0A0"   # 中性灰 (皮层下结构)
    )
    
    # 确保 network 列是因子且顺序与颜色一致
    data$network <- factor(data$network, levels = names(nature_colors))
    
    # 修改后的绘图代码
    p <- p + geom_node_point(size = node.size, aes(colour = data$network)) + 
      scale_colour_manual(values = nature_colors)
  }
  
  ## add labs
  if(directed == T && labels==T){
    p <- p + geom_node_text(aes(label = data$ROI.Name, filter = as.logical(facet)),
                            size=label.size, repel=TRUE,
                            nudge_x = node.size+2, nudge_y = node.size)
  }
  
  if (directed == F && labels==T){
    p <- p + geom_node_text(aes(label = data$ROI.Name),
                            size=label.size, repel=TRUE,
                            nudge_x = node.size+2, nudge_y = node.size)
  }
  #remove gridlines
  p <- p + theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y =element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          
    )
  #legend
  if (show.legend==F){p <- p + theme(legend.position="none")}
  if (show.legend==T){p <- p + scale_colour_manual(values = nature_colors)}
  p
}