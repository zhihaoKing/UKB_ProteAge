setwd("/Users/zhihaojin/Documents/neuroscience_PhD_in_FDU/PhD_project/bioinformatics/UKBiobank_version_2")
## http://www.sthda.com/english/wiki/cox-proportional-hazards-model
## https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html

library(tidyverse)
library(ggplot2)
library(cowplot)     
library(ggforestplot)
library(psych)
library(broom)        
library(gtsummary)
library(htmlTable)
library(tableone)   
library(survival)    
library(Hmisc)
library(rms)
library(splines) 
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(haven)
library(patchwork)
library(grid)
library(gridExtra)
library(splines)
library(ggeffects)
library(ggpubr)
library(Metrics)
library(Hmisc)
library(reshape2)
library(RColorBrewer)
library(ggsci)
library(ggview)
library(smplot2)
library(ggpmisc)
library(cowplot)
library(doMC)
library(ClusterGVis)
library(limma)
library(tidyverse)
library(readr)
library(magrittr)
library(dplyr)
library(BioAge)
library(haven)
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)
library(gridExtra)
library(psych)
library(splines) # for natural spline model to calculate BA residuals
library(ggeffects)
library(ggpubr)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggview)
library(pROC)
library(ggvenn)
library(ggVennDiagram)
library(pheatmap)
library(mice)
rm(list=ls());gc()
# female proteage #992224
# female KDM #E3625D
# female phenoage #F0C284
# male proteage #3E4F94
# male KDM #3E90BF
# male phenoage #58B6E9

vvisCluster <- function (object = NULL, ht.col.list = list(col_range = c(-2, 
    0, 2), col_color = c("#3E4F94", "white", "#992224")), border = TRUE, 
    plot.type = c("line", "heatmap", "both"), ms.col = c("#0099CC", 
        "grey90", "#CC3333"), line.size = 0.1, line.col = "grey90", 
    add.mline = TRUE, mline.size = 2, mline.col = "#CC3333", 
    ncol = 4, ctAnno.col = NULL, set.md = "median", textbox.pos = c(0.5, 
        0.8), textbox.size = 8, panel.arg = c(2, 0.25, 4, "grey90", 
        NA), ggplot.panel.arg = c(2, 0.25, 4, "grey90", NA), 
    annoTerm.data = NULL, annoTerm.mside = "right", termAnno.arg = c("grey95", 
        "grey50"), add.bar = FALSE, bar.width = 8, textbar.pos = c(0.8, 
        0.8), go.col = NULL, go.size = NULL, by.go = "anno_link", 
    annoKegg.data = NULL, annoKegg.mside = "right", keggAnno.arg = c("grey95", 
        "grey50"), add.kegg.bar = FALSE, kegg.col = NULL, kegg.size = NULL, 
    by.kegg = "anno_link", word_wrap = TRUE, add_new_line = TRUE, 
    add.box = FALSE, boxcol = NULL, box.arg = c(0.1, "grey50"), 
    add.point = FALSE, point.arg = c(19, "orange", "orange", 
        1), add.line = TRUE, line.side = "right", markGenes = NULL, 
    markGenes.side = "right", genes.gp = c("italic", 10, NA), 
    term.text.limit = c(10, 18), mulGroup = NULL, lgd.label = NULL, 
    show_row_names = FALSE, subgroup.anno = NULL, annnoblock.text = TRUE, 
    annnoblock.gp = c("white", 8), add.sampleanno = TRUE, sample.group = NULL, 
    sample.col = NULL, sample.order = NULL, cluster.order = NULL, 
    sample.cell.order = NULL, HeatmapAnnotation = NULL, column.split = NULL, 
    cluster_columns = FALSE, pseudotime_col = NULL, gglist = NULL, 
    ...) 
{
    ComplexHeatmap::ht_opt(message = FALSE)
    if (is.null(ht.col.list[["col_range"]])) {
        col_range = c(-2, 0, 2)
    }
    else {
        col_range = ht.col.list[["col_range"]]
    }
    if (is.null(ht.col.list[["col_color"]])) {
        col_color = c("#3E4F94", "white", "#992224")
    }
    else {
        col_color = ht.col.list[["col_color"]]
    }
    col_fun = circlize::colorRamp2(col_range, col_color)
    plot.type <- match.arg(plot.type)
    if (plot.type == "line") {
        data <- data.frame(object$long.res)
        if (!is.null(sample.order)) {
            data$cell_type <- factor(data$cell_type, levels = sample.order)
        }
        line <- ggplot2::ggplot(data, ggplot2::aes(x = cell_type, 
            y = norm_value))
        if (object$type == "mfuzz") {
            line <- line + ggplot2::geom_line(ggplot2::aes(color = membership, 
                group = gene), size = line.size) + ggplot2::scale_color_gradient2(low = ms.col[1], 
                mid = ms.col[2], high = ms.col[3], midpoint = 0.5)
        }
        else {
            line <- line + ggplot2::geom_line(ggplot2::aes(group = gene), 
                color = line.col, size = line.size)
        }
        if (add.mline == TRUE) {
            if (object$type == "wgcna") {
                linec <- unique(data$modulecol)
                names(linec) <- linec
                line <- line + ggplot2::geom_line(stat = "summary", 
                  fun = "median", size = mline.size, ggplot2::aes(group = 1, 
                    color = modulecol)) + ggplot2::scale_color_manual(values = linec)
            }
            else {
                line <- line + ggplot2::geom_line(stat = "summary", 
                  fun = "median", colour = mline.col, size = mline.size, 
                  ggplot2::aes(group = 1))
            }
        }
        else {
            line <- line
        }
        line1 <- line + ggplot2::theme_classic(base_size = 14) + 
            ggplot2::ylab("Normalized expression") + ggplot2::xlab("") + 
            ggplot2::theme(axis.ticks.length = ggplot2::unit(0.1, 
                "cm"), axis.text.x = ggplot2::element_text(angle = 45, 
                hjust = 1, color = "black"), strip.background = ggplot2::element_blank()) + 
            ggplot2::facet_wrap(~cluster_name, ncol = ncol, scales = "free")
        return(line1)
    }
    else {
        data <- data.frame(object$wide.res, check.names = FALSE) %>% 
            dplyr::arrange(as.numeric(as.character(cluster)))
        if (object$type == "mfuzz") {
            mat <- data %>% dplyr::arrange(as.numeric(as.character(cluster))) %>% 
                dplyr::select(-gene, -cluster, -membership)
        }
        else if (object$type == "wgcna") {
            mat <- data %>% dplyr::arrange(as.numeric(as.character(cluster))) %>% 
                dplyr::select(-gene, -cluster, -modulecol)
        }
        else if (object$type == "scRNAdata") {
            mat <- data %>% dplyr::arrange(as.numeric(as.character(cluster))) %>% 
                dplyr::select(-gene, -cluster)
        }
        else if (object$type == "monocle") {
            mat <- data %>% dplyr::arrange(as.numeric(as.character(cluster))) %>% 
                dplyr::select(-gene, -cluster)
        }
        else {
            mat <- data %>% dplyr::arrange(as.numeric(as.character(cluster))) %>% 
                dplyr::select(-gene, -cluster)
        }
        rownames(mat) <- data$gene
        if (object$geneMode == "all" | ncol(mat) > 20) {
            use_raster = TRUE
        }
        else {
            use_raster = FALSE
        }
        if (!is.null(sample.order)) {
            mat <- mat[, sample.order]
        }
        cl.info <- data.frame(table(data$cluster)) %>% dplyr::mutate(Var1 = as.numeric(as.character(Var1))) %>% 
            dplyr::arrange(Var1)
        cluster.num <- nrow(cl.info)
        subgroup <- lapply(1:nrow(cl.info), function(x) {
            nm <- rep(as.character(cl.info$Var1[x]), cl.info$Freq[x])
            paste("C", nm, sep = "")
        }) %>% unlist()
        if (!is.null(cluster.order)) {
            subgroup <- factor(subgroup, levels = paste("C", 
                cluster.order, sep = ""))
            cluster_row_slices = FALSE
        }
        else {
            cluster_row_slices = TRUE
        }
        if (object$geneMode == "all" & object$type == "scRNAdata") {
            celltype <- sapply(strsplit(colnames(mat), split = "\\|"), 
                "[", 2)
            cell.num.info <- table(celltype)[unique(celltype)]
            if (is.null(sample.cell.order)) {
                column_split = factor(rep(names(cell.num.info), 
                  cell.num.info), levels = unique(celltype))
            }
            else {
                column_split = factor(rep(names(cell.num.info), 
                  cell.num.info), levels = sample.cell.order)
            }
            if (is.null(sample.col)) {
                block.col = 1:length(cell.num.info)
            }
            else {
                block.col = sample.col
            }
        }
        else {
            if (is.null(sample.group)) {
                sample.info = colnames(mat)
                if (is.null(HeatmapAnnotation)) {
                  if (object$geneType == "branched") {
                    if (ncol(mat) == 200) {
                      column_split = rep(c("branch1", "branch2"), 
                        each = 100)
                    }
                    else {
                      column_split = rep(levels(object$pseudotime), 
                        rep(100, ncol(mat)/100))
                    }
                  }
                  else {
                    column_split = NULL
                  }
                }
                else {
                  column_split = column.split
                }
            }
            else {
                sample.info = sample.group
                column_split = sample.group
            }
            if (object$type != "monocle") {
                sample.info <- factor(sample.info, levels = unique(sample.info))
                if (is.null(sample.col)) {
                  scol <- circlize::rand_color(n = length(sample.info))
                  names(scol) <- sample.info
                }
                else {
                  scol <- sample.col
                  names(scol) <- sample.info
                }
            }
            else {
                sample.info <- factor(object$pseudotime, levels = unique(object$pseudotime))
                if (is.null(pseudotime_col)) {
                  if (object$geneType == "branched") {
                    if (length(unique(object$pseudotime)) == 
                      3) {
                      pseudotime_col <- c("red", "grey80", "blue")
                    }
                    else {
                      pseudotime_col <- circlize::rand_color(n = length(unique(object$pseudotime)))
                    }
                  }
                  else {
                    pseudotime_col <- c("blue", "red")
                  }
                }
                else {
                  pseudotime_col <- pseudotime_col
                }
                if (is.null(sample.col)) {
                  if (object$type != "monocle") {
                    scol <- circlize::rand_color(n = length(sample.info))
                    names(scol) <- sample.info
                  }
                  else {
                    if (object$geneType == "branched") {
                      if (length(unique(object$pseudotime)) == 
                        3) {
                        scol <- rep(pseudotime_col, table(object$pseudotime)[unique(object$pseudotime)])
                        names(scol) <- sample.info
                      }
                      else {
                        scol <- rep(pseudotime_col, table(object$pseudotime)[unique(object$pseudotime)])
                        names(scol) <- sample.info
                      }
                    }
                    else {
                      scol <- (grDevices::colorRampPalette(pseudotime_col))(100)
                      names(scol) <- sample.info
                    }
                  }
                }
                else {
                  scol <- sample.col
                  names(scol) <- sample.info
                }
            }
        }
        if (add.sampleanno == TRUE) {
            if (object$geneMode == "all" & object$type == "scRNAdata") {
                topanno = ComplexHeatmap::HeatmapAnnotation(cluster = ComplexHeatmap::anno_block(gp = grid::gpar(fill = block.col), 
                  labels = NULL), show_annotation_name = FALSE)
            }
            else {
                if (is.null(HeatmapAnnotation)) {
                  topanno = ComplexHeatmap::HeatmapAnnotation(sample = sample.info, 
                    col = list(sample = scol), gp = grid::gpar(col = ifelse(object$type == 
                      "monocle", NA, "white")), show_legend = ifelse(object$type == 
                      "monocle", FALSE, TRUE), show_annotation_name = FALSE)
                }
                else {
                  topanno = HeatmapAnnotation
                }
            }
        }
        else {
            topanno = NULL
        }
        if (is.null(ctAnno.col)) {
            colanno <- jjAnno::useMyCol("stallion", n = cluster.num)
        }
        else {
            colanno <- ctAnno.col
        }
        names(colanno) <- 1:cluster.num
        align_to = split(1:nrow(mat), subgroup)
        anno.block <- ComplexHeatmap::anno_block(align_to = align_to, 
            panel_fun = function(index, nm) {
                npos = as.numeric(unlist(strsplit(nm, split = "C"))[2])
                grid::grid.rect(gp = grid::gpar(fill = colanno[npos], 
                  col = NA))
                if (annnoblock.text == TRUE) {
                  grid::grid.text(label = paste("n:", length(index), 
                    sep = ""), rot = 90, gp = grid::gpar(col = annnoblock.gp[1], 
                    fontsize = as.numeric(annnoblock.gp[2])))
                }
            }, which = "row")
        if (!is.null(markGenes)) {
            rowGene <- rownames(mat)
            annoGene <- markGenes
            gene.col <- data %>% dplyr::select(gene, cluster) %>% 
                dplyr::filter(gene %in% annoGene)
            gene.col <- purrr::map_df(1:cluster.num, function(x) {
                tmp <- gene.col %>% dplyr::filter(as.numeric(cluster) == 
                  x) %>% dplyr::mutate(col = colanno[x])
            })
            gene.col <- gene.col[match(annoGene, gene.col$gene), 
                ]
            if (is.na(genes.gp[3])) {
                gcol = gene.col$col
            }
            else {
                gcol = genes.gp[3]
            }
            index <- match(annoGene, rowGene)
            geneMark = ComplexHeatmap::anno_mark(at = index, 
                labels = annoGene, which = "row", side = markGenes.side, 
                labels_gp = grid::gpar(fontface = genes.gp[1], 
                  fontsize = as.numeric(genes.gp[2]), col = gcol))
        }
        else {
            geneMark = NULL
        }
        right_annotation = ComplexHeatmap::rowAnnotation(gene = geneMark, 
            cluster = anno.block)
        if (object$type == "monocle" | object$geneMode == "all" | 
            ncol(mat) > 20) {
            show_column_names = FALSE
        }
        else {
            show_column_names = TRUE
        }
        if (object$geneType == "non-branched") {
            rg <- range(as.numeric(as.character(sample.info)))
            col_fun2 = circlize::colorRamp2(c(rg[1], rg[2]), 
                pseudotime_col)
            lgd = ComplexHeatmap::Legend(col_fun = col_fun2, 
                title = "pseudotime")
            lgd_list = list(lgd)
        }
        else if (object$geneType == "branched") {
            if (length(levels(sample.info)) == 3) {
                lgd = ComplexHeatmap::Legend(labels = levels(sample.info), 
                  legend_gp = grid::gpar(fill = pseudotime_col), 
                  title = "branch")
            }
            else {
                lgd = ComplexHeatmap::Legend(labels = levels(sample.info), 
                  legend_gp = grid::gpar(fill = pseudotime_col), 
                  title = "branch")
            }
            lgd_list = list(lgd)
        }
        else {
            lgd_list <- NULL
        }
        if (plot.type == "heatmap") {
            htf <- ComplexHeatmap::Heatmap(as.matrix(mat), name = "Z-score", 
                cluster_columns = cluster_columns, show_row_names = show_row_names, 
                border = border, column_split = column_split, 
                row_split = subgroup, cluster_row_slices = cluster_row_slices, 
                column_names_side = "top", show_column_names = show_column_names, 
                top_annotation = topanno, right_annotation = right_annotation, 
                col = col_fun, use_raster = use_raster, ...)
            ComplexHeatmap::draw(htf, merge_legend = TRUE, annotation_legend_list = lgd_list)
        }
        else {
            rg = range(mat)
            if (!is.null(gglist)) {
                anno_ggplot2 = ComplexHeatmap::anno_zoom(align_to = align_to, 
                  which = "row", panel_fun = function(index, 
                    nm) {
                    g <- gglist[[nm]]
                    g <- grid::grid.grabExpr(print(g))
                    grid::pushViewport(grid::viewport())
                    grid::grid.rect()
                    grid::grid.draw(g)
                    grid::popViewport()
                  }, size = grid::unit(as.numeric(ggplot.panel.arg[1]), 
                    "cm"), gap = grid::unit(as.numeric(ggplot.panel.arg[2]), 
                    "cm"), width = grid::unit(as.numeric(ggplot.panel.arg[3]), 
                    "cm"), side = "right", link_gp = grid::gpar(fill = ggplot.panel.arg[4], 
                    col = ggplot.panel.arg[5]))
            }
            else {
                anno_ggplot2 = NULL
            }
            panel_fun = function(index, nm) {
                if (add.box == TRUE & add.line != TRUE) {
                  xscale = c(-0.1, 1.1)
                }
                else {
                  xscale = c(0, 1)
                }
                grid::pushViewport(grid::viewport(xscale = xscale, 
                  yscale = c(0, 1)))
                grid::grid.rect()
                if (object$geneMode == "all" & object$type == 
                  "scRNAdata") {
                  mulGroup <- cell.num.info
                  grid::grid.lines(x = c(0, 1), y = rep(0.5, 
                    2), gp = grid::gpar(col = "black", lty = "dashed"))
                  cu <- cumsum(mulGroup)
                  seqn <- data.frame(st = c(1, cu[1:(length(cu) - 
                    1)] + 1), sp = c(cu[1], cu[2:length(cu)]))
                }
                else {
                  if (is.null(mulGroup)) {
                    mulGroup <- ncol(mat)
                    seqn <- data.frame(st = 1, sp = ncol(mat))
                  }
                  else {
                    mulGroup <- mulGroup
                    grid::grid.lines(x = c(0, 1), y = rep(0.5, 
                      2), gp = grid::gpar(col = "black", lty = "dashed"))
                    cu <- cumsum(mulGroup)
                    seqn <- data.frame(st = c(1, cu[1:(length(cu) - 
                      1)] + 1), sp = c(cu[1], cu[2:length(cu)]))
                  }
                }
                if (object$geneMode == "all" & object$type == 
                  "scRNAdata") {
                  cell.ave <- purrr::map_dfr(1:nrow(seqn), function(x) {
                    tmp <- seqn[x, ]
                    tmpmat <- mat[index, c(tmp$st:tmp$sp)]
                    rg <- range(mat[index, ])
                    if (set.md == "mean") {
                      mdia <- mean(rowMeans(tmpmat))
                    }
                    else if (set.md == "median") {
                      mdia <- stats::median(apply(tmpmat, 1, 
                        stats::median))
                    }
                    else {
                      message("supply mean/median !")
                    }
                    res <- data.frame(x = x, val = mdia)
                    return(res)
                  })
                  if (add.line == TRUE) {
                    grid::grid.lines(x = scales::rescale(cell.ave$x, 
                      to = c(0, 1)), y = scales::rescale(cell.ave$val, 
                      to = c(0.1, 0.9)), gp = grid::gpar(lwd = 3, 
                      col = mline.col))
                  }
                }
                else {
                  lapply(1:nrow(seqn), function(x) {
                    tmp <- seqn[x, ]
                    tmpmat <- mat[index, c(tmp$st:tmp$sp)]
                    if (set.md == "mean") {
                      mdia <- colMeans(tmpmat)
                    }
                    else if (set.md == "median") {
                      mdia <- apply(tmpmat, 2, stats::median)
                    }
                    else {
                      message("supply mean/median !")
                    }
                    pos = scales::rescale(1:ncol(tmpmat), to = c(0, 
                      1))
                    if (is.null(boxcol)) {
                      boxcol <- rep("grey90", ncol(tmpmat))
                    }
                    else {
                      boxcol <- boxcol
                    }
                    if (add.box == TRUE) {
                      lapply(1:ncol(tmpmat), function(x) {
                        ComplexHeatmap::grid.boxplot(scales::rescale(tmpmat[, 
                          x], to = c(0, 1), from = c(rg[1] - 
                          0.5, rg[2] + 0.5)), pos = pos[x], direction = "vertical", 
                          box_width = as.numeric(box.arg[1]), 
                          outline = FALSE, gp = grid::gpar(col = box.arg[2], 
                            fill = boxcol[x]))
                      })
                    }
                    if (add.point == TRUE) {
                      grid::grid.points(x = scales::rescale(1:ncol(tmpmat), 
                        to = c(0, 1)), y = scales::rescale(mdia, 
                        to = c(0, 1), from = c(rg[1] - 0.5, rg[2] + 
                          0.5)), pch = as.numeric(point.arg[1]), 
                        gp = grid::gpar(fill = point.arg[2], 
                          col = point.arg[3]), size = grid::unit(as.numeric(point.arg[4]), 
                          "char"))
                    }
                    if (add.line == TRUE) {
                      grid::grid.lines(x = scales::rescale(1:ncol(tmpmat), 
                        to = c(0, 1)), y = scales::rescale(mdia, 
                        to = c(0, 1), from = c(rg[1] - 0.5, rg[2] + 
                          0.5)), gp = grid::gpar(lwd = 3, col = mline.col[x]))
                    }
                  })
                }
                grid.textbox <- utils::getFromNamespace("grid.textbox", 
                  "ComplexHeatmap")
                text <- paste("Number of proteins:", nrow(mat[index, ]), 
                  sep = " ")
                grid.textbox(text, x = textbox.pos[1], y = textbox.pos[2], 
                  gp = grid::gpar(fontsize = textbox.size, fontface = "italic", 
                    ...))
                grid::popViewport()
            }
            if (!is.null(subgroup.anno)) {
                align_to = split(1:nrow(mat), subgroup)
                align_to = align_to[subgroup.anno]
            }
            else {
                align_to = subgroup
            }
            anno = ComplexHeatmap::anno_link(align_to = align_to, 
                which = "row", panel_fun = panel_fun, size = grid::unit(as.numeric(panel.arg[1]), 
                  "cm"), gap = grid::unit(as.numeric(panel.arg[2]), 
                  "cm"), width = grid::unit(as.numeric(panel.arg[3]), 
                  "cm"), side = line.side, link_gp = grid::gpar(fill = panel.arg[4], 
                  col = panel.arg[5]))
            if (!is.null(annoTerm.data)) {
                termanno <- annoTerm.data
                if (ncol(termanno) == 2) {
                  colnames(termanno) <- c("id", "term")
                }
                else if (ncol(termanno) == 3) {
                  colnames(termanno) <- c("id", "term", "pval")
                }
                else if (ncol(termanno) == 4) {
                  colnames(termanno) <- c("id", "term", "pval", 
                    "ratio")
                }
                else {
                  message("No more than 4 columns!")
                }
                if (is.null(go.col)) {
                  gocol <- circlize::rand_color(n = nrow(termanno))
                }
                else {
                  gocol <- go.col
                }
                if (is.null(go.size)) {
                  gosize <- rep(12, nrow(termanno))
                }
                else {
                  if (go.size == "pval") {
                    termanno.tmp <- purrr::map_df(unique(termanno$id), 
                      function(x) {
                        tmp <- termanno %>% dplyr::filter(id == 
                          x) %>% dplyr::mutate(size = scales::rescale(-log10(pval), 
                          to = term.text.limit))
                      })
                    gosize <- termanno.tmp$size
                  }
                  else {
                    gosize <- go.size
                  }
                }
                termanno <- termanno %>% dplyr::ungroup() %>% 
                  dplyr::mutate(col = gocol, fontsize = gosize)
                term.list <- lapply(1:length(unique(termanno$id)), 
                  function(x) {
                    tmp = termanno[which(termanno$id == unique(termanno$id)[x]), 
                      ]
                    df <- data.frame(text = tmp$term, col = tmp$col, 
                      fontsize = tmp$fontsize)
                    return(df)
                  })
                names(term.list) <- unique(termanno$id)
                if (!is.null(subgroup.anno)) {
                  align_to2 = split(seq_along(subgroup), subgroup)
                  align_to2 = align_to2[subgroup.anno]
                  term.list = term.list[subgroup.anno]
                }
                else {
                  align_to2 = subgroup
                  term.list = term.list
                }
                textbox = ComplexHeatmap::anno_textbox(align_to2, 
                  term.list, word_wrap = word_wrap, add_new_line = add_new_line, 
                  side = annoTerm.mside, background_gp = grid::gpar(fill = termAnno.arg[1], 
                    col = termAnno.arg[2]), by = by.go)
                if (ncol(termanno) - 2 > 2) {
                  anno_gobar <- function(data = NULL, bar.width = 0.1, 
                    align_to = NULL, panel.arg = panel.arg, ...) {
                    if (ncol(data) - 2 == 3) {
                      data <- data %>% dplyr::mutate(bary = -log10(pval))
                    }
                    else {
                      data <- data %>% dplyr::mutate(bary = ratio)
                    }
                    ComplexHeatmap::anno_zoom(align_to = align_to, 
                      which = "row", panel_fun = function(index, 
                        nm) {
                        grid::pushViewport(grid::viewport(xscale = c(0, 
                          1), yscale = c(0, 1)))
                        grid::grid.rect()
                        tmp <- data %>% dplyr::filter(id == nm)
                        grid::grid.segments(x0 = rep(0, nrow(tmp)), 
                          x1 = scales::rescale(rev(tmp$bary), 
                            to = c(0.1, 0.9)), y0 = scales::rescale(1:nrow(tmp), 
                            to = c(0.1, 0.9)), y1 = scales::rescale(1:nrow(tmp), 
                            to = c(0.1, 0.9)), gp = grid::gpar(lwd = bar.width, 
                            col = rev(tmp$col), lineend = "butt"))
                        grid.textbox <- utils::getFromNamespace("grid.textbox", 
                          "ComplexHeatmap")
                        text <- nm
                        grid.textbox(text, x = textbar.pos[1], 
                          y = textbar.pos[2], gp = grid::gpar(fontsize = textbox.size, 
                            fontface = "italic", col = unique(tmp$col), 
                            ...))
                        grid::popViewport()
                      }, size = grid::unit(as.numeric(panel.arg[1]), 
                        "cm"), gap = grid::unit(as.numeric(panel.arg[2]), 
                        "cm"), width = grid::unit(as.numeric(panel.arg[3]), 
                        "cm"), side = "right", link_gp = grid::gpar(fill = termAnno.arg[1], 
                        col = termAnno.arg[2]), ...)
                  }
                  baranno = anno_gobar(data = termanno, align_to = align_to2, 
                    panel.arg = panel.arg, bar.width = bar.width)
                }
                if (add.bar == TRUE) {
                  baranno
                }
                else {
                  baranno = NULL
                }
            }
            else {
                textbox = NULL
                baranno = NULL
            }
            if (!is.null(annoKegg.data)) {
                termanno <- annoKegg.data
                if (ncol(termanno) == 2) {
                  colnames(termanno) <- c("id", "term")
                }
                else if (ncol(termanno) == 3) {
                  colnames(termanno) <- c("id", "term", "pval")
                }
                else if (ncol(termanno) == 4) {
                  colnames(termanno) <- c("id", "term", "pval", 
                    "ratio")
                }
                else {
                  message("No more than 4 columns!")
                }
                if (is.null(kegg.col)) {
                  gocol <- circlize::rand_color(n = nrow(termanno))
                }
                else {
                  gocol <- kegg.col
                }
                if (is.null(kegg.size)) {
                  gosize <- rep(12, nrow(termanno))
                }
                else {
                  if (kegg.size == "pval") {
                    termanno.tmp <- purrr::map_df(unique(termanno$id), 
                      function(x) {
                        tmp <- termanno %>% dplyr::filter(id == 
                          x) %>% dplyr::mutate(size = scales::rescale(-log10(pval), 
                          to = term.text.limit))
                      })
                    gosize <- termanno.tmp$size
                  }
                  else {
                    gosize <- kegg.size
                  }
                }
                termanno <- termanno %>% dplyr::ungroup() %>% 
                  dplyr::mutate(col = gocol, fontsize = gosize)
                term.list <- lapply(1:length(unique(termanno$id)), 
                  function(x) {
                    tmp = termanno[which(termanno$id == unique(termanno$id)[x]), 
                      ]
                    df <- data.frame(text = tmp$term, col = tmp$col, 
                      fontsize = tmp$fontsize)
                    return(df)
                  })
                names(term.list) <- unique(termanno$id)
                if (!is.null(subgroup.anno)) {
                  align_to2 = split(seq_along(subgroup), subgroup)
                  align_to2 = align_to2[subgroup.anno]
                  term.list = term.list[subgroup.anno]
                }
                else {
                  align_to2 = subgroup
                  term.list = term.list
                }
                textbox.kegg = ComplexHeatmap::anno_textbox(align_to2, 
                  term.list, word_wrap = word_wrap, add_new_line = add_new_line, 
                  side = annoKegg.mside, background_gp = grid::gpar(fill = keggAnno.arg[1], 
                    col = keggAnno.arg[2]), by = by.kegg)
                if (ncol(termanno) - 2 > 2) {
                  anno_keggbar <- function(data = NULL, bar.width = 0.1, 
                    align_to = NULL, panel.arg = panel.arg, ...) {
                    if (ncol(data) - 2 == 3) {
                      data <- data %>% dplyr::mutate(bary = -log10(pval))
                    }
                    else {
                      data <- data %>% dplyr::mutate(bary = ratio)
                    }
                    ComplexHeatmap::anno_zoom(align_to = align_to, 
                      which = "row", panel_fun = function(index, 
                        nm) {
                        grid::pushViewport(grid::viewport(xscale = c(0, 
                          1), yscale = c(0, 1)))
                        grid::grid.rect()
                        tmp <- data %>% dplyr::filter(id == nm)
                        grid::grid.segments(x0 = rep(0, nrow(tmp)), 
                          x1 = scales::rescale(rev(tmp$bary), 
                            to = c(0.1, 0.9)), y0 = scales::rescale(1:nrow(tmp), 
                            to = c(0.1, 0.9)), y1 = scales::rescale(1:nrow(tmp), 
                            to = c(0.1, 0.9)), gp = grid::gpar(lwd = bar.width, 
                            col = rev(tmp$col), lineend = "butt"))
                        grid.textbox <- utils::getFromNamespace("grid.textbox", 
                          "ComplexHeatmap")
                        text <- nm
                        grid.textbox(text, x = textbar.pos[1], 
                          y = textbar.pos[2], gp = grid::gpar(fontsize = textbox.size, 
                            fontface = "italic", col = unique(tmp$col), 
                            ...))
                        grid::popViewport()
                      }, size = grid::unit(as.numeric(panel.arg[1]), 
                        "cm"), gap = grid::unit(as.numeric(panel.arg[2]), 
                        "cm"), width = grid::unit(as.numeric(panel.arg[3]), 
                        "cm"), side = "right", link_gp = grid::gpar(fill = keggAnno.arg[1], 
                        col = keggAnno.arg[2]), ...)
                  }
                  baranno.kegg = anno_keggbar(data = termanno, 
                    align_to = align_to2, panel.arg = panel.arg, 
                    bar.width = bar.width)
                }
                if (add.kegg.bar == TRUE) {
                  baranno.kegg
                }
                else {
                  baranno.kegg = NULL
                }
            }
            else {
                textbox.kegg = NULL
                baranno.kegg = NULL
            }
            if (line.side == "right") {
                if (markGenes.side == "right") {
                  right_annotation2 = ComplexHeatmap::rowAnnotation(gene = geneMark, 
                    cluster = anno.block, line = anno, anno_ggplot2 = anno_ggplot2, 
                    textbox = textbox, bar = baranno, textbox.kegg = textbox.kegg, 
                    baranno.kegg = baranno.kegg)
                  left_annotation = NULL
                }
                else {
                  right_annotation2 = ComplexHeatmap::rowAnnotation(cluster = anno.block, 
                    line = anno, anno_ggplot2 = anno_ggplot2, 
                    textbox = textbox, bar = baranno, textbox.kegg = textbox.kegg, 
                    baranno.kegg = baranno.kegg)
                  left_annotation = ComplexHeatmap::rowAnnotation(gene = geneMark)
                }
            }
            else {
                if (markGenes.side == "right") {
                  right_annotation2 = ComplexHeatmap::rowAnnotation(gene = geneMark, 
                    cluster = anno.block, anno_ggplot2 = anno_ggplot2, 
                    textbox = textbox, bar = baranno, textbox.kegg = textbox.kegg, 
                    baranno.kegg = baranno.kegg)
                  left_annotation = ComplexHeatmap::rowAnnotation(line = anno)
                }
                else {
                  right_annotation2 = ComplexHeatmap::rowAnnotation(cluster = anno.block, 
                    anno_ggplot2 = anno_ggplot2, textbox = textbox, 
                    bar = baranno, textbox.kegg = textbox.kegg, 
                    baranno.kegg = baranno.kegg)
                  left_annotation = ComplexHeatmap::rowAnnotation(line = anno, 
                    gene = geneMark)
                }
            }
            if (object$type == "monocle" | object$geneMode == 
                "all" | ncol(mat) > 20) {
                show_column_names = FALSE
            }
            else {
                show_column_names = TRUE
            }
            htf <- ComplexHeatmap::Heatmap(as.matrix(mat), name = "Z-score", 
                cluster_columns = cluster_columns, show_row_names = show_row_names, 
                border = border, column_split = column_split, 
                top_annotation = topanno, right_annotation = right_annotation2, 
                left_annotation = left_annotation, column_names_side = "top", 
                show_column_names = show_column_names, row_split = subgroup, 
                cluster_row_slices = cluster_row_slices, col = col_fun, 
                use_raster = use_raster, ...)
            if (is.null(mulGroup)) {
                ComplexHeatmap::draw(htf, merge_legend = TRUE, 
                  annotation_legend_list = lgd_list)
            }
            else {
                if (is.null(lgd.label)) {
                  lgd.label <- paste("group", 1:length(mulGroup), 
                    sep = "")
                }
                else {
                  lgd.label <- lgd.label
                }
                lgd_list2 = ComplexHeatmap::Legend(labels = lgd.label, 
                  type = "lines", legend_gp = grid::gpar(col = mline.col, 
                    lty = 1))
                if (!is.null(lgd_list)) {
                  lgd_list_com <- ComplexHeatmap::packLegend(lgd_list, 
                    lgd_list2)
                }
                else {
                  lgd_list_com <- lgd_list2
                }
                ComplexHeatmap::draw(htf, annotation_legend_list = lgd_list_com, 
                  merge_legend = TRUE)
            }
        }
    }
}

### load data
p <- read.csv("data/protein_data.csv")
colnames(p)[-1] <- toupper(colnames(p)[-1])
f <- read.csv("data/protein_cohort_sex_age.csv") 
#f$X <- NULL

d <- merge(f,p,by = "eid") %>% dplyr ::rename(Sex = p31,Age = p21022)
d1 <- d %>% filter(Sex == "Male")  
d2 <- d %>% filter(Sex == "Female")


save(d1,d2,vvisCluster,file = "result/proteins_changed_with_age/data.RData")
rm(list=ls());gc()
load("result/proteins_changed_with_age/data.RData")

#### calculated correlation 
# male
for(i in 4:ncol(d1)) {
  d1[,i][is.na(d1[,i])] <- mean(d1[,i], na.rm=TRUE)
}
y1 <- as.numeric(d1[,"Age"])
d11 <- d1[,-c(1:3)]
colnames1 <- colnames(d11)
cor_data_df1 <- data.frame(colnames1)
for (i in 1:2923){
test1 <- cor.test(as.numeric(d11[,i]),y1,type="pearson")
cor_data_df1[i,2] <- test1$estimate
cor_data_df1[i,3] <- test1$p.value
}
names(cor_data_df1) <- c("symbol","correlation","pvalue")
cor_data_df1$FDR <- p.adjust(cor_data_df1$pvalue,method = "BH")




# female
for(i in 4:ncol(d2)) {
  d2[,i][is.na(d2[,i])] <- mean(d2[,i], na.rm=TRUE)
}
y2 <- as.numeric(d2[,"Age"])
d22 <- d2[,-c(1:3)]
colnames2 <- colnames(d22)
cor_data_df2 <- data.frame(colnames2)
for (i in 1:2923){
test2 <- cor.test(as.numeric(d22[,i]),y2,type="pearson")
cor_data_df2[i,2] <- test2$estimate
cor_data_df2[i,3] <- test2$p.value
}
names(cor_data_df2) <- c("symbol","correlation","pvalue")
cor_data_df2$FDR <- p.adjust(cor_data_df2$pvalue,method = "BH")




##### male cohort 
d1$Sex <- NULL
rownames(d1) <- d1$eid
d1$eid <- NULL
g1 <- list()
age <- c(40:70)
for(i in seq(age)){
  g1[[i]] <- d1 %>% filter(Age == age[i]) %>% as.tibble()
}
mean1 <- list()
for(i in seq(age)){
  mean1[[i]] <- apply(g1[[i]][,-1],2,function(x){mean(na.omit(x))})
}
names(mean1) <- age
dat1 <- data.frame(mean1[[1]],mean1[[2]],mean1[[3]],mean1[[4]],mean1[[5]],mean1[[6]],mean1[[7]],
                  mean1[[8]],mean1[[9]],mean1[[10]],mean1[[11]],mean1[[12]],mean1[[13]],mean1[[14]],
                  mean1[[15]],mean1[[16]],mean1[[17]],mean1[[18]],mean1[[19]],mean1[[20]],mean1[[21]],
                  mean1[[22]],mean1[[23]],mean1[[24]],mean1[[25]],mean1[[26]],mean1[[27]],mean1[[28]],
                  mean1[[29]],mean1[[30]],mean1[[31]])
colnames(dat1) <- age
dat1[,1] <- NULL
dat1 <- na.omit(dat1)

### 5 cluster
cm1 <- clusterData(exp = dat1,
                  cluster.method = "mfuzz",
                  cluster.num = 5,
                  seed = 987)
res1 <- cm1[["wide.res"]]
write.csv(res1,"result/proteins_changed_with_age/male_cohort.csv")
### up-regulated protein
res1_up <- res1 %>% filter(cluster == 2)
#write.csv(res1_up,"submit/supplyment/cluster_4_protein_male_cohort.csv")
res1_cor1_up <- merge(cor_data_df1,res1_up,by.x = "symbol",by.y = "gene") %>% filter(FDR < 0.05) %>% arrange(desc(correlation))
cor1_top_up <- res1_cor1_up[c(1:20),1]
#write.csv(res1_up,"result/select_protein/change_with_age/male_cohort_up_protein.csv")
g1_f_up <- mapIds(org.Hs.eg.db,res1_up$gene,'ENTREZID','SYMBOL') 

## kegg
kegg1_f_up <- enrichKEGG(g1_f_up,
                        organism = "hsa",
                        keyType = "kegg",
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        #minGSSize = 10,
                        #maxGSSize = 500,
                       qvalueCutoff = 0.1)
kegg1_f_up@result$GeneRatio <- sapply(kegg1_f_up@result$GeneRatio, function(x) eval(parse(text=x)))
kegg1_f_up@result$GeneRatio <- round(kegg1_f_up@result$GeneRatio,2)
d1_f_up <- kegg1_f_up@result
write.csv(d1_f_up %>% arrange(desc(Count)),"result/proteins_changed_with_age/kegg_go/kegg analysis of up-regulated proteins in male cohort.csv")
d1_f_up$logP <- -log10(d1_f_up$p.adjust)
d1_f_up <- d1_f_up %>% arrange(desc(Count))
d1_f_up$Description <- factor(d1_f_up$Description, levels = d1_f_up$Description)
d1_f_up <- d1_f_up[c(1:10),]
d1_f_up$group <- "Up-regulaed in Male Cohort"
p1_f_up <- ggplot(d1_f_up,aes(x = Count,y = Description,size = logP)) +
                  geom_point(color = "#F0C284") +
                  theme_classic(base_size = 14) +
                  theme(axis.text.y = element_text(size = 10,face = "bold"),
                        axis.text.x = element_text(size = 10,face = "bold"),
                        axis.title = element_text(size = 12,face = "bold"),
                        legend.title = element_text(face = "bold",size = 10),
                        plot.title = element_text(size = 16,face = "bold"),
                        legend.background = element_rect(linetype="solid",colour ="black"),
                        legend.key.size = unit(0.01, units = "cm"),
                        legend.position = c(0.7, 0.7)) +
                  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
                  #scale_x_continuous(limits = c(20,110),breaks = seq(20,110,10)) +
                  labs(y = "", x = "Count",
                       title = "KEGG analysis", 
                       size = "-logP.adj")
ggview(w = 4.7,h = 4.5)
ggsave("1st submission NC/Figures JZH 20241226/Figure 6/Figure 6E.pdf",w = 4.7,h = 4.5,dpi = 300)
#write.csv(d1_f_up,"result/select_protein/change_with_age/kegg_male_cohort_up_protein.csv")

## go
go1_f_up <- enrichGO(gene       = g1_f_up,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENTREZID',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE)
go1_f_up@result$GeneRatio <- sapply(go1_f_up@result$GeneRatio, function(x) eval(parse(text=x)))
go1_f_up@result$GeneRatio <- round(go1_f_up@result$GeneRatio,2)                    
d1_f_up <- go1_f_up@result
write.csv(d1_f_up %>% arrange(desc(Count)),"result/proteins_changed_with_age/kegg_go/go analysis of up-regulated proteins in male cohort.csv")
d1_f_up$logP <- -log10(d1_f_up$p.adjust)

### rank by count
d1_f_up <- d1_f_up %>% arrange(desc(Count))
d1_f_up$Description <- factor(d1_f_up$Description, levels = d1_f_up$Description)
d1_f_up <- d1_f_up[c(1:10),]
d1_f_up$group <- "Up-regulaed in Male Cohort"

p1_f_up <- ggplot(d1_f_up,aes(x = Count,y = Description,size = logP)) +
                  geom_point(color = "#F0C284") +
                  theme_classic(base_size = 14) +
                  theme(axis.text.y = element_text(size = 10,face = "bold"),
                        axis.text.x = element_text(size = 10,face = "bold"),
                        axis.title = element_text(size = 12,face = "bold"),
                        legend.title = element_text(face = "bold",size = 10),
                        plot.title = element_text(size = 16,face = "bold"),
                        legend.background = element_rect(linetype="solid",colour ="black"),
                        legend.key.size = unit(0.01, units = "cm"),
                        legend.position = c(0.8, 0.8)) +
                  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
                  labs(y = "", x = "Count",
                       title = "GO analysis", 
                       size = "-logP.adj")
ggview(w = 5.1,h = 4.5)
ggsave("1st submission NC/Supplemental figures JZH 20241226/supplemental figure 2 JZH 20241226/supplemental figure 2c JZH 20241226.pdf",w = 5.1,h = 4.5,dpi = 300)



### down-regulated protein
res1_down <- res1 %>% filter(cluster == 5)
#write.csv(res1_down,"submit/supplyment/cluster_1_protein_male_cohort.csv")
res1_cor1_down <- merge(cor_data_df1,res1_down,by.x = "symbol",by.y = "gene") %>% filter(FDR < 0.05) %>% arrange(correlation)
cor1_top_down <- res1_cor1_down[c(1:20),1]
#write.csv(res1_down,"result/select_protein/change_with_age/male_cohort_down_protein.csv")
g1_f_down <- mapIds(org.Hs.eg.db,res1_down$gene,'ENTREZID','SYMBOL') 

## kegg
kegg1_f_down <- enrichKEGG(g1_f_down,
                        organism = "hsa",
                        keyType = "kegg",
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        #minGSSize = 10,
                        #maxGSSize = 500,
                        qvalueCutoff = 0.1)
kegg1_f_down@result$GeneRatio <- sapply(kegg1_f_down@result$GeneRatio, function(x) eval(parse(text=x)))
kegg1_f_down@result$GeneRatio <- round(kegg1_f_down@result$GeneRatio,2)    
d1_f_down <- kegg1_f_down@result
write.csv(d1_f_down %>% arrange(desc(Count)),"result/proteins_changed_with_age/kegg_go/kegg analysis of down-regulated proteins in male cohort.csv")
d1_f_down$logP <- -log10(d1_f_down$p.adjust)
#write.csv(d1_f_down,"result/select_protein/change_with_age/kegg_male_cohort_down_protein.csv")
d1_f_down <- d1_f_down %>% arrange(desc(Count))
d1_f_down$Description <- factor(d1_f_down$Description, levels = d1_f_down$Description)
d1_f_down <- d1_f_down[c(1:10),]
d1_f_down$group <- "Down-regulaed in Male Cohort"
p1_f_down <- ggplot(d1_f_down,aes(x = Count,y = Description,size = logP)) +
                  geom_point(color = "#58B6E9") +
                  theme_classic(base_size = 14) +
                  theme(axis.text.y = element_text(size = 10,face = "bold"),
                        axis.text.x = element_text(size = 10,face = "bold"),
                        axis.title = element_text(size = 12,face = "bold"),
                        legend.title = element_text(face = "bold",size = 10),
                        plot.title = element_text(size = 16,face = "bold"),
                        legend.background = element_rect(linetype="solid",colour ="black"),
                        legend.key.size = unit(0.01, units = "cm"),
                        legend.position = c(0.7, 0.7)) +
                  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
                  #scale_x_continuous(limits = c(20,110),breaks = seq(20,110,10)) +
                  labs(y = "", x = "Count",
                       title = "KEGG analysis", 
                       size = "-logP.adj")
ggview(w = 4.7,h = 4.5)
ggsave("1st submission NC/Figures JZH 20241226/Figure 6/Figure 6F.pdf",w = 4.7,h = 4.5,dpi = 300)


## go
go1_f_down <- enrichGO(gene       = g1_f_down,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENTREZID',
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
go1_f_down@result$GeneRatio <- sapply(go1_f_down@result$GeneRatio, function(x) eval(parse(text=x)))
go1_f_down@result$GeneRatio <- round(go1_f_down@result$GeneRatio,2) 
d1_f_down <- go1_f_down@result
write.csv(d1_f_down %>% arrange(desc(Count)),"result/proteins_changed_with_age/kegg_go/go analysis of down-regulated proteins in male cohort.csv")
d1_f_down$logP <- -log10(d1_f_down$p.adjust)
#write.csv(d1_f_down,"result/select_protein/change_with_age/go_male_cohort_down_protein.csv")
### rank by Count
d1_f_down <- d1_f_down %>% arrange(desc(Count))
#d1_f_down[9,3] <- "PI3K/PKB signal transduction"
d1_f_down$Description <- factor(d1_f_down$Description, levels = d1_f_down$Description)
d1_f_down <- d1_f_down[c(1:10),]
d1_f_down$group <- "Down-regulaed in Male Cohort"
p1_f_down <- ggplot(d1_f_down,aes(x = Count,y = Description,size = logP)) +
                  geom_point(color = "#58B6E9") +
                  theme_classic(base_size = 14) +
                  theme(axis.text.y = element_text(size = 10,face = "bold"),
                        axis.text.x = element_text(size = 10,face = "bold"),
                        axis.title = element_text(size = 12,face = "bold"),
                        legend.title = element_text(face = "bold",size = 10),
                        plot.title = element_text(size = 16,face = "bold"),
                        legend.background = element_rect(linetype="solid",colour ="black"),
                        legend.key.size = unit(0.01, units = "cm"),
                        legend.position = c(0.83, 0.5)) +
                  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
                  #scale_x_continuous(limits = c(20,110),breaks = seq(20,110,10)) +
                  labs(y = "", x = "Count",
                       title = "GO analysis", 
                       size = "-logP.adj")
ggview(w = 5.1,h = 4.5)
ggsave("1st submission NC/Supplemental figures JZH 20241226/supplemental figure 2 JZH 20241226/supplemental figure 2d JZH 20241226.pdf",w = 5.1,h = 4.5,dpi = 300)

### rank by logP
#d1_f_down <- d1_f_down %>% arrange(desc(logP))
#d1_f_down$Description <- factor(d1_f_down$Description, levels = d1_f_down$Description)
#d1_f_down <- d1_f_down[c(1:10),]
#d1_f_down$group <- "Up-regulaed in Female Cohort"
#p1_f_down <- ggplot(d1_f_down) +
#  geom_bar(stat = "identity",aes(x = Description,y = logP),fill = "#58B6E9") +
#  coord_flip()  +
#  theme_classic(base_size = 14) +
#  theme(axis.text.y = element_text(size = 10,face = "bold"),
#        axis.text.x = element_text(size = 10,face = "bold"),
#        axis.title = element_text(size = 12,face = "bold"),
#        plot.title = element_text(size = 16,hjust = 1,face = "bold")) +
#  labs(y = "-logP.adj", x = "",title = "Down-regulated proteins in Male Cohort")
#ggview(w = 6.5,h = 3)
#ggsave("submit/figure/figure4/male_down_go_9_13_10_17_rank_by_logP.pdf",w = 6.5,h = 3,dpi = 300)

markGenes1 = c(cor1_top_down,cor1_top_up)
#png("result/proteins_changed_with_age/male_cohort_5_cluster.png",w = 8.5,h = 8, units = 'in', res = 600)
pdf("1st submission NC/Figures JZH 20241226/Figure 6/Figure 6B.pdf",w = 8.5,h = 8)
vvisCluster(object = cm1,
           plot.type = "both",
           column_names_rot = 45,
           markGenes = c(cor1_top_down,cor1_top_up),
           genes.gp = c("bold", 8, NA))
dev.off()  






##### female cohort          
d2$Sex <- NULL
rownames(d2) <- d2$eid
d2$eid <- NULL
g2 <- list()
age <- c(39:70)
for(i in seq(age)){
  g2[[i]] <- d2 %>% filter(Age == age[i]) %>% as.tibble()
}
mean2 <- list()
for(i in seq(age)){
  mean2[[i]] <- apply(g2[[i]][,-1],2,function(x){mean(na.omit(x))})
}
names(mean2) <- age
dat2 <- data.frame(mean2[[1]],mean2[[2]],mean2[[3]],mean2[[4]],mean2[[5]],mean2[[6]],mean2[[7]],
                  mean2[[8]],mean2[[9]],mean2[[10]],mean2[[11]],mean2[[12]],mean2[[13]],mean2[[14]],
                  mean2[[15]],mean2[[16]],mean2[[17]],mean2[[18]],mean2[[19]],mean2[[20]],mean2[[21]],
                  mean2[[22]],mean2[[23]],mean2[[24]],mean2[[25]],mean2[[26]],mean2[[27]],mean2[[28]],
                  mean2[[29]],mean2[[30]],mean2[[31]],mean2[[32]])
colnames(dat2) <- age
dat2[,1] <- NULL
dat2 <- na.omit(dat2)

### 5 cluster
cm2 <- clusterData(exp = dat2,
                  cluster.method = "mfuzz",
                  cluster.num = 5,
                  seed = 987)

vvisCluster(object = cm2,
           plot.type = "both",
           column_names_rot = 45)           


res2 <- cm2[["wide.res"]]
write.csv(res2,"result/proteins_changed_with_age/female_cohort.csv")

#### up-regulated protein
res2_up <- res2 %>% filter(cluster == 2)
#write.csv(res2_up,"submit/supplyment/cluster_1_protein_female_cohort.csv")
res2_cor2_up <- merge(cor_data_df2,res2_up,by.x = "symbol",by.y = "gene") %>% filter(FDR < 0.05) %>% arrange(desc(correlation))
cor2_top_up <- res2_cor2_up[c(1:20),1]
#write.csv(res2_up,"result/select_protein/change_with_age/female_cohort_up_protein.csv")
g2_f_up <- mapIds(org.Hs.eg.db,res2_up$gene,'ENTREZID','SYMBOL') 

## kegg
kegg2_f_up <- enrichKEGG(g2_f_up,
                        organism = "hsa",
                        keyType = "kegg",
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        #minGSSize = 10,
                        #maxGSSize = 500,
                        qvalueCutoff = 0.1)
kegg2_f_up@result$GeneRatio <- sapply(kegg2_f_up@result$GeneRatio, function(x) eval(parse(text=x)))
kegg2_f_up@result$GeneRatio <- round(kegg2_f_up@result$GeneRatio,2) 
d2_f_up <- kegg2_f_up@result
write.csv(d2_f_up %>% arrange(desc(Count)),"result/proteins_changed_with_age/kegg_go/kegg analysis of up-regulated proteins in female cohort.csv")
d2_f_up$logP <- -log10(d2_f_up$p.adjust)
d2_f_up <- d2_f_up %>% arrange(desc(Count))
d2_f_up$Description <- factor(d2_f_up$Description, levels = d2_f_up$Description)
d2_f_up <- d2_f_up[c(1:10),]
d2_f_up$group <- "Up-regulaed in Female Cohort"
p2_f_up <- ggplot(d2_f_up,aes(x = Count,y = Description,size = logP)) +
                  geom_point(color = "#E3625D") +
                  theme_classic(base_size = 14) +
                  theme(axis.text.y = element_text(size = 10,face = "bold"),
                        axis.text.x = element_text(size = 10,face = "bold"),
                        axis.title = element_text(size = 12,face = "bold"),
                        legend.title = element_text(face = "bold",size = 10),
                        plot.title = element_text(size = 16,face = "bold"),
                        legend.background = element_rect(linetype="solid",colour ="black"),
                        legend.key.size = unit(0.01, units = "cm"),
                        legend.position = c(0.7, 0.7)) +
                  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
                  #scale_x_continuous(limits = c(20,110),breaks = seq(20,110,10)) +
                  labs(y = "", x = "Count",
                       title = "KEGG analysis", 
                       size = "-logP.adj")
ggview(w = 4.7,h = 4.5)
ggsave("1st submission NC/Figures JZH 20241226/Figure 6/Figure 6C.pdf",w = 4.7,h = 4.5,dpi = 300)

## go 
go2_f_up <- enrichGO(gene       = g2_f_up,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENTREZID',
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
go2_f_up@result$GeneRatio <- sapply(go2_f_up@result$GeneRatio, function(x) eval(parse(text=x)))
go2_f_up@result$GeneRatio <- round(go2_f_up@result$GeneRatio,2) 
d2_f_up <- go2_f_up@result
write.csv(d2_f_up %>% arrange(desc(Count)),"result/proteins_changed_with_age/kegg_go/go analysis of up-regulated proteins in female cohort.csv")
d2_f_up$logP <- -log10(d2_f_up$p.adjust)
#write.csv(d2_f_up,"result/select_protein/change_with_age/go_female_cohort_up_protein.csv")
### rank by Count
d2_f_up <- d2_f_up %>% arrange(desc(Count))
d2_f_up$Description <- factor(d2_f_up$Description, levels = d2_f_up$Description)
d2_f_up <- d2_f_up[c(1:10),]
d2_f_up$group <- "Up-regulaed in Female Cohort"
p2_f_up <- ggplot(d2_f_up,aes(x = Count,y = Description,size = logP)) +
                  geom_point(color = "#E3625D") +
                  theme_classic(base_size = 14) +
                  theme(axis.text.y = element_text(size = 10,face = "bold"),
                        axis.text.x = element_text(size = 10,face = "bold"),
                        axis.title = element_text(size = 12,face = "bold"),
                        legend.title = element_text(face = "bold",size = 10),
                        plot.title = element_text(size = 16,face = "bold"),
                        legend.background = element_rect(linetype="solid",colour ="black"),
                        legend.key.size = unit(0.01, units = "cm"),
                        legend.position = c(0.7, 0.7)) +
                  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
                  #scale_x_continuous(limits = c(20,110),breaks = seq(20,110,10)) +
                  labs(y = "", x = "Count",
                       title = "GO analysis", 
                       size = "-logP.adj")
ggview(w = 5.1,h = 4.5)
ggsave("1st submission NC/Supplemental figures JZH 20241226/supplemental figure 2 JZH 20241226/supplemental figure 2a JZH 20241226.pdf",w = 5.1,h = 4.5,dpi = 300)


### down-regulated protein
res2_down <- res2 %>% filter(cluster == 5)
#write.csv(res2_down,"submit/supplyment/cluster_5_protein_female_cohort.csv")
res2_cor2_down <- merge(cor_data_df2,res2_down,by.x = "symbol",by.y = "gene") %>% filter(FDR < 0.05) %>% arrange(correlation)
cor2_top_down <- res2_cor2_down[c(1:20),1]
#write.csv(res2_down,"result/select_protein/change_with_age/female_cohort_down_protein.csv")
g2_f_down <- mapIds(org.Hs.eg.db,res2_down$gene,'ENTREZID','SYMBOL') 

## kegg
kegg2_f_down <- enrichKEGG(g2_f_down,
                        organism = "hsa",
                        keyType = "kegg",
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        #minGSSize = 10,
                        #maxGSSize = 500,
                        qvalueCutoff = 0.1)
kegg2_f_down@result$GeneRatio <- sapply(kegg2_f_down@result$GeneRatio, function(x) eval(parse(text=x)))
kegg2_f_down@result$GeneRatio <- round(kegg2_f_down@result$GeneRatio,2) 
d2_f_down <- kegg2_f_down@result
write.csv(d2_f_down %>% arrange(desc(Count)),"result/proteins_changed_with_age/kegg_go/kegg analysis of down-regulated proteins in female cohort.csv")
d2_f_down$logP <- -log10(d2_f_down$p.adjust)
#write.csv(d2_f_down,"result/select_protein/change_with_age/kegg_female_cohort_down_protein.csv")
d2_f_down <- d2_f_down %>% arrange(desc(Count))
d2_f_down$Description <- factor(d2_f_down$Description, levels = d2_f_down$Description)
d2_f_down <- d2_f_down[c(1:10),]
d2_f_down$group <- "Down-regulaed in Female Cohort"
p2_f_down <- ggplot(d2_f_down,aes(x = Count,y = Description,size = logP)) +
                  geom_point(color = "#3E90BF") +
                  theme_classic(base_size = 14) +
                  theme(axis.text.y = element_text(size = 10,face = "bold"),
                        axis.text.x = element_text(size = 10,face = "bold"),
                        axis.title = element_text(size = 12,face = "bold"),
                        legend.title = element_text(face = "bold",size = 10),
                        plot.title = element_text(size = 16,face = "bold"),
                        legend.background = element_rect(linetype="solid",colour ="black"),
                        legend.key.size = unit(0.01, units = "cm"),
                        legend.position = c(0.7, 0.7)) +
                  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
                  #scale_x_continuous(limits = c(20,110),breaks = seq(20,110,10)) +
                  labs(y = "", x = "Count",
                       title = "KEGG analysis", 
                       size = "-logP.adj")
ggview(w = 4.7,h = 4)
ggsave("1st submission NC/Figures JZH 20241226/Figure 6/Figure 6D.pdf",w = 4.7,h = 4.5,dpi = 300)


## go
go2_f_down <- enrichGO(gene       = g2_f_down,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENTREZID',
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)
go2_f_down@result$GeneRatio <- sapply(go2_f_down@result$GeneRatio, function(x) eval(parse(text=x)))
go2_f_down@result$GeneRatio <- round(go2_f_down@result$GeneRatio,2) 
d2_f_down <- go2_f_down@result
write.csv(d2_f_down %>% arrange(desc(Count)),"result/proteins_changed_with_age/kegg_go/go analysis of down-regulated proteins in female cohort.csv")
d2_f_down$logP <- -log10(d2_f_down$p.adjust)
#write.csv(d2_f_down,"result/select_protein/change_with_age/go_female_cohort_down_protein.csv")
### rank by count
d2_f_down <- d2_f_down %>% arrange(desc(Count))
d2_f_down$Description <- factor(d2_f_down$Description, levels = d2_f_down$Description)
d2_f_down <- d2_f_down[c(1:10),]
d2_f_down$group <- "Down-regulaed in Female Cohort"

p2_f_down <- ggplot(d2_f_down,aes(x = Count,y = Description,size = logP)) +
                  geom_point(color = "#3E90BF") +
                  theme_classic(base_size = 14) +
                  theme(axis.text.y = element_text(size = 10,face = "bold"),
                        axis.text.x = element_text(size = 10,face = "bold"),
                        axis.title = element_text(size = 12,face = "bold"),
                        legend.title = element_text(face = "bold",size = 10),
                        plot.title = element_text(size = 16,face = "bold"),
                        legend.background = element_rect(linetype="solid",colour ="black"),
                        legend.key.size = unit(0.01, units = "cm"),
                        legend.position = c(0.8, 0.5)) +
                  scale_y_discrete(labels = function(x) str_wrap(x, width = 35)) +
                  #scale_x_continuous(limits = c(20,110),breaks = seq(20,110,10)) +
                  labs(y = "", x = "Count",
                       title = "GO analysis", 
                       size = "-logP.adj")
ggview(w = 5.1,h = 4.5)
ggsave("1st submission NC/Supplemental figures JZH 20241226/supplemental figure 2 JZH 20241226/supplemental figure 2b JZH 20241226.pdf",w = 5.1,h = 4.5,dpi = 300)

### rank by logP
#d2_f_down <- d2_f_down %>% arrange(desc(logP))
#d2_f_down$Description <- factor(d2_f_down$Description, levels = d2_f_down$Description)
#d2_f_down <- d2_f_down[c(1:10),]
#d2_f_down$group <- "Up-regulaed in Female Cohort"
#g2_f_down <- ggplot(d2_f_down) +
#  geom_bar(stat = "identity",aes(x = Description,y = logP),fill = "#3E90BF") +
#  coord_flip()  +
#  theme_classic(base_size = 14) +
#  theme(axis.text.y = element_text(size = 10,face = "bold"),
#        axis.text.x = element_text(size = 10,face = "bold"),
#        axis.title = element_text(size = 12,face = "bold"),
#        plot.title = element_text(size = 16,hjust = 1,face = "bold")) +
#  labs(y = "-logP.adj", x = "",title = "Down-regulated proteins in Female Cohort")
#ggview(w = 6.5,h = 2.5)
#ggsave("submit/figure/figure4/female_down_go_09_13_10_20_rank_by_logP.pdf",w = 6.5,h = 3,dpi = 300)

markGenes2 = c(cor2_top_down,cor2_top_up)
#png("result/proteins_changed_with_age/female_cohort_5_cluster.png",w = 8.5,h = 8, units = 'in', res = 600)
pdf("1st submission NC/Figures JZH 20241226/Figure 6/Figure 6A.pdf",w = 8.5,h = 8)
vvisCluster(object = cm2,
           plot.type = "both",
           column_names_rot = 45,
           markGenes = c(cor2_top_down,cor2_top_up),
           genes.gp = c("bold", 8, NA))
dev.off()  


save(cm1,cm2,dat1,dat2,file = "result/proteins_changed_with_age/cluster_result.RData")

########################## select specific proteins 
rm(list=ls())
load("result/proteins_changed_with_age/cluster_result.RData")

######## female cohort
cm_f <- cm2
cm_m <- cm1
dat_f <- dat2
dat_m <- dat1
rm(cm1)
rm(cm2)
rm(dat1)
rm(dat2)

col_color = c("#3E4F94", "white", "#992224")
col_range = c(-2, 0, 2)
col_fun = circlize::colorRamp2(col_range, col_color)

###### female cohort
res_f <- cm_f[["wide.res"]]
rownames(res_f) <- res_f$gene
res_f$gene <- NULL

res_f_up <- filter(res_f,cluster == 3)[,-c(32,33)]
#png("result/proteins_changed_with_age/cluster_3_protein_female_cohort.png",w = 6,h = 18, units = 'in', res = 600)
pdf("result/proteins_changed_with_age/cluster_3_protein_female_cohort.pdf",w = 6,h = 18)
ComplexHeatmap::Heatmap(as.matrix(res_f_up),
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        name = "Z-score",
                        col = col_fun,
                        column_title = "Up-regulated proteins with age in female cohort",
                        row_names_gp = gpar(fontsize = 4),
                        column_names_gp = gpar(fontsize = 8))
dev.off()

res_f_down <- filter(res_f,cluster == 5)[,-c(32,33)]
#png("result/proteins_changed_with_age/cluster_5_protein_female_cohort.png",w = 6,h = 20, units = 'in', res = 600)
pdf("result/proteins_changed_with_age/cluster_5_protein_female_cohort.pdf",w = 6,h = 20)
ComplexHeatmap::Heatmap(as.matrix(res_f_down),
                        cluster_rows = T,
                        cluster_columns = FALSE,
                        show_row_dend = F,
                        name = "Z-score",
                        col = col_fun,
                        column_title = "Down-regulated proteins with age in female cohort",
                        row_names_gp = gpar(fontsize = 4),
                        column_names_gp = gpar(fontsize = 8))
dev.off()

###### male cohort
res_m <- cm_m[["wide.res"]]
rownames(res_m) <- res_m$gene
res_m$gene <- NULL

res_m_up <- filter(res_m,cluster == 2)[,-c(31,32)]
#png("result/proteins_changed_with_age/cluster_4_protein_male_cohort.png",w = 6,h = 25, units = 'in', res = 600)
pdf("result/proteins_changed_with_age/cluster_4_protein_male_cohort.pdf",w = 6,h = 25)
ComplexHeatmap::Heatmap(as.matrix(res_m_up),
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        name = "Z-score",
                        col = col_fun,
                        column_title = "Up-regulated proteins with age in male cohort",
                        row_names_gp = gpar(fontsize = 2),
                        column_names_gp = gpar(fontsize = 8))
dev.off()

res_m_down <- filter(res_m,cluster == 5)[,-c(31,32)]
png("result/proteins_changed_with_age/cluster_1_protein_male_cohort.png",w = 6,h = 20, units = 'in', res = 600)
#pdf("result/proteins_changed_with_age/cluster_1_protein_male_cohort.pdf",w = 6,h = 20)
ComplexHeatmap::Heatmap(as.matrix(res_m_down),
                        cluster_rows = T,
                        cluster_columns = FALSE,
                        show_row_dend = F,
                        name = "Z-score",
                        col = col_fun,
                        column_title = "Down-regulated proteins with age in male cohort",
                        row_names_gp = gpar(fontsize = 3),
                        column_names_gp = gpar(fontsize = 8))
dev.off()















