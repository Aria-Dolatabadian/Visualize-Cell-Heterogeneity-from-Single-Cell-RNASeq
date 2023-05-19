#https://jokergoo.github.io/ComplexHeatmap-reference/book/more-examples.html#the-measles-vaccine-heatmap



expr = read.table("mouse_scRNAseq_corrected.txt", sep = "\t", header = TRUE)
expr = expr[!duplicated(expr[[1]]), ]
rownames(expr) = expr[[1]]
expr = expr[-1]
expr = as.matrix(expr)


expr = expr[apply(expr, 1, function(x) sum(x > 0)/length(x) > 0.5), , drop = FALSE]

get_correlated_variable_genes = function(mat, n = nrow(mat), cor_cutoff = 0, n_cutoff = 0) {
    ind = order(apply(mat, 1, function(x) {
            q = quantile(x, c(0.1, 0.9))
            x = x[x > q[1] & x < q[2]]
            var(x)/mean(x)
        }), decreasing = TRUE)[1:n]
    mat2 = mat[ind, , drop = FALSE]
    dt = cor(t(mat2), method = "spearman")
    diag(dt) = 0
    dt[abs(dt) < cor_cutoff] = 0
    dt[dt < 0] = -1
    dt[dt > 0] = 1

    i = colSums(abs(dt)) > n_cutoff

    mat3 = mat2[i, ,drop = FALSE]
    return(mat3)
}

mat = get_correlated_variable_genes(expr, cor_cutoff = 0.5, n_cutoff = 20)
mat2 = t(apply(mat, 1, function(x) {
    q10 = quantile(x, 0.1)
    q90 = quantile(x, 0.9)
    x[x < q10] = q10
    x[x > q90] = q90
    scale(x)
}))
colnames(mat2) = colnames(mat)

cc = readRDS("mouse_cell_cycle_gene.rds")
ccl = rownames(mat) %in% cc
cc_gene = rownames(mat)[ccl]

rp = readRDS("mouse_ribonucleoprotein.rds")
rpl = rownames(mat) %in% rp

base_mean = rowMeans(mat)

library(GetoptLong)
ht_list = Heatmap(mat2, col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")), 
    name = "scaled_expr", column_title = qq("relative expression for @{nrow(mat)} genes"),
    show_column_names = FALSE, width = unit(8, "cm"),
    heatmap_legend_param = list(title = "Scaled expr")) +
    Heatmap(base_mean, name = "base_expr", width = unit(5, "mm"),
        heatmap_legend_param = list(title = "Base expr")) +
    Heatmap(rpl + 0, name = "ribonucleoprotein", col = c("0" = "white", "1" = "purple"), 
        show_heatmap_legend = FALSE, width = unit(5, "mm")) +
    Heatmap(ccl + 0, name = "cell_cycle", col = c("0" = "white", "1" = "red"), 
        show_heatmap_legend = FALSE, width = unit(5, "mm")) +
    rowAnnotation(link = anno_mark(at = which(ccl & base_mean > quantile(base_mean, 0.25)), 
        labels = rownames(mat)[ccl & base_mean > quantile(base_mean, 0.25)], 
        labels_gp = gpar(fontsize = 10), padding = unit(1, "mm"))) +
    Heatmap(cor(t(mat2)), name = "cor", 
        col = colorRamp2(c(-1, 0, 1), c("green", "white", "red")), 
        show_row_names = FALSE, show_column_names = FALSE, row_dend_side = "right", 
        show_column_dend = FALSE, column_title = "pairwise correlation between genes",
        heatmap_legend_param = list(title = "Correlation"))
ht_list = draw(ht_list, main_heatmap = "cor")
decorate_column_dend("scaled_expr", {
    tree = column_dend(ht_list)$scaled_expr
    ind = cutree(as.hclust(tree), k = 2)[order.dendrogram(tree)]

    first_index = function(l) which(l)[1]
    last_index = function(l) { x = which(l); x[length(x)] }
    x1 = c(first_index(ind == 1), first_index(ind == 2)) - 1
    x2 = c(last_index(ind == 1), last_index(ind == 2))
    grid.rect(x = x1/length(ind), width = (x2 - x1)/length(ind), just = "left",
        default.units = "npc", gp = gpar(fill = c("#FF000040", "#00FF0040"), col = NA))
})


