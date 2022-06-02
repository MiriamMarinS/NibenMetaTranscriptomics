# Circos y heatmap.

#Data with quintiles information.
data <- read.table("C:/Users/milil/Downloads/Nextcloud/Documents/Miriam/Javi/Expression/Score/raw_data/data_quintile_round_inNiben261.txt", sep="\t", header=TRUE, row.names = 1)

library(circlize)

dataChrscytobands <- read.table("C:/Users/milil/Downloads/Nextcloud/Documents/Miriam/Javi/Expression/Score/raw_data/cytobands_chr.txt", sep="\t", header=TRUE)
dataChrscytobands <- subset(dataChrscytobands, select=-gieStain)
circos.initializeWithIdeogram(dataChrscytobands)

# bed: chr, start, end, value.
z <- read.table("C:/Users/milil/Downloads/Nextcloud/Documents/Miriam/Javi/Expression/Score/results/findIR_complete_v2.txt", sep="\t", header=TRUE, row.names = 1)
z_new <- z[,c("Chr", "Start", "End", "Suitability", "Mean.Length.IR", "requisite", "Window", "Interval", "N.genes.assigned.to.expression.category", "N.genes.total", "Mean...GC", "Mean...N", "N.IR")]
z_new <- subset(z_new, !is.na(z_new[,"Mean.Length.IR"]))
quantile <- quantile(z_new$Suitability, 0.75)[[1]]
z_new <- subset(z_new, Suitability > quantile)
col <- ifelse(z_new$requisite == 1, "black", "lightseagreen")
z_heatmap <- subset(z_new, select=-c(Mean.Length.IR, requisite, Window, Interval, N.genes.assigned.to.expression.category, N.genes.total, Mean...GC, Mean...N, N.IR))
max <- max(z_heatmap$Suitability)
min <- min(z_heatmap$Suitability)
mean <- mean(z_heatmap$Suitability)
col_fun = colorRamp2(c(min, mean, max), c("#CDF0FF", "#37C6FF", "#0072A0"))
circos.genomicHeatmap(z_heatmap, col = col_fun, side = "inside", border = "black", heatmap_height = 0.05, line_col = col)

z_rainfall <- z_new
z_rainfall <- subset(z_rainfall, select=-c(Suitability, requisite, Window, Interval, N.genes.assigned.to.expression.category, N.genes.total, Mean...GC, Mean...N, N.IR))
circos.genomicRainfall(z_rainfall, pch = 16, cex = 0.45, col = "lightseagreen", track.height = 0.06, normalize_to_width = FALSE, ylim = NULL)

circos.genomicDensity(subset(data, quintile == 5), col = c("#FF000095"), track.height = 0.06)
circos.genomicDensity(subset(data, quintile == 4), col = c("#FF000045"), track.height = 0.06)
circos.genomicDensity(subset(data, quintile == 3), col = c("#FF000030"), track.height = 0.06)
circos.genomicDensity(subset(data, quintile == 2), col = c("#0000FF45"), track.height = 0.06)
circos.genomicDensity(subset(data, quintile == 1), col = c("#0000FF95"), track.height = 0.06)
circos.genomicDensity(subset(data, quintile == 0), col = c("grey"), track.height = 0.06)

library(graphics)
lgd_expr = legend(x = 1, y = 1, title = "Suitability", c(min, round(mean, 2), max), fill = c("#CDF0FF", "#37C6FF", "#0072A0"), cex = 0.8)
#lgd_gene_type = legend(x = 1, y = -0.6, title = "Expression category", c("VHE", "HE", "ME", "LE", "VLE", "NE"), 
#                       fill = c("#FF000095", "#FF000045", "#FF000030", "#0000FF45", "#0000FF95", "grey"), cex = 1)

