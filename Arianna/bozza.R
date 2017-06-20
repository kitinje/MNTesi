y$samples$cell_line = c( rep("HCC1599",12), rep("MB157",12) )
y$samples$treatment = factor( rep(c(rep("A",3), rep("AD", 3), rep("C", 3), rep("D", 3)  ),2) )

pcaY <- prcomp(t(y.cpm.log2), scale. = TRUE)

pcaDf <- data.frame(pcaY$x, "cell_line" = y$samples$cell_line, "treatment" = y$samples$treatment)
plot(pcaY, type = "l") # Prinicipal components explain variance:

# 2D STATIC PCA-PLOT
plot.pca_2d <- ggplot(pcaDf,aes(x=PC1, y=PC2)) +
  theme_classic() +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_vline(xintercept = 0, color = "gray70") +
  geom_point(aes(color = cell_line, shape=treatment), alpha = 0.55, size = 3) +
  xlab("PC1") +
  ylab("PC2") + 
  ggtitle("PCA from cell-line data")
plot.pca_2d_numbered = plot.pca_2d + geom_text_repel(aes(y = PC2 + 0.25, label = c(seq(1:dim(y)[2]) )))
plot.pca_2d_numbered

# 2D INTERACTIVE PCA-PLOT
plot.pca_2d_int <- plot_ly(pcaDf, x = pcaDf[,1] , y = pcaDf[,2],  text = pcaDf$group, mode = "markers", color = pcaDf$cell_line, symbol = y$samples$treatment,  marker = list(size = 11), type = "scatter" ) 
plot.pca_2d_int <- layout(plot.pca_2d_int,title = "Principal Component Analysis", 
                          xaxis = list(title = "PC 1"),
                          yaxis = list(title = "PC 2"))
plot.pca_2d_int

# 3D INTERACTIVE PCA-PLOT
plot.pca_3d_int <- plot_ly(pcaDf, x = pcaDf[,1] , y = pcaDf[,2], z = pcaDf[,3],  text = pcaDf$group,
                           mode = "markers", color = pcaDf$cell_line, marker = list(size = 7), type="scatter3d", symbol = y$samples$treatment) 
plot.pca_3d_int <- layout(plot.pca_3d_int,  title = "Principal Component Analysis")
plot.pca_3d_int


mb157 <- y[,grep(x = colnames(y),  pattern = "MB157")]

mb157 <- mb157[,grep(x = colnames(mb157),  pattern = "AD", invert = TRUE)]

#create design matrix
# A <- factor(c(1,1,1,1,1,1,0,0,0,0,0,0))
# D <- factor(c(0,0,0,1,1,1,0,0,0,1,1,1))
# design = model.matrix(~ 0 + A + D)

# design[10:12,1] = 0 
# colnames(design)=c("C","A", "D")

A <- factor(c(1,1,1,0,0,0,0,0,0))
D <- factor(c(0,0,0,0,0,0,1,1,1))
design = model.matrix(~ 0 + A + D)
design[7:9,1] = 0 
colnames(design)=c("C","A", "D")

mb157 <- estimateDisp(mb157, design, robust=TRUE)
plotBCV(mb157) #plot NB dispersion
fit_int <- glmQLFit(mb157, design, robust=TRUE)

AvsC=makeContrasts(A-C, levels=design)
DvsC=makeContrasts(D-C, levels=design)

mb157.qlf.DvsC <- glmQLFTest(fit_int,  contrast = DvsC)
topTags(mb157.qlf.DvsC)
