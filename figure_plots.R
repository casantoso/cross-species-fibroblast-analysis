emb <- Embeddings(zf_merged, "umap")
xr  <- range(emb[,1])
yr  <- range(emb[,2])
pad_x <- 0.03 * diff(xr)
pad_y <- 0.03 * diff(yr)

lims_x <- c(xr[1]-pad_x, xr[2]+pad_x)
lims_y <- c(yr[1]-pad_y, yr[2]+pad_y)

human_filtered
human_markers

mouse_filtered
mouse_markers

macaque_integrated
macaque_markers

macaque_colon_SCT
macaque_colon_markers

chicken_merged
chicken_markers

pig_filtered
pig_markers

frog_merged
frog_markers

zf_merged
zf_markers


# human 
#   frac --> 0.06
#   ECM --> 0.15
#   
# Mouse
#   frac --> 0.08
#   ECM --> 0.2
# 
# Macaque_GSE196792
#   frac --> 0.08
#   ECM --> 0.2 
# 
# Macaque_NHPCA (colon)
#   frac --> 0.08
#   ECM --> 0.2
# 
# Chicken
#   frac --> 0.2
#   ECM --> 0.3
# 
# Pig_GSE293713
#   frac --> 0.12
#   ECM --> 0.2
# 
# Frog
#   frac --> 0.12
#   ECM --> 0.12
# 
# Zebrafish
#   frac --> 0.05
#   ECM --> 0.1



#fraction umap w/o text
Frac_img <- FeaturePlot(zf_merged, features = "ECM_fraction",  reduction = 'umap')+
  scale_colour_gradientn(
    colours = c("blue", "green", "yellow", "red"),
    limits  = c(0, 0.05)) +
  ggtitle(NULL) + 
  theme(
    axis.title  = element_blank(),
    axis.text   = element_blank(), 
    axis.ticks  = element_blank()
  ) + 
  guides(
    color = guide_colorbar(
      title = NULL,       # no title
      label = FALSE,      # remove text labels
      ticks = FALSE       # hide tick marks (optional)
      # barheight = unit(60, "pt"), barwidth = unit(8, "pt")  # size, if you want
    )) +
  coord_cartesian(xlim = lims_x, ylim = lims_y) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))


#fraction umap with text
Frac_txt_img <- FeaturePlot(zf_merged, features = "ECM_fraction",  reduction = 'umap')+
  scale_colour_gradientn(
    colours = c("blue", "green", "yellow", "red"),
    limits  = c(0, 0.05)) +
  ggtitle(NULL) + 
  theme(
    axis.title  = element_blank(),
    axis.text   = element_blank(), 
    axis.ticks  = element_blank()
  ) + 
  guides(
    color = guide_colorbar(
      title = NULL,       # no title
      ticks = FALSE       # hide tick marks (optional)
    )) +
  coord_cartesian(xlim = lims_x, ylim = lims_y) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))


#ucell umap w/o text
UCell_img <- FeaturePlot(zf_merged, features = "ECM_UCell",  reduction = 'umap')+
  scale_colour_gradientn(
    colours = c("blue", "green", "yellow", "red"),
    limits  = c(0, 0.1)) +
  ggtitle(NULL) + 
  theme(
    axis.title  = element_blank(),
    axis.text   = element_blank(), 
    axis.ticks  = element_blank()
  ) + 
  guides(
    color = guide_colorbar(
      title = NULL,       
      label = FALSE,     
      ticks = FALSE     
    )) +
  coord_cartesian(xlim = lims_x, ylim = lims_y) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))

#ucell umap with text
UCell_txt_img <- FeaturePlot(zf_merged, features = "ECM_UCell",  reduction = 'umap')+
  scale_colour_gradientn(
    colours = c("blue", "green", "yellow", "red"),
    limits  = c(0, 0.1)) +
  ggtitle(NULL) + 
  theme(
    axis.title  = element_blank(),
    axis.text   = element_blank(), 
    axis.ticks  = element_blank()
  ) + 
  guides(
    color = guide_colorbar(
      title = NULL,       
      ticks = FALSE     
    )) +
  coord_cartesian(xlim = lims_x, ylim = lims_y) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))


#cell type umap
cellType_img <- DimPlot(zf_merged, group.by = "cell_type", cols = celltype_colors) +
  ggtitle(NULL) + 
  theme(
    axis.title  = element_blank(),
    axis.text   = element_blank(), 
    axis.ticks  = element_blank()
  ) +                     
  guides(
    colour = guide_legend(
      title = NULL,                  
      label = FALSE,                 
      override.aes = list(size = 5)   
    )) +
  coord_cartesian(xlim = lims_x, ylim = lims_y) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))


ggsave("/Users/chiarasantoso/Desktop/purm/scRNAseq/figures/zebrafish/ECM_Frac.jpeg", Frac_img, width = 6, height = 5, units = "in", dpi = 100, bg = "white")
ggsave("/Users/chiarasantoso/Desktop/purm/scRNAseq/figures/zebrafish/ECM_UCell.jpeg", UCell_img, width = 6, height = 5, units = "in", dpi = 100, bg = "white")
ggsave("/Users/chiarasantoso/Desktop/purm/scRNAseq/figures/zebrafish/ECM_Frac_wtext.jpeg", Frac_txt_img, width = 6, height = 5, units = "in", dpi = 100, bg = "white")
ggsave("/Users/chiarasantoso/Desktop/purm/scRNAseq/figures/zebrafish/ECM_UCell_wtext.jpeg", UCell_txt_img, width = 6, height = 5, units = "in", dpi = 100, bg = "white")
ggsave("/Users/chiarasantoso/Desktop/purm/scRNAseq/figures/zebrafish/celltype.jpeg", cellType_img, width = 6, height = 5, units = "in", dpi = 100, bg = "white")

celltype_colors <- c(
  "Fibroblasts" = "#F8766D",
  "SMCs" = "#00BE67",
  "Pericytes/SMCs" = "#00BE67",
  "Glial cells" = "#ABA300",
  "Plasma cells" = "#00B8E7",
  "Endothelial" = "#C77CFF",
  "Pericytes" = "#8494FF",
  "B cells" = "#E68613",
  "T cells" = "#FF61CC",
  "Epithelial" = "#7570B3",
  "Goblet cells" = "#1B9E77",
  "Macrophages" = "#7CAE00",
  "ICCs" = "#00A9FF",
  "Neurons" = "#CD9600", 
  "Immune cells" = "#FF68A1",
  "Erythrocytes" = "#E6AB02",
  "Astrocytes" = "#D95F02", 
  "Enterocytes" = "#00C19A", 
  "Monocytes" = "#00BFC4", 
  "Schwann cells" = "#A6761D",
  "Mast cells" = "#666666",
  "Leukocytes" = "#FFD300", 
  "Acinar cells" = "#7CAE00", 
  "Epidermis" = "#00B8E7", 
  "Hepatocytes" = "#D95F02",
  "Muscle cells" = "#A6761D", 
  "Enteroendocrine" = "#00C19A"
)

"Myeloid cells" = "#B39DDB"
"Enteroendocrine" = "#00C19A" 
"Dendritic cells" = "#FFD300"

Fibroblast_markers_CAPS <- c("BMP4" , "BMP5", "CCL19", "CCL8", "COL14A1" , "COL1A1","COL3A1","COL4A5","COL4A6", "CPM", "CXCL13", "C3", 
                        "DCN", "FGFR2", "GREM1", "LUM", "MMP1", "OGN", "PDGFRA", "SFRP2", "SOX6", "SPARC", "VCAM1", "DPT", "MGP", "LPAR1")
Fibroblast_markers <- c("Bmp4" , "Bmp5", "Ccl19", "Ccl8", "Col14a1" , "Col1a1","Col3a1","Col4a5","Col4a6", "Cpm", "Cxcl13", "C3", 
                       "Dcn", "Fgfr2", "Grem1", "Lum", "Mmp1", "Ogn", "Pdgfra", "Sfp2", "Sox6", "Sparc", "Vcam1", "Dpt", "Mgp", "Lpar1")
Fibroblast_markers_frog <- c("bmp4.L", "bmp4.S", "bmp5.L", "bmp5.S",  "col14a1.L", "col14a1.S", "col1a1.L" ,
                       "col1a1.S","col3a1.L", "col3a1.S" ,"col4a5.L", "col4a5.S","col4a6.L ","col4a6.S",
                       "cpm.L", "cpm.S", "cxcl13.L","cxcl13.S",  "c3.L", "c3.S", "dcn.L","dcn.S" ,
                       "fgfr2.L", "fgfr2.S", "grem1.L", "lum.L","lum.S",  "mmp1.L","mmp1.S",  "ogn.L", "ogn.S",
                       "pdgfra.L", "pdgfra.S","sfrp2.L", "sfrp2.S",  "sox6.L", "sox6.S",  "sparc.L", "sparc.S",  
                       "vcam1.L","vcam1.S",  "dpt.L","dpt.S" ,"mgp.L","mgp.S",  "lpar1.L", "lpar1.S")
Fibroblast_markers_zf <- c("bmp4" , "bmp5", "col14a1a", "col1a1a", "cpm" , "dcn","fgfr2","grem1a","lum", "mmp13a", "ogna", "ognb", 
                       "pdgfra", "sfrp2", "sfrp2", "sparc", "dpt", "mgp", "c3", "lpar1")


additional
chicken <- c("COL1A2")
zebrafish <- c("pmp22a", "aqp1a.1", "podxl", "cavin2b", "cavin1b", "rhag")



#fibroblast plots 
for (marker in zebrafish) {
  if (marker %in% zf_markers$gene) {
    p <- FeaturePlot(zf_merged, features = marker, order = TRUE) +
      scale_colour_gradientn(colours = c("blue", "green", "yellow", "red")) +
      ggtitle(NULL) + 
      theme(
        axis.title  = element_blank(),
        axis.text   = element_blank(), 
        axis.ticks  = element_blank()
      ) + 
      guides(
        color = guide_colorbar(
          title = NULL,       
          label = FALSE,     
          ticks = FALSE     
        )) +
      coord_cartesian(xlim = lims_x, ylim = lims_y) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) 
    ggsave(filename = paste0("/Users/chiarasantoso/Desktop/purm/scRNAseq/figures/zebrafish/feature_plots/", marker, "_FeaturePlot.jpeg"),
           plot = p, width = 6, height = 5, dpi = 100)
  } else {
    message(paste("Gene", marker, "not found in the dataset. Skipping."))
  }
}

Idents(zf_merged) <- "seurat_clusters"

for (marker in zebrafish) {
  if (marker %in% zf_markers$gene) {
    p <- VlnPlot(zf_merged, features = marker, assay = "RNA", layer = "data") +
      ggtitle(NULL) +
      guides(fill = guide_legend(title = NULL, label = FALSE)) +  
      scale_x_discrete(breaks = NULL, labels = NULL) +
      theme(
        axis.title  = element_blank(),
        axis.text   = element_blank(),
        axis.ticks  = element_blank()
      )
    ggsave(filename = paste0("/Users/chiarasantoso/Desktop/purm/scRNAseq/figures/zebrafish/VlnPlots/", marker, "_VlnPlot.jpeg"),
           plot = p, width = 6, height = 5, units = "in", dpi = 100)
  } else {
    message(paste("Gene", marker, "not found in the dataset. Skipping."))
  }
}


