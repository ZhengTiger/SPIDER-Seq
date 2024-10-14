## Figure 6E

```{r}
seu <- Adult.IT.PT.barcode
seu <- SetupForWGCNA(
  seu,
  features = all_gene,
  wgcna_name = "Adult.IT.PT.barcode"
)

seu <- MetacellsByGroups(
  seurat_obj = seu,
  group.by = c("SubType"),
  reduction = 'pca',
  k = 25,
  max_shared = 10,
  ident.group = 'SubType'
)
seu <- NormalizeMetacells(seu)

seu <- SetDatExpr(
  seu,
  group_name = unique(seu$SubType),
  group.by='SubType',
  assay = 'RNA',
  slot = 'data'
)

seu <- TestSoftPowers(
  seu,
  networkType = 'signed'
)
plot_list <- PlotSoftPowers(seu)
wrap_plots(plot_list, ncol=2)

seu <- ConstructNetwork(
  seu,
  soft_power=10,
  tom_name = 'L_R',
  minModuleSize = 10
)

seu <- ScaleData(seu, features=VariableFeatures(seu))
seu <- ModuleEigengenes(seu)
seu <- ModuleConnectivity(
  seu,
  group_name = unique(seu$SubType),
  group.by='SubType'
)
seu <- ResetModuleNames(
  seu,
  new_name = "M"
)
```


```{r}
modules <- GetModules(seu)
mods <- levels(modules$module)
mods <- mods[mods != 'grey']

HubGeneNetworkPlot(
  seu,
  n_hubs = 5, n_other=100,
  edge_prop = 0.75,
  mods = mods
)
```



```{r fig.width=4, fig.height=4}
set.seed(20240703)
seu <- RunModuleUMAP(
  seu,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=25, # neighbors parameter for UMAP
  min_dist=0.5
)

ModuleUMAPPlot(
  seu,
  edge.alpha=0.5,
  sample_edges=TRUE,
  edge_prop=0.5, # proportion of edges to sample (20% here)
  label_hubs=3 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE,
  vertex.label.cex = 0.01
)
```


```{r}
save.image(file="../data/csv/transmitter_and_receptor/hdwgcna.RData")
```






## AGs CAMs

```{r}
seu.hdwgcna.AGMs_CAMs <- Adult.IT.PT.barcode
seu.hdwgcna.AGMs_CAMs <- SetupForWGCNA(
  seu.hdwgcna.AGMs_CAMs,
  features = c(all_gene),
  wgcna_name = "AGMs_CAMs"
)

seu.hdwgcna.AGMs_CAMs <- MetacellsByGroups(
  seurat_obj = seu.hdwgcna.AGMs_CAMs,
  group.by = c("SubType"),
  reduction = 'pca',
  k = 25,
  max_shared = 10,
  ident.group = 'SubType'
)
seu.hdwgcna.AGMs_CAMs <- NormalizeMetacells(seu.hdwgcna.AGMs_CAMs)

seu.hdwgcna.AGMs_CAMs <- SetDatExpr(
  seu.hdwgcna.AGMs_CAMs,
  group_name = unique(seu.hdwgcna.AGMs_CAMs$SubType),
  group.by='SubType',
  assay = 'RNA',
  slot = 'data'
)

seu.hdwgcna.AGMs_CAMs <- TestSoftPowers(
  seu.hdwgcna.AGMs_CAMs,
  networkType = 'signed'
)
plot_list <- PlotSoftPowers(seu.hdwgcna.AGMs_CAMs)
wrap_plots(plot_list, ncol=2)

seu.hdwgcna.AGMs_CAMs <- ConstructNetwork(
  seu.hdwgcna.AGMs_CAMs,
  soft_power=16,
  tom_name = 'AGMs_CAMs',
  minModuleSize = 20
)

seu.hdwgcna.AGMs_CAMs <- ScaleData(seu.hdwgcna.AGMs_CAMs,
                                   features=VariableFeatures(seu.hdwgcna.AGMs_CAMs))
seu.hdwgcna.AGMs_CAMs <- ModuleEigengenes(seu.hdwgcna.AGMs_CAMs)
seu.hdwgcna.AGMs_CAMs <- ModuleConnectivity(
  seu.hdwgcna.AGMs_CAMs,
  group_name = unique(seu.hdwgcna.AGMs_CAMs$SubType),
  group.by='SubType'
)
seu.hdwgcna.AGMs_CAMs <- ResetModuleNames(
  seu.hdwgcna.AGMs_CAMs,
  new_name = "M"
)
```


```{r}
PlotDendrogram(seu.hdwgcna.AGMs_CAMs, main='INH hdWGCNA Dendrogram')
```





