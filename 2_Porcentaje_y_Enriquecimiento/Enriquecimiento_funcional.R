  """
  Name: Rafael Barrero Rodríguez
  Date: 2020-03-11
  Description: En este script recogemos y realizamos el enriquecimiento funcional
  a partir de los distintos ficheros BED con los intervalos genómicos donde se
  encuentra nuestro estado. 
  """
  
  root_folder <- paste0("/home/rafael/Master_UAM/Transcriptomica_RegulacionGenomica_Epigenomica/",
                        "3_Regulacion_Genomica_Epigenomica/Trabajo_Polycomb/anotation/",
                        "enriquecimiento_funcional/")
  setwd(root_folder)
  
  
  library (clusterProfiler); packageDescription ("clusterProfiler", fields = "Version") 
  library(topGO)
  library(org.Hs.eg.db)
  
  library("FGNet")
  library("AnnotationDbi")
  library("KEGGprofile")	# If not installed, we’ll use pathview instead (hopefully it will install)
  library("pathview")		# If not installed, then go to KEGG Mapper (web page https://www.genome.jp/kegg/tool/map_pathway2.html)
  library("gProfileR")	# If not installed (web page instead https://biit.cs.ut.ee/gprofiler/)
  
  
  
  #########################################
  # ALL SEGMENTS (segments_collapsed)
  #########################################
  
  # Empezamos realizando el análisis con todos los segmentos. 
  
  
  gene_table <- read.table("segments_collapsed/gene_table.tsv", header = TRUE, sep = "\t")
  gene_table <- dplyr::distinct(gene_table)
  
  """
  Lo primero que haremos será ver los distintos tipos de genes en los que se
  encuentra (o en los que aparece asociado nuestro estado)
  """
  
  gene_type <- as.data.frame(table(gene_table$Gene.type), stringsAsFactors = FALSE)
  
  # Pie Chart with Percentages
  slices <- gene_type$Freq[gene_type$Freq > 1000]
  lbls <- gene_type$Var1[gene_type$Freq > 1000]
  
  # Estamos desconsiderando elementos que tienen menos de 1000 ocurrencias
  slices <- c(slices, sum(gene_type$Freq)-sum(slices))
  lbls <- c(lbls, "Others")
  
  pct <- round(slices/sum(gene_type$Freq)*100)
  lbls <- paste(lbls, pct) # add percents to labels
  lbls <- paste(lbls,"%",sep="") # ad % to labels
  pie(slices,labels = lbls, col=rainbow(length(lbls)),
      main="Gene Type Distribution")
  
  
  """
  A continuación haremos el enriquecimiento funcional usando solo los genes del
  tipo protein_coding
  """
  
  # Tomamos los gene_symbols
  gene_symbols <- as.character(gene_table$Gene.name[gene_table$Gene.type == "protein_coding"])
  gene_symbols <- levels(factor(gene_symbols[!is.na(gene_symbols)]))
  
  # Tomamos ENTREZ ID
  gene_entrez_id <- as.character(gene_table$EntrezGene.ID[gene_table$Gene.type == "protein_coding"])
  gene_entrez_id <- levels(factor(gene_entrez_id[!is.na(gene_entrez_id)]))
  
  # Tomamos ENSEMBL ID
  gene_ensembl_id <- as.character(gene_table$Gene.stable.ID[gene_table$Gene.type == "protein_coding"])
  gene_ensembl_id <- levels(factor(gene_ensembl_id[!is.na(gene_ensembl_id)]))
  
  # Usamos enrichGO para enriquecer en GO
  
  ego1 <- enrichGO(gene          = gene_symbols,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
  results <- ego1@result
  
  dotplot(ego1, showCategory=15)
  # enrichMap(ego1, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
  # cnetplot(ego1)
  # par(cex = 0.65)
  # plotGOgraph(ego1, firstSigNodes = 5)
  
  """
  Los términos GO que aparecen más enriquecidos son los asociados con la activación
  de la respuesta inmune mediada por neutrófilos. Recordemos que nuestro estado
  epigenético se asociaba a genes que estaban preparados para expresarse cuando
  fuera necesario. Por ello, estos resultados son coherentes, al menos parcialmente,
  teniendo en cuenta que nuestras células de partida eran células del sistema
  inmune, en concreto, monocitos. Puede ser que muchos de los genes implicados
  en la respuesta por neutrófilos también participen en la activación de monocitos.
  """
  
  # A contniación haremos enriquecimiento de KEGG PATHWAYS
  
  geneLabels <- unlist(as.list(org.Hs.egSYMBOL))
  
  geneList <- geneLabels[which(geneLabels %in% gene_symbols)]
  
  kegg_enrich <- enrichKEGG(gene = names(geneList),
                            organism = "hsa",
                            keyType = "kegg",
                            pvalueCutoff = 0.05,
                            pAdjustMethod = "BH",
                            universe,
                            minGSSize = 10,
                            maxGSSize = 500,
                            qvalueCutoff = 0.2,
                            use_internal_data = FALSE)
  kegg_enrich_df <- kegg_enrich@result
  dotplot(kegg_enrich, showCategory=10)

dir.create("./segments_collapsed/KEGG_Profile")
setwd("./segments_collapsed/KEGG_Profile")

list_variantspergene <- setNames(rep(NA, length(gene_ensembl_id)), gene_ensembl_id) 
#procesos_kegg <-levels(ENSGenes$pathway)   # 6 pathways
procesos_kegg <- kegg_enrich_df$ID[c(2,3,5,8,9)]
for (keggNumber in procesos_kegg) { 
  plotKegg(paste0("hsa",keggNumber), geneExpr=list_variantspergene, geneIDtype="ENSEMBL") 
}


setwd(root_folder)

write.table(results[1:15, -c(8,9)], file = "go_results.tsv", row.names = FALSE,
            col.names = TRUE, quote = FALSE, sep = "\t")

write.table(kegg_enrich_df[1:15, -c(8,9)], file = "kegg_results.tsv", row.names = FALSE,
            col.names = TRUE, quote = FALSE, sep = "\t")
