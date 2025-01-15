devtools::install_github("elisabettasciacca/DEGGs", build_vignettes = TRUE)
data("BRCA_metadata") 
data("BRCA_normCounts")

subnetworks_object <- generate_subnetworks(normalised_counts = BRCA_normCounts,  
                                           metadata = BRCA_metadata,   subgroup_variable = "SUBTYPE",  
                                           subgroups = c("BRCA_Her2",   "BRCA_LumA"),   entrezIDs = TRUE,  
                                           convert_to_gene_symbols = TRUE,   cores = 2)
View_interactive_subnetwork(subnetworks_object)

table <- extract_sig_deggs(subnetworks_object)

print_regressions(gene_A = "NOTCH2", gene_B = "DTX4", deggs_object = subnetworks_object, legend_position = "bottomright")


# My data
SraRunTable <- as.data.frame(SraRunTable)
rownames(SraRunTable) <- SraRunTable$`Sample Name`
SraRunTable$`Sample Name` <-NULL

df1 <- as.data.frame(df1)
rownames(df1) <- df1$genes
df1$genes <-NULL
colnames(df1)

colnames(df1) <- rownames(SraRunTable)
subnetworks_object_new <- generate_subnetworks(normalised_counts = df1,  
                                           metadata = SraRunTable,   subgroup_variable = "SUBTYPE",  
                                           subgroups = c("Healthy",   "Rheumatoid arthritis"),   entrezIDs = FALSE,  
                                           convert_to_gene_symbols = FALSE,   cores = 2)
View_interactive_subnetwork(subnetworks_object_new)
table2 <- extract_sig_deggs(subnetworks_object_new)
write.csv(table2, "table2.csv")

print_regressions(gene_A = "ROCK1", gene_B = "PPP1CA", deggs_object = subnetworks_object_new, legend_position = "bottomright")

