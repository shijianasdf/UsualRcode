KEGG_enrichment <- function (genes){
	#@param genes:基因集合,必须是entrezID
	library(KEGGprofile)
    pho_KEGGresult <- find_enriched_pathway(genes, species = "hsa") 
	return(pho_KEGGresult)
}
