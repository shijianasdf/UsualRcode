KEGG_enrichment <- function (genes){
	#@param genes:���򼯺�,������entrezID
	library(KEGGprofile)
    pho_KEGGresult <- find_enriched_pathway(genes, species = "hsa") 
	return(pho_KEGGresult)
}