#' Function to prepare genetic predictors given a list of seed SNPs together with the significance level (e.g. GWAS reported p-values)
#'
#' \code{oPierSNPsAdv} is supposed to prepare genetic predictors given a list of seed SNPs together with the significance level (e.g. GWAS reported p-values). Internally it calls \code{\link{oPierSNPs}} to prepare the distance predictor, the eQTL predictors (if required) and the HiC predictors (if required). It returns a list of class "pNode" objects.
#'
#' @param data a named input vector containing the sinificance level for nodes (dbSNP). For this named vector, the element names are dbSNP ID (or in the format such as 'chr16:28525386'), the element values for the significance level (measured as p-value or fdr). Alternatively, it can be a matrix or data frame with two columns: 1st column for dbSNP, 2nd column for the significance level
#' @param include.LD additional SNPs in LD with Lead SNPs are also included. By default, it is 'NA' to disable this option. Otherwise, LD SNPs will be included based on one or more of 5 super-populations from 1000 Genomics Project data (phase 3). They are "AFR", "AMR", "EAS", "EUR", and "SAS". Explanations for population code can be found at \url{http://www.1000genomes.org/faq/which-populations-are-part-your-study}
#' @param LD.customised a user-input matrix or data frame with 3 columns: 1st column for Lead SNPs, 2nd column for LD SNPs, and 3rd for LD r2 value. It is designed to allow the user analysing their pre-calculated LD info. This customisation (if provided) has the high priority over built-in LD SNPs
#' @param LD.r2 the LD r2 value. By default, it is 0.8, meaning that SNPs in LD (r2>=0.8) with input SNPs will be considered as LD SNPs. It can be any value from 0.8 to 1
#' @param significance.threshold the given significance threshold. By default, it is set to NULL, meaning there is no constraint on the significance level when transforming the significance level of SNPs into scores. If given, those SNPs below this are considered significant and thus scored positively. Instead, those above this are considered insigificant and thus receive no score
#' @param score.cap the maximum score being capped. By default, it is set to 10. If NULL, no capping is applied
#' @param distance.max the maximum distance between genes and SNPs. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby SNPs per gene
#' @param decay.kernel a character specifying a decay kernel function. It can be one of 'slow' for slow decay, 'linear' for linear decay, and 'rapid' for rapid decay. If no distance weight is used, please select 'constant'
#' @param decay.exponent an integer specifying a decay exponent. By default, it sets to 2
#' @param GR.SNP the genomic regions of SNPs. By default, it is 'dbSNP_Common', that is, Common SNPs from dbSNP (version 151) plus GWAS SNPs and their LD SNPs (hg19). Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to dbSNP IDs. Then, tell "GR.SNP" with your RData file name (with or without extension), plus specify your file RData path in "placeholder". Note: you can also load your customised GR object directly
#' @param GR.Gene the genomic regions of genes. By default, it is 'UCSC_knownGene', that is, UCSC known genes (together with genomic locations) based on human genome assembly hg19. Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to Gene Symbols. Then, tell "GR.Gene" with your RData file name (with or without extension), plus specify your file RData path in "placeholder". Note: you can also load your customised GR object directly
#' @param include.TAD TAD boundary regions are also included. By default, it is 'none' to disable this option. Otherwise, inclusion of a TAD dataset to pre-filter SNP-nGene pairs (i.e. only those within a TAD region will be kept). TAD datasets can be one of "GM12878"  (lymphoblast), "IMR90" (fibroblast), "MSC" (mesenchymal stem cell) ,"TRO" (trophoblasts-like cell), "H1" (embryonic stem cell), "MES" (mesendoderm) and "NPC" (neural progenitor cell). Explanations can be found at \doi{10.1016/j.celrep.2016.10.061}
#' @param include.QTL the eQTL supported currently. By default, it is 'NA' to disable this option. Pre-built eQTL datasets are detailed in \code{\link{oDefineQTL}}
#' @param QTL.customised a user-input matrix or data frame with 4 columns: 1st column for SNPs/eQTLs, 2nd column for Genes, 3rd for eQTL mapping significance level (p-values or FDR), and 4th for contexts (required even though only one context is input). Alternatively, it can be a file containing these 4 columns. It is designed to allow the user analysing their eQTL data. This customisation (if provided) will populate built-in eQTL data
#' @param include.RGB genes linked to input SNPs are also included. By default, it is 'NA' to disable this option. Otherwise, those genes linked to SNPs will be included according to Promoter Capture HiC (PCHiC) datasets. Pre-built HiC datasets are detailed in \code{\link{oDefineRGB}}
#' @param cdf.function a character specifying a Cumulative Distribution Function (cdf). It can be one of 'exponential' based on exponential cdf, 'empirical' for empirical cdf
#' @param scoring.scheme the method used to calculate seed gene scores under a set of SNPs. It can be one of "sum" for adding up, "max" for the maximum, and "sequential" for the sequential weighting. The sequential weighting is done via: \eqn{\sum_{i=1}{\frac{R_{i}}{i}}}, where \eqn{R_{i}} is the \eqn{i^{th}} rank (in a descreasing order)
#' @param network the built-in network. Currently two sources of network information are supported: the STRING database (version 10) and the Pathway Commons database (version 7). STRING is a meta-integration of undirect interactions from the functional aspect, while Pathways Commons mainly contains both undirect and direct interactions from the physical/pathway aspect. Both have scores to control the confidence of interactions. Therefore, the user can choose the different quality of the interactions. In STRING, "STRING_highest" indicates interactions with highest confidence (confidence scores>=900), "STRING_high" for interactions with high confidence (confidence scores>=700), "STRING_medium" for interactions with medium confidence (confidence scores>=400), and "STRING_low" for interactions with low confidence (confidence scores>=150). For undirect/physical interactions from Pathways Commons, "PCommonsUN_high" indicates undirect interactions with high confidence (supported with the PubMed references plus at least 2 different sources), "PCommonsUN_medium" for undirect interactions with medium confidence (supported with the PubMed references). For direct (pathway-merged) interactions from Pathways Commons, "PCommonsDN_high" indicates direct interactions with high confidence (supported with the PubMed references plus at least 2 different sources), and "PCommonsUN_medium" for direct interactions with medium confidence (supported with the PubMed references). In addition to pooled version of pathways from all data sources, the user can also choose the pathway-merged network from individual sources, that is, "PCommonsDN_Reactome" for those from Reactome, "PCommonsDN_KEGG" for those from KEGG, "PCommonsDN_HumanCyc" for those from HumanCyc, "PCommonsDN_PID" for those froom PID, "PCommonsDN_PANTHER" for those from PANTHER, "PCommonsDN_ReconX" for those from ReconX, "PCommonsDN_TRANSFAC" for those from TRANSFAC, "PCommonsDN_PhosphoSite" for those from PhosphoSite, and "PCommonsDN_CTD" for those from CTD. For direct (pathway-merged) interactions sourced from KEGG, it can be 'KEGG' for all, 'KEGG_metabolism' for pathways grouped into 'Metabolism', 'KEGG_genetic' for 'Genetic Information Processing' pathways, 'KEGG_environmental' for 'Environmental Information Processing' pathways, 'KEGG_cellular' for 'Cellular Processes' pathways, 'KEGG_organismal' for 'Organismal Systems' pathways, and 'KEGG_disease' for 'Human Diseases' pathways. 'REACTOME' for protein-protein interactions derived from Reactome pathways
#' @param STRING.only the further restriction of STRING by interaction type. If NA, no such restriction. Otherwide, it can be one or more of "neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score". Useful options are c("experimental_score","database_score"): only experimental data (extracted from BIND, DIP, GRID, HPRD, IntAct, MINT, and PID) and curated data (extracted from Biocarta, BioCyc, GO, KEGG, and Reactome) are used
#' @param weighted logical to indicate whether edge weights should be considered. By default, it sets to false. If true, it only works for the network from the STRING database
#' @param network.customised an object of class "igraph". By default, it is NULL. It is designed to allow the user analysing their customised network data that are not listed in the above argument 'network'. This customisation (if provided) has the high priority over built-in network. If the user provides the "igraph" object with the "weight" edge attribute, RWR will assume to walk on the weighted network
#' @param seeds.inclusive logical to indicate whether non-network seed genes are included for prioritisation. If TRUE (by default), these genes will be added to the netowrk
#' @param normalise the way to normalise the adjacency matrix of the input graph. It can be 'laplacian' for laplacian normalisation, 'row' for row-wise normalisation, 'column' for column-wise normalisation, or 'none'
#' @param restart the restart probability used for Random Walk with Restart (RWR). The restart probability takes the value from 0 to 1, controlling the range from the starting nodes/seeds that the walker will explore. The higher the value, the more likely the walker is to visit the nodes centered on the starting nodes. At the extreme when the restart probability is zero, the walker moves freely to the neighbors at each step without restarting from seeds, i.e., following a random walk (RW)
#' @param normalise.affinity.matrix the way to normalise the output affinity matrix. It can be 'none' for no normalisation, 'quantile' for quantile normalisation to ensure that columns (if multiple) of the output affinity matrix have the same quantiles
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param verbose.details logical to indicate whether the detailed messages from being-called functions will be displayed in the screen. By default, it sets to FALSE enabling messages 
#' @param placeholder the characters to tell the location of built-in RData files. See \code{\link{oRDS}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{oRDS}} for details
#' @return
#' A list of class "pNode" objects, each object having a list with following components:
#' \itemize{
#'  \item{\code{priority}: a matrix of nNode X 6 containing node priority information, where nNode is the number of nodes in the input graph, and the 6 columns are "name" (node names), "node" (1 for network genes, 0 for non-network seed genes), "seed" (1 for seeds, 0 for non-seeds), "weight" (weight values),  "priority" (the priority scores that are rescaled to the range [0,1]), "rank" (ranks of the priority scores), "description" (node description)}
#'  \item{\code{g}: an input "igraph" object}
#'  \item{\code{SNP}: a data frame of nSNP X 4 containing input SNPs and/or LD SNPs info, where nSNP is the number of input SNPs and/or LD SNPs, and the 4 columns are "SNP" (dbSNP), "Score" (the SNP score), "Pval" (the SNP p-value), "Flag" (indicative of Lead SNPs or LD SNPs)}
#'  \item{\code{Gene2SNP}: a data frame of nPair X 3 containing Gene-SNP pair info, where nPair is the number of Gene-SNP pairs, and the 3 columns are "Gene" (seed genes), "SNP" (dbSNP), "Score" (an SNP's genetic influential score on a seed gene)}
#'  \item{\code{nGenes}: if not NULL, it is a data frame containing nGene-SNP pair info}
#'  \item{\code{eGenes}: if not NULL, it is a data frame containing eGene-SNP pair info per context}
#'  \item{\code{cGenes}: if not NULL, it is a data frame containing cGene-SNP pair info per context}
#' }
#' @note This function calls \code{\link{oPierSNPs}} in a loop way generating the distance predictor, the eQTL predictors (if required) and the HiC predictors (if required).
#' @export
#' @seealso \code{\link{oPierSNPsAdv}}
#' @include oPierSNPsAdv.r
#' @examples
#' \dontrun{
#' # a) provide the SNPs with the significance info
#' ImmunoBase <- oRDS('ImmunoBase', placeholder=placeholder)
#' gr <- ImmunoBase$AS$variants
#' AS <- as.data.frame(gr) %>% select(Variant, Pvalue)
#'
#' # b) perform priority analysis
#' ls_pNode <- oPierSNPsAdv(data=AS, include.TAD='GM12878', include.QTL="JKng_mono", include.RGB='Monocytes', network="PCommonsUN_medium", restart=0.7, placeholder=placeholder)
#' #ls_pNode <- oPierSNPsAdv(data=AS, include.TAD='GM12878', include.QTL="JKng_mono", include.RGB='Monocytes', network="PCommonsUN_medium", restart=0.7, placeholder=placeholder, QTL.customised='QTL.customised.Artery.txt')
#' }

oPierSNPsAdv <- function(data, include.LD=NA, LD.customised=NULL, LD.r2=0.8, significance.threshold=5e-8, score.cap=100, distance.max=20000, decay.kernel=c("constant","slow","linear","rapid"), decay.exponent=2, GR.SNP="dbSNP_Common", GR.Gene="UCSC_knownGene", include.TAD=c("none","GM12878","IMR90","MSC","TRO","H1","MES","NPC"), include.QTL=NA, QTL.customised=NULL, include.RGB=NA, cdf.function=c("empirical","exponential"), scoring.scheme=c("max","sum","sequential"), network=c("STRING_highest","STRING_high","STRING_medium","STRING_low","PCommonsUN_high","PCommonsUN_medium","PCommonsDN_high","PCommonsDN_medium","PCommonsDN_Reactome","PCommonsDN_KEGG","PCommonsDN_HumanCyc","PCommonsDN_PID","PCommonsDN_PANTHER","PCommonsDN_ReconX","PCommonsDN_TRANSFAC","PCommonsDN_PhosphoSite","PCommonsDN_CTD", "KEGG","KEGG_metabolism","KEGG_genetic","KEGG_environmental","KEGG_cellular","KEGG_organismal","KEGG_disease","REACTOME"), STRING.only=c(NA,"neighborhood_score","fusion_score","cooccurence_score","coexpression_score","experimental_score","database_score","textmining_score")[1], weighted=FALSE, network.customised=NULL, seeds.inclusive=TRUE, normalise=c("laplacian","row","column","none"), restart=0.7, normalise.affinity.matrix=c("none","quantile"), verbose=TRUE, verbose.details=FALSE, placeholder=NULL, guid=NULL)
{

    startT <- Sys.time()
    if(verbose){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    decay.kernel <- match.arg(decay.kernel)
    cdf.function <- match.arg(cdf.function)
    scoring.scheme <- match.arg(scoring.scheme)
    network <- match.arg(network)
    normalise <- match.arg(normalise)
    normalise.affinity.matrix <- match.arg(normalise.affinity.matrix)
    
    ## force verbose.details to be FALSE if verbose is FALSE
    if(verbose==FALSE){
    	verbose.details <- FALSE
    }
    ####################################################################################
	if(verbose){
		message(sprintf("Preparing the distance predictor (%s) ...", as.character(Sys.time())), appendLF=TRUE)
	}
	relative.importance <- c(1,0,0)
    pNode_distance <- oPierSNPs(data=data, include.LD=include.LD, LD.customised=LD.customised, LD.r2=LD.r2, significance.threshold=significance.threshold, score.cap=score.cap, distance.max=distance.max, decay.kernel=decay.kernel, decay.exponent=decay.exponent, GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.TAD=include.TAD, include.QTL=NA, QTL.customised=NULL, include.RGB=NA, cdf.function=cdf.function, relative.importance=relative.importance, scoring.scheme=scoring.scheme, network=network, weighted=weighted, network.customised=network.customised, seeds.inclusive=seeds.inclusive, normalise=normalise, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, verbose=verbose.details, placeholder=placeholder, guid=guid)
    ls_pNode_distance <- list(pNode_distance)
    names(ls_pNode_distance) <- paste('nGene_', distance.max, '_', decay.kernel, sep='')
    
    ####################################################################################
    
    ls_pNode_eQTL <- NULL
    
    include.QTLs <- include.QTL[!is.na(include.QTL)]
    if(length(include.QTLs)>0){
		names(include.QTLs) <- include.QTLs
		ls_pNode_eQTL <- pbapply::pblapply(include.QTLs, function(x){
			if(verbose){
				message(sprintf("\nPreparing the eQTL predictor '%s' (%s) ...", x, as.character(Sys.time())), appendLF=TRUE)
			}
			relative.importance <- c(0,1,0)
			pNode <- oPierSNPs(data=data, include.LD=include.LD, LD.customised=LD.customised, LD.r2=LD.r2, significance.threshold=significance.threshold, score.cap=score.cap, distance.max=distance.max, decay.kernel=decay.kernel, decay.exponent=decay.exponent, GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.QTL=x, QTL.customised=NULL, include.RGB=NA, cdf.function=cdf.function, relative.importance=relative.importance, scoring.scheme=scoring.scheme, network=network, weighted=weighted, network.customised=network.customised, seeds.inclusive=seeds.inclusive, normalise=normalise, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, verbose=verbose.details, placeholder=placeholder, guid=guid)
			if(verbose & is.null(pNode)){
				message(sprintf("\tNote: this predictor '%s' is NULL", x))
			}
			return(pNode)
		})
		names(ls_pNode_eQTL) <- paste('eGene_', names(ls_pNode_eQTL), sep='')
    }
    
    ################################
    ################################
    ls_pNode_eQTL_customised <- NULL
    df_SGS_customised <- NULL
    if(!is.null(QTL.customised)){
    
		if(is.vector(QTL.customised)){
			# assume a file
			df <- utils::read.delim(file=QTL.customised, header=TRUE, row.names=NULL, stringsAsFactors=FALSE)
		}else if(is.matrix(QTL.customised) | is.data.frame(QTL.customised)){
			df <- QTL.customised %>% as.data.frame()
		}
		
		if(!is.null(df)){
			colnames(df) <- c("SNP", "Gene", "Sig", "Context")
			SGS_customised <- df
			#SGS_customised <- cbind(df, Context=rep('Customised',nrow(df)))
			
			############################
			# remove Gene if NA
			# remove SNP if NA
			df_SGS_customised <- SGS_customised[!is.na(SGS_customised[,1]) & !is.na(SGS_customised[,2]),]
			############################
		}
    }
    if(!is.null(df_SGS_customised)){
		ls_df <- split(x=df_SGS_customised, f=df_SGS_customised$Context)
		ls_pNode_eQTL_customised <- pbapply::pblapply(1:length(ls_df), function(i){
			if(verbose){
				message(sprintf("\nPreparing the customised eQTL predictor '%s' (%s) ...", names(ls_df)[i], as.character(Sys.time())), appendLF=TRUE)
			}
			relative.importance <- c(0,1,0)
			pNode <- oPierSNPs(data=data, include.LD=include.LD, LD.customised=LD.customised, LD.r2=LD.r2, significance.threshold=significance.threshold, score.cap=score.cap, distance.max=distance.max, decay.kernel=decay.kernel, decay.exponent=decay.exponent, GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.QTL=NA, QTL.customised=ls_df[[i]], include.RGB=NA, cdf.function=cdf.function, relative.importance=relative.importance, scoring.scheme=scoring.scheme, network=network, weighted=weighted, network.customised=network.customised, seeds.inclusive=seeds.inclusive, normalise=normalise, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, verbose=verbose.details, placeholder=placeholder, guid=guid)
			if(verbose & is.null(pNode)){
				message(sprintf("\tNote: this predictor '%s' is NULL", names(ls_df)[i]), appendLF=TRUE)
			}
			return(pNode)
		})
		names(ls_pNode_eQTL_customised) <- paste('eGene_', names(ls_df), sep='')
    }
    ls_pNode_eQTL <- c(ls_pNode_eQTL, ls_pNode_eQTL_customised)
    ################################
    ################################
    
    include.RGBs <- include.RGB[!is.na(include.RGB)]
    if(length(include.RGBs)>0){
		names(include.RGBs) <- include.RGBs
		ls_pNode_HiC <- pbapply::pblapply(include.RGBs, function(x){
			if(verbose){
				message(sprintf("\nPreparing the HiC predictor '%s' (%s) ...", x, as.character(Sys.time())), appendLF=TRUE)
			}
			relative.importance <- c(0,0,1)
			pNode <- oPierSNPs(data=data, include.LD=include.LD, LD.customised=LD.customised, LD.r2=LD.r2, significance.threshold=significance.threshold, score.cap=score.cap, distance.max=distance.max, decay.kernel=decay.kernel, decay.exponent=decay.exponent, GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.QTL=NA, QTL.customised=NULL, include.RGB=x, cdf.function=cdf.function, relative.importance=relative.importance, scoring.scheme=scoring.scheme, network=network, weighted=weighted, network.customised=network.customised, seeds.inclusive=seeds.inclusive, normalise=normalise, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, verbose=verbose.details, placeholder=placeholder, guid=guid)
			if(verbose & is.null(pNode)){
				message(sprintf("\tNote: this predictor '%s' has NULL", x), appendLF=TRUE)
			}
			return(pNode)
		})
		names(ls_pNode_HiC) <- paste('cGene_', names(ls_pNode_HiC), sep='')
	}else{
		ls_pNode_HiC <- NULL
	}
    
    ##########################################################################################
    ## prioritisation equally
    #relative.importance <- c(1/3,1/3,1/3)
    #pNode_all <- oPierSNPs(data=data, include.LD=include.LD, LD.customised=LD.customised, LD.r2=LD.r2, significance.threshold=significance.threshold, score.cap=score.cap, distance.max=distance.max, decay.kernel=decay.kernel, decay.exponent=decay.exponent, GR.SNP=GR.SNP, GR.Gene=GR.Gene, include.QTL=include.QTLs, QTL.customised=NULL, include.RGB=include.RGBs, cdf.function=cdf.function, relative.importance=relative.importance, scoring.scheme=scoring.scheme, network=network, weighted=weighted, network.customised=network.customised, seeds.inclusive=seeds.inclusive, normalise=normalise, restart=restart, normalise.affinity.matrix=normalise.affinity.matrix, verbose=verbose, placeholder=placeholder, guid=guid)
    ##########################################################################################
    ls_pNode <- c(ls_pNode_distance, ls_pNode_eQTL, ls_pNode_HiC)
    
    ####################################################################################
    endT <- Sys.time()
    if(verbose){
        message(paste(c("\nFinish at ",as.character(endT)), collapse=""), appendLF=TRUE)
    
		runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
		message(paste(c("Runtime in total (oPierSNPsAdv): ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }
    
    invisible(ls_pNode)
}
