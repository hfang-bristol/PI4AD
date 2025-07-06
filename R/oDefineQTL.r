#' Function to extract eQTL-gene pairs given a list of SNPs or a customised eQTL mapping data
#'
#' \code{oDefineQTL} is supposed to extract eQTL-gene pairs given a list of SNPs or a customised eQTL mapping data.
#'
#' @param data NULL or an input vector containing SNPs. If NULL, all SNPs will be considered. If a input vector containing SNPs, SNPs should be provided as dbSNP ID (ie starting with rs). Alternatively, they can be in the format of 'chrN:xxx', where N is either 1-22 or X, xxx is number; for example, 'chr16:28525386'
#' @param include.QTL genes modulated by eQTL (also Lead SNPs or in LD with Lead SNPs) are also included. By default, it is 'NA' to disable this option. Otherwise, those genes modulated by eQTL will be included. Pre-built eQTL datasets are detailed in the section 'Note'
#' @param QTL.customised a user-input matrix or data frame with 4 columns: 1st column for SNPs/eQTLs, 2nd column for Genes, 3rd for eQTL mapping significance level (p-values or FDR), and 4th for contexts (required even though only one context is input). Alternatively, it can be a file containing these 4 columns. It is designed to allow the user analysing their eQTL data. This customisation (if provided) will populate built-in eQTL data; mysql -e "use pi; SELECT rs_id_dbSNP147_GRCh37p13,gene_name,pval_nominal,Tissue FROM GTEx_V7_pair WHERE rs_id_dbSNP147_GRCh37p13!='.';" > /var/www/bigdata/QTL.customised.txt
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param placeholder the characters to tell the location of built-in RDS files. See \code{\link{oRDS}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{oRDS}} for details
#' @return
#' a data frame with following columns:
#' \itemize{
#'  \item{\code{SNP}: eQTLs}
#'  \item{\code{Gene}: eQTL-containing genes}
#'  \item{\code{Sig}: the eQTL mapping significant level}
#'  \item{\code{Context}: the context in which eQTL data was generated}
#' }
#' @note Pre-built eQTL datasets are described below according to the data sources.\cr
#' 1. eQTL_xMEdb
#' \itemize{
#'  \item{\code{eQTL_xMEdb_CD14}: cis and trans-eQTLs in the resting/CD14+ state.}
#'  \item{\code{eQTL_xMEdb_LPS2}: cis and trans-eQTLs in the activating state induced by 2-hour LPS.}
#'  \item{\code{eQTL_xMEdb_LPS24}: cis and trans-eQTLs in the activating state induced by 24-hour LPS.}
#'  \item{\code{eQTL_xMEdb_IFN}: cis and trans-eQTLs in the activating state induced by 24-hour interferon-gamma.}
#'  \item{\code{eQTL_xMEdb_Bcell}: cis and trans-eQTLs in B cells.}
#'  \item{\code{eQTL_xMEdb_Blood}: cis and trans-eQTLs in the blood.}
#'  \item{\code{eQTL_xMEdb_CD4}: cis and trans-eQTLs in the CD4 cells.}
#'  \item{\code{eQTL_xMEdb_CD8}: cis and trans-eQTLs in the CD8 cells.}
#'  \item{\code{eQTL_xMEdb_Monocyte}: cis and trans-eQTLs in the monocytes.}
#'  \item{\code{eQTL_xMEdb_Neutrophil}: cis and trans-eQTLs in the neutrophils.}
#'  \item{\code{eQTL_xMEdb_NK}: cis and trans-eQTLs in the NK cells.}
#' }
#' 2. eQTL_DICE (PMID:30449622)
#' \itemize{
#'  \item{\code{eQTL_DICE_Bcell_naive}: cis-eQTLs in the naive B cells.}
#'  \item{\code{eQTL_DICE_CD4_naive}: cis-eQTLs in the naive CD4 cells.}
#'  \item{\code{eQTL_DICE_CD4_stimu}: cis-eQTLs in the stimulated CD4 cells.}
#'  \item{\code{eQTL_DICE_CD4_Tfh}: cis-eQTLs in the CD4 Tfh cells.}
#'  \item{\code{eQTL_DICE_CD4_Th1}: cis-eQTLs in the CD4 Th1 cells.}
#'  \item{\code{eQTL_DICE_CD4_Th17}: cis-eQTLs in the CD4 Th17 cells.}
#'  \item{\code{eQTL_DICE_CD4_Th2}: cis-eQTLs in the CD4 Th2 cells.}
#'  \item{\code{eQTL_DICE_CD8_naive}: cis-eQTLs in the naive CD8 cells.}
#'  \item{\code{eQTL_DICE_CD8_stimu}: cis-eQTLs in the stimulated CD8 cells.}
#'  \item{\code{eQTL_DICE_Macrophages_M2}: cis-eQTLs in the macrophage M2 cells.}
#'  \item{\code{eQTL_DICE_Monocytes}: cis-eQTLs in the monocyte cells.}
#'  \item{\code{eQTL_DICE_NK}: cis-eQTLs in the NK cells.}
#'  \item{\code{eQTL_DICE_Treg_memo}: cis-eQTLs in the memory T regulatory cells.}
#'  \item{\code{eQTL_DICE_Treg_naive}: cis-eQTLs in the naive T regulatory cells.}
#' }
#' @export
#' @seealso \code{\link{oDefineQTL}}
#' @include oDefineQTL.r
#' @examples
#' \dontrun{
#' # a) provide the SNPs with the significance info
#' gr <- ImmunoBase$AS$variants
#' data <- gr$Variant
#'
#' # b) define eQTL genes
#' df_SGS <- oDefineQTL(data, include.QTL="eQTL_xMEdb_Bcell", placeholder=placeholder)
#' }

oDefineQTL <- function(data=NULL, include.QTL=NA, QTL.customised=NULL, verbose=TRUE, placeholder=NULL, guid=NULL)
{

    ######################################################
    # Link to targets based on eQTL
    ######################################################
    
    default.include.QTL <- c("eQTL_xMEdb_Bcell","eQTL_xMEdb_Blood","eQTL_xMEdb_CD14","eQTL_xMEdb_CD4","eQTL_xMEdb_CD8","eQTL_xMEdb_IFN","eQTL_xMEdb_LPS2","eQTL_xMEdb_LPS24","eQTL_xMEdb_Monocyte","eQTL_xMEdb_Neutrophil","eQTL_xMEdb_NK", "eQTL_DICE_Bcell_naive","eQTL_DICE_CD4_naive","eQTL_DICE_CD4_stimu","eQTL_DICE_CD4_Tfh","eQTL_DICE_CD4_Th1","eQTL_DICE_CD4_Th17","eQTL_DICE_CD4_Th1h17","eQTL_DICE_CD4_Th2","eQTL_DICE_CD8_naive","eQTL_DICE_CD8_stimu","eQTL_DICE_Macrophages_M2","eQTL_DICE_Monocytes","eQTL_DICE_NK","eQTL_DICE_Treg_memo","eQTL_DICE_Treg_naive","eQTL_GTExV8_Adipose_Subcutaneous","eQTL_GTExV8_Adipose_Visceral_Omentum","eQTL_GTExV8_Adrenal_Gland","eQTL_GTExV8_Artery_Aorta","eQTL_GTExV8_Artery_Coronary","eQTL_GTExV8_Artery_Tibial","eQTL_GTExV8_Brain_Amygdala","eQTL_GTExV8_Brain_Anterior_cingulate_cortex_BA24","eQTL_GTExV8_Brain_Caudate_basal_ganglia","eQTL_GTExV8_Brain_Cerebellar_Hemisphere","eQTL_GTExV8_Brain_Cerebellum","eQTL_GTExV8_Brain_Cortex","eQTL_GTExV8_Brain_Frontal_Cortex_BA9","eQTL_GTExV8_Brain_Hippocampus","eQTL_GTExV8_Brain_Hypothalamus","eQTL_GTExV8_Brain_Nucleus_accumbens_basal_ganglia","eQTL_GTExV8_Brain_Putamen_basal_ganglia","eQTL_GTExV8_Brain_Spinal_cord_cervical_c1","eQTL_GTExV8_Brain_Substantia_nigra","eQTL_GTExV8_Breast_Mammary_Tissue","eQTL_GTExV8_Cells_Cultured_fibroblasts","eQTL_GTExV8_Cells_EBVtransformed_lymphocytes","eQTL_GTExV8_Colon_Sigmoid","eQTL_GTExV8_Colon_Transverse","eQTL_GTExV8_Esophagus_Gastroesophageal_Junction","eQTL_GTExV8_Esophagus_Mucosa","eQTL_GTExV8_Esophagus_Muscularis","eQTL_GTExV8_Heart_Atrial_Appendage","eQTL_GTExV8_Heart_Left_Ventricle","eQTL_GTExV8_Kidney_Cortex","eQTL_GTExV8_Liver","eQTL_GTExV8_Lung","eQTL_GTExV8_Minor_Salivary_Gland","eQTL_GTExV8_Muscle_Skeletal","eQTL_GTExV8_Nerve_Tibial","eQTL_GTExV8_Ovary","eQTL_GTExV8_Pancreas","eQTL_GTExV8_Pituitary","eQTL_GTExV8_Prostate","eQTL_GTExV8_Skin_Not_Sun_Exposed_Suprapubic","eQTL_GTExV8_Skin_Sun_Exposed_Lower_leg","eQTL_GTExV8_Small_Intestine_Terminal_Ileum","eQTL_GTExV8_Spleen","eQTL_GTExV8_Stomach","eQTL_GTExV8_Testis","eQTL_GTExV8_Thyroid","eQTL_GTExV8_Uterus","eQTL_GTExV8_Vagina","eQTL_GTExV8_Whole_Blood", "eQTL_eQTLGen", "eQTL_eQTLGen_PMID34475573", "pQTL_Plasma", "pQTL_Plasma_PMID29875488", "pQTL_Plasma_PMID34857953", "pQTL_Plasma_PMID37794186", "eQTL_LCL")
    
    default.include.QTL <- c(default.include.QTL, "eQTL_eQTLCatalogue_Alasoo_2018_macrophage_IFNg","eQTL_eQTLCatalogue_Alasoo_2018_macrophage_IFNgSalmonella","eQTL_eQTLCatalogue_Alasoo_2018_macrophage_Salmonella","eQTL_eQTLCatalogue_Alasoo_2018_macrophage_naive","eQTL_eQTLCatalogue_BLUEPRINT_Tcell","eQTL_eQTLCatalogue_BLUEPRINT_monocyte","eQTL_eQTLCatalogue_BLUEPRINT_neutrophil","eQTL_eQTLCatalogue_BrainSeq_brain","eQTL_eQTLCatalogue_CEDAR_microarray_Bcell_CD19","eQTL_eQTLCatalogue_CEDAR_microarray_Tcell_CD4","eQTL_eQTLCatalogue_CEDAR_microarray_Tcell_CD8","eQTL_eQTLCatalogue_CEDAR_microarray_ileum","eQTL_eQTLCatalogue_CEDAR_microarray_monocyte_CD14","eQTL_eQTLCatalogue_CEDAR_microarray_neutrophil_CD15","eQTL_eQTLCatalogue_CEDAR_microarray_platelet","eQTL_eQTLCatalogue_CEDAR_microarray_rectum","eQTL_eQTLCatalogue_CEDAR_microarray_transverse_colon","eQTL_eQTLCatalogue_FUSION_adipose_naive","eQTL_eQTLCatalogue_FUSION_muscle_naive","eQTL_eQTLCatalogue_Fairfax_2012_microarray_Bcell_CD19","eQTL_eQTLCatalogue_Fairfax_2014_microarray_monocyte_IFN24","eQTL_eQTLCatalogue_Fairfax_2014_microarray_monocyte_LPS2","eQTL_eQTLCatalogue_Fairfax_2014_microarray_monocyte_LPS24","eQTL_eQTLCatalogue_Fairfax_2014_microarray_monocyte_naive","eQTL_eQTLCatalogue_GENCORD_LCL","eQTL_eQTLCatalogue_GENCORD_Tcell","eQTL_eQTLCatalogue_GENCORD_fibroblast","eQTL_eQTLCatalogue_GEUVADIS_LCL","eQTL_eQTLCatalogue_GTEx_LCL")
    
    default.include.QTL <- c(default.include.QTL,"eQTL_eQTLCatalogue_GTEx_V8_Adipose_Subcutaneous","eQTL_eQTLCatalogue_GTEx_V8_Adipose_Visceral_Omentum","eQTL_eQTLCatalogue_GTEx_V8_Adrenal_Gland","eQTL_eQTLCatalogue_GTEx_V8_Artery_Aorta","eQTL_eQTLCatalogue_GTEx_V8_Artery_Coronary","eQTL_eQTLCatalogue_GTEx_V8_Artery_Tibial","eQTL_eQTLCatalogue_GTEx_V8_Brain_Amygdala","eQTL_eQTLCatalogue_GTEx_V8_Brain_Anterior_cingulate_cortex_BA24","eQTL_eQTLCatalogue_GTEx_V8_Brain_Caudate_basal_ganglia","eQTL_eQTLCatalogue_GTEx_V8_Brain_Cerebellar_Hemisphere","eQTL_eQTLCatalogue_GTEx_V8_Brain_Cerebellum","eQTL_eQTLCatalogue_GTEx_V8_Brain_Cortex","eQTL_eQTLCatalogue_GTEx_V8_Brain_Frontal_Cortex_BA9","eQTL_eQTLCatalogue_GTEx_V8_Brain_Hippocampus","eQTL_eQTLCatalogue_GTEx_V8_Brain_Hypothalamus","eQTL_eQTLCatalogue_GTEx_V8_Brain_Nucleus_accumbens_basal_ganglia","eQTL_eQTLCatalogue_GTEx_V8_Brain_Putamen_basal_ganglia","eQTL_eQTLCatalogue_GTEx_V8_Brain_Spinal_cord_cervical_c1","eQTL_eQTLCatalogue_GTEx_V8_Brain_Substantia_nigra","eQTL_eQTLCatalogue_GTEx_V8_Breast_Mammary_Tissue","eQTL_eQTLCatalogue_GTEx_V8_Cells_Cultured_fibroblasts","eQTL_eQTLCatalogue_GTEx_V8_Cells_EBVtransformed_lymphocytes","eQTL_eQTLCatalogue_GTEx_V8_Colon_Sigmoid","eQTL_eQTLCatalogue_GTEx_V8_Colon_Transverse","eQTL_eQTLCatalogue_GTEx_V8_Esophagus_Gastroesophageal_Junction","eQTL_eQTLCatalogue_GTEx_V8_Esophagus_Mucosa","eQTL_eQTLCatalogue_GTEx_V8_Esophagus_Muscularis","eQTL_eQTLCatalogue_GTEx_V8_Heart_Atrial_Appendage","eQTL_eQTLCatalogue_GTEx_V8_Heart_Left_Ventricle","eQTL_eQTLCatalogue_GTEx_V8_Kidney_Cortex","eQTL_eQTLCatalogue_GTEx_V8_Liver","eQTL_eQTLCatalogue_GTEx_V8_Lung","eQTL_eQTLCatalogue_GTEx_V8_Minor_Salivary_Gland","eQTL_eQTLCatalogue_GTEx_V8_Muscle_Skeletal","eQTL_eQTLCatalogue_GTEx_V8_Nerve_Tibial","eQTL_eQTLCatalogue_GTEx_V8_Ovary","eQTL_eQTLCatalogue_GTEx_V8_Pancreas","eQTL_eQTLCatalogue_GTEx_V8_Pituitary","eQTL_eQTLCatalogue_GTEx_V8_Prostate","eQTL_eQTLCatalogue_GTEx_V8_Skin_Not_Sun_Exposed_Suprapubic","eQTL_eQTLCatalogue_GTEx_V8_Skin_Sun_Exposed_Lower_leg","eQTL_eQTLCatalogue_GTEx_V8_Small_Intestine_Terminal_Ileum","eQTL_eQTLCatalogue_GTEx_V8_Spleen","eQTL_eQTLCatalogue_GTEx_V8_Stomach","eQTL_eQTLCatalogue_GTEx_V8_Testis","eQTL_eQTLCatalogue_GTEx_V8_Thyroid","eQTL_eQTLCatalogue_GTEx_V8_Uterus","eQTL_eQTLCatalogue_GTEx_V8_Vagina","eQTL_eQTLCatalogue_GTEx_V8_Whole_Blood")

    #default.include.QTL <- c(default.include.QTL,"eQTL_eQTLCatalogue_GTEx_adipose_subcutaneous","eQTL_eQTLCatalogue_GTEx_adipose_visceral","eQTL_eQTLCatalogue_GTEx_adrenal_gland","eQTL_eQTLCatalogue_GTEx_artery_aorta","eQTL_eQTLCatalogue_GTEx_artery_coronary","eQTL_eQTLCatalogue_GTEx_artery_tibial","eQTL_eQTLCatalogue_GTEx_blood","eQTL_eQTLCatalogue_GTEx_brain_amygdala","eQTL_eQTLCatalogue_GTEx_brain_anterior_cingulate_cortex","eQTL_eQTLCatalogue_GTEx_brain_caudate","eQTL_eQTLCatalogue_GTEx_brain_cerebellar_hemisphere","eQTL_eQTLCatalogue_GTEx_brain_cerebellum","eQTL_eQTLCatalogue_GTEx_brain_cortex","eQTL_eQTLCatalogue_GTEx_brain_frontal_cortex","eQTL_eQTLCatalogue_GTEx_brain_hippocampus","eQTL_eQTLCatalogue_GTEx_brain_hypothalamus","eQTL_eQTLCatalogue_GTEx_brain_nucleus_accumbens","eQTL_eQTLCatalogue_GTEx_brain_putamen","eQTL_eQTLCatalogue_GTEx_brain_spinal_cord","eQTL_eQTLCatalogue_GTEx_brain_substantia_nigra","eQTL_eQTLCatalogue_GTEx_breast","eQTL_eQTLCatalogue_GTEx_colon_sigmoid","eQTL_eQTLCatalogue_GTEx_colon_transverse","eQTL_eQTLCatalogue_GTEx_esophagus_gej","eQTL_eQTLCatalogue_GTEx_esophagus_mucosa","eQTL_eQTLCatalogue_GTEx_esophagus_muscularis","eQTL_eQTLCatalogue_GTEx_fibroblast","eQTL_eQTLCatalogue_GTEx_heart_atrial_appendage","eQTL_eQTLCatalogue_GTEx_heart_left_ventricle","eQTL_eQTLCatalogue_GTEx_kidney_cortex","eQTL_eQTLCatalogue_GTEx_liver","eQTL_eQTLCatalogue_GTEx_lung","eQTL_eQTLCatalogue_GTEx_minor_salivary_gland","eQTL_eQTLCatalogue_GTEx_muscle","eQTL_eQTLCatalogue_GTEx_nerve_tibial","eQTL_eQTLCatalogue_GTEx_ovary","eQTL_eQTLCatalogue_GTEx_pancreas","eQTL_eQTLCatalogue_GTEx_pituitary","eQTL_eQTLCatalogue_GTEx_prostate","eQTL_eQTLCatalogue_GTEx_skin_not_sun_exposed","eQTL_eQTLCatalogue_GTEx_skin_sun_exposed","eQTL_eQTLCatalogue_GTEx_small_intestine","eQTL_eQTLCatalogue_GTEx_spleen","eQTL_eQTLCatalogue_GTEx_stomach","eQTL_eQTLCatalogue_GTEx_testis","eQTL_eQTLCatalogue_GTEx_thyroid","eQTL_eQTLCatalogue_GTEx_uterus","eQTL_eQTLCatalogue_GTEx_vagina")
    
	default.include.QTL <- c(default.include.QTL,   "eQTL_eQTLCatalogue_HipSci_iPSC","eQTL_eQTLCatalogue_Kasela_2017_microarray_Tcell_CD4","eQTL_eQTLCatalogue_Kasela_2017_microarray_Tcell_CD8","eQTL_eQTLCatalogue_Lepik_2017_blood","eQTL_eQTLCatalogue_Naranbhai_2015_microarray_neutrophil_CD16","eQTL_eQTLCatalogue_Nedelec_2016_macrophage_Listeria","eQTL_eQTLCatalogue_Nedelec_2016_macrophage_Salmonella","eQTL_eQTLCatalogue_Nedelec_2016_macrophage_naive","eQTL_eQTLCatalogue_Quach_2016_monocyte_IAV","eQTL_eQTLCatalogue_Quach_2016_monocyte_LPS","eQTL_eQTLCatalogue_Quach_2016_monocyte_Pam3CSK4","eQTL_eQTLCatalogue_Quach_2016_monocyte_R848","eQTL_eQTLCatalogue_Quach_2016_monocyte_naive","eQTL_eQTLCatalogue_ROSMAP_brain_naive","eQTL_eQTLCatalogue_Schmiedel_2018_Bcell_naive","eQTL_eQTLCatalogue_Schmiedel_2018_CD4_Tcell_antiCD3CD28","eQTL_eQTLCatalogue_Schmiedel_2018_CD4_Tcell_naive","eQTL_eQTLCatalogue_Schmiedel_2018_CD8_Tcell_antiCD3CD28","eQTL_eQTLCatalogue_Schmiedel_2018_CD8_Tcell_naive","eQTL_eQTLCatalogue_Schmiedel_2018_NKcell_naive","eQTL_eQTLCatalogue_Schmiedel_2018_Tfh_memory","eQTL_eQTLCatalogue_Schmiedel_2018_Th117_memory","eQTL_eQTLCatalogue_Schmiedel_2018_Th17_memory","eQTL_eQTLCatalogue_Schmiedel_2018_Th1_memory","eQTL_eQTLCatalogue_Schmiedel_2018_Th2_memory","eQTL_eQTLCatalogue_Schmiedel_2018_Treg_memory","eQTL_eQTLCatalogue_Schmiedel_2018_Treg_naive","eQTL_eQTLCatalogue_Schmiedel_2018_monocyte_CD16_naive","eQTL_eQTLCatalogue_Schmiedel_2018_monocyte_naive","eQTL_eQTLCatalogue_Schwartzentruber_2018_sensory_neuron","eQTL_eQTLCatalogue_TwinsUK_LCL","eQTL_eQTLCatalogue_TwinsUK_blood","eQTL_eQTLCatalogue_TwinsUK_fat","eQTL_eQTLCatalogue_TwinsUK_skin","eQTL_eQTLCatalogue_van_de_Bunt_2015_pancreatic_islet")
    
    ######################################################
    # eQTLCatalogue
    ###########################
	default.include.QTL <- c(default.include.QTL, "eQTL_eQTLCatalogue_PMID22446964_Fairfax_2012_Bcell_CD19","eQTL_eQTLCatalogue_PMID23755361_GENCORD_LCL","eQTL_eQTLCatalogue_PMID23755361_GENCORD_Tcell","eQTL_eQTLCatalogue_PMID23755361_GENCORD_fibroblast","eQTL_eQTLCatalogue_PMID24037378_GEUVADIS_LCL","eQTL_eQTLCatalogue_PMID24604202_Fairfax_2014_monocyte_IFN24","eQTL_eQTLCatalogue_PMID24604202_Fairfax_2014_monocyte_LPS2","eQTL_eQTLCatalogue_PMID24604202_Fairfax_2014_monocyte_LPS24","eQTL_eQTLCatalogue_PMID24604202_Fairfax_2014_monocyte_naive","eQTL_eQTLCatalogue_PMID25436857_TwinsUK_LCL","eQTL_eQTLCatalogue_PMID25436857_TwinsUK_blood","eQTL_eQTLCatalogue_PMID25436857_TwinsUK_fat","eQTL_eQTLCatalogue_PMID25436857_TwinsUK_skin","eQTL_eQTLCatalogue_PMID26151758_Naranbhai_2015_neutrophil_CD16","eQTL_eQTLCatalogue_PMID26624892_van_de_Bunt_2015_pancreatic_islet","eQTL_eQTLCatalogue_PMID27768888_Quach_2016_monocyte_IAV","eQTL_eQTLCatalogue_PMID27768888_Quach_2016_monocyte_LPS","eQTL_eQTLCatalogue_PMID27768888_Quach_2016_monocyte_Pam3CSK4","eQTL_eQTLCatalogue_PMID27768888_Quach_2016_monocyte_R848","eQTL_eQTLCatalogue_PMID27768888_Quach_2016_monocyte_naive","eQTL_eQTLCatalogue_PMID27768889_Nedelec_2016_macrophage_Listeria","eQTL_eQTLCatalogue_PMID27768889_Nedelec_2016_macrophage_Salmonella","eQTL_eQTLCatalogue_PMID27768889_Nedelec_2016_macrophage_naive","eQTL_eQTLCatalogue_PMID27863251_BLUEPRINT_Tcell","eQTL_eQTLCatalogue_PMID27863251_BLUEPRINT_monocyte","eQTL_eQTLCatalogue_PMID27863251_BLUEPRINT_neutrophil","eQTL_eQTLCatalogue_PMID28248954_Kasela_2017_Tcell_CD4","eQTL_eQTLCatalogue_PMID28248954_Kasela_2017_Tcell_CD8","eQTL_eQTLCatalogue_PMID28388432_PhLiPS_HLC","eQTL_eQTLCatalogue_PMID28388432_PhLiPS_iPSC","eQTL_eQTLCatalogue_PMID28410642_iPSCORE_iPSC","eQTL_eQTLCatalogue_PMID28489815_HipSci_iPSC","eQTL_eQTLCatalogue_PMID28814792_KimHellmuth_2017_monocyte_IVT_6h","eQTL_eQTLCatalogue_PMID28814792_KimHellmuth_2017_monocyte_IVT_90min","eQTL_eQTLCatalogue_PMID28814792_KimHellmuth_2017_monocyte_LPS_6h","eQTL_eQTLCatalogue_PMID28814792_KimHellmuth_2017_monocyte_LPS_90min","eQTL_eQTLCatalogue_PMID28814792_KimHellmuth_2017_monocyte_MDP_6h","eQTL_eQTLCatalogue_PMID28814792_KimHellmuth_2017_monocyte_MDP_90min","eQTL_eQTLCatalogue_PMID28814792_KimHellmuth_2017_monocyte_naive","eQTL_eQTLCatalogue_PMID28869584_ROSMAP_brain_naive","eQTL_eQTLCatalogue_PMID28922377_Lepik_2017_blood","eQTL_eQTLCatalogue_PMID29229984_Schwartzentruber_2018_sensory_neuron","eQTL_eQTLCatalogue_PMID29379200_Alasoo_2018_macrophage_IFNg","eQTL_eQTLCatalogue_PMID29379200_Alasoo_2018_macrophage_IFNgSalmonella","eQTL_eQTLCatalogue_PMID29379200_Alasoo_2018_macrophage_Salmonella","eQTL_eQTLCatalogue_PMID29379200_Alasoo_2018_macrophage_naive","eQTL_eQTLCatalogue_PMID29930244_CEDAR_Bcell_CD19","eQTL_eQTLCatalogue_PMID29930244_CEDAR_Tcell_CD4","eQTL_eQTLCatalogue_PMID29930244_CEDAR_Tcell_CD8","eQTL_eQTLCatalogue_PMID29930244_CEDAR_ileum","eQTL_eQTLCatalogue_PMID29930244_CEDAR_monocyte_CD14","eQTL_eQTLCatalogue_PMID29930244_CEDAR_neutrophil_CD15","eQTL_eQTLCatalogue_PMID29930244_CEDAR_platelet","eQTL_eQTLCatalogue_PMID29930244_CEDAR_rectum","eQTL_eQTLCatalogue_PMID29930244_CEDAR_transverse_colon","eQTL_eQTLCatalogue_PMID30050107_BrainSeq_brain")
	
	default.include.QTL <- c(default.include.QTL, "eQTL_eQTLCatalogue_PMID30449622_Schmiedel_2018","eQTL_eQTLCatalogue_PMID30449622_Schmiedel_2018_Bcell_naive","eQTL_eQTLCatalogue_PMID30449622_Schmiedel_2018_CD4_Tcell_antiCD3CD28","eQTL_eQTLCatalogue_PMID30449622_Schmiedel_2018_CD4_Tcell_naive","eQTL_eQTLCatalogue_PMID30449622_Schmiedel_2018_CD8_Tcell_antiCD3CD28","eQTL_eQTLCatalogue_PMID30449622_Schmiedel_2018_CD8_Tcell_naive","eQTL_eQTLCatalogue_PMID30449622_Schmiedel_2018_NKcell_naive","eQTL_eQTLCatalogue_PMID30449622_Schmiedel_2018_Tfh_memory","eQTL_eQTLCatalogue_PMID30449622_Schmiedel_2018_Th117_memory","eQTL_eQTLCatalogue_PMID30449622_Schmiedel_2018_Th17_memory","eQTL_eQTLCatalogue_PMID30449622_Schmiedel_2018_Th1_memory","eQTL_eQTLCatalogue_PMID30449622_Schmiedel_2018_Th2_memory","eQTL_eQTLCatalogue_PMID30449622_Schmiedel_2018_Treg_memory","eQTL_eQTLCatalogue_PMID30449622_Schmiedel_2018_Treg_naive","eQTL_eQTLCatalogue_PMID30449622_Schmiedel_2018_monocyte_CD16_naive","eQTL_eQTLCatalogue_PMID30449622_Schmiedel_2018_monocyte_naive","eQTL_eQTLCatalogue_PMID30596636_Peng_2018_placenta_naive","eQTL_eQTLCatalogue_PMID31076557_FUSION_adipose_naive","eQTL_eQTLCatalogue_PMID31076557_FUSION_muscle_naive","eQTL_eQTLCatalogue_PMID31551426_CommonMind_DLPFC_naive","eQTL_eQTLCatalogue_PMID31626773_Walker_2019_Neocortex","eQTL_eQTLCatalogue_PMID32098967_Braineac2_putamen","eQTL_eQTLCatalogue_PMID32098967_Braineac2_substantia_nigra","eQTL_eQTLCatalogue_PMID32787775_CAP_LCL_naive","eQTL_eQTLCatalogue_PMID32787775_CAP_LCL_statin")
	
	default.include.QTL <- c(default.include.QTL, "eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Adipose_Subcutaneous","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Adipose_Visceral_Omentum","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Adrenal_Gland","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Artery_Aorta","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Artery_Coronary","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Artery_Tibial","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Brain_Amygdala","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Brain_Anterior_cingulate_cortex_BA24","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Brain_Caudate_basal_ganglia","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Brain_Cerebellar_Hemisphere","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Brain_Cerebellum","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Brain_Cortex","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Brain_Frontal_Cortex_BA9","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Brain_Hippocampus","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Brain_Hypothalamus","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Brain_Nucleus_accumbens_basal_ganglia","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Brain_Putamen_basal_ganglia","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Brain_Spinal_cord_cervical_c1","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Brain_Substantia_nigra","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Breast_Mammary_Tissue","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Cells_Cultured_fibroblasts","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Cells_EBVtransformed_lymphocytes","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Colon_Sigmoid","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Colon_Transverse","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Esophagus_Gastroesophageal_Junction","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Esophagus_Mucosa","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Esophagus_Muscularis","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Heart_Atrial_Appendage","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Heart_Left_Ventricle","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Kidney_Cortex","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Liver","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Lung","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Minor_Salivary_Gland","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Muscle_Skeletal","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Nerve_Tibial","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Ovary","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Pancreas","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Pituitary","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Prostate","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Skin_Not_Sun_Exposed_Suprapubic","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Skin_Sun_Exposed_Lower_leg","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Small_Intestine_Terminal_Ileum","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Spleen","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Stomach","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Testis","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Thyroid","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Uterus","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Vagina","eQTL_eQTLCatalogue_PMID32913098_GTEx_V8_Whole_Blood")
	
	#default.include.QTL <- c(default.include.QTL, "eQTL_eQTLCatalogue_PMID32913098_GTEx_LCL","eQTL_eQTLCatalogue_PMID32913098_GTEx_adipose_subcutaneous","eQTL_eQTLCatalogue_PMID32913098_GTEx_adipose_visceral","eQTL_eQTLCatalogue_PMID32913098_GTEx_adrenal_gland","eQTL_eQTLCatalogue_PMID32913098_GTEx_artery_aorta","eQTL_eQTLCatalogue_PMID32913098_GTEx_artery_coronary","eQTL_eQTLCatalogue_PMID32913098_GTEx_artery_tibial","eQTL_eQTLCatalogue_PMID32913098_GTEx_blood","eQTL_eQTLCatalogue_PMID32913098_GTEx_brain_amygdala","eQTL_eQTLCatalogue_PMID32913098_GTEx_brain_anterior_cingulate_cortex","eQTL_eQTLCatalogue_PMID32913098_GTEx_brain_caudate","eQTL_eQTLCatalogue_PMID32913098_GTEx_brain_cerebellar_hemisphere","eQTL_eQTLCatalogue_PMID32913098_GTEx_brain_cerebellum","eQTL_eQTLCatalogue_PMID32913098_GTEx_brain_cortex","eQTL_eQTLCatalogue_PMID32913098_GTEx_brain_frontal_cortex","eQTL_eQTLCatalogue_PMID32913098_GTEx_brain_hippocampus","eQTL_eQTLCatalogue_PMID32913098_GTEx_brain_hypothalamus","eQTL_eQTLCatalogue_PMID32913098_GTEx_brain_nucleus_accumbens","eQTL_eQTLCatalogue_PMID32913098_GTEx_brain_putamen","eQTL_eQTLCatalogue_PMID32913098_GTEx_brain_spinal_cord","eQTL_eQTLCatalogue_PMID32913098_GTEx_brain_substantia_nigra","eQTL_eQTLCatalogue_PMID32913098_GTEx_breast","eQTL_eQTLCatalogue_PMID32913098_GTEx_colon_sigmoid","eQTL_eQTLCatalogue_PMID32913098_GTEx_colon_transverse","eQTL_eQTLCatalogue_PMID32913098_GTEx_esophagus_gej","eQTL_eQTLCatalogue_PMID32913098_GTEx_esophagus_mucosa","eQTL_eQTLCatalogue_PMID32913098_GTEx_esophagus_muscularis","eQTL_eQTLCatalogue_PMID32913098_GTEx_fibroblast","eQTL_eQTLCatalogue_PMID32913098_GTEx_heart_atrial_appendage","eQTL_eQTLCatalogue_PMID32913098_GTEx_heart_left_ventricle","eQTL_eQTLCatalogue_PMID32913098_GTEx_kidney_cortex","eQTL_eQTLCatalogue_PMID32913098_GTEx_liver","eQTL_eQTLCatalogue_PMID32913098_GTEx_lung","eQTL_eQTLCatalogue_PMID32913098_GTEx_minor_salivary_gland","eQTL_eQTLCatalogue_PMID32913098_GTEx_muscle","eQTL_eQTLCatalogue_PMID32913098_GTEx_nerve_tibial","eQTL_eQTLCatalogue_PMID32913098_GTEx_ovary","eQTL_eQTLCatalogue_PMID32913098_GTEx_pancreas","eQTL_eQTLCatalogue_PMID32913098_GTEx_pituitary","eQTL_eQTLCatalogue_PMID32913098_GTEx_prostate","eQTL_eQTLCatalogue_PMID32913098_GTEx_skin_not_sun_exposed","eQTL_eQTLCatalogue_PMID32913098_GTEx_skin_sun_exposed","eQTL_eQTLCatalogue_PMID32913098_GTEx_small_intestine","eQTL_eQTLCatalogue_PMID32913098_GTEx_spleen","eQTL_eQTLCatalogue_PMID32913098_GTEx_stomach","eQTL_eQTLCatalogue_PMID32913098_GTEx_testis","eQTL_eQTLCatalogue_PMID32913098_GTEx_thyroid","eQTL_eQTLCatalogue_PMID32913098_GTEx_uterus","eQTL_eQTLCatalogue_PMID32913098_GTEx_vagina")
	
	default.include.QTL <- c(default.include.QTL, "eQTL_eQTLCatalogue_PMID33637762_Steinberg_2020_high_grade_cartilage_naive","eQTL_eQTLCatalogue_PMID33637762_Steinberg_2020_low_grade_cartilage_naive","eQTL_eQTLCatalogue_PMID33637762_Steinberg_2020_synovium_naive","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_Astro_Naive_D52","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_Astro_ROT_D52","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_DA_Naive_D30","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_DA_Naive_D52","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_DA_ROT_D52","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_Epen1_Naive_D30","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_Epen1_Naive_D52","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_Epen1_ROT_D52","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_Epen2_Naive_D52","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_Epen2_ROT_D52","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_FPP_Naive_D11","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_FPP_Naive_D30","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_FPP_Naive_D52","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_FPP_ROT_D52","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_NB_Naive_D11","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_P_FPP_Naive_D11","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_P_FPP_Naive_D30","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_P_FPP_Naive_D52","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_P_FPP_ROT_D52","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_P_Sert_Naive_D52","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_P_Sert_ROT_D52","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_Sert_Naive_D30","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_Sert_Naive_D52","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_Sert_ROT_D52","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_U_Neur1_Naive_D30","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_U_Neur1_Naive_D52","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_U_Neur1_ROT_D52","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_U_Neur2_Naive_D30","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_U_Neur3_Naive_D52","eQTL_eQTLCatalogue_PMID33664506_Jerber_2021_U_Neur3_ROT_D52","eQTL_eQTLCatalogue_PMID34083789_Young_2019_microglia_naive","eQTL_eQTLCatalogue_PMID34416157_Aygun_2021_Neuron","eQTL_eQTLCatalogue_PMID34416157_Aygun_2021_Progenitor","eQTL_eQTLCatalogue_PMID34644572_PISA_pancreatic_islet","eQTL_eQTLCatalogue_PMID34822289_Randolph_2021_B_NI","eQTL_eQTLCatalogue_PMID34822289_Randolph_2021_B_flu","eQTL_eQTLCatalogue_PMID34822289_Randolph_2021_CD4_T_NI","eQTL_eQTLCatalogue_PMID34822289_Randolph_2021_CD4_T_flu","eQTL_eQTLCatalogue_PMID34822289_Randolph_2021_CD8_T_NI","eQTL_eQTLCatalogue_PMID34822289_Randolph_2021_CD8_T_flu","eQTL_eQTLCatalogue_PMID34822289_Randolph_2021_NK_NI","eQTL_eQTLCatalogue_PMID34822289_Randolph_2021_NK_flu","eQTL_eQTLCatalogue_PMID34822289_Randolph_2021_highly_infected_flu","eQTL_eQTLCatalogue_PMID34822289_Randolph_2021_infected_monocytes_flu","eQTL_eQTLCatalogue_PMID34822289_Randolph_2021_monocytes_NI","eQTL_eQTLCatalogue_PMID34822289_Randolph_2021_monocytes_flu")
	
	default.include.QTL <- c(default.include.QTL,"eQTL_eQTLCatalogue_PMID35389779_OneK1K_B_intermediate","eQTL_eQTLCatalogue_PMID35389779_OneK1K_B_memory","eQTL_eQTLCatalogue_PMID35389779_OneK1K_B_naive","eQTL_eQTLCatalogue_PMID35389779_OneK1K_CD14_Mono","eQTL_eQTLCatalogue_PMID35389779_OneK1K_CD16_Mono","eQTL_eQTLCatalogue_PMID35389779_OneK1K_CD4_CTL","eQTL_eQTLCatalogue_PMID35389779_OneK1K_CD4_Naive","eQTL_eQTLCatalogue_PMID35389779_OneK1K_CD4_TCM","eQTL_eQTLCatalogue_PMID35389779_OneK1K_CD4_TEM","eQTL_eQTLCatalogue_PMID35389779_OneK1K_CD8_Naive","eQTL_eQTLCatalogue_PMID35389779_OneK1K_CD8_TCM","eQTL_eQTLCatalogue_PMID35389779_OneK1K_CD8_TEM","eQTL_eQTLCatalogue_PMID35389779_OneK1K_HSPC","eQTL_eQTLCatalogue_PMID35389779_OneK1K_MAIT","eQTL_eQTLCatalogue_PMID35389779_OneK1K_NK","eQTL_eQTLCatalogue_PMID35389779_OneK1K_NK_CD56bright","eQTL_eQTLCatalogue_PMID35389779_OneK1K_NK_Proliferating","eQTL_eQTLCatalogue_PMID35389779_OneK1K_Plasmablast","eQTL_eQTLCatalogue_PMID35389779_OneK1K_Platelet","eQTL_eQTLCatalogue_PMID35389779_OneK1K_Treg","eQTL_eQTLCatalogue_PMID35389779_OneK1K_cDC2","eQTL_eQTLCatalogue_PMID35389779_OneK1K_dnT","eQTL_eQTLCatalogue_PMID35389779_OneK1K_gdT","eQTL_eQTLCatalogue_PMID35389779_OneK1K_pDC","eQTL_eQTLCatalogue_PMID35389781_Perez_2022_B","eQTL_eQTLCatalogue_PMID35389781_Perez_2022_NK","eQTL_eQTLCatalogue_PMID35389781_Perez_2022_Prolif","eQTL_eQTLCatalogue_PMID35389781_Perez_2022_T4","eQTL_eQTLCatalogue_PMID35389781_Perez_2022_T8","eQTL_eQTLCatalogue_PMID35389781_Perez_2022_cDC","eQTL_eQTLCatalogue_PMID35389781_Perez_2022_cM","eQTL_eQTLCatalogue_PMID35389781_Perez_2022_ncM","eQTL_eQTLCatalogue_PMID35389781_Perez_2022_pDC","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD48_PD1TIGIT","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_CCR4","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_CCR4ICOS_central","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_CCR4_central","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_CCR5_cytotoxic","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_CD161_Th1","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_CD161_Th2","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_CD161_cytotoxic","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_CD27","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_CD27CD161","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_CD38ICOS_central","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_HLADR","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_RORC_Treg","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_Th1","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_Th17","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_Th171","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_Th2","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_Treg","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_activated","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_central","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_cytotoxic","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD4_lncRNA","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD8_CXCR3","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD8_GZMB","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD8_GZMK","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD8_activated","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_CD8_central","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_Vd1","eQTL_eQTLCatalogue_PMID35545678_Nathan_2022_Vd2","eQTL_eQTLCatalogue_PMID35591976_BossiniCastillo_2019_Treg_naive")
	
	default.include.QTL <- c(default.include.QTL,"eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_HSP_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TCM1_Lowly_Active_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TCM1_Lowly_Active_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TCM2_Lowly_Active_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TCM2_Lowly_Active_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TCM_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TCM_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TCM_STIM_5D","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TCM_UNS_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TEMRA_Lowly_Active_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TEMRA_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TEMRA_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TEMRA_STIM_5D","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TEMRA_UNS_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TEM_HLA_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TEM_HLA_STIM_5D","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TEM_Lowly_Active_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TEM_Lowly_Active_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TEM_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TEM_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TEM_STIM_5D","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TEM_UNS_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TM_ERstress_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TM_cycling_STIM_5D","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN1_Lowly_Active_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN1_Lowly_Active_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN1_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN1_UNS_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN2_Lowly_Active_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN2_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN2_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN2_UNS_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN3_Lowly_Active_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN3_UNS_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TNTCM_CXCR4_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN_A_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN_A_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN_B_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN_C_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN_HSP_STIM_5D","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN_IFN_Lowly_Active_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN_IFN_Lowly_Active_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN_IFN_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN_IFN_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN_IFN_STIM_5D","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN_NFKB_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN_NFKB_STIM_5D","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN_STIM_5D","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN_cycling_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_TN_cycling_STIM_5D","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_T_ERstress_STIM_5D","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_T_NFKB_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_combined_CD4_Memory_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_combined_CD4_Memory_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_combined_CD4_Memory_STIM_5D","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_combined_CD4_Memory_UNS_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_combined_CD4_Naive_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_combined_CD4_Naive_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_combined_CD4_Naive_STIM_5D","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_combined_CD4_Naive_UNS_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_nTreg_STIM_16H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_nTreg_STIM_40H","eQTL_eQTLCatalogue_PMID35618845_Cytoimmgen_nTreg_UNS_16H")
	
	default.include.QTL <- c(default.include.QTL,"eQTL_eQTLCatalogue_PMID35835762_Gilchrist_2021_NKcell_naive")
	
	#default.include.QTL <- c(default.include.QTL, )
    ###########################
    ######################################################
    
	ind <- match(default.include.QTL, include.QTL)
	include.QTL <- default.include.QTL[!is.na(ind)]
    
    df_SGS <- NULL
    if(length(include.QTL) > 0){
		###########################	
		# built-in eQTL
		###########################	
    	if(0){
			# only load once 'xQTLdb.xMEdb'
			if(any(grepl("eQTL_xMEdb_", include.QTL, perl=TRUE))){
				xQTLdb.xMEdb <- oRDS('xQTLdb.xMEdb', placeholder=placeholder, guid=guid, verbose=verbose)
			}
			# only load once 'xQTLdb.DICE'
			if(any(grepl("eQTL_DICE_", include.QTL, perl=TRUE))){
				eQTL_DICE <- oRDS('xQTLdb.DICE', placeholder=placeholder, guid=guid, verbose=verbose)
			}
		}
		
		res_list <- lapply(include.QTL, function(x){

			if(verbose){
				now <- Sys.time()
				message(sprintf("Processing %s ...", x), appendLF=TRUE)
			}
			
			if(any(grep("eQTL_xMEdb_", x, perl=TRUE))){
				if(0){
					ind <- match(xQTLdb.xMEdb$Contexts, x)
					df <- xQTLdb.xMEdb[!is.na(ind), ] %>% as.data.frame()
				}else{
					x <- x %>% str_replace('eQTL_xMEdb_','xQTLdb.eQTL_xMEdb_')
					df <- oRDS(x, placeholder=placeholder, guid=guid, verbose=verbose) %>% as.data.frame()
				}
			
			}else if(any(grepl("eQTL_DICE_", x, perl=TRUE))){
				if(0){
					ind <- match(eQTL_DICE$Contexts, x)
					df <- eQTL_DICE[!is.na(ind), ] %>% as.data.frame()
				}else{
					x <- x %>% str_replace('eQTL_DICE_','xQTLdb.eQTL_DICE_')
					df <- oRDS(x, placeholder=placeholder, guid=guid, verbose=verbose) %>% as.data.frame()
				}

			}else if(any(grepl("eQTL_GTExV8_", x, perl=TRUE))){
				x <- x %>% str_replace('eQTL_GTExV8_','xQTLdb.eQTL_GTExV8_')
				df <- oRDS(x, placeholder=placeholder, guid=guid, verbose=verbose) %>% as.data.frame()

			#}else if(x=="eQTL_eQTLGen"){
			#	df <- oRDS("xQTLdb.eQTL_eQTLGen", placeholder=placeholder, guid=guid, verbose=verbose) %>% as.data.frame()
			
			}else if(any(grepl("eQTL_eQTLGen", x, perl=TRUE))){
				x <- x %>% str_replace('eQTL_eQTLGen','xQTLdb.eQTL_eQTLGen')
				df <- oRDS(x, placeholder=placeholder, guid=guid, verbose=verbose) %>% as.data.frame()
			
			}else if(any(grepl("eQTL_eQTLCatalogue_", x, perl=TRUE))){
				x <- x %>% str_replace('eQTL_eQTLCatalogue_','xQTLdb.eQTL_eQTLCatalogue_')
				df <- oRDS(x, placeholder=placeholder, guid=guid, verbose=verbose) %>% as.data.frame()
			
			#}else if(x=="pQTL_Plasma"){
			#	df <- oRDS("xQTLdb.pQTL_Plasma", placeholder=placeholder, guid=guid, verbose=verbose) %>% as.data.frame()
			
			}else if(any(grepl("pQTL_Plasma", x, perl=TRUE))){
				x <- x %>% str_replace('pQTL_Plasma','xQTLdb.pQTL_Plasma')
				df <- oRDS(x, placeholder=placeholder, guid=guid, verbose=verbose) %>% as.data.frame()
			
			}else if(x=="eQTL_LCL"){
				df <- oRDS("xQTLdb.eQTL_LCL", placeholder=placeholder, guid=guid, verbose=verbose) %>% as.data.frame()
			
			}else{
				df <- NULL
			}
			
			return(df)
		})
		## get data frame (SNP Gene FDR)
		SGS <- do.call(rbind, res_list)
	
		############################
		# remove Gene if NA
		# remove SNP if NA
		df_SGS <- SGS[!is.na(SGS[,1]) & !is.na(SGS[,2]),]
		############################
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("%d eGenes are built-in", length(unique(df_SGS[,2]))), appendLF=TRUE)
		}
	
	}

	###########################	
	# customised eQTL
	###########################
	df_SGS_customised <- NULL
	if(!is.null(QTL.customised)){
			
		if(is.vector(QTL.customised)){
			# assume a file
			df <- utils::read.delim(file=QTL.customised, header=TRUE, row.names=NULL, stringsAsFactors=FALSE)
		}else if(is.matrix(QTL.customised) | is.data.frame(QTL.customised)){
			df <- QTL.customised
		}
		
		if(!is.null(df)){
			colnames(df) <- c("SNP", "Gene", "Sig", "Context")
			SGS_customised <- df
			#SGS_customised <- cbind(df, Context=rep('Customised',nrow(df)), stringsAsFactors=FALSE)
			
			############################
			# remove Gene if NA
			# remove SNP if NA
			df_SGS_customised <- SGS_customised[!is.na(SGS_customised[,1]) & !is.na(SGS_customised[,2]),]
			############################
			
			if(verbose){
				now <- Sys.time()
				message(sprintf("%d eGenes are customised for %d contexts", length(unique(df_SGS_customised[,2])), length(unique(df_SGS_customised[,4]))), appendLF=TRUE)
			}
		}
	}
	
	#########################################
	df_SGS <- do.call(rbind, list(df_SGS, df_SGS_customised))
	#########################################
		
	if(!is.null(df_SGS)){
		############################
		# remove Gene if ''
		# remove SNP if ''
		df_SGS <- df_SGS[df_SGS[,1]!='' & df_SGS[,2]!='',]
		############################
	}
	
	###########################################
	if(!is.null(data)){
		## replace '_' with ':'
		data <- gsub("_", ":", data, perl=TRUE)
		## replace 'imm:' with 'chr'
		data <- gsub("imm:", "chr", data, perl=TRUE)
	
		data <- unique(data)
	
		## eQTL weight for input SNPs
		ind <- match(df_SGS[,1], data)
		df_SGS <- data.frame(df_SGS[!is.na(ind),])
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("A total of %d input SNPs with %d eGenes", length(data), length(unique(df_SGS[,2]))), appendLF=TRUE)
		}

	}
    
	
    invisible(df_SGS)
}
