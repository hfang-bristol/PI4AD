#' Function to visualise effect-by-removal using a upset plot
#'
#' \code{oVisAttack} is supposed to visualise effect-by-removal using a upset plot. The top is the kite plot, and visualised below is the combination matrix for nodes removed. It returns an object of class "ggplot".
#'
#' @param data a data frame. It contains all these columns including  'frac.disconnected', 'nodes.removed' and 'i' (number of nodes removed)
#' @param color the point color
#' @param shape the point shape
#' @param size the point size
#' @param label.height.unit the unit specifying the per-row height of the combination matrix. By default, it is NULL. If specified (such as 8), it will be used to decide 'combmatrix.label.height' and 'combmatrix.label.text'
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{oVisAttack}}
#' @include oVisAttack.r
#' @examples
#' \dontrun{
#' obj %>% transmute(value=nO, member=overlap) -> data
#' #data %>% arrange(-value) %>% mutate(member=fct_inorder(member)) -> data
#' data %>% oVisAttack() -> gp
#' gp + geom_point(aes(fill=as.factor(i)), shape=22, size=2, color="transparent") + geom_line(aes(group=i,color=as.factor(i)))
#' 
#' # customised levels for members
#' data %>% separate_rows(member,sep=',') %>% arrange(member) %>% pull(member) %>% unique() -> levels
#' data %>% oVisAttack(member.levels='customised',levels.customised=levels)
#' }

oVisAttack <- function(data, color="cyan4", shape=18, size=2, label.height.unit=7)
{
    
	i <- frac.disconnected <- nodes.removed <- NULL
    
    value <- member <- NULL
    
    df <- NULL
	if(is(data,'data.frame')){
		if(all(c('i','frac.disconnected','nodes.removed') %in% colnames(data))){
			df <- data %>% dplyr::distinct() %>% dplyr::arrange(i,frac.disconnected) %>% dplyr::mutate(x=forcats::fct_inorder(nodes.removed)) %>% dplyr::transmute(value=frac.disconnected, member=nodes.removed)
		}
	}

	##########################
	##########################
	if(nrow(df)==0){
		return(NULL)
	}
	##########################
	##########################
	
	levels.customised <- df %>% tidyr::separate_rows(member,sep=",") %>% pull(member) %>% unique()
	
	gp <- df %>% oUpsetAdv(member.levels="customised", shape=shape, size=size, color=color, levels.customised=levels.customised, label.height.unit=label.height.unit)
	gp <- gp + scale_y_continuous(limits=NULL) + geom_point(aes(fill=as.factor(i)), shape=22, size=2, color="transparent") + geom_line(aes(group=i,color=as.factor(i))) + guides(color="none")
	gp <- gp + ylab("Effect-by-removal\n(% nodes disconnected upon removal)") + scale_fill_discrete(guide=guide_legend("# nodes removed", title.position="top", nrow=1)) + theme(legend.position="top", axis.title.y=element_text(size=10), axis.text.y=element_text(size=8))

	invisible(gp)
}

