#' Dataset of subsampling sims of logit_DRR102664 Arimizu E. coli isolate 
#' 
#' This data contains the subsampling results at various coverages of E.coli isolate DRR102664 for eaeA gene presence
#' 
#' 
#' @format A data frame 1200 x 7: 
#'  \describe{
#'          \item{ID}{E.coli sample ID}
#'          \item{Depth}{sequencing depth}
#'          \item{genome length}{total length of the genome, rounded}
#'          \item{coverage}{mean coverage of the genome} 
#'          \item{GC_content}{GC content}
#'          \item{genomelength2}{unrounded genome length}
#'          \item{presence}{eaeA gene presence or absence}
#'           }
#'           
#'           
#'  @source {This data was taken from Arimizu et al. 2019 https://genome.cshlp.org/content/29/9/1495}
#'  
#'
#'  @examples
#'  data(logit_data)
"logit_data"