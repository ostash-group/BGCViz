#' BGC list from antiSMASH for S.coelicolor
#'
#' Dataset containing list of clusters from antiSMASH5, found in S.coelicolor
#'
#' @format A data frame with 27 rows and 6 variables:
#' \describe{
#'   \item{Cluster}{Cluster #, numeric}
#'   \item{Start}{Start location of a cluster, in bp}
#'   \item{Stop}{Stop location of a cluster, in bp}
#'   \item{Type}{Type of a cluster}
#'   \item{chromosome}{chromosome, on which cluster is, filled automatically ipon file upload}
#'   \item{Type2}{Type of a cluster, filled automatically upon file upload. Same as Type}
#' }
#' @source \url{https://raw.githubusercontent.com/pavlohrab/BGCViz-datasets/main/example_data/sco_antismash.csv}
"anti_data"
#' Gene list from ARTS for S.coelicolor
#'
#' Dataset containing list of core and resistance genes from ARTS, found in S.coelicolor. 
#' It is generated automatically upon zip folder from ARTS upload
#'
#' @format A data frame with 161 rows and 13 variables:
#' \describe{
#'   \item{Cluster}{Gene #, numeric}
#'   \item{Start}{Start location of a gene, in bp}
#'   \item{Stop}{Stop location of a gene, in bp}
#'   \item{Type}{Type of a gene}
#'   \item{Hit}{Unique core gene hit identifier, based on count of genes}
#'   \item{Type2}{Type of a gene, same as Type}
#'   \item{Core}{Type of core model from ARTS}
#'   \item{Description}{Annotation of a gene}
#'   \item{Count}{Count of a gene hit, to keep track of multiple hits}
#'   \item{ID}{Gene #, same as Cluster}
#'   \item{Evalue}{E-value of a resistance gene hit}
#'   \item{Bitscore}{Bitscore of a resistance gene hit}
#'   \item{Model}{Type of Model of gene hit. If core gene, then the model is Core}
#' }
#' @source \url{https://raw.githubusercontent.com/pavlohrab/BGCViz-datasets/main/example_data/sco_arts.csv}
"arts_data"
#' BGC list from DeepBGC for S.coelicolor
#'
#' Dataset containing list of clusters from DeepBGC, found in S.coelicolor
#' For full description see DeepBGC documentation!
#'
#' @format A data frame with 172 rows and 28 variables:
#' @source \url{https://raw.githubusercontent.com/pavlohrab/BGCViz-datasets/main/example_data/sco_deep.tsv}
"deep_data"
#' BGC list from GECCO for S.coelicolor
#'
#' Dataset containing list of clusters from GECCO, found in S.coelicolor
#'
#' @format A data frame with 27 rows and 6 variables:
#' \describe{
#'   \item{Cluster}{Cluster #, numeric}
#'   \item{Start}{Start location of a cluster, in bp}
#'   \item{Stop}{Stop location of a cluster, in bp}
#'   \item{Type}{Type of a cluster}
#'   \item{chromosome}{chromosome, on which cluster is, filled automatically ipon file upload}
#'   \item{Type2}{Type of a cluster, filled automatically upon file upload. Same as Type}
#'   \item{ID}{Cluster #. Same as Cluster}
#'   \item{average_p}{Average probability of a cluster}
#'   \item{max_p}{Maximum probability of a cluster}
#'   \item{bgc_id}{ID of a cluster, filled by GECCO}
#'   \item{sequence_id}{ID of annotated sequence}
#'   \item{type}{Type of a cluster, filled by GECCO}
#'   \item{proteins}{list of identified proteins in a cluster, filled by GECCO}
#'   \item{domains}{list of identified domains in a cluster, filled by GECCO}
#'   \item{pks}{Score that cluster is pks, filled by GECCO}
#'   \item{other}{Score that cluster is other, filled by GECCO}
#'   \item{nrps}{Score that cluster is nrps, filled by GECCO}
#'   \item{alkaloid}{Score that cluster is alkaloid, filled by GECCO}
#'   \item{terpene}{Score that cluster is terpene, filled by GECCO}
#'   \item{sacharide}{Score that cluster is sacharide, filled by GECCO}
#'   \item{ripp}{Score that cluster is ripp, filled by GECCO}
#'   \item{num_prot}{Number of proteins in a cluster}
#'   \item{num_domains}{Number of domains in a cluster}
#' }
#' @source \url{https://raw.githubusercontent.com/pavlohrab/BGCViz-datasets/main/example_data/sco_gecco.tsv}
"gecco_data"
#' BGC list from PRISM for S.coelicolor
#'
#' Dataset containing list of clusters from PRISM, found in S.coelicolor
#' Generated automatically from json results file
#' @format A data frame with 19 rows and 4 variables:
#' \describe{
#'   \item{Cluster}{Cluster #, numeric}
#'   \item{Start}{Start location of a cluster, in bp}
#'   \item{Stop}{Stop location of a cluster, in bp}
#'   \item{Type}{Type of a cluster}
#' }
#' @source \url{https://raw.githubusercontent.com/pavlohrab/BGCViz-datasets/main/example_data/sco_prism.json}
"prism_data"
#' Gene list from PRISM for S.coelicolor
#'
#' Dataset containing list of regulatory and resistance from PRISM, found in S.coelicolor.
#' Generated automatically from json results file
#'
#' @format A data frame with 125 rows and 9 variables:
#' \describe{
#'   \item{Cluster}{Gene #, numeric}
#'   \item{Start}{Start location of a gene, in bp}
#'   \item{Stop}{Stop location of a gene, in bp}
#'   \item{Type}{Type of a gene}
#'   \item{Full_name}{Annotation of a gene hit}
#'   \item{ID}{Gene #, same as Cluster}
#'   \item{Type2}{Type of a gene. Same as Type}
#'   \item{Score}{Score of a gene hite}
#'   \item{Name}{Name of model, or gene for a hit}
#' }
#' @source \url{https://raw.githubusercontent.com/pavlohrab/BGCViz-datasets/main/example_data/sco_prism.json}
"prism_supp_data"
#' RRE list from RREFinder for S.coelicolor
#'
#' Dataset containing list of RRE elements from RREFinder, found in S.coelicolor
#'
#' @format A data frame with 27 rows and 6 variables:
#' \describe{
#'   \item{Cluster}{RRE #, numeric}
#'   \item{Start}{Start location of a RRE, in bp}
#'   \item{Stop}{Stop location of a RRE, in bp}
#'   \item{Type}{Type of a RRE}
#'   \item{chromosome}{chromosome, on which RRE is, filled automatically ipon file upload}
#'   \item{Type2}{Type of a RRE, filled automatically upon file upload. Same as Type}
#'   \item{Sequence}{Annotated sequence ID}
#'   \item{Locus_tag}{Locus tag of a RRE element}
#'   \item{BGC.ID}{Id, of a BGC}
#'   \item{BGC.product}{Product of BGC}
#'   \item{Domain.name}{Domain name of RRE element}
#'   \item{E.value}{E-value of a hit}
#'   \item{Bitscore}{Bitscore of a hit}
#'   \item{End}{End of RRE element, bp}
#'   \item{ID}{RRE ID, same as Cluster}
#' }
#' @source \url{https://raw.githubusercontent.com/pavlohrab/BGCViz-datasets/main/example_data/sco_rre.txt}
"rre_data"
#' BGC list from SEMPI for S.coelicolor
#'
#' Dataset containing list of clusters from SEMPI, found in S.coelicolor
#'
#' @format A data frame with 33 rows and 5 variables:
#' \describe{
#'   \item{Cluster}{Cluster #, numeric}
#'   \item{Start}{Start location of a cluster, in bp}
#'   \item{Stop}{Stop location of a cluster, in bp}
#'   \item{Type}{Type of a cluster}
#'   \item{Type2}{Type of a cluster, filled automatically upon file upload. Same as Type}
#' }
#' @source \url{https://raw.githubusercontent.com/pavlohrab/BGCViz-datasets/main/example_data/sco_sempi.csv}
"sempi_data"