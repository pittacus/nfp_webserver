if (!require('devtools'))  install.packages("devtools") 
if (!require('NFP')) devtools::install_github("yiluheihei/NFP")
require(NFP)
require(igraph)
require(KEGGgraph)
require(AnnotationDbi)

data(kegg_refnet)
net=igraph.from.graphNEL(parseKGML2Graph(getKGMLurl("04630",organism = 'hsa')))

refnet=unlist2(kegg_refnet@network,recursive=FALSE)
refname=unlist2(kegg_refnet@name)
refgroup=names(refnet)

ref.net=list(
  networks=refname,
  annotation=data.frame(
    name=refname, description=refname, category=refgroup, rank=1:length(refnet)
    ),
  description="Reference networks from KEGG Pathway Database, updated by 2015/10/26",
  name="KEGG Signalling Pathway",
  create_time = "2015/10/1",
  create_user = "system",
  public = TRUE
  )

option=read.table(text=
"algorithm, parameters, default_value
  NFP, algorithm, 'APCLUSTER' #('APCLUSTER', 'GHOST','SPINAL','ISORANK')
  NFP, sim, 'geneontology' #('sequence','geneontology')
  NFP,nperm, 10 #(10~1000,integer)
  isorank, K, 10 #(10~30,integer)
  isorank, thresh, 1e-4#(1e-5 ~ 1e-3)
  isorank, alpha, 0.6#(0.2~1)
  isorank, beta,0.5 #(0.25~0.75)
  isorank, maxveclen, 1000000 #(200k~5000k,integer)
  ghost,hops, 4 #(2~6,integer)
  ghost,nneighbors,-1 #(-1~10,integer)
  ghost,searchiter,10 #(5~15,integer)
  ghost,ratio,8 #(4~10,integer)
  spinal,alpha,0.7 #(0~1)
  apcluster,q, NA#(0~1)
  apcluster,maxits,1000 #(100~10000,integer)
  apcluster,convits,100 #(50~150,integer)
  apcluster, lam, 0.9 #(0.5~1)
  apcluster, nonoise, FALSE #(TRUE, FALSE)
"
, header=TRUE,sep=",",strip.white=TRUE,stringsAsFactors=FALSE)

#calc_sim_score <- function(net, refnet,
#                           option = dataframe, homedir){

#	if(option[option$algorith==="NFP" & option$parameters=="algorithm",3] == "APCLUSTER"){
		nperm=as.integer(option[option$algorith=="NFP" & option$parameters=="nperm",3])
		nfp=NFP::calc_sim_score(net,kegg_refnet,nperm=nperm,plot=FALSE)
#    		return(list(networkfingerprint = nfp@standardize_score,
#              		netnode_cluster_index = nfp@cluster,
#              		refnetnode_cluster_index = refnetnode_cluster_index,
#              		score_of_cluster = score_of_cluster))
#  	}
#}


