if (!require('devtools'))  install.packages("devtools") 
if (!require('NFP')) devtools::install_github("yiluheihei/NFP")
require(NFP)
require(igraph)
require(KEGGgraph)
require(plyr)
require(dplyr)
require(AnnotationDbi)

data(kegg_refnet) # This is original kegg_refnet in a list format

refnetwork=unlist2(kegg_refnet@network,recursive=FALSE)
refname = unlist(kegg_refnet@name, recursive=FALSE, use.names=FALSE)
refname = gsub(" - Homo sapiens (human)","",refname,fixed = TRUE)
refname_underscore = gsub(" ","_",refname,fixed = TRUE)
refgroup=names(refnetwork)

# Definition example using kegg
kegg_refnet134 = list(
  network=refnetwork,
  annotation = data.frame( name = refname_underscore, 
                           description=refname,
                           category=refgroup, 
                           rank=1:length(refname), 
                           stringsAsFactors=FALSE
                           ),
  description = "Reference networks from KEGG Pathway Database, updated by 2015/10/26",
  name = "KEGG Signalling Pathway",
  create_time = "2015/10/1",
  create_user = "system",
  update_time = "2015/10/28",
  update_user = "system",
  public = TRUE,
  organism = "Homo sapiens"
  )

option= read.table(text=
 "algorithm, parameters, default_value
  NFP, algorithm, 'APCLUSTER' #('APCLUSTER', 'GHOST','SPINAL','ISORANK')
  NFP, sim, 'geneontology' #('sequence','geneontology')
  NFP,nperm, 2 #(10~1000,integer)
  ISORANK, K, 10 #(10~30,integer)
  ISORANK, thresh, 1e-4#(1e-5 ~ 1e-3)
  ISORANK, alpha, 0.6#(0.2~1)
  ISORANK, beta,0.5 #(0.25~0.75)
  ISORANK, maxveclen, 1000000 #(200k~5000k,integer)
  GHOST, hops, 4 #(2~6,integer)
  GHOST, nneighbors,-1 #(-1~10,integer)
  GHOST, searchiter,10 #(5~15,integer)
  GHOST, ratio,8 #(4~10,integer)
  SPINAL, alpha,0.7 #(0~1)
  APCLUSTER, q, NA#(0~1)
  APCLUSTER, maxits,1000 #(100~10000,integer)
  APCLUSTER, convits,100 #(50~150,integer)
  APCLUSTER, lam, 0.9 #(0.5~1)
  APCLUSTER, nonoise, FALSE #(TRUE, FALSE)
"
, header=TRUE,sep=",",strip.white=TRUE, stringsAsFactors=FALSE)


calc_sim_score_apcluster <- function(net, refnet, option, homedir){
  
  # split refnet by category
#   NFPref<-new("NFPRefnet",network = split(refnet$network,refnet$annotation$category), 
#               name = split(refnet$annotation$name,refnet$annotation$category) , 
#               group = unique(refnet$annotation$category),organism = "hsa")
#   rank2 = split(refnet$annotation$rank,refnet$annotation$category)
#   
# extract option
#   nperm=as.integer(option[option$algorith=="NFP" & option$parameters=="nperm",3])

#   nfp=NFP::calc_sim_score(net,NFPref,nperm=nperm,plot=FALSE)
#   
#   # cluster info convert to entrez id format
#   node_cluster_index <- function(clusters){
#     node_cluster_index <- c()
#     node_cluster_names <- c()
#     for (i in 1:length(clusters)){
#       node_cluster_index[clusters[[i]]] <- i
#       node_cluster_names[clusters[[i]]] <- names(clusters[[i]])
#     }
#     names(node_cluster_index) <- node_cluster_names
#     return(node_cluster_index)
#   }
#   node_cluster_index_all <- llply(nfp@cluster, node_cluster_index)
#   
#   netnode_cluster_index <- llply(node_cluster_index_all, function(x){
#     x[names(x) %in% gsub("hsa:","",names(V(net)))]
#     })
#   refnetnode_cluster_index <- llply(1:length(node_cluster_index_all), function(i){
#     node_cluster_index_all[[i]][names(node_cluster_index_all[[1]]) %in% 
#                                   gsub("hsa:","",names(V(refnet$network[[i]])))]
#   })
#   names(refnetnode_cluster_index) <- names(node_cluster_index_all)
#   
  # generate mock data 2015-10-29
  netfingerprint <- rnorm(length(refnet$network),mean = 0,sd=2)
  
  netnode_cluster_index <- llply(1:length(refnet$network),
                                 function(i) sample(1:5,length(V(net)),replace=TRUE))
  for (i in 1:length(refnet$network)){
    names(netnode_cluster_index[[i]]) <- names(V(net))
  }

  refnetnode_cluster_index <- llply(1:length(refnet$network),
                               function(i) sample(1:5,length(V(refnet$network[[i]])),replace=TRUE))
  for (i in 1:length(refnet$network)){
    names(refnetnode_cluster_index[[i]]) <- names(V(refnet$network[[i]]))
  }                               
  
  score_of_cluster <-llply(1:length(refnet$network),
                           function(i) runif(5,min =0,max =1))

  return(list(networkfingerprint = netfingerprint,
              netnode_cluster_index = netnode_cluster_index,
              refnetnode_cluster_index = refnetnode_cluster_index,
              score_of_cluster = score_of_cluster)
         )
}

calc_sim_score_ghost <- function(net, refnet, option, homedir){
  
}


calc_sim_score_spinal <- function(net, refnet, option, homedir){
  
}


calc_sim_score_isorank <- function(net, refnet, option, homedir){
  
}

calc_sim_score_main <- function(net, refnet,option, homedir){

	algorithm=option[option$algorithm=="NFP" & option$parameters=="algorithm",3]
  if(algorithm == "APCLUSTER")
	  ret=calc_sim_score_apcluster(net, refnet,option, homedir)
  else if (algorithm == "GHOST")
    ret=calc_sim_score_ghost(net, refnet,option, homedir)
	else if (algorithm == "SPINAL")
	  ret=calc_sim_score_spinal(net, refnet,option, homedir)
	else if (algorithm == "ISORANK")
	  ret=calc_sim_score_isorank(net, refnet,option, homedir)
	else return("Undefined algorithm")
	
  return(ret)	
}

refnet <- kegg_refnet134
# save(kegg_refnet134,file ="kegg_refnet134.rdata")
# This is original kegg_refnet in a list format
load("kegg_refnet134.rdata") 
load("nci_signaling_refnet.rdata")

# Net processing step remove hsa: from original igraph.
T1DM <- igraph.from.graphNEL(parseKGML2Graph(getKGMLurl("04930",organism = 'hsa')))                           
# MAPK <- igraph.from.graphNEL(parseKGML2Graph(getKGMLurl("04010",organism = 'hsa')))                           
# PanCancer <- igraph.from.graphNEL(parseKGML2Graph(getKGMLurl("05212",organism = 'hsa')))
V(T1DM)$name <- gsub("hsa:","",V(T1DM)$name)
          
test <- calc_sim_score_main(T1DM,nci_signaling_refnet,option, "~")

# Take the following two network as example of visualization
plot(nci_signaling_refnet$network[[2]])
plot(nci_signaling_refnet$network[[23]])
net1 <- nci_signaling_refnet$network[[2]]
net2 <- nci_signaling_refnet$network[[23]]
# pairwise alignment: net1 is ready as igraph; 
# net2 shoud be prepare as a refnet list
net2_refnet = list(
  network = list(net2),
  annotation = data.frame( name = names(net2), 
                           stringsAsFactors=FALSE
  ),
  public = FALSE,
  organism = "Homo sapiens"
)
test2 <- calc_sim_score_main(net1,net2_refnet,option, "~")
