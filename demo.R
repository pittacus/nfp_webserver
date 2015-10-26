if (!require('devtools'))  install.packages("devtools") 
if (!require('NFP')) devtools::install_github("yiluheihei/NFP")
require(NFP)
require(igraph)
require(KEGGgraph)
require(plyr)

data(kegg_refnet)
net=igraph.from.graphNEL(parseKGML2Graph(getKGMLurl("04630",organism = 'hsa')))

refnetwork=unlist2(kegg_refnet@network,recursive=FALSE)
refname=unlist(kegg_refnet@name, recursive=FALSE, use.names=FALSE)
refgroup=names(refnetwork)

refnet=list(
  network=refnetwork,
  annotation=data.frame(
    name=refname, description=refname, category=refgroup, rank=1:length(refname), stringsAsFactors=FALSE
    ),
  description="Reference networks from KEGG Pathway Database, updated by 2015/10/26",
  name="KEGG Signalling Pathway",
  create_time = "2015/10/1",
  update_time = "2015/10/1",
  create_user = "system",
  public = TRUE,
  organism = "hsa"
  )

option=read.table(text=
"algorithm, parameters, default_value
  NFP, algorithm, 'APCLUSTER' #('APCLUSTER', 'GHOST','SPINAL','ISORANK')
  NFP, sim, 'geneontology' #('sequence','geneontology')
  NFP,nperm, 10 #(10~1000,integer)
  ISORANK, K, 10 #(10~30,integer)
  ISORANK, thresh, 1e-4#(1e-5 ~ 1e-3)
  ISORANK, alpha, 0.6#(0.2~1)
  ISORANK, beta,0.5 #(0.25~0.75)
  ISORANK, maxveclen, 1000000 #(200k~5000k,integer)
  GHOST,hops, 4 #(2~6,integer)
  GHOST,nneighbors,-1 #(-1~10,integer)
  GHOST,searchiter,10 #(5~15,integer)
  GHOST,ratio,8 #(4~10,integer)
  SPINAL,alpha,0.7 #(0~1)
  APCLUSTER,q, NA#(0~1)
  APCLUSTER,maxits,1000 #(100~10000,integer)
  APCLUSTER,convits,100 #(50~150,integer)
  APCLUSTER, lam, 0.9 #(0.5~1)
  APCLUSTER, nonoise, FALSE #(TRUE, FALSE)
"
, header=TRUE,sep=",",strip.white=TRUE,stringsAsFactors=FALSE)


calc_sim_score_apcluster <- function(net, refnet, option, homedir){
  

  NFPref<-new("NFPRefnet",network = split(refnet$network,refnet$annotation$category), 
              name = split(refnet$annotation$name,refnet$annotation$category) , 
              group = unique(refnet$annotation$category),organism = "hsa")
  
  nperm=as.integer(option[option$algorith=="NFP" & option$parameters=="nperm",3])
  nfp=NFP::calc_sim_score(net,NFPref,nperm=nperm,plot=FALSE)
  
  node_cluster_index <- function(clusters){
    node_cluster_index <- c()
    node_cluster_names <- c()
    for (i in 1:length(clusters)){
      node_cluster_index[clusters[[i]]] <- i
      node_cluster_names[clusters[[i]]] <- names(clusters[[i]])
    }
    names(node_cluster_index) <- node_cluster_names
    return(node_cluster_index)
  }
  node_cluster_index_all <- llply(nfp@cluster, node_cluster_index)
  netnode_cluster_index <- llply(node_cluster_index_all, function(x){
    x[names(x) %in% gsub("hsa:","",names(V(net)))]
    })
  refnetnode_cluster_index <- llply(1:length(node_cluster_index_all), function(i){
    node_cluster_index_all[[i]][names(node_cluster_index_all[[1]]) %in% gsub("hsa:","",names(V(refnet$network[[i]])))]
  })
  names(refnetnode_cluster_index) <- names(node_cluster_index_all)
  
  
  return(list(networkfingerprint = nfp@standardized_score,
	  netnode_cluster_index = netnode_cluster_index,
    refnetnode_cluster_index = refnetnode_cluster_index,
    score_of_cluster = rnorm(length(nfp@cluster))
	  ))
}

calc_sim_score_ghost <- function(net, refnet, option, homedir){
  
}


calc_sim_score_spinal <- function(net, refnet, option, homedir){
  
}


calc_sim_score_isorank <- function(net, refnet, option, homedir){
  
}

calc_sim_score <- function(net, refnet,option, homedir){

	algorithm=option[option$algorith=="NFP" & option$parameters=="algorithm",3]
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

calc_sim_score(net,refnet,option, "~")
