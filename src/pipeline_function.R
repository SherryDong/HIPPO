library(RColorBrewer)
###########################################
#' Prepare dataset for HIPPO calculation. 
#'
#' \code{HIPPO.prepare_dataset} Read in read-region content matrix, 
#' reference information files, and output required dataset for calculation
#'
#' @param mat data.frame, original read-region content matrix.
#' @param mod_region data.frame, region position file in bed format. 
#' @param mod_info data.frame, position transfer file. 
#' @param int_pos character, intrested positions in genome position. 
#' 
#' @return A list contain: all_mat, group_info, group_name and pre_define. 
#'
#' @examples
#' \dontrun{
#'  load('demo/HIPPO_demo_20201204.RData')
#'  res <- HIPPO.prepare_dataset(mat=demo_mat,mod_region=mod_region,
#'         mod_info=mod_info,int_pos=int_pos)
#'  all_mat <- res$all_mat
#'  group_info <- res$group_info
#'  group_name <- res$group_name
#'  pre_define <- res$pre_define
#' }
#' @export
HIPPO.prepare_dataset <- function(mat=NULL,mod_region=NULL,mod_info=NULL,int_pos=NULL){
  int_pos <- as.numeric(unlist(strsplit(int_pos,',')))
  int_region <- sprintf('TID_%s',int_pos)
  r1$colname <- r1$V4
  r1$colname <- gsub('-','.',r1$colname)
  rownames(r1) <- r1$colname
  r1 <- r1[which(r1$colname %in% colnames(d1)),]
  d2 <- d1[,r1$colname]
  int_region_len <- r1[int_region,3]-r1[int_region,2]
  ref_nt <- do.call(rbind,lapply(1:length(int_region),function(x){
    xp <- int_pos[x]
    xl <- int_region_len[x]
    x1 <- f1[which(f1$genome_pos %in% xp:(xp+xl-1)),]
    x2 <- c(int_region[x],paste(x1[,4],collapse=''),paste(x1[,5],collapse=''),x1[1,3])
    names(x2) <- c('region',colnames(f1)[c(4,5,3)])
    x2
  }))
  ref_nt <- as.data.frame(ref_nt,stringsAsFactors=FALSE)
  all_mat <- d2[,int_region]
  group_info <- ref_nt[,2:(ncol(ref_nt)-1)]
  group_name <- ref_nt[,ncol(ref_nt)]
  pre_define <- colorRampPalette(brewer.pal(8,'Spectral'))(length(group_info))
  names(pre_define) <- names(group_info)
  return(list(all_mat=all_mat,group_info=group_info,group_name=group_name,pre_define=pre_define))
}
###########################################
## funcitons for HIPPO.plot_readsComponent(); HIPPO.plot_adjacentCombination()
get.class.color <- function(x,use_color=NULL,pre_define=NULL) {
  if(is.null(pre_define)==FALSE & is.null(names(pre_define))==TRUE){
    message('No class name for the color vector, please check and re-try !');return(FALSE);
  }
  if(is.null(use_color)==TRUE){
    use_color <- brewer.pal(9, 'Set1')
  }
  if (base::length(base::intersect(x, names(pre_define))) == 0) {
    w1 <- base::length(base::unique(x))
    if(w1 < length(use_color)){
      cc2 <- use_color[1:w1]
    }else{
      cc2 <- grDevices::colorRampPalette(use_color)(base::length(base::unique(x)))
    }
    names(cc2) <- base::unique(x)
    cc2 <- cc2[x]
  } else{
    x1 <- base::unique(x)
    x2 <- base::setdiff(x1, names(pre_define))
    cc1 <- NULL
    w1 <- base::length(x2)
    if (w1 > 0) {
      if(w1 < length(use_color)){
        cc1 <- use_color[1:w1]
      }else{
        cc1 <- grDevices::colorRampPalette(use_color)(w1)
      }
      names(cc1) <- x2
    }
    cc2 <- c(pre_define, cc1)
    cc2 <- cc2[x]
  }
  return(cc2)
}
plot_new <- function(xlim=c(),ylim=c()){
  plot(1,xlim=xlim,ylim=ylim,
       xaxt='n',yaxt='n',xlab='',ylab='',bty='n',col='white')
}
curve_sankey <- function (x0, x1, y0, y1,nsteps = 1000){
  xx <- seq(-pi/2, pi/2, length.out = nsteps)
  yy <- y0 + (y1 - y0) * (sin(xx) + 1)/2
  xx <- seq(x0, x1, length.out = nsteps)
  list(x=xx,y=yy)
}
sankey_polygon <- function(Fx0, Fx1, Fy0, Fy1,Sx0, Sx1, Sy0, Sy1,
                           nsteps = 100,col='red',border=NA,...){
  curve_pos1 <- curve_sankey(Fx0,Fx1,Fy0,Fy1)
  curve_pos2 <- curve_sankey(Sx0,Sx1,Sy0,Sy1)
  polygon(x=c(curve_pos1$x,rev(curve_pos2$x)),
          y=c(curve_pos1$y,rev(curve_pos2$y)),
          col=col,border=border,xpd=TRUE,...)
}

###########################################
#' Plot read positional-components
#'
#' \code{HIPPO.plot_readsComponent} plot reads component for each position.
#'
#' @param all_mat data.frame, read-region content matrix with interested positions.
#' @param remove_char character, character for none matching sequence content, default is ".". 
#' @param group_info data.frame, sequence context at each position for true and pseudogene. 
#' @param group_name character, genome position for each interested positions. 
#'
#' @examples
#' \dontrun{
#' HIPPO.plot_readsComponent(all_mat,group_info=group_info,
#'                           remove_char='.',
#'                           group_name=group_name)
#' }
#' @export
HIPPO.plot_readsComponent <- function(all_mat,group_info=NULL,
                        remove_char='.',
                        group_name=NULL){
  n <- ncol(all_mat); 
  all_r <- paste(rep(remove_char,n),collapse = '_')
  tmp1 <- apply(all_mat,1,function(x)paste(x,collapse = '_'))
  tmp2 <- sort(table(tmp1),decreasing = T)
  tmp2 <- tmp2[which(names(tmp2)!=all_r)]
  m <- length(group_info[[1]])
  par(mar=c(2,2,3,2))
  plot_new(xlim=c(0,m),ylim=c(-0.5,1.25))
  text(x=1:m,y=0,group_info[[1]],xpd=TRUE,cex=0.9)
  text(x=1:m,y=1,group_info[[2]],xpd=TRUE,cex=0.9)
  text(x=0,y=0,names(group_info)[1],pos=2,xpd=TRUE,cex=0.7)
  text(x=0,y=1,names(group_info)[2],pos=2,xpd=TRUE,cex=0.7)
  text(x=1:m,y=1.2,group_name,xpd=TRUE,adj=0,srt=90,cex=0.5)
  cc <- rev(colorRampPalette(c('black','red'))(length(tmp2)))
  dyy <- 0.8/length(tmp2)
  for(j in 1:length(tmp2)){
    x1 <- unlist(strsplit(names(tmp2)[j],'_'))
    w0 <- which(x1!='.')
    if(length(w0)>1){
      for(i in 1:(length(w0)-1)){
        w1 <- which(group_info[w0[i],] == x1[w0[i]])
        w2 <- which(group_info[w0[i+1],] == x1[w0[i+1]])
        if(length(w1)>0 & length(w2)>0){
          segments(x0=w0[i],x1=w0[i+1],
                   y0=w1-1-dyy*j,y1=w2-1-dyy*j,
                   col=adjustcolor(cc[j],0.5),
                   lwd=3*tmp2[j]/max(tmp2)+1,xpd=TRUE,
                   lty=ifelse(w0[i+1]-w0[i]==1,1,2))
          points(w0[i],w1-1-dyy*j,xpd=TRUE,cex=0.5)
          points(w0[i+1],w2-1-dyy*j,xpd=TRUE,cex=0.5)
        }
      }
    }
  }
}

###########################################
#' Plot read positional-components
#'
#' \code{HIPPO.plot_adjacentCombination} plot reads component for each position.
#'
#' @param all_mat data.frame, read-region content matrix with interested positions.
#' @param top_n, numeric, number of top sequence context at each position. Default is 2.
#' @param each_width, numeric, bar width for each position. Default is 0.25. 
#' @param remove_char character, character for none matching sequence content, default is ".". 
#' @param count_link_thre numeric, threshold for the percentage of each intersected sequence context to the total number of valid sequence.   
#' Default is 0.1.
#' @param part_thre numeric, threshold for the percentage of each intersected sequence context to the total number of start-side valid sequence. 
#' Default is 0.1.
#' @param group_info data.frame, sequence context at each position for true and pseudogene. 
#' @param only_use_group logical, if true, sequence not match group_info will not be plotted. Default is FALSE.
#' @param group_name character, genome position for each interested positions. 
#' @param strategy character, choose from 'percentage' and 'count'. 
#' If "percentage", the height for each position is normalized to the total number of valid sequence. 
#' If "count", the height for each position is the true read number for each sequence context. 
#' Default is 'percentage'. 
#'
#' @examples
#' \dontrun{
#' HIPPO.plot_adjacentCombination(all_mat,remove_char='.',
#'                                group_info=group_info,
#'                                group_name=group_name,only_use_group=FALSE,
#'                                count_link_thre=0,part_thre=0,
#'                                top_n=2)
#' }
#' @export
HIPPO.plot_adjacentCombination <- function(all_mat,top_n=3,each_width=0.25,
                       remove_char='.',
                       count_link_thre=0.1,part_thre=0.1,
                       group_info=NULL,
                       only_use_group=TRUE,
                       group_name=NULL,strategy='percentage'){
  # strategy = 'count'
  n <- ncol(all_mat)
  if(is.null(group_info)==FALSE){
    par(mar=c(5,3,ifelse(is.null(group_name)==TRUE,2,5),5))
  }else{
    par(mar=c(5,3,ifelse(is.null(group_name)==TRUE,2,5),2))
  }
  plot_new(xlim=c(1.1,n+0.25),ylim=c(0,1))
  hap_start <- list(); hap_end <- list(); all_top <- list();
  all_valid <- list(); all_topc <- list()
  if(strategy=='count'){
    w1 <- apply(all_mat,2,function(x)length(which(x!=remove_char)))
    max_w1 <- max(w1); scale_w1 <- w1/max_w1;
  }else{
    scale_w1 <- rep(1,n)
  }
  for(i in 1:n){
    w1 <- which(all_mat[,i]!=remove_char)
    all_valid[[i]] <- w1
    tmp1 <- sort(table(all_mat[w1,i]),decreasing = T)
    tmp1 <- tmp1/sum(tmp1)
    if(is.null(group_info)==FALSE){
      w1 <- unlist(lapply(group_info,function(x)x[[i]]))
      if(only_use_group==TRUE){
        tmp1 <- tmp1[c(w1)]
        tmp1c <- names(group_info)
      }else{
        tmp1 <- tmp1[c(w1,setdiff(names(tmp1),w1))]
        tmp1c <- c(names(group_info),setdiff(names(tmp1),w1))
      }
    }
    if(length(tmp1)>top_n){
      tmp1 <- c(tmp1[1:top_n],other=sum(tmp1)-sum(tmp1[1:top_n]))
      top_k <- top_n+1
    }else{
      top_k <- length(tmp1)
    }
    tmp1 <- tmp1*scale_w1[i]
    tmp2 <- cumsum(rev(tmp1))
    tmp31 <- c(0,tmp2[1:(top_k-1)]); tmp32 <- c(tmp2[1:top_k])
    hap_start[[i]] <- rev(tmp31);
    hap_end[[i]]   <- rev(tmp32);
    all_top[[i]] <- names(tmp1);
    cc <- get.class.color(tmp1c,use_color=brewer.pal(8,'Pastel2'),
                                  pre_define = pre_define)
    names(cc) <- names(tmp1)
    all_topc[[i]] <- cc;
    for(j in 1:top_k){
      rect(xleft=i,xright = i+each_width,
           ybottom = tmp31[j],ytop=tmp32[j],col=cc[names(tmp32[j])],xpd=TRUE)
    }
    legend(x=i,y=-0.05,names(tmp1),cc,
           cex=0.5,xpd=T,border = NA,bty='n',yjust = 1)
  }
  if(is.null(group_name)==FALSE){
    text(1:n+each_width/2,1.05,group_name,srt=90,adj = 0,xpd=TRUE,cex=0.7)
  }
  mat1 <- as.matrix(all_mat)
  for(i in 1:(n-1)){
    start_left <- rep(0,length(all_top[[i]]))
    end_right <- rep(0,length(all_top[[i+1]]))
    for(j1 in 1:length(all_top[[i]])){
      for(j2 in 1:length(all_top[[i+1]])){
        inter_count <- length(which(mat1[,i]==all_top[[i]][j1] & 
                                      mat1[,i+1]==all_top[[i+1]][j2]))
        count_link <- inter_count/length(intersect(all_valid[[i]],all_valid[[i+1]]))
        per1_part <- inter_count/length(which(mat1[,i]==all_top[[i]][j1]))
        per1 <- inter_count/length(all_valid[[i]])*scale_w1[i]
        per2 <- inter_count/length(all_valid[[i+1]])*scale_w1[i+1]
        if(count_link > count_link_thre & per1_part>part_thre){
          Fy0 <- hap_end[[i]][j1]-start_left[j1]
          Fy1 <- hap_end[[i+1]][j2]-end_right[j2]
          Sy0 <- Fy0-per1
          Sy1 <- Fy1-per2
          sankey_polygon(i+each_width,i+1,Fy0,Fy1,i+each_width,i+1,Sy0,Sy1,
                         col = adjustcolor(all_topc[[i]][all_top[[i]][j1]],0.3),
                         border = 1,lwd=0.4)
          start_left[j1] <- start_left[j1]+per1
          end_right[j2] <- end_right[j2]+per2
        }
      }
    }
  }
  if(is.null(group_info)==FALSE){
    legend(x=n+0.5,y=0.5,names(group_info),
           fill=cc[1:length(names(group_info))],xpd=T,cex=0.75,
           border = NA,bty='n')
  }
  if(strategy == 'count'){
    ss <- round(seq(1,max_w1,length.out = 11))
    axis(side=2,at=c(0:10)/10,labels = ss,
         las=2,cex.axis=0.7,cex=0.7)
  }else{
    axis(side=2,at=c(0:10)/10,labels = sprintf('%s%s',c(0:10)*10,'%'),
         las=2,cex.axis=0.7,cex=0.7)
  }
}

###########################################
#' Summarize imputation results. 
#'
#' \code{HIPPO.summ_readsComponent} plot reads component for each position.
#'
#' @param all_mat data.frame, read-region content matrix with interested positions.
#' @param group_info data.frame, sequence context at each position for true and pseudogene. 
#' @param group_name character, genome position for each interested positions. 
#'
#' @examples
#' \dontrun{
#'res_readsComponent <- HIPPO.summ_readsComponent(all_mat,group_info,group_name)
#' @export
HIPPO.summ_readsComponent <- function(all_mat,group_info,group_name){
  d4 <- sapply(1:ncol(all_mat),function(i){
    x1 <- all_mat[,i]
    w1 <- which(x1==group_info[i,1])
    w2 <- which(x1==group_info[i,2])
    x1[w1] <- sprintf("%s:%s",colnames(group_info)[1],x1[w1])
    x1[w2] <- sprintf("%s:%s",colnames(group_info)[2],x1[w2])
    x1
  })
  tmp1 <- apply(d4,1,function(x)paste(x,collapse = '_'))
  tmp2 <- sort(table(tmp1),decreasing = T)
  #tmp2 <- tmp2[which(tmp2>=nrow(d3)*0.001)]
  all_res <- list()
  for(j in 1:length(tmp2)){
    x1 <- unlist(strsplit(names(tmp2)[j],'_'))
    w0 <- which(x1!='.')
    if(length(w0)>1){
      tmp11 <- apply(d4[,w0],1,function(x)paste(x,collapse = '_'))
      r_count <- length(which(tmp11==paste(x1[w0],collapse = '_')))
      all_res[[j]] <- c(x1,r_count)
    }
  }
  all_res <- do.call(rbind,all_res)
  all_res <- all_res[order(as.numeric(all_res[,ncol(all_res)]),decreasing = T),]
  all_res <- rbind(cbind(t(group_info),'.'),all_res)
  all_res <- data.frame(all_res,stringsAsFactors = F)
  colnames(all_res)[-ncol(all_res)] <- group_name
  colnames(all_res)[ncol(all_res)] <- 'Support_Read_Number'
  rownames(all_res) <- c(colnames(group_info),
                         sprintf('Top%s',1:(nrow(all_res)-ncol(group_info)))) 
  all_res
}

###########################################
## funcitons for HIPPO.impute_Haplotype()
get_valid <- function(all_mat,x,remove_character='.'){
  all_mat[which(all_mat[,x]!=remove_character),]
}
get_from <- function(group_info,each_t,x){
  x1 <- unlist(lapply(group_info,function(xx)xx[each_t]))
  x2 <- names(group_info)[which(x1==x)]
  if(length(which(x1==x))==0) x2 <- 'other'
  x2
}
compare_twoSeq <- function(x1,x2,min_overlap=2,remove_character='.'){ # x1-ref, x2-obs
  x1 <- unlist(strsplit(x1,'_'))
  x2 <- unlist(strsplit(x2,'_'))
  u1 <- which(x1!=remove_character); u2 <- which(x2!=remove_character)
  x1 <- sprintf('P%s:%s',1:length(x1),x1)[u1]
  x2 <- sprintf('P%s:%s',1:length(x2),x2)[u2]
  w0 <- intersect(x1,x2)
  w12 <- setdiff(x1,x2)
  w21 <- setdiff(x2,x1)
  if(length(w21)==0 & length(w0)>=min_overlap){
    return(1)
  }else{
    return(0)
  }
}
compare_twoSeq_count <- function(x1,x2,min_overlap=2,remove_character='.'){ # x1-ref, x2-obs
  x1 <- unlist(strsplit(x1,'_'))
  x2 <- unlist(strsplit(x2,'_'))
  u1 <- which(x1!=remove_character); u2 <- which(x2!=remove_character)
  x1 <- sprintf('P%s:%s',1:length(x1),x1)[u1]
  x2 <- sprintf('P%s:%s',1:length(x2),x2)[u2]
  w0 <- intersect(x1,x2)
  w12 <- setdiff(x1,x2)
  w21 <- setdiff(x2,x1)
  if(length(w21)==0 & length(w0)>=min_overlap){
    return(length(w0))
  }else{
    return(0)
  }
}



###########################################
#' Impute for haplotype
#'
#' \code{HIPPO.impute_Haplotype} perform imputation for haplotypes. 
#'
#' @param all_mat data.frame, read-region content matrix with interested positions.
#' @param group_info data.frame, sequence context at each position for true and pseudogene. 
#' @param top_each, numeric, number of top haplotypes to leave at each step. Default is 3.
#' @param use_site, character, if not NULL, only impute for the input sites. Default is NULL. 
#' @param min_overlap, numeric, minimun matching_score for a read assigned to a haplotype. Default is 2.  
#' @param check_thre, numeric, threshold for the remaining possibility when only consider top_each haplotypes at each step. 
#' If the remaining possibility is higher than check_thre, we will include more haplotypes until the remaining is possibility less than check_thre. 
#' Default is 0.05.
#' @param adjust_per logical, whether to perform frequencey iteration step to adjust haplotype frquency.Default is TRUE. 
#' @param perm_k numeric, maximum number of permutation, default is 5000. 
#' @param perm_strategy character, choose from 'sample' and 'frequency'. 
#' If "sample", the read’s haplotype belonging is random-sampled with matching_score*haplotype frequency as probability weights.
#' If "frequency", read’s haplotype belonging is distributed according to the weight value.
#' Default is 'sample'. 
#' @param min_frequency numeric, minimun frequency for the haplotype to output. Default is 0.01. 
#' @param remove_char character, character for none matching sequence content, default is ".". 
#'
#' @examples
#' \dontrun{
#' res_trace <- HIPPO.impute_Haplotype(all_mat=all_mat,
#'                                     group_info=group_info,
#'                                     top_each=3)
#' }
#' @export
HIPPO.impute_Haplotype <- function(all_mat,group_info,top_each=3,use_site=NULL,
                         min_overlap=2,check_thre=0.05,
                         adjust_per=TRUE,perm_k=5000,perm_strategy='sample',min_frequency=0.01,
                         remove_char='.'){
  w1 <- apply(all_mat,1,function(x)length(setdiff(unique(x),remove_char)))
  all_mat <- all_mat[which(w1>0),]
  if(is.null(use_site)==FALSE){
    all_mat <- all_mat[,use_site]
    group_info <- lapply(group_info,function(x)x[use_site])
  }
  ## each_t: for each Site
  ## S1, S2, ... Smax_t
  max_t <- length(group_info[[1]]) ## each_t
  result_t <- list()
  check_prob <- c()
  ##
  each_t <- 1
  use_mat <- get_valid(all_mat,each_t)
  tmp1 <- table(use_mat[,each_t])/nrow(use_mat)
  # use top
  o1 <- order(tmp1,decreasing = T)[1:min(top_each,length(tmp1))]
  check_prob[[each_t]] <- sum(tmp1)-sum(tmp1[o1])
  tmp1 <- tmp1[o1]
  # normalize
  tmp1 <- tmp1/sum(tmp1)
  res <- tmp1
  result_t[[each_t]] <- list(seq=lapply(names(res),c),
                             trace=lapply(names(res),function(x)get_from(group_info,each_t,x)),
                             record=lapply(1:length(res),function(x)x),
                             prob=as.numeric(res),
                             check_prob=check_prob[[each_t]])
  ##
  for(each_t in 2:max_t){
    use_mat <- get_valid(all_mat,each_t)
    use_res <- result_t[[each_t-1]]
    tmp1 <- use_res$seq
    tmp2 <- sort(table(use_mat[,each_t]),decreasing = T)
    tmp2 <- names(tmp2)
    res <- list(seq=list(),trace=list(),record=list(),prob=c(),check_prob=c())
    count <- 0
    for(i in 1:length(tmp1)){ ## 1~each_t-1
      this_count <- c()
      for(j in 1:length(tmp2)){ ## each_t
        count <- count + 1
        this_count <- c(this_count,count)
        ref <- c(use_res$seq[[i]],tmp2[j])
        obs <- apply(use_mat[,1:each_t,drop=F],1,function(x){
          compare_twoSeq(ref,x,min_overlap=min(min_overlap,each_t))
        })
        if(length(which(obs==1))>0){
          obs_pv  <- length(which(obs==1))/length(obs)
        }else{
          obs_pv <- 0
        }
        post_pv <- obs_pv
        res$seq[[count]] <- ref
        res$prob <- c(res$prob,post_pv)
        res$trace[[count]] <- c(use_res$trace[[i]],
                                get_from(group_info,each_t,tmp2[j]))
        res$record[[count]] <- use_res$record[[i]]
      }
      if(sum(res$prob[this_count])>0){
        ###### normalize for each 1-(each_t-1)
        res$prob[this_count] <- res$prob[this_count]/sum(res$prob[this_count]) # normalize
        res$prob[this_count] <- res$prob[this_count]*use_res$prob[i] # maximum --> use_res$prob[i]
        ######
      }
    }
		res$prob <- res$prob/sum(res$prob)
    ori_o <- order(res$prob,decreasing = T)
    prob_o <- 1-cumsum(res$prob[ori_o])
		if(length(which(prob_o<=check_thre))>0){
    	w1 <- min(which(prob_o<=check_thre))
		}else{
			w1 <- length(prob_o)
		}
    o1 <- ori_o[1:max(w1,min(top_each,length(res$prob)))]
    check_prob[[each_t]] <- sum(res$prob)-sum(res$prob[o1])
    res$seq <- res$seq[o1]
    res$prob <- res$prob[o1]
    res$trace <- res$trace[o1]
    res$record <- res$record[o1]
    res$check_prob <- check_prob[[each_t]]
    # normalize
    res$prob <- res$prob/sum(res$prob)
    for(i in 1:length(res$record)){
      res$record[[i]] <- c(res$record[[i]],i)
    }
    result_t[[each_t]] <- res
  }
  ## 迭代
  if(adjust_per == TRUE){
    ##### re-calculate to refine probability
    ## record read used for each haplotype
    final_t <- result_t[[max_t]]
    read2hap <- matrix(0,nrow(all_mat),ncol=length(final_t$seq)); 
    rownames(read2hap) <- rownames(all_mat)
    for(i in 1:length(final_t$seq)){
      ref <- final_t$seq[[i]]
      obs <- apply(all_mat,1,function(x){
        compare_twoSeq_count(ref,x,min_overlap=min(min_overlap,each_t))
      })
      obs1 <- obs[which(obs>0)]
      read2hap[names(obs1),i] <- obs1
    }
    ##
    read2hap <- read2hap[which(rowSums(read2hap)>0),]
    new_prob <- colSums(read2hap)/nrow(read2hap); 
    new_prob <- new_prob/sum(new_prob)
    all_prob <- list()
    for(k in 1:perm_k){
      read2hap_new <- t(apply(read2hap,1,function(x){
        if(perm_strategy=='sample'){
          x1 <- sample(1:length(x),size=1,prob=x*new_prob)
          x[x1] <- 1; x[setdiff(1:length(x),x1)] <- 0; return(x)
        }
        if(perm_strategy=='frequency'){
          x <- x*new_prob; x <- x/sum(x);
        }
      }))
      new_prob1 <- colSums(read2hap_new)/nrow(read2hap_new)
      s1 <- sum(abs(new_prob1-new_prob))
      if(sum(abs(new_prob1-new_prob))<1e-30){
        break
      }else{
        new_prob <- new_prob1
      }
      all_prob[[k]] <- c(new_prob,s1)
    }
    all_prob <- do.call(rbind,all_prob) ## just for test
    new_prob[which(new_prob<min_frequency)] <- 0; ## remove lowest frequency
    new_prob <- new_prob/sum(new_prob); 
    w1 <- which(new_prob>0); w1 <- w1[order(new_prob[w1],decreasing = T)]
    final_t$prob <- new_prob[w1]
    final_t$record <- final_t$record[w1]
    final_t$trace <- final_t$trace[w1]
    for(i in 1:length(w1)){
      final_t$record[[i]][max_t] <- i
    }
    final_t$seq <- final_t$seq[w1]
    result_t[[max_t]] <- final_t
  }
  ## output
  return(result_t)
}


###########################################
#' Plot for imputation results
#'
#' \code{HIPPO.plot_imputeHaplotype} plot the imputation results. 
#'
#' @param result_t list, output object from HIPPO.impute_Haplotype(). 
#' @param prob_thre, numeric, threshold for the haplotype frequency to show. Default is 0.01.
#' @param group_name character, genome position for each interested positions. 
#' @param consider_site_per logical, whether to display each site according to the possibility at each step. 
#' Default is FALSE.
#'
#' @examples
#' \dontrun{
#'       HIPPO.plot_imputeHaplotype(res_trace,prob_thre=0.01,
#'                           group_name=group_name)
#' }
#' @export
HIPPO.plot_imputeHaplotype <- function(result_t,prob_thre=0.01,group_name=NULL,
                              consider_site_per=FALSE){
  par(mar=c(2,2,5,8))
  if(consider_site_per==TRUE){
    use_result_t <- lapply(result_t, function(x){
      w1 <- which(x$prob>=prob_thre)
      list(seq=x$seq[w1],trace=x$trace[w1],prob=x$prob[w1],record=x$record[w1])
    })
    max_x <- length(use_result_t)
    max_y <- max(unlist(lapply(use_result_t,function(x)length(x$prob))))
    all_m <- unique(unlist(lapply(use_result_t,function(x)unlist(x$trace))))
    cc <- get.class.color(all_m,use_color=brewer.pal(8,'Pastel2'),
                          pre_define = pre_define)
    # points
    plot_new(xlim = c(1,max_x),ylim=c(1,max_y))
    for(i in 2:max_x){
      each_r <- use_result_t[[i]]
      each_t <- each_r$trace
      each_d <- each_r$record
      for(j in 1:length(each_t)){
        for(j1 in 1:(i-1)){
          use_t <- each_t[[j]][j1:(j1+1)]
          segments(x0=j1,x1=j1+1,
                   y0=1+max_y-each_d[[j]][j1],
                   y1=1+max_y-each_d[[j]][j1+1],
                   lwd=3+5*each_r$prob[j],
                   col=adjustcolor('dark grey',0.5))
        }
      }
    }
    ##
    for(i in 1:max_x){
      each_r <- use_result_t[[i]]
      each_b <- each_r$prob
      each_p <- unlist(lapply(each_r$trace,function(x)x[i]))
      each_n <- unlist(lapply(each_r$seq,function(x)x[i]))
      yy <- 1+max_y-(1:length(each_p))
      points(x=rep(i,length.out=length(yy)),y=yy,
             col=cc[each_p],pch=16,cex=each_b+2,xpd=TRUE)
      text(x=i,y=yy,each_n,xpd=TRUE,adj=0.5)
    }
    legend(x=max_x+max_x/5,y=max_y/2,all_m,
           fill=cc[all_m],xpd=T,cex=0.75,
           border = NA,bty='n',xjust=0,yjust=0.5)
    text(x=1:max_x,y=max_y+max_y/10,group_name,srt=90,adj=0,cex=0.7,xpd=TRUE)
    text(x=i,y=yy,sprintf('%s%s',round(each_b*10000)/100,'%'),xpd=TRUE,pos=4)
  }else{
    max_x <- length(result_t)
    final_t <- result_t[[max_x]]
    w1 <- which(final_t$prob>=prob_thre)
    final_t_mat <- do.call(rbind,final_t$seq)[w1,]
    final_s_mat <- do.call(rbind,final_t$trace)[w1,]
    max_x <- ncol(final_t_mat)
    max_y <- nrow(final_t_mat)
    all_m <- unique(unlist(final_t$trace[w1]))
    cc <- get.class.color(all_m,use_color=brewer.pal(8,'Pastel2'),
                          pre_define = pre_define)
    # points
    plot_new(xlim = c(1,max_x),ylim=c(1,max_y))
    for(i in 1:max_x){
      each_p <- final_s_mat[,i]
      each_n <- final_t_mat[,i]
      yy <- 1+max_y-c(1:length(each_p))
      points(x=rep(i,length.out=length(yy)),y=yy,
             col=cc[each_p],pch=16,cex=2,xpd=TRUE)
      text(x=i,y=yy,each_n,xpd=TRUE,adj=0.5)
    }   
    each_b <- final_t$prob[w1]
    for(i in 1:max_y){
      segments(x0=1,x1=max_x,
             y0=i,y1=i,
             lwd=3+5*each_b[i],
             col=adjustcolor('dark grey',0.5),xpd=T)
    }
    legend(x=max_x+max_x/5,y=max_y/2,all_m,
           fill=cc[all_m],xpd=T,cex=0.75,
           border = NA,bty='n',xjust=0,yjust=0.5)
    text(x=1:max_x,y=max_y+max_y/10,group_name,srt=90,adj=0,cex=0.7,xpd=TRUE)
    text(x=max_x,y=max_y:1,
         sprintf('%s%s',round(each_b*10000)/100,'%'),xpd=TRUE,pos=4)
  }
  ##
}

###########################################
#' Calculate allele frequency for each site. 
#'
#' \code{HIPPO.cal_siteProb} calculate allele frequency for each site.
#'
#' @param final_t list, one list object from HIPPO.impute_Haplotype(). 
#'
#' @examples
#' \dontrun{
#'   res_finalHaplotype <- res_trace[[max(length(res_trace))]]
#'   site_prob <- HIPPO.cal_siteProb(res_finalHaplotype)
#' }
#' @export
HIPPO.cal_siteProb <- function(final_t){
  n <- length(final_t$seq[[1]])
  prob <- final_t$prob
  s1 <- do.call(rbind,final_t$seq)
  t1 <- do.call(rbind,final_t$trace)
  all_r <- list()
  for(i in 1:n){
    r2 <- aggregate(prob,list(t1[,i]),sum)
    rownames(r2) <- r2$Group.1
    all_r[[i]] <- r2
  }
  all_n <- unique(unlist(final_t$trace))
  all_r <- lapply(all_r,function(x){
    x1 <- x$x; names(x1) <- x$Group.1; x1[all_n]
  })
  all_r <- do.call(cbind,all_r)
  rownames(all_r) <- all_n
  all_r <- t(t(all_r)/colSums(all_r,na.rm=T))
  all_r
}

###########################################
#' Output statistics to text files. 
#'
#' \code{HIPPO.output_result} output result statistics to text files. 
#' 
#' @param res_readsComponent list, object from HIPPO.summ_readsComponent(). 
#' @param res_finalHaplotype list, one list object from HIPPO.impute_Haplotype(). 
#' @param group_name character, genome position for each interested positions. 
#' @param output_file character, output file name. 
#'
#' @examples
#' \dontrun{
#' res_readsComponent <- HIPPO.summ_readsComponent(all_mat,group_info,group_name)
#' res_finalHaplotype <- res_trace[[max(length(res_trace))]]
#' HIPPO.output_result(res_readsComponent=res_readsComponent,
#'                     res_finalHaplotype=res_finalHaplotype,
#'                     group_name=group_name,
#'                     output_file=output_file)
#' }
#' @export
HIPPO.output_result <- function(res_readsComponent,
                                 res_finalHaplotype,group_name,
                                 output_file){
  ## output reads component results
  write.table(res_readsComponent,file=gsub(".txt$","_readsComponent.txt",output_file),
              row.names = T,col.names = T,quote = F,sep='\t')
  ## output site probability
  final_res_trace <- do.call(rbind,lapply(1:length(res_finalHaplotype$seq),function(i){
    cbind(rbind(res_finalHaplotype$seq[[i]],res_finalHaplotype$trace[[i]]),
          signif(res_finalHaplotype$prob[[i]],3))
  }))
  colnames(final_res_trace)[ncol(final_res_trace)] <- 'Support_Percentage'
  colnames(final_res_trace)[1:(ncol(final_res_trace)-1)] <- group_name
  write.table(final_res_trace,file=gsub(".txt$","_haplotype.txt",output_file),
              row.names = T,col.names = T,quote = F,sep='\t')
  ## output site probability
  site_prob <- HIPPO.cal_siteProb(res_finalHaplotype)
  colnames(site_prob) <- colnames(final_res_trace)[1:ncol(site_prob)]
  write.table(site_prob,file=gsub(".txt$","_siteFreq.txt",output_file),
              row.names = T,col.names = T,quote = F,sep='\t')
  return(TRUE)
}
