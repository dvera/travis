bedPoissonTest <- function( test_r1, test_r2, cntl_r1, cntl_r2 , windows , adjustp="fdr" , pvalue=0.05 , threads=getOption("threads",1L), lower=FALSE ){

  test_r1_lines <- filelines(test_r1)
  test_r2_lines <- filelines(test_r2)
  cntl_r1_lines <- filelines(cntl_r1)
  cntl_r2_lines <- filelines(cntl_r2)

  refcount <- min(c(test_r1_lines,test_r2_lines,cntl_r1_lines,cntl_r2_lines))

  test_r1_cov <- ( bedtoolsCoverage (test_r1,windows,scalar=1) )
  test_r2_cov <- ( bedtoolsCoverage (test_r2,windows,scalar=1) )
  cntl_r1_cov <- ( bedtoolsCoverage (cntl_r1,windows,scalar=1) )
  cntl_r2_cov <- ( bedtoolsCoverage (cntl_r2,windows,scalar=1) )


  bgl <- bgRead(c(test_r1_cov,test_r2_cov,cntl_r1_cov,cntl_r2_cov),threads=threads)

  test_rep_diff <- abs( (bgl[,1] * refcount/test_r1_lines) - (bgl[,2] * refcount/test_r2_lines) )
  cntl_rep_diff <- abs( (bgl[,3] * refcount/cntl_r1_lines) - (bgl[,4] * refcount/cntl_r2_lines) )
  maxm_rep_diff <- apply(data.frame(test_rep_diff,cntl_rep_diff),1,max)
  #maxm_rep_diff <- test_rep_diff
  #c_gt_t <- which(cntl_rep_diff>test_rep_diff)
  #maxm_rep_diff[c_gt_t] <- cntl_rep_diff[c_gt_t]

  test_ctl_diff1 <- ( (bgl[,1] * refcount/test_r1_lines) - (bgl[,3] * refcount/cntl_r1_lines) )
  test_ctl_diff2 <- ( (bgl[,1] * refcount/test_r1_lines) - (bgl[,4] * refcount/cntl_r2_lines) )
  test_ctl_diff3 <- ( (bgl[,2] * refcount/test_r2_lines) - (bgl[,3] * refcount/cntl_r1_lines) )
  test_ctl_diff4 <- ( (bgl[,2] * refcount/test_r2_lines) - (bgl[,4] * refcount/cntl_r2_lines) )

  all_test_diffs <- data.frame(test_ctl_diff1,test_ctl_diff2,test_ctl_diff3,test_ctl_diff4)
  signs <- sign(rowMeans(all_test_diffs))
  minm_tst_diff <- apply(abs(all_test_diffs,1,min))


  pvals <- ppois(minm_tst_diff,maxm_rep_diff,lower=lower)

  if( is.null(adjustp) | is.na(adjustp) ){
    padjust=TRUE
    pvals2 <- pvals
  } else{
    pvals2 <- p.adjust(pvals, method=adjustp)
  }

  pos_segs <- which(pvals2<=pvalue & signs==1)
  neg_segs <- which(pvals2<=pvalue & signs==-1)

  tmpbg <- read_tsv(test_r1_cov,col_names=F)

  pvals3=(-signs*log10(pvals2))
  pvals3[is.infinite(pvals3)]<-0
  pvals3[is.na(pvals3)]<-0


}
