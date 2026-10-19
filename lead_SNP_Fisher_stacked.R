invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
out_dir <- "outputs/Top_SNP_regions_GCA029852925_20260921"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pheno_file <- "C:/Users/15280/OneDrive/Desktop/Data availability/Genetic mapping identifies regional candidates and genotype–phenotype associations/1.Phenotype.xlsx"
geno_file <- "C:/Users/15280/OneDrive/Desktop/Data availability/Genetic mapping identifies regional candidates and genotype–phenotype associations/2.SNP_GENOTYPE1224_012.xlsx"
pheno <- as.data.frame(readxl::read_excel(pheno_file, sheet="Sheet1"))
head <- readxl::read_excel(geno_file,sheet="Sheet2",range="A1:DA1",col_names=FALSE,col_types="text",.name_repair="minimal")
ids <- as.character(unlist(head[1,6:105], use.names=FALSE))
pheno$ID <- trimws(as.character(pheno$ID))
stopifnot(length(ids)==100L, !anyDuplicated(ids), !anyDuplicated(pheno$ID), setequal(ids,pheno$ID))
pheno <- pheno[match(ids,pheno$ID),]
stopifnot(identical(ids,pheno$ID),all(pheno$`A-Type` %in% c("A-","aa")),all(pheno$`E-Type` %in% c("E-","ee")))
ald <- as.integer(pheno$`A-Type`=="A-"); ester <- as.integer(pheno$`E-Type`=="E-")
phen4 <- ifelse(ald==1 & ester==1,"AE",ifelse(ald==1,"A",ifelse(ester==1,"E","Null")))
meta <- data.frame(SNP=c("SNP12656","SNP39462","SNP24243"),Index=c(12656L,39462L,24243L),
                  Chr=c("Chr17","Chr15","Chr9"),POS=c(4190608,11539227,18714086),
                  SourceContig=c("CTG_23","CTG_49","CTG_31"),SourcePOS=c(17150070,11411344,17490591),
                  REF=c("G","A","T"),ALT=c("T","C","G"),Swap=c(FALSE,TRUE,FALSE))
aligned <- data.frame(ID=ids,Ald=ald,Ester=ester,Phenotype4=phen4)
tabs <- list(); summary <- list(); long <- list(); old <- list()
for(i in 1:3){
  m <- meta[i,]; row <- m$Index+1L
  raw <- as.character(unlist(readxl::read_excel(geno_file,sheet="Sheet2",range=sprintf("A%d:DA%d",row,row),
                        col_names=FALSE,col_types="text",.name_repair="minimal"),use.names=FALSE))
  stopifnot(as.integer(raw[1])==m$Index,raw[2]==m$SourceContig,as.numeric(raw[3])==m$SourcePOS)
  # The current genotype workbook already contains numeric ALT dosages
  # (0 = reference homozygote, 1 = heterozygote, 2 = alternate homozygote).
  g <- suppressWarnings(as.integer(raw[6:105]))
  stopifnot(length(g)==100L, !anyNA(g), all(g %in% 0:2))
  oldg <- g
  if(m$Swap)g <- 2L-g
  aligned[[paste0(m$SNP,"_source")]] <- oldg;aligned[[paste0(m$SNP,"_target")]] <- g
  y <- if(i==1)factor(ald,levels=0:1) else if(i==2)factor(ester,levels=0:1) else factor(phen4,levels=c("A","AE","E","Null"))
  tab <- table(Genotype=factor(g,levels=0:2),Phenotype=y)
  tabold <- table(Genotype=factor(oldg,levels=0:2),Phenotype=y)
  f <- fisher.test(tab,hybrid=FALSE,simulate.p.value=FALSE,workspace=20000000L)
  f0 <- fisher.test(tabold,hybrid=FALSE,simulate.p.value=FALSE,workspace=20000000L)
  stopifnot(abs(f$p.value-f0$p.value)<1e-12,sum(tab)==100L)
  tabs[[i]] <- tabold
  letters_gt <- c(paste0(m$REF,m$REF),paste0(m$REF,m$ALT),paste0(m$ALT,m$ALT))
  original_letters <- if(m$Swap)rev(letters_gt) else letters_gt
  primary_table <- data.frame(SNP=m$SNP,Original_dosage=0:2,Chromosome_genotype=original_letters,N=rowSums(tabold),as.data.frame.matrix(tabold),check.names=FALSE)
  write.csv(primary_table,file.path(out_dir,paste0(m$SNP,"_contingency.csv")),row.names=FALSE)
  tabout <- data.frame(SNP=m$SNP,Target_dosage=0:2,Genotype=letters_gt,N=rowSums(tab),as.data.frame.matrix(tab),check.names=FALSE)
  write.csv(tabout,file.path(out_dir,paste0(m$SNP,"_contingency_target.csv")),row.names=FALSE)
  write.csv(data.frame(Source_dosage=0:2,as.data.frame.matrix(tabold),check.names=FALSE),file.path(out_dir,paste0(m$SNP,"_contingency_source.csv")),row.names=FALSE)
  counts <- as.data.frame(tabold);counts$N_genotype <- rowSums(tabold)[as.integer(counts$Genotype)]
  counts$Percent <- counts$Freq/counts$N_genotype*100;counts$SNP <- m$SNP
  counts$Actual_genotype <- original_letters[as.integer(counts$Genotype)]
  counts$Encoding <- "Original source 0/1/2; letters on chromosome-reference strand"
  long[[i]] <- counts
  summary[[i]] <- data.frame(Analysis=c("Ald","Ester","Four_class")[i],SNP=m$SNP,Chromosome=m$Chr,Position_bp=m$POS,
                             Target_REF=m$REF,Target_ALT=m$ALT,Genotype_0_2_swapped=m$Swap,N=100,Table=paste(dim(tab),collapse="x"),
                             Fisher_exact_P=f$p.value,Simulated=FALSE,Minimum_genotype_n=min(rowSums(tab)))
}
summ <- do.call(rbind,summary);summ$Holm_P_three_selected_tests <- p.adjust(summ$Fisher_exact_P,method="holm")
write.csv(summ,file.path(out_dir,"Fisher_summary.csv"),row.names=FALSE)
write.csv(do.call(rbind,long),file.path(out_dir,"Stacked_bar_data.csv"),row.names=FALSE)
write.csv(aligned,file.path(out_dir,"Top_SNP_100_samples.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out_dir,"R_sessionInfo_Fisher.txt"))

draw <- function(){
  layout(matrix(1:3,1,3));par(mar=c(4.3,4.0,5.0,0.6),oma=c(2.1,0,0,0),family="Arial",cex=1,
       col.axis="#203D57",col.lab="#203D57",fg="#203D57",mgp=c(2.3,0.7,0),xaxs="i",yaxs="i")
  for(i in 1:3){
    m <- meta[i,];tab <- tabs[[i]];pct <- t(tab/rowSums(tab)*100)
    cols <- if(i==1)c("#DDE3E8","#367DAE") else if(i==2)c("#DDE3E8","#D68E33") else c("#367DAE","#6B579D","#D68E33","#DDE3E8")
    labels <- if(i<3)c("Absent (0)","Present (1)") else c("A","AE","E","Null")
    gtlabel <- c(paste0(m$REF,m$REF),paste0(m$REF,m$ALT),paste0(m$ALT,m$ALT))
    if(m$Swap)gtlabel <- rev(gtlabel)
    bp <- barplot(pct,col=cols,border="white",space=.55,width=.85,ylim=c(0,111),axes=FALSE,
                  names.arg=paste0(0:2,"\n",gtlabel),cex.names=.92,ylab="Phenotype within genotype (%)")
    axis(2,at=seq(0,100,20),labels=seq(0,100,20),las=1,cex.axis=.85,lwd=.8)
    axis(1,at=range(bp)+c(-.7,.7),labels=FALSE,tick=FALSE,lwd=.8)
    for(j in 1:3){
      text(bp[j],104,paste0("n = ",rowSums(tab)[j]),cex=.78)
      bottoms <- c(0,head(cumsum(pct[,j]),-1));centers <- bottoms+pct[,j]/2
      for(k in seq_len(nrow(pct)))if(pct[k,j]>=7){
        color <- if((i<3 && k==1)||(i==3 && k==4))"#203D57" else "white"
        text(bp[j],centers[k],sprintf("%.1f%%",floor(pct[k,j]*10+.5)/10),cex=.79,col=color,font=2)
      }
    }
    mtext(paste0(c("Ald","Ester","4-class phenotype")[i]," | ",m$SNP),side=3,line=3.8,cex=1.00,font=2)
    mtext(sprintf("%s:%s",m$Chr,format(m$POS,big.mark=",",scientific=FALSE)),side=3,line=2.6,cex=.82)
    mtext(sprintf("Fisher exact P = %.3g",summ$Fisher_exact_P[i]),side=3,line=1.55,cex=.90)
    legend("top",legend=labels,fill=cols,border=NA,horiz=TRUE,bty="n",cex=.72,xpd=NA,inset=c(0,-.10))
    mtext("Genotype (original 0/1/2 coding)",side=1,line=3.0,cex=.86)
  }
  mtext("N = 100; original genotype codes retained. Allele labels use the chromosome-reference strand. Fisher tests are supplementary.",
        side=1,outer=TRUE,line=.45,cex=.71,col="#586A7B")
}
for(dev in c("preview","png","pdf","svg")){
  if(dev %in% c("preview","png"))ragg::agg_png(file.path(out_dir,if(dev=="preview")"Fisher_stacked_preview.png" else "Fisher_genotype_phenotype.png"),width=12.6,height=4.85,units="in",res=if(dev=="preview")130 else 450)
  if(dev=="pdf")cairo_pdf(file.path(out_dir,"Fisher_genotype_phenotype.pdf"),width=12.6,height=4.85,family="Arial")
  if(dev=="svg")svglite::svglite(file.path(out_dir,"Fisher_genotype_phenotype.svg"),width=12.6,height=4.85)
  draw();dev.off()
}
print(summ)

