require(tidyverse); require(magrittr);

splitVCFbyChr<-function(Chr,vcfIn,filters=NULL,outPath,outSuffix){
  system(paste0("vcftools --gzvcf ",vcfIn," ",
                "--chr ",Chr," ",filters," ",
                "--recode --stdout | bgzip -c -@ 24 > ",
                outPath,"chr",Chr,"_",outSuffix,".vcf.gz")) }

runBeagle5<-function(targetVCF,refVCF,mapFile,outName,
                     nthreads,maxmem="500g",impute=TRUE,ne=100000,samplesToExclude=NULL){
  system(paste0("java -Xms2g -Xmx",maxmem," -jar /programs/beagle/beagle.jar ",
                "gt=",targetVCF," ",
                "map=",mapFile," ",
                "ref=",refVCF," ",
                "out=",outName," ",
                "nthreads=",nthreads," impute=",impute," ne=",ne,
                ifelse(!is.null(samplesToExclude),paste0(" excludesamples=",samplesToExclude),""))) }

postImputeFilter<-function(inPath=NULL,inName,outPath=NULL,outName,DR2thresh=0.75,HWEthresh=1e-20,MAFthresh=0.005){
  require(magrittr); require(dplyr)
  # Extract imputation quality scores (DR2 and AF) from VCF
  system(paste0("vcftools --gzvcf ",inPath,inName,".vcf.gz --get-INFO DR2 --get-INFO AF --out ",outPath,inName))
  system(paste0("vcftools --gzvcf ",inPath,inName,".vcf.gz --hardy --out ",outPath,inName))

  # Read scores into R
  INFO<-read.table(paste0(outPath,inName,".INFO"),
                   stringsAsFactors = F, header = T)
  hwe<-read.table(paste0(outPath,inName,".hwe"),
                  stringsAsFactors = F, header = T)
  stats2filterOn<-left_join(INFO,hwe %>% rename(CHROM=CHR))
  # Compute MAF from AF and make sure numeric
  stats2filterOn %<>%
    dplyr::mutate(DR2=as.numeric(DR2),
                  AF=as.numeric(AF)) %>%
    dplyr::filter(!is.na(DR2),
                  !is.na(AF)) %>%
    dplyr::mutate(MAF=ifelse(AF>0.5,1-AF,AF))
  # Identify sites passing filter
  sitesPassingFilters<-stats2filterOn %>%
    dplyr::filter(DR2>=DR2thresh,
                  P_HWE>HWEthresh,
                  MAF>MAFthresh) %>%
    dplyr::select(CHROM,POS)
  print(paste0(nrow(sitesPassingFilters)," sites passing filter"))

  # Write a list of positions passing filter to disk
  write.table(sitesPassingFilters,
              file = paste0(outPath,inName,".sitesPassing"),
              row.names = F, col.names = F, quote = F)
  # Apply filter to vcf file with vcftools
  system(paste0("vcftools --gzvcf ",inPath,inName,".vcf.gz"," ",
                "--positions ",outPath,inName,".sitesPassing"," ",
                "--recode --stdout | bgzip -c -@ 24 > ",
                outPath,outName,".vcf.gz"))
  print(paste0("Filtering Complete: ",outName))
}

mergeVCFs<-function(inPath=NULL,inVCF1,inVCF2,outPath=NULL,outName){
  system(paste0("tabix -f -p vcf ",inPath,inVCF1,".vcf.gz"))
  system(paste0("tabix -f -p vcf ",inPath,inVCF2,".vcf.gz"))
  system(paste0("bcftools merge ",
                "--output ",outPath,outName,".vcf.gz ",
                "--merge snps --output-type z --threads 6 ",
                inPath,inVCF1,".vcf.gz"," ",
                inPath,inVCF2,".vcf.gz"))
}

filter_positions<-function(inPath=NULL,inVCF,positionFile,outPath=NULL,outName){
  system(paste0("vcftools --gzvcf ",inPath,inVCF," ",
                "--positions ",inPath,positionFile," ",
                "--recode --stdout | bgzip -c -@ 24 > ",
                outPath,outName,".vcf.gz"))
}

convertVCFtoDosage<-function(pathIn,pathOut,vcfName){
  #' @vcfName expects file to be *.vcf.gz, expects name _not_ to include extension .vcf.gz
  # Make binary plink
  system(paste0("export PATH=/programs/plink-1.9-x86_64-beta3.30:$PATH;",
                "plink --vcf ",pathIn,vcfName,".vcf.gz ",
                "--make-bed --const-fid --keep-allele-order ",
                "--out ",pathOut,vcfName))
  # Write first 5 columns of VCF file to txt file *.sitesWithAlleles
  system(paste0("zcat ",pathIn,vcfName,".vcf.gz ",
                "| cut -f1-5 > ",pathOut,vcfName,".sitesWithAlleles"))
  # Read *.sitesWithAlleles, take columns required by plink and write back to disk
  read.table(paste0(pathOut,vcfName,".sitesWithAlleles"),
             stringsAsFactors = F, header = F, sep = c("\t")) %>%
    select(V3,V5) %>%
    write.table(.,file=paste0(pathOut,vcfName,".alleleToCount"),
                row.names = F, sep = c("\t"), quote = F, col.names = F)
  # Recode to dosage
  system(paste0("export PATH=/programs/plink-1.9-x86_64-beta3.30:$PATH;",
                "plink --bfile ",pathOut,vcfName," --keep-allele-order --recode A ",
                "--recode-allele ",pathOut,vcfName,".alleleToCount ",
                "--out ",pathOut,vcfName))
}

createGenomewideDosage<-function(pathIn,chroms,nameOfchromWiseDosageFiles){
  snps<-future_map(chroms,~read.table(paste0(pathIn,"chr",.,nameOfchromWiseDosageFiles,".raw"), stringsAsFactor=F, header = T) %>%
                     dplyr::select(-FID,-PAT,-MAT,-SEX,-PHENOTYPE) %>%
                     column_to_rownames(var = "IID") %>%
                     as.matrix()) %>%
    reduce(.,cbind)
  saveRDS(snps,file = paste0(pathIn,"DosageMatrix",nameOfchromWiseDosageFiles,".rds"))
}

runBeagle4pt1GL<-function(targetVCF,refVCF,mapFile,outName,
                          nthreads,maxmem="500g",impute=TRUE,ne=100000,samplesToExclude=NULL,niter=10){
  system(paste0("java -Xms2g -Xmx",maxmem," -jar /programs/beagle41/beagle41.jar ",
                "gl=",targetVCF," ",
                "map=",mapFile," ",
                "ref=",refVCF," ",
                "out=",outName," ",
                "niterations=",niter," ",
                "nthreads=",nthreads," impute=",impute," ne=",ne,
                ifelse(!is.null(samplesToExclude),paste0(" excludesamples=",samplesToExclude),""))) }

postImputeFilterBeagle4pt1<-function(inPath=NULL,inName,outPath=NULL,outName,AR2thresh=0.75,HWEthresh=1e-20,MAFthresh=0.005){
  require(magrittr); require(dplyr)
  # Extract imputation quality scores (DR2 and AF) from VCF
  system(paste0("vcftools --gzvcf ",inPath,inName,".vcf.gz --get-INFO AR2 --get-INFO AF --out ",outPath,inName))
  system(paste0("vcftools --gzvcf ",inPath,inName,".vcf.gz --hardy --out ",outPath,inName))

  # Read scores into R
  INFO<-read.table(paste0(outPath,inName,".INFO"),
                   stringsAsFactors = F, header = T)
  hwe<-read.table(paste0(outPath,inName,".hwe"),
                  stringsAsFactors = F, header = T)
  stats2filterOn<-left_join(INFO,hwe %>% rename(CHROM=CHR))
  # Compute MAF from AF and make sure numeric
  stats2filterOn %<>%
    dplyr::mutate(AR2=as.numeric(AR2),
                  AF=as.numeric(AF)) %>%
    dplyr::filter(!is.na(AR2),
                  !is.na(AF)) %>%
    dplyr::mutate(MAF=ifelse(AF>0.5,1-AF,AF))
  # Identify sites passing filter
  sitesPassingFilters<-stats2filterOn %>%
    dplyr::filter(AR2>=AR2thresh,
                  P_HWE>HWEthresh,
                  MAF>MAFthresh) %>%
    dplyr::select(CHROM,POS)
  print(paste0(nrow(sitesPassingFilters)," sites passing filter"))

  # Write a list of positions passing filter to disk
  write.table(sitesPassingFilters,
              file = paste0(outPath,inName,".sitesPassing"),
              row.names = F, col.names = F, quote = F)
  # Apply filter to vcf file with vcftools
  system(paste0("vcftools --gzvcf ",inPath,inName,".vcf.gz"," ",
                "--positions ",outPath,inName,".sitesPassing"," ",
                "--recode --stdout | bgzip -c -@ 24 > ",
                outPath,outName,".vcf.gz"))
  print(paste0("Filtering Complete: ",outName))
}

# targetVCF=paste0(targetVCFpath,"chr",1,"_DCas20_5360.vcf.gz")
# refVCF=paste0(refVCFpath,"chr",1,"_ImputationReferencePanel_EMBRAPA_Phased_102619.vcf.gz")
# mapFile=paste0(mapPath,"chr",1,"_cassava_cM_pred.v6_91019.map")
# outName=paste0(outPath,"chr",1,"_DCas20_5360_REFimputed")
# nthreads=112
# maxmem="500g"
# impute=TRUE
# ne=100000
