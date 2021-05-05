# Clean TP data
readDBdata<-function(phenotypeFile,metadataFile=NULL){
  indata<-read.csv(phenotypeFile,
                   na.strings = c("#VALUE!",NA,".",""," ","-","\""),
                   stringsAsFactors = F)
  if(!is.null(metadataFile)){
    meta<-read.csv(metadataFile,
                   na.strings = c("#VALUE!",NA,".",""," ","-","\""),
                   stringsAsFactors = F) %>%
      rename(programName=breedingProgramName,
             programDescription=breedingProgramDescription,
             programDbId=breedingProgramDbId)
    indata<-left_join(indata,meta) }
  indata %<>%
    filter(observationLevel=="plot")
  return(indata) }

makeTrialTypeVar<-function(indata){
  # So far, this function is not very general
  # Handles IITA and NRCRI trial names as of September 2020.
  # Can customize this or add lines to grab TrialTypes for each breeding program
  if(indata$programName=="IITA"){
    outdata<-indata %>%
      mutate(TrialType=ifelse(grepl("CE|clonal|13NEXTgenC1",studyName,ignore.case = T),"CET",NA),
             TrialType=ifelse(grepl("EC",studyName,ignore.case = T),"ExpCET",TrialType),
             TrialType=ifelse(grepl("PYT",studyName,ignore.case = T),"PYT",TrialType),
             TrialType=ifelse(grepl("AYT",studyName,ignore.case = T),"AYT",TrialType),
             TrialType=ifelse(grepl("UYT",studyName,ignore.case = T),"UYT",TrialType),
             TrialType=ifelse(grepl("geneticgain|gg|genetic gain",studyName,ignore.case = T),"GeneticGain",TrialType),
             TrialType=ifelse(grepl("Cassava",studyName,ignore.case = T) & grepl("/",studyName),"GeneticGain",TrialType),
             TrialType=ifelse((grepl("clonal evaluation trial",!grepl("genetic gain",studyDescription,ignore.case = T),
                                     ignore.case = T)),"CET",TrialType),
             TrialType=ifelse(grepl("preliminary yield trial",studyDescription,ignore.case = T),"PYT",TrialType),
             TrialType=ifelse(grepl("Crossingblock|\\.CB\\.|cross",studyName) & is.na(TrialType),"CrossingBlock",TrialType),
             TrialType=ifelse(grepl("NCRP",studyName) & is.na(TrialType),"NCRP",TrialType),
             TrialType=ifelse(grepl("conservation",studyName) & is.na(TrialType),"Conservation",TrialType),
             TrialType=ifelse(grepl("seedling|\\.SN",studyName),"SN",TrialType)) }
  if(indata$programName=="NRCRI"){
    outdata<-indata %>%
      mutate(TrialType=ifelse(grepl("TP1",studyName,ignore.case = T),"TP1",NA),
             TrialType=ifelse(grepl("TP2",studyName,ignore.case = T),"TP2",TrialType),
             TrialType=ifelse(grepl("C1a",studyName,ignore.case = T),"C1a",TrialType),
             TrialType=ifelse(grepl("C1b",studyName,ignore.case = T),"C1b",TrialType),
             TrialType=ifelse(grepl("C2a",studyName,ignore.case = T),"C2a",TrialType),
             TrialType=ifelse(grepl("C2b",studyName,ignore.case = T),"C2b",TrialType),
             TrialType=ifelse(grepl("NCRP",studyName) & is.na(TrialType),"NCRP",TrialType),
             TrialType=ifelse(grepl("15nextgen60gs-cbUM|crossnblk|crossingblock",studyName,ignore.case = T) &
                                !grepl("CET",studyName),
                              "CrossingBlock",TrialType),
             TrialType=ifelse(grepl("seedling",studyName,ignore.case = T),NA,TrialType)) }

  if(indata$programName=="TARI"){
    outdata<-indata %>%
      mutate(TrialType=ifelse(grepl("Advanced Yield|AYT", trialType, ignore.case = T),"AYT",NA),
             TrialType=ifelse(grepl("Clonal|CET", trialType, ignore.case = T),"CET",TrialType),
             TrialType=ifelse(grepl("Preliminary|PYT", trialType, ignore.case = T),"PYT",TrialType),
             TrialType=ifelse(grepl("Regional", trialType, ignore.case = T),"RegionalTrial",TrialType),
             TrialType=ifelse(grepl("Uniform|UYT", trialType, ignore.case = T),"UYT",TrialType),
             TrialType=ifelse(grepl("Variety Release", trialType, ignore.case = T),"VarietyTrial",TrialType),
             TrialType=ifelse(grepl("CROSSING", trialType, ignore.case = T),"CrossingBlock",TrialType),
             TrialType=ifelse(grepl("GWAS", trialType, ignore.case = T),"GWAS",TrialType)) }

  return(outdata) }


#' @param  traitabbrevs data.frame with 2 cols (TraitAbbrev and TraitName). TraitName should match exactly to cassava ontology names
#' @param  indata data.frame read from cassavabase download
#' @param  customColsToKeep char. vec. of any custom cols you added and want to keep
renameAndSelectCols<-function(traitabbrevs,indata,
                              customColsToKeep=NULL){
  outdata<-indata %>%
    select(any_of(c("studyYear","programName","locationName","studyName","studyDesign",
                    "plotWidth","plotLength","fieldSize","plantingDate","harvestDate",
                    "germplasmName","observationUnitDbId",
                    "replicate","blockNumber","plotNumber","rowNumber","colNumber","entryType",
                    "trialType","plantsPerPlot","numberBlocks","numberReps")),
           any_of(customColsToKeep),
           any_of(traitabbrevs$TraitName)) %>% ungroup() %>%
    mutate(across(any_of(traitabbrevs$TraitName), as.numeric)) %>% ungroup() %>%
    pivot_longer(cols = any_of(traitabbrevs$TraitName),
                 names_to = "TraitName",
                 values_to = "Value") %>%
    left_join(.,traitabbrevs) %>%
    select(-TraitName) %>%
    pivot_wider(names_from = TraitAbbrev,
                values_from = "Value")
  return(outdata) }

# Curate by trial
nestByTrials<-function(indata){
  nested_indata<-indata %>%
    # Create some explicitly nested variables including loc and year to nest with the trial data
    mutate(yearInLoc=paste0(programName,"_",locationName,"_",studyYear),
           trialInLocYr=paste0(yearInLoc,"_",studyName),
           repInTrial=paste0(trialInLocYr,"_",replicate),
           blockInRep=paste0(repInTrial,"_",blockNumber)) %>%
    nest(TrialData=-c(programName,locationName,studyYear,TrialType,studyName))
  return(nested_indata)
}

detectExptDesigns<-function(indata){
  # nestByTrials
  nestedDBdata<-indata %>%
    # Create some explicitly nested variables including loc and year to nest with the trial data
    mutate(yearInLoc=paste0(programName,"_",locationName,"_",studyYear),
           trialInLocYr=paste0(yearInLoc,"_",studyName),
           repInTrial=paste0(trialInLocYr,"_",replicate),
           blockInRep=paste0(repInTrial,"_",blockNumber)) %>%
    nest(TrialData=-c(programName,locationName,studyYear,TrialType,studyName))

  # Define complete blocks
  nestedDBdata %>%
    mutate(Nobs=map_dbl(TrialData,~nrow(.)),
           MaxNOHAV=map_dbl(TrialData,~unique(.$MaxNOHAV)),
           Nrep=map_dbl(TrialData,~length(unique(.$replicate))),
           Nblock=map_dbl(TrialData,~length(unique(.$blockInRep))),
           Nclone=map_dbl(TrialData,~length(unique(.$germplasmName))),
           # median number of obs per clone
           medObsPerClone=map_dbl(TrialData,~count(.,germplasmName) %$% round(median(n),1)),
           # median number of obs per replicate
           medObsPerRep=map_dbl(TrialData,~count(.,replicate) %$% round(median(n),1)),
           # Define complete block effects based on the "replicate" variable
           CompleteBlocks=ifelse(Nrep>1 & medObsPerClone==Nrep & Nobs!=Nrep,TRUE,FALSE),
           # Additional trials with imperfect complete blocks
           CompleteBlocks=ifelse(Nrep>1 & medObsPerClone!=Nrep & medObsPerClone>1 & Nobs!=Nrep,TRUE,CompleteBlocks)) -> x
  x %>%
    # Some complete blocks may only be represented by the "blockNumber" column
    mutate(medBlocksPerClone=map_dbl(TrialData,~select(.,blockInRep,germplasmName) %>%
                                       # median number of blockInRep per clone
                                       distinct %>%
                                       count(germplasmName) %$%
                                       round(median(n))),
           # If CompleteBlocks==FALSE (complete blocks not detected based on replicate)
           # and if more than half the clones are represented in more than one block based on the blockInRep variable
           # Copy the blockInRep values into the repInTrial column
           # Recompute Nrep
           # and declare CompleteBlocks==TRUE
           TrialData=ifelse(medBlocksPerClone>1 & CompleteBlocks==FALSE,map(TrialData,~mutate(.,repInTrial=blockInRep)),TrialData),
           Nrep=map_dbl(TrialData,~length(unique(.$repInTrial))),
           CompleteBlocks=ifelse(medBlocksPerClone>1 & CompleteBlocks==FALSE,TRUE,CompleteBlocks)) -> y
  # Define incomplete blocks
  y %>%
    mutate(repsEqualBlocks=map_lgl(TrialData,~all(.$replicate==.$blockNumber)),
           NrepEqualNblock=ifelse(Nrep==Nblock,TRUE,FALSE),
           medObsPerBlockInRep=map_dbl(TrialData,~count(.,blockInRep) %$% round(median(n),1))) -> z
  # Define complete blocked trials with nested sub-blocks
  z %<>%
    mutate(IncompleteBlocks=ifelse(CompleteBlocks==TRUE & Nobs!=Nblock & Nblock>1 & medObsPerBlockInRep>1 & NrepEqualNblock==FALSE,TRUE,FALSE))
  # Define clearly unreplicated (CompleteBlocks==FALSE & Nrep==1) trials with nested sub-blocks
  z %<>%
    mutate(IncompleteBlocks=ifelse(CompleteBlocks==FALSE & Nobs!=Nblock & Nblock>1 & medObsPerBlockInRep>1 & Nrep==1,TRUE,IncompleteBlocks))
  # Define additional trials with incomplete blocks (blockInRep) where CompleteBlocks==FALSE but Nrep>1 and Nrep==Block
  z %<>%
    mutate(IncompleteBlocks=ifelse(CompleteBlocks==FALSE & IncompleteBlocks==FALSE &
                                     Nobs!=Nblock & Nblock>1 &  Nobs!=Nrep &
                                     medObsPerBlockInRep>1 & Nrep>1 & NrepEqualNblock==TRUE,TRUE,IncompleteBlocks))
  # Last few cases (2 trials actually) where Nrep>1 and Nblock>1 and Nrep!=Nblock but CompleteBlocks==FALSE
  z %<>%
    mutate(IncompleteBlocks=ifelse(CompleteBlocks==FALSE & IncompleteBlocks==FALSE &
                                     Nobs!=Nblock & Nobs!=Nrep &
                                     medObsPerBlockInRep>1 & Nrep>1,TRUE,IncompleteBlocks))
  z %<>%
    dplyr::select(-MaxNOHAV) %>%
    unnest(TrialData)
  return(z)
}

nestDesignsDetectedByTraits<-function(indata,traits){
  indata %<>%
    select(programName,locationName,studyYear,TrialType,studyName,
           CompleteBlocks,IncompleteBlocks,
           yearInLoc,trialInLocYr,repInTrial,blockInRep,observationUnitDbId,
           germplasmName,FullSampleName,GID,all_of(traits),PropNOHAV) %>%
    mutate(IncompleteBlocks=ifelse(IncompleteBlocks==TRUE,"Yes","No"),
           CompleteBlocks=ifelse(CompleteBlocks==TRUE,"Yes","No")) %>%
    pivot_longer(cols = all_of(traits), names_to = "Trait", values_to = "Value") %>%
    filter(!is.na(Value),
           !is.na(GID)) %>%
    nest(MultiTrialTraitData=c(-Trait))
  return(indata)
}

nestTrialsByTrait<-function(indata,traits){
  nested_trialdata<-dbdata %>%
    select(-MaxNOHAV) %>%
    unnest(TrialData) %>%
    pivot_longer(cols = any_of(traits),
                 names_to = "Trait",
                 values_to = "TraitValue") %>%
    nest(TraitByTrialData=-c(Trait,studyYear,programName,locationName,studyName,TrialType))
  return(nested_trialdata)
}

calcPropMissing<-function(TraitValues){ length(which(is.na(TraitValues))) / length(TraitValues) }

curateTrialOneTrait<-function(Trait,TraitByTrialData,GID="GID"){
  require(lme4)

  modelFormula<-paste0("TraitValue ~ (1|",GID,")")
  modelFormula<-ifelse(all(TraitByTrialData$CompleteBlocks),
                       paste0(modelFormula,"+(1|repInTrial)"),modelFormula)
  modelFormula<-ifelse(all(TraitByTrialData$IncompleteBlocks),
                       paste0(modelFormula,"+(1|blockInRep)"),modelFormula)
  modelFormula<-ifelse(grepl("logRTNO",Trait) | grepl("logFYLD",Trait) | grepl("logTOPYLD",Trait),
                       paste0(modelFormula,"+PropNOHAV"),modelFormula)

  propMiss<-calcPropMissing(TraitByTrialData$TraitValue)
  fit_model<-possibly(function(modelFormula,TraitByTrialData){
    model_out<-lmer(as.formula(modelFormula),data=TraitByTrialData)
    if(!is.na(model_out)){
      outliers<-which(abs(rstudent(model_out))>=3.3)
      if(length(outliers)>0){
        model_out<-lmer(as.formula(modelFormula),data=TraitByTrialData,
                        subset=abs(rstudent(model_out))<3.3)
      }
    }
    return(list(model_out=model_out,outliers=outliers)) },
    otherwise = NA)
  model_out<-fit_model(modelFormula,TraitByTrialData)
  if(is.na(model_out)){
    out <-tibble(H2=NA,VarComps=list(NULL),BLUPs=list(NULL),Model=modelFormula,Noutliers=NA,Outliers=NA,propMiss=propMiss)
  } else {
    varcomps<-as.data.frame(VarCorr(model_out[["model_out"]]))[,c("grp","vcov")] %>%
      spread(grp,vcov)
    Vg<-varcomps$GID
    H2<-Vg/(Vg+varcomps$Residual)
    BLUP<-ranef(model_out[["model_out"]], condVar=TRUE)[[GID]]
    PEV <- c(attr(BLUP, "postVar"))
    blups<-tibble(GID=rownames(BLUP),BLUP=BLUP$`(Intercept)`,PEV=PEV) %>%
      mutate(REL=1-(PEV/Vg),
             drgBLUP=BLUP/REL,
             WT=(1-H2)/((0.1 + (1-REL)/REL)*H2))
    out <- tibble(H2=H2,
                  VarComps=list(varcomps),
                  BLUPs=list(blups),
                  Model=modelFormula,
                  Noutliers=length(model_out[["outliers"]]),
                  Outliers=list(model_out[["outliers"]]),
                  propMiss=propMiss) }
  return(out)
}

curateTrialsByTrait<-function(nestedTrialData,traits){
  outdata<-nestedTrialData %>%
    mutate(modelOutput=map2(Trait,TraitByTrialData,~curateTrialOneTrait(Trait = .x,TraitByTrialData = .y))) %>%
    dplyr::select(-TraitByTrialData) %>%
    unnest(modelOutput)
  return(outdata)
}

# Get BLUPs
nestForMultiTrialAnalysis<-function(curatedTrialData){
  nested_trialdata<-curatedTrialData %>%
    # remove trait-trial models that failed
    filter(!is.na(H2)) %>%
    # remove some per-trial summaries we don't want at this stage
    select(-H2,-VarComps,-Model,-Noutliers,-propMiss) %>%
    unnest(BLUPs) %>%
    nest(MultiTrialTraitData=c(-Trait))
  return(nested_trialdata)
}

fitMultiTrialModel<-function(curatedTrialData,GID="GID"){
  require(lme4)
  modelFormula<-paste0("drgBLUP ~ (1|",GID,")")
  fit_model<-possibly(function(modelFormula,curatedTrialData){
    model_out<-lmer(as.formula(modelFormula),
                    data=curatedTrialData,
                    weights = WT)
    return(model_out) },
    otherwise = NA)
  model_out<-fit_model(modelFormula,curatedTrialData)
  summary(model_out)
  if(is.na(model_out)){
    out <-tibble(H2=NA,VarComps=list(NULL),BLUPs=list(NULL),Model=modelFormula)
  } else {
    varcomps<-as.data.frame(VarCorr(model_out))[,c("grp","vcov")] %>%
      spread(grp,vcov)
    Vg<-varcomps$GID
    H2<-Vg/(Vg+varcomps$Residual)
    BLUP<-ranef(model_out, condVar=TRUE)[[GID]]
    PEV <- c(attr(BLUP, "postVar"))
    blups<-tibble(GID=rownames(BLUP),BLUP=BLUP$`(Intercept)`,PEV=PEV) %>%
      mutate(REL=1-(PEV/Vg),
             drgBLUP=BLUP/REL,
             WT=(1-H2)/((0.1 + (1-REL)/REL)*H2))
    out <- tibble(H2=H2,
                  VarComps=list(varcomps),
                  BLUPs=list(blups),
                  Model=modelFormula) }
  return(out)
}

# Cross-validation

maf_filter<-function(snps,thresh){
  freq<-colMeans(snps, na.rm=T)/2; maf<-freq;
  maf[which(maf > 0.5)]<-1-maf[which(maf > 0.5)]
  snps1<-snps[,which(maf>thresh)];
  return(snps1) }

#' kinship function
#'
#' Function to create additive and dominance genomic relationship matrices from biallelic dosages.
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID.
#' @param type string, "add" or "dom". type="add" gives same as rrBLUP::A.mat(), i.e. Van Raden, Method 1. type="dom" gives classical parameterization according to Vitezica et al. 2013.
#'
#' @return square symmetic genomic relationship matrix
#' @export
#'
#' @examples
#' K<-kinship(M,"add")
kinship<-function(M,type){
  M<-round(M)
  freq <- colMeans(M,na.rm=T)/2
  P <- matrix(rep(freq,nrow(M)),byrow=T,ncol=ncol(M))
  if(type=="add"){
    Z <- M-2*P
    varD<-sum(2*freq*(1-freq))
    K <- tcrossprod(Z)/ varD
    return(K)
  }
  if(type=="dom"){
    W<-M;
    W[which(W==1)]<-2*P[which(W==1)];
    W[which(W==2)]<-(4*P[which(W==2)]-2);
    W <- W-2*(P^2)
    varD<-sum((2*freq*(1-freq))^2)
    D <- tcrossprod(W) / varD
    return(D)
  }
}


#' @param byGroup logical, if TRUE, assumes a column named "Group" is present which unique classifies each GID into some genetic grouping.
#' @param modelType string, A, AD or ADE representing model with Additive-only, Add. plus Dominance, and Add. plus Dom. plus. AxD Epistasis (AD), respectively.
#' @param grms list of GRMs where each element is named either A, D, or, AD. Matrices supplied must match required by A, AD and ADE models. For ADE grms=list(A=A,D=D,AD=AD)...
#' @param augmentTP option to supply an additional set of training data, which will be added to each training model but never included in the test set.
#' @param TrainTestData data.frame with de-regressed BLUPs, BLUPs and weights (WT) for training and test. If byGroup==TRUE, a column with Group as the header uniquely classifying GIDs into genetic groups, is expected.
runCrossVal<-function(TrainTestData,modelType,grms,nrepeats,nfolds,ncores=1,
                          byGroup=FALSE,augmentTP=NULL,gid="GID",...){
  require(sommer); require(rsample)
  # Set-up replicated cross-validation folds
  # splitting by clone (if clone in training dataset, it can't be in testing)
  if(byGroup){
    cvsamples<-tibble(GroupName=unique(TrainTestData$Group))
  } else { cvsamples<-tibble(GroupName="None") }
  cvsamples<-cvsamples %>%
    mutate(Splits=map(GroupName,function(GroupName){
      if(GroupName!="None"){
        thisgroup<-TrainTestData %>%
          filter(Group==GroupName) } else { thisgroup<-TrainTestData }
      out<-tibble(repeats=1:nrepeats,
                  splits=rerun(nrepeats,group_vfold_cv(thisgroup, group = gid, v = nfolds))) %>%
        unnest(splits)
      return(out)
    })) %>%
    unnest(Splits)

  ## Internal function
  ## fits prediction model and calcs. accuracy for each train-test split

  fitModel<-possibly(function(splits,modelType,augmentTP,TrainTestData,GroupName,grms){
    starttime<-proc.time()[3]
    # Set-up training set
    trainingdata<-training(splits)
    ## Make sure, if there is an augmentTP, no GIDs in test-sets
    if(!is.null(augmentTP)){
      ## remove any test-set members from augment TP before adding to training data
      training_augment<-augmentTP %>% filter(!(!!sym(gid) %in% testing(splits)[[gid]]))
      trainingdata<-bind_rows(trainingdata,training_augment) }
    if(GroupName!="None"){ trainingdata<-bind_rows(trainingdata,
                                                   TrainTestData %>%
                                                     filter(Group!=GroupName,
                                                            !(!!sym(gid) %in% testing(splits)[[gid]]))) }
    # Subset kinship matrices
    traintestgids<-union(trainingdata[[gid]],testing(splits)[[gid]])
    A1<-grms[["A"]][traintestgids,traintestgids]
    trainingdata[[paste0(gid,"a")]]<-factor(trainingdata[[gid]],levels=rownames(A1))
    if(modelType %in% c("AD","ADE")){
      D1<-grms[["D"]][traintestgids,traintestgids]
      trainingdata[[paste0(gid,"d")]]<-factor(trainingdata[[gid]],levels=rownames(D1))
      if(modelType=="ADE"){
        #AA1<-grms[["AA"]][traintestgids,traintestgids]
        AD1<-grms[["AD"]][traintestgids,traintestgids]
        diag(AD1)<-diag(AD1)+1e-06
        #DD1<-grms[["DD"]][traintestgids,traintestgids]
        #trainingdata[[paste0(gid,"aa")]]<-factor(trainingdata[[gid]],levels=rownames(AA1))
        trainingdata[[paste0(gid,"ad")]]<-factor(trainingdata[[gid]],levels=rownames(AD1))
        #trainingdata[[paste0(gid,"dd")]]<-factor(trainingdata[[gid]],levels=rownames(DD1))
      }
    }
    # Set-up random model statements
    randFormula<-paste0("~vs(",gid,"a,Gu=A1)")
    if(modelType %in% c("AD","ADE")){
      randFormula<-paste0(randFormula,"+vs(",gid,"d,Gu=D1)")
      if(modelType=="ADE"){
        randFormula<-paste0(randFormula,"+vs(",gid,"ad,Gu=AD1)")
        #"+vs(",gid,"aa,Gu=AA1)",
        #"+vs(",gid,"ad,Gu=AD1)")
        #"+vs(",gid,"dd,Gu=DD1)")
      }
    }
    # Fit genomic prediction model
    fit <- mmer(fixed = drgBLUP ~1,
                random = as.formula(randFormula),
                weights = WT,
                data=trainingdata)
    # Gather the BLUPs
    gblups<-tibble(GID=as.character(names(fit$U[[paste0("u:",gid,"a")]]$drgBLUP)),
                   GEBV=as.numeric(fit$U[[paste0("u:",gid,"a")]]$drgBLUP))
    if(modelType %in% c("AD","ADE")){
      gblups %<>% mutate(GEDD=as.numeric(fit$U[[paste0("u:",gid,"d")]]$drgBLUP))
      if(modelType=="ADE"){
        gblups %<>% mutate(#GEEDaa=as.numeric(fit$U[[paste0("u:",gid,"aa")]]$drgBLUP),
          GEEDad=as.numeric(fit$U[[paste0("u:",gid,"ad")]]$drgBLUP))
        #GEEDdd=as.numeric(fit$U[[paste0("u:",gid,"dd")]]$drgBLUP))
      }
    }
    # Calc GETGVs
    ## Note that for modelType=="A", GEBV==GETGV
    gblups %<>%
      mutate(GETGV=rowSums(.[,grepl("GE",colnames(.))]))
    # Test set validation data
    validationData<-TrainTestData %>%
      dplyr::select(gid,BLUP) %>%
      filter(GID %in% testing(splits)[[gid]])
    # Measure accuracy in test set
    ## cor(GEBV,BLUP)
    ## cor(GETGV,BLUP)
    accuracy<-gblups %>%
      mutate(GETGV=rowSums(.[,grepl("GE",colnames(.))])) %>%
      filter(GID %in% testing(splits)[[gid]]) %>%
      left_join(validationData) %>%
      summarize(accGEBV=cor(GEBV,BLUP, use = 'complete.obs'),
                accGETGV=cor(GETGV,BLUP, use = 'complete.obs'))
    computeTime<-proc.time()[3]-starttime
    accuracy %<>% mutate(computeTime=computeTime)
    return(accuracy)
  },otherwise = NA)
  ## Run models across all train-test splits
  ## Parallelize
  require(furrr); plan(multiprocess); options(mc.cores=ncores);
  cvsamples<-cvsamples %>%
    mutate(accuracy=future_map2(splits,GroupName,
                                ~fitModel(splits=.x,GroupName=.y,
                                          modelType=modelType,augmentTP=NULL,TrainTestData=TrainTestData,grms=grms),
                                .progress = FALSE)) %>%
    unnest(accuracy)
  return(cvsamples)
}

#' @param blups nested data.frame with list-column "TrainingData" containing BLUPs
#' @param modelType string, A, AD or ADE representing model with Additive-only, Add. plus Dominance, and Add. plus Dom. plus. Epistasis (AA+AD+DD), respectively.
#' @param grms list of GRMs. Any genotypes in the GRMs get predicted with, or without phenotypes. Each element is named either A, D, AA, AD, DD. Matrices supplied must match required by A, AD and ADE models. For ADE grms=list(A=A,D=D,AA=AA,AD=AD,DD=DD).
runGenomicPredictions<-function(blups,modelType,grms,ncores=1,gid="GID",...){
  require(sommer);
  runOnePred<-possibly(function(trainingdata,modelType,grms){
    trainingdata[[paste0(gid,"a")]]<-factor(trainingdata[[gid]],levels=rownames(grms[["A"]]))
    if(modelType %in% c("AD","ADE")){ trainingdata[[paste0(gid,"d")]]<-factor(trainingdata[[gid]],levels=rownames(grms[["D"]]))
    if(modelType=="ADE"){
      #trainingdata[[paste0(gid,"aa")]]<-factor(trainingdata[[gid]],levels=rownames(grms[["AA"]]))
      trainingdata[[paste0(gid,"ad")]]<-factor(trainingdata[[gid]],levels=rownames(grms[["AD"]]))
      diag(grms[["AD"]])<-diag(grms[["AD"]])+1e-06
      #trainingdata[[paste0(gid,"dd")]]<-factor(trainingdata[[gid]],levels=rownames(grms[["DD"]]))
    }
    }
    # Set-up random model statements
    randFormula<-paste0("~vs(",gid,"a,Gu=A)")
    if(modelType %in% c("AD","ADE")){
      randFormula<-paste0(randFormula,"+vs(",gid,"d,Gu=D)")
      if(modelType=="ADE"){
        randFormula<-paste0(randFormula,"+vs(",gid,"ad,Gu=AD)")
        # "+vs(",gid,"aa,Gu=AA)",
        # "+vs(",gid,"dd,Gu=DD)")
      }
    }
    # Fit genomic prediction model
    fit <- mmer(fixed = drgBLUP ~1,
                random = as.formula(randFormula),
                weights = WT,
                data=trainingdata,getPEV = TRUE)

    # Gather the BLUPs
    gblups<-tibble(GID=as.character(names(fit$U[[paste0("u:",gid,"a")]]$drgBLUP)),
                   GEBV=as.numeric(fit$U[[paste0("u:",gid,"a")]]$drgBLUP))
    pev<-diag((fit$PevU[["u:GIDa"]]$drgBLUP))
    pev<-tibble(GID=names(pev),
                PEVa=pev)
    gblups %<>% left_join(pev)
    if(modelType %in% c("AD","ADE")){
      gblups %<>% mutate(GEDD=as.numeric(fit$U[[paste0("u:",gid,"d")]]$drgBLUP))
      pev<-diag((fit$PevU[["u:GIDd"]]$drgBLUP))
      pev<-tibble(GID=names(pev),
                  PEVd=pev)
      gblups %<>% left_join(pev)

      if(modelType=="ADE"){
        gblups %<>% mutate(GEEDad=as.numeric(fit$U[[paste0("u:",gid,"ad")]]$drgBLUP))
        pev<-diag((fit$PevU[["u:GIDad"]]$drgBLUP))
        pev<-tibble(GID=names(pev),
                    PEVad=pev)
        gblups %<>% left_join(pev)
      }
    }
    # Calc GETGVs
    ## Note that for modelType=="A", GEBV==GETGV
    gblups %<>%
      mutate(GETGV=rowSums(.[,grepl("GE",colnames(.))]))
    varcomps<-summary(fit)$varcomp
    out<-tibble(gblups=list(gblups),varcomps=list(varcomps))
    return(out)
  },otherwise = NA)
  ## Run models across all train-test splits
  ## Parallelize
  require(furrr); plan(multiprocess); options(mc.cores=ncores);
  predictions<-blups %>%
    mutate(genomicPredOut=future_map(TrainingData,~runOnePred(trainingdata=.,
                                                              modelType=modelType,grms=grms)))
  return(predictions)
}

# runCrossVal(TrainTestData=traintestdf,
#             modelType="A",
#             grms=list(A=A),
#             byGroup=FALSE,augmentTP=NULL,
#             nrepeats=5,nfolds=4,ncores=1,gid="GID")
#
# runCrossVal(TrainTestData=traintestdf,
#             modelType="A",
#             grms=list(A=A),
#             byGroup=FALSE,augmentTP=augmentDF,
#             nrepeats=5,nfolds=4,ncores=1,gid="GID")
#
# runCrossVal(TrainTestData=traintestdf,
#             modelType="A",
#             grms=list(A=A),
#             byGroup=TRUE,augmentTP=augmentDF,
#             nrepeats=5,nfolds=4,ncores=1,gid="GID")
