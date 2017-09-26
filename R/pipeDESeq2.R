#' @title A pipeline for DESeq2
#'
#' @description
#' \code{pipeDESeq2} Run a pipeline of DESeq2 functions for differential
#' gene expression.
#'
#' @param counts A count matrix, row names should be gene IDs
#' @param info An experimental design matrix
#' @param formula If specified, perform a likelihood ratio test. A vector of
#' character strings that can be coerced into a formula. Must be of the
#' same length as reduced. Specifies base formulas for LRT.
#' @param reduced A character vector of the same length as formula, which
#' can be coerced to a formula that represents a sub model to formula.
#' All reduced formulae must be sub models of the respective
#' formula.
#' @param testNames A character string of names for the LRT, if specified,
#' must be the same length as formula and reduced.
#' @param contrasts A list of character vectors that can be input into
#' DESeq2::results contrasts argument.
#' @param contrastFormula The formula specifying the design matrix for the
#' contrasts. Keep in mind that this formula specifies how the model is fit
#' and can affect inference of the effects of contrasts.
#' @param contrastNames Names for the contrasts. Must be the same length
#' as contrasts
#' @param returnNormData If NULL (default), no transformed data is returned.
#' Other options are "rlog" and "vst" which return rlog and variance stabilized data
#' respecitvel.
#' @param geneIDs The names of genes. If NA, use row names from counts matrix
#' @param verbose Logical, return progress updates?
#' @param quiet Should messages from DESeq2 be suppressed?
#' @param fitType Passed to DESeq. Default is "parametric"
#' @param betaPrior Passed to DESeq. Default is F
#' @param modelMatrixType Passed to DESeq. Default is "standard"
#' @param parallel Passed to DESeq. Default is F
#' @param minReplicatesForReplace Passed to DESeq. Default is 7
#' @param ... additional arguments to pass to DESeq2::results
#'
#' @details This function runs the following pipeline:
##' \itemize{
##'  \item{1. }{DESeq's pipeline using the specified test}
##'  \item{2. }{Extraction of results (from DESeq2::results)}
##'  \item{3. }{Renaming of columns and combining results across tests}
##' }
##'
#' @return a list with 3 elements (if test is not run, elements are NULL):
##' \itemize{
##'  \item{1. }{"LRT_results": The results from the likelihood ratio tests}
##'  \item{2. }{"CONT_results": The results from the contrasts}
##'  \item{3. }{"normData": normalized data}
##' }
##'
#' @examples
#' \dontrun{
#' library(deTools)
#' data(counts)
#' data(info)

#' stats<-pipeDESeq2(counts=counts, info=info,
#'    formula = c(" ~ trt"),
#'    reduced = c(" ~ 1"),
#'    testNames = c("trt"))
#' stats<-pipeDESeq2(counts=counts, info=info,
#'    contrastFormula = " ~ trt + geno",
#'    contrasts = list(c("geno","HAL2","FIL2")),
#'    contrastNames = "HvF")
#'
#' stats<-pipeDESeq2(counts=counts, info=info,)
#' stats<-pipeDESeq2(counts=counts, info=info,
#'    formula = c(" ~ trt + geno + trt:geno"," ~ trt + geno", " ~ trt + geno"),
#'    reduced = c(" ~ trt + geno"," ~ trt"," ~ geno"),
#'    testNames = c("GxE","geno","trt"),
#'    contrastFormula = " ~ trt + geno",
#'    contrasts = list(c("trt","Wet","Dry"),c("geno","HAL2","FIL2")),
#'    contrastNames = c("WvD","HvF"),returnNormData = "vst")
#' }
#'
#' @import DESeq2
#' @import edgeR
#' @importFrom stats as.formula model.matrix p.adjust
#' @export
pipeDESeq2<-function(counts, info,
                     formula = NULL, reduced = NULL, testNames = NULL,
                     contrasts = NULL, contrastFormula = NULL, contrastNames = NULL,
                     returnNormData = "none",
                     geneIDs=NA, verbose=TRUE,
                     fitType = "parametric", quiet = TRUE, betaPrior = F,
                     modelMatrixType = "standard", parallel = F,
                     minReplicatesForReplace = 7, ...){
  
  ######################
  # Set up the entire analysis
  
  if(!requireNamespace("DESeq2", quietly = TRUE)){
    stop("install the DESeq2 package to run this function\n")
  }else{
    require("DESeq2")
  }
  
  se<-SummarizedExperiment(assays = data.matrix(counts),
                           colData = DataFrame(info))
  if(is.na(geneIDs)){
    geneIDs<-rownames(se)
  }
  
  ######################
  # Set up LRT
  if(!is.null(formula)){
    if(length(formula)!=length(reduced)){
      stop("must have the same number of full(formula) and reduced models \n")
    }else{
      ntests<-length(formula)
    }
    if(!is.null(testNames) & length(testNames) != ntests){
      warning("if testnames are specified, must be the same length as formulae ...
              dropping testnames")
    }
    if(is.null(testNames)){
      testNames<-sapply(1:ntests, function(x){
        f<-gsub("~","",formula[x], fixed = T)
        r<-gsub("~","",reduced[x], fixed = T)
        for(i in c("+","-")) f<-gsub(i,"",f, fixed = T)
        for(i in c("+","-")) r<-gsub(i,"",r, fixed = T)
        tn<-gsub(r,"",f,fixed = T)
        tn<-gsub("*",".",tn,fixed = T)
        return(tn)
      })
    }
    
    if(verbose)  cat("running Likelihood ratio tests for results for:\n")
    
    runs<-lapply(1:ntests, function(x){
      if(verbose) cat(formula[x]," vs. ",reduced[x],"\n")
      if(quiet){
        suppressMessages(dds<- DESeqDataSet(se = se, design = as.formula(formula[x])))
        des<-DESeq(dds, test="LRT", full = as.formula(formula[x]),
                   reduced= as.formula(reduced[x]), quiet = T)
      }else{
        dds<- DESeqDataSet(se = se, design = as.formula(formula[x]))
        des<-DESeq(dds, test="LRT", full = as.formula(formula[x]),
                   reduced= as.formula(reduced[x]), quiet = F)
      }
      
      resname<-resultsNames(des)
      cat("possible results to output:", paste(resname, collapse = ", "))
      if(grepl(".",testNames[x], fixed = T)){
        res2get<-resname[grepl(".",resname, fixed = T) &
                           grepl(strsplit(testNames[x],".", fixed = T)[[1]][1],resname, fixed = T) & 
                           grepl(strsplit(testNames[x],".", fixed = T)[[1]][2],resname, fixed = T)]
      }else{
        res2get<-resname[grep(testNames[x],resname)]
      }
      cat(" (testing", res2get,")\n")
      res<-data.frame(gene=geneIDs, results(des, name = res2get),
                      stringsAsFactors = F)
      res$padj<-p.adjust(res$pvalue, method = "fdr")
      colnames(res)[-1]<-paste0(colnames(res)[-1],"_",testNames[x])
      return(res)
    })
    
    runs.comb<-Reduce(function(x, y)
      merge(x, y, by = "gene", all = TRUE), runs)
    
    }else{
      runs.comb<-NULL
    }
  
  
  ######################
  # Set up contrasts test
  if(!is.null(contrasts)){
    ncontrasts<-length(contrasts)
    if(is.null(contrastFormula)){
      stop("specify a formula for the design matrix")
    }
    if(verbose)  cat("running contrast tests for:",contrastFormula,"\n")
    if(is.null(contrastNames)){
      contrastNames<-1:ncontrasts
    }
    
    if(quiet){
      suppressMessages(dds<- DESeqDataSet(se = se, design = as.formula(contrastFormula)))
      des<-DESeq(dds, test="Wald",  quiet = T)
    }else{
      dds<- DESeqDataSet(se = se, design = as.formula(contrastFormula))
      des<-DESeq(dds, test="Wald",  quiet = F)
    }
    
    
    conts<-lapply(1:ncontrasts, function(x){
      if(verbose) cat(contrasts[[x]],"\n")
      
      res<-data.frame(gene=geneIDs,
                      results(des, contrast  = contrasts[[x]], ...),
                      stringsAsFactors = F)
      res$padj<-p.adjust(res$pvalue, method = "fdr")
      colnames(res)[-1]<-paste0(colnames(res)[-1],"_",contrastNames[x])
      return(res)
    })
    
    cont.comb<-Reduce(function(x, y)
      merge(x, y, by = "gene", all = TRUE), conts)
  }else{
    cont.comb = NULL
  }
  
  if(returnNormData == "rlog"){
    if(verbose) cat("conducting rlog tranformation\n")
    normDat <- rlog(counts)
  }else{
    if(returnNormData == "vst"){
      if(verbose) cat("conducting variance stabilizing tranformation\n")
      normDat <- varianceStabilizingTransformation(counts)
    }else{
      normDat = NULL
    }
  }
  return(list(LRT_results = runs.comb,
              CONT_results = cont.comb,
              normData = normDat))
  }
