#	Welcome to R!
#	Statements preceded by "#" are comments,
#	"<-" commands assign the variable to the left with the value to the right
#	"()" contain the targets of functions
#	"{}" contain conditional statements or loops
#	"~" specifies a relationship in a model formula

#	first, clean and tidy workspace
rm(list = ls())
graphics.off()
#################################################################
#	Useful Packages												#
#################################################################
#NO GUARANTEE THIS WORKS IN Rstudio!
chooseCRANmirror(graphics=FALSE, ind = which(getCRANmirrors()$Country == 'Germany')) # 'SouthAfrica')[2])
install.packages('multcomp')#this one is proving elusive
install.packages('nparcomp')#this one is quite new
library(multcomp)#multiple comparisons for parametric models
library(nparcomp)#multiple comparisons for nonparametric models



#################################################################
#	Useful Functions											#
#################################################################
#lettering algorithms based on
#Piepho, H. P. (2004). An algorithm for a letter-based representation of all-pairwise comparisons. Journal of Computational and Graphical Statistics, 13(2), 456-466.
#All functions shamelessly pirated from https://rdrr.io/cran/multcomp/src/R/cld.R



#	"Insert absorb" subfunction from cld in the multcomp package	#
insert_absorb <- function( x, Letters=c(letters, LETTERS), separator=".", decreasing = FALSE, 
                           comps = NULL, lvl_order){

  obj_x <- deparse(substitute(x))
  if (is.null(comps)) {
      namx <- names(x)
      namx <- gsub(" ", "", names(x))
      if(length(namx) != length(x))
          stop("Names required for ", obj_x)
      split_names <- strsplit(namx, "-")
      stopifnot( sapply(split_names, length) == 2 )
      comps <- t(as.matrix(as.data.frame(split_names)))
  } 
  rownames(comps) <- names(x)
  lvls <- lvl_order
  n <- length(lvls)
  lmat <- array(TRUE, dim=c(n,1), dimnames=list(lvls, NULL) )

  if( sum(x) == 0 ){                                                        # no differences
    ltrs <- rep(get_letters(1, Letters=Letters, separator=separator), length(lvls) )
    names(ltrs) <- lvls
    colnames(lmat) <- ltrs[1]
    msl <- ltrs
    ret <- list(Letters=ltrs, monospacedLetters=msl, LetterMatrix=lmat)
    class(ret) <- "multcompLetters"
    return(ret)
  }
  else{
    signifs <- comps[x,,drop=FALSE]
    
    absorb <- function(m){
      for(j in 1:(ncol(m)-1)){
        for(k in (j+1):ncol(m)){
          if( all(m[which(m[,k]),k] & m[which(m[,k]),j]) ){                 # column k fully contained in column j
            m <- m[,-k, drop=FALSE]
            return(absorb(m))
          }
          else if( all(m[which(m[,j]),k] & m[which(m[,j]),j]) ){            # column j fully contained in column k
            m <- m[,-j, drop=FALSE]
            return(absorb(m))
          }
        }
      }
      return(m)
    }
    for( i in 1:nrow(signifs) ){                                            # insert
      tmpcomp <- signifs[i,]
      wassert <- which(lmat[tmpcomp[1],] & lmat[tmpcomp[2],])               # which columns wrongly assert nonsignificance
      if(any(wassert)){
        tmpcols <- lmat[,wassert,drop=FALSE]
        tmpcols[tmpcomp[2],] <- FALSE
        lmat[tmpcomp[1],wassert] <- FALSE
        lmat <- cbind(lmat, tmpcols)
        colnames(lmat) <- get_letters( ncol(lmat), Letters=Letters,
                                       separator=separator)
        if(ncol(lmat) > 1){                                                 # absorb columns if possible
          lmat <- absorb(lmat)
          colnames(lmat) <- get_letters( ncol(lmat),  Letters=Letters,
                                         separator=separator )
        }
      }
    }
  }
  lmat <- lmat[,order(apply(lmat, 2, sum))]
  lmat <- sweepLetters(lmat)                                                                  # 1st
  lmat <- lmat[,names(sort(apply(lmat,2, function(x) return(min(which(x))))))]                # reorder columns
  colnames(lmat) <- get_letters( ncol(lmat),  Letters=Letters,
                                 separator=separator)
  lmat <- lmat[,order(apply(lmat, 2, sum))]                                                   # 2nd sweep
  lmat <- sweepLetters(lmat)
  lmat <- lmat[,names(sort(apply(lmat,2, function(x) return(min(which(x)))), 
                           decreasing = decreasing))]                # reorder columns
  colnames(lmat) <- get_letters( ncol(lmat),  Letters=Letters,
                                 separator=separator)
  ltrs <- apply(lmat,1,function(x) return(paste(names(x)[which(x)], sep="", collapse="") ) )
  msl <- matrix(ncol=ncol(lmat), nrow=nrow(lmat))                                             # prepare monospaced letters
  for( i in 1:nrow(lmat) ){
    msl[i,which(lmat[i,])] <- colnames(lmat)[which(lmat[i,])]
    absent <- which(!lmat[i,])
    if( length(absent) < 2 ){
      if( length(absent) == 0 )
        next
      else{
        msl[i,absent] <- paste( rep(" ", nchar(colnames(lmat)[absent])), collapse="" )
      }
    }
    else{
      msl[i,absent] <- unlist( lapply( sapply( nchar(colnames(lmat)[absent]),
                                               function(x) return(rep( " ",x)) ),
                                       paste, collapse="") )
    }
  }
  msl <- apply(msl, 1, paste, collapse="")
  names(msl) <- rownames(lmat)
  ret <- list( Letters=ltrs, monospacedLetters=msl, LetterMatrix=lmat, 
               aLetters = Letters, aseparator = separator )
  class(ret) <- "multcompLetters"
  return(ret)
}#insert_absorb

#	"Sweep letters" subfunction from cld in the multcomp package#
sweepLetters <- function(mat, start.col=1, Letters=c(letters, LETTERS), separator="."){

  stopifnot( all(start.col %in% 1:ncol(mat)) )
  locked <- matrix(rep(0,ncol(mat)*nrow(mat)), ncol=ncol(mat))          # 1 indicates that another letter dependes on this entry
  cols <- 1:ncol(mat)
  cols <- cols[c( start.col, cols[-start.col] )]
  if( any(is.na(cols) ) )
    cols <- cols[-which(is.na(cols))]

  for( i in cols){
    tmp <- matrix(rep(0,ncol(mat)*nrow(mat)), ncol=ncol(mat))
    tmp[which(mat[,i]),] <- mat[which(mat[,i]),]                        # get items of those rows which are TRUE in col "i"
    one <- which(tmp[,i]==1)

    if( all(apply(tmp[,-i,drop=FALSE], 1, function(x) return( any(x==1) ))) ){     # there is at least one row "l" where mat[l,i] is the only item which is TRUE i.e. no item can be removed in this column
      next
    }
    for( j in one ){                                                    # over all 1's
      if( locked[j,i] == 1 ){                                           # item is locked
        next
      }
      chck <- 0
      lck <- list()
      for( k in one ){
        if( j==k ){
          next
        }
        else{                                                           # pair j-k
          rows <- tmp[c(j,k),]
          dbl <- rows[1,] & rows[2,]
          hit <- which(dbl)
          hit <- hit[-which(hit==i)]
          dbl <- rows[1,-i,drop=FALSE] & rows[2,-i,drop=FALSE]
          if( any(dbl) ){
            chck <- chck + 1
            lck[[chck]] <- list(c(j,hit[length(hit)]), c(k,hit[length(hit)]))      # record items which have to be locked, use last column if multiple hits
          }
        }
      }
      if( (chck == (length(one)-1)) && chck != 0 ){                     # item is redundant
        for( k in 1:length(lck) ){                                      # lock items
          locked[ lck[[k]][[1]][1], lck[[k]][[1]][2] ] <- 1
          locked[ lck[[k]][[2]][1], lck[[k]][[2]][2] ] <- 1
        }
        mat[j,i] <- FALSE                                               # delete redundant entry
      }
    }
    if(all(mat[,i]==FALSE)){                                           # delete column where each entry is FALSE and restart
      mat <- mat[,-i,drop=FALSE]
      colnames(mat) <- get_letters( ncol(mat), Letters=Letters, separator=separator)
      return(sweepLetters(mat, Letters=Letters, separator=separator))
    }
  }
  onlyF <- apply(mat, 2, function(x) return(all(!x)))
  if( any(onlyF) ){                                                     # There are columns with just FALSE entries
    mat <- mat[,-which(onlyF),drop=FALSE]
    colnames(mat) <- get_letters( ncol(mat), Letters=Letters, separator=separator)
  }
  return( mat )
}#sweepLetters



#	"Get letters" subfunction from cld in the multcomp package#
get_letters <- function( n, Letters=c(letters, LETTERS), separator="." ){

  n.complete <- floor(n / length(Letters))        # number of complete sets of Letters
  n.partial <- n %% length(Letters)               # number of additional Letters
  lett <- character()
  separ=""
  if( n.complete > 0 ){
    for( i in 1:n.complete ){
      lett <- c(lett, paste(separ, Letters, sep="") )
      separ <- paste( separ, separator, sep="" )
    }
  }
  if(n.partial > 0 )
    lett <- c(lett, paste(separ, Letters[1:n.partial], sep="") )
  return(lett)
}#get_letters


#####################################################################
#	Input Variables													#
#####################################################################

#I'll invent come data, but you can also read in tab sep data using:
# dta <- read.table(paste0(Sys.getenv('HOME'), '/Datafolder/', 'data.txt'), sep = '\t', header = T)
nn <- 20
#data here will be log-normal
#rnorm produces normal distribution with a mean at 0 and standard deviation of 1
#by giving different conditions different means, we can determine what differences there should be
#by adding a constant and logging the data, we ensure that this is no longer normally distributed
yvariable <- log10(10+c(rnorm(nn), #Matthew is very average
				rnorm(nn,mean =2.1), #Mark is above average
				rnorm(nn,mean =0.1), #Luke is like Matthew
				rnorm(nn,mean =0.7 #John is a bit like everyone
				))) 
condition <- as.factor(c(rep('Matthew', nn),
							rep('Mark', nn),
							rep('Luke', nn),
							rep('John', nn)))
dta <- data.frame(condition = condition, yvariable = yvariable)
#inspect data
stripchart(yvariable~condition, data = dta, method = 'jitter', vertical = T, pch = 20)# a lot of overlapping variances that would be hard to deal with by parametric methods

#####################################################################
#	Method for Parametric Models (for reference)					#
#####################################################################
#Courtesy of Altaf Hussain https://www.researchgate.net/post/How_to_denote_letters_to_mark_significant_differences_in_a_bar_chart_plot
#thankyou helpful strangers!
summary(mymodel<-aov(yvariable ~condition)) # if significant, do the following:
tuk <- glht(mymodel, linfct = mcp(condition = "Tukey"))
summary(tuk) # standard display
tuk.cld <- cld(tuk)
opar <- par(mai=c(1,1,1.5,1))
plot(tuk.cld) # letter-based display on a box and whisker plot

#####################################################################
#	Method for Nonparametric Models									#
#####################################################################
print(kwmodel<-kruskal.test(yvariable ~condition)) # if significant, do the following:
#supposedly the non-parametric version of glht (Gen.Linear Hypothesis Tests)
#"nonparametric relative contrasts"
#Konietschke, F., Brunner, E., Hothorn, L.A. (2008). Nonparametric Relative Contrast Effects: Asymptotic Theory and Small Sample Approximations.
# Munzel. U., Hothorn, L.A. (2001). A unified Approach to Simultaneous Rank Tests Procedures in the Unbalanced One-way Layout. Biometric Journal, 43, 553-569
np.tuk <- nparcomp(yvariable ~condition, data = dta)
summary(np.tuk)#some significant comparisons some not


#	let's look a little closer at how these comparisons were done
#this is the matrix of comparisons, we'll need this for later
np.tuk$Contrast#-1 is reference level, 1 is comparison, 0 ignored
#this is the comparison with estimated difference in mean rank
np.tuk$Analysis#it has a "$p.Value" column (note the capital V!)
with(np.tuk$Analysis,cbind(Comparison,p.Value = round(np.tuk$Analysis$p.Value,2)))

#	now set up input for "insert_absorb" algorithm that assigns letters
#input based on cld.summary.glht from cld in the multcomp package
# x <- factor
#by default, comparisons are conducted in alphabetical order
#see ?relevel and ?factor for methods to change this order
lvl_order <- levels(condition)#levels(ret$x)[order(tapply(as.numeric(ret$y)[1:length(ret$x)], ret$x, mean))]
#contrast matrix of comparisons
K <- np.tuk$Contrast#contrMat(table(x), type = "Tukey")
# matrix of comparisons with comparison and ref. level
comps <- cbind(apply(K, 1, function(k) levels(condition)[k == 1]),
               apply(K, 1, function(k) levels(condition)[k == -1]))
#identity of p values that are significant at 0.05
signif <- np.tuk$Analysis$p.Value <0.05#(object$test$pvalues < level)


#	run the insert & absorb letters algorithm with these inputs#
mcletters <- insert_absorb(signif,#significant differences
					decreasing = F, #letters will increase
                          comps = comps,#comparison matrix to search
                          lvl_order = lvl_order)#, ...)

#Everybody in the house make some noise for MC Letters!
print(mcletters$Letters)