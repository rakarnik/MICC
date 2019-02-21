###########################################################################################
########### MICCMain.R
########### functions:
###########            InputMatrixFormatted: process input PET clusters into MICC main model input
###########            MICCMainLearn: fit model parameters
###########            FDRcompute: compute FDR
###########			   MICCoutput: implement MICC model and get output 
###########################################################################################
####=============================================================================

InputMatrixFormatted <- function(data) {
	x <- data
	cAB <- x[,7]
	cA <- x[,8]
	cB <- x[,9]
	intra <- as.character(x[,1])==as.character(x[,4])
	inter <- as.character(x[,1])!=as.character(x[,4])
	distance <- (x[,5]+x[,6]-x[,3]-x[,2]) / 2 / 1000
	distance[inter] <- Inf
	data_formatted <- list( cAB=cAB, cA=cA, cB=cB, distance=distance )
}

MICCMainLearn <- function(data_formatted, params.init=NULL, reltol=1e-5, abstol=1e-3, step=200, restart=5, MinConfident=5) {
	Par <- EMIter( data_formatted, params.init=params.init, reltol=reltol, abstol=abstol, step=step, restart=restart, MinConfident=MinConfident )
	Par
}

FDRcompute <- function( data_formatted, params, PostProb ) {
	fdr <- FDRestimate( data_formatted, params, PostProb )
	fdr
}

MICCoutput <- function( data, outfilename, params.init=NULL, reltol=1e-5, abstol=1e-3, step=200, restart=5, MinConfident=5 ) {
	data_formatted <- InputMatrixFormatted(data)
	Par <- MICCMainLearn( data_formatted, params.init=params.init, reltol=reltol, abstol=abstol, step=step, restart=restart, MinConfident=MinConfident )
	params <- Par$params
	PostProb <- Par$PostProb
	fdr <- FDRcompute( data_formatted, params, PostProb[,1] )
	output.colnames <- c("chr.", "start", "end", "chr.", "start", "end", "cAB", "cA", "cB", "-log10(1-PostProb)", "fdr")
	y <- cbind(data, -log10(1-PostProb[,1]), fdr)
	colnames(y) <- output.colnames
	write.table( y, file=outfilename, sep="\t", row.names=F, quote=F )
}
	
	