randomTrials <- function(count = 100, varCount = 2, range = matrix(c(0,0,5,3),2,2)) {
        trialStack <- NULL
        for (i in 1:count) {
		trial <- NULL
		for (j in 1:varCount)
			trial <- c( trial, runif(1,range[j,1],range[j,2]) )
		trialStack <- rbind(trialStack,trial)
	}
        trialStack
}

inputFunction1 <- function(values) {
	4*values[1]**2 + 4*values[2]**2
}

inputFunction2 <- function(values) {
	(values[1]-5)**2 + (values[2]-5)**2
}

constraints <- function(values, epsilon=0.0001) {
	max ( (abs(values[2] - values[1]**2) - epsilon) , 0 )
}

fitness <- function(values, func=inputFunction, constraints=TRUE, constrFunc=constraints, penalty=1) {
	if (constraints==FALSE)
		return(inputFunction(values))
	else 
		return( inputFunction(values) + penalty*constraints(values) )
}

getFitenssStack <- function(trialStacks, func= fitness, inputFunc=inputFunction, constraints=TRUE, constrFunc=constraints, penalty=1) {
        fitnessStack <- NULL
        for (i in 1:nrow(trialStacks))
                fitnessStack <- c( fitnessStack, func(trialStacks[i,], inputFunc,constraints=TRUE , constrFunc, penalty) )
        fitnessStack
}

naturalSelection <- function( trialStacks , type = "tournament", func = fitness, inputFunc1=inputFunction1, inputFunc2=inputFunction2 ,constraints=TRUE ,constrFunc=constraints, penalty=1 ) {
        if ( type=="tournament" ) {
		solStack <- cbind(apply(trialStacks, 1, inputFunc1), apply(trialStacks, 1, inputFunc2))
		trialSolStack <- cbind(trialStacks, solStack)
		trialSolStack <- getMultObjTable (trialSolStack)
		cols <- ncol(trialSolStack)
		rows <- nrow(trialSolStack)
		ranks <- unique(trialSolStack[,cols])
		crowdDist <- NULL
		for (i in ranks) {
			rankRows <- which(trialSolStack[,cols]==i)
			crowdDist <- c( crowdDist , getCrowdDistance(trialSolStack[rankRows,c(cols-2,cols-1)]) )
		}
		trialSolStack <- cbind(trialSolStack,crowdDist)
                winners <- NULL
                for (round in 1:rows) {
                        pair <- sample(rows,2)
                        if ( trialSolStack[ pair[1], cols ] < trialSolStack[ pair[2] , cols ] )
                                winners <- rbind(winners,trialSolStack[pair[1],1:ncol(trialStacks)])
			else if ( (trialSolStack[ pair[1], cols ] == trialSolStack[ pair[2] , cols ]) )
				if ( (trialSolStack[ pair[1], cols+1 ] > trialSolStack[ pair[2] , cols+1 ]) )
                                	winners <- rbind(winners,trialSolStack[pair[1],1:ncol(trialStacks)])
				else
                                	winners <- rbind(winners,trialSolStack[pair[2],1:ncol(trialStacks)])
                        else
                                winners <- rbind(winners,trialSolStack[pair[2],1:ncol(trialStacks)])
                }
                return (winners)
        }
}

crossOver <- function(trialStacks, crossProb = 0.8 ,type = "wholeArithmetic") {
        if (type=="wholeArithmetic") {
                len <- nrow(trialStacks)
                digits <- ncol(trialStacks)
		alpha <- runif(1)
                for (round in 1:len) {
                        pair <- sample(len,2)
                        if (runif(1) < crossProb) {
                                temp <- trialStacks[ pair[1], ]*alpha + trialStacks[ pair[2], ]*(1-alpha)
                                trialStacks[ pair[1], ]  <- trialStacks[ pair[2], ]*alpha + trialStacks[ pair[1], ]*(1-alpha)
                                trialStacks[ pair[2], ]  <- temp
                        }
                }
                return(trialStacks)
        }
}

mutation <- function(trialStacks, mutate_prob=(1/nrow(trialStacks)) , type = "realMutation" , timestep, totalTimestep, beta=1) {
        if (type=="realMutation") {
                len <- nrow(trialStacks)
                digits <- ncol(trialStacks)
		mins <- apply(trialStacks,2,min)
		maxs <- apply(trialStacks,2,max)
                for (i in 1:len)
                        for (j in 1:digits)
                                if (runif(1) < mutate_prob)
					if (runif(1) < 0.5)
                                        	trialStacks[i,j] <- trialStacks[i,j] + ( maxs[j] - trialStacks[i,j] )*(1-runif(1)^(1-timestep/totalTimestep))^beta
					else
                                        	trialStacks[i,j] <- trialStacks[i,j] - ( trialStacks[i,j] - mins[j] )*(1-runif(1)^(1-timestep/totalTimestep))^beta
                return(trialStacks)
        }
}

GeneOpt <- function( inputFunc1=inputFunction1, inputFunc2=inputFunction2,  fitnessFunc =fitness, iterations=50, crossProb = 0.8, mutate_prob=0.001, count=100, varCount = 2, range = matrix(c(-1,-1,1,1),2,2), typeNS = "tournament", constraints=FALSE, constrFunc=constraints, penalty=1 , typeCross = "wholeArithmetic", typeMutate = "realMutation", beta=1, display=FALSE ) {

        trialStacks <- randomTrials(count, varCount, range)
        for ( i in 1:iterations) {
		trialStacks <- naturalSelection(trialStacks,typeNS,func,inputFunc1,inputFunc2 ,constraints,constrFunc,penalty) 
                trialStacks <- crossOver( trialStacks, crossProb = crossProb ,type = typeCross )
                trialStacks <- mutation ( trialStacks, mutate_prob=mutate_prob , type = typeMutate , i, iterations, beta )
		if (display==TRUE) {
			solStack <- cbind(apply(trialStacks, 1, inputFunc1), apply(trialStacks, 1, inputFunc2))
			trialSolStack <- cbind(trialStacks, solStack)
			trialSolStack <- getMultObjTable (trialSolStack)
			cols <- ncol(trialSolStack) 
			if (i==1) {
				minX <- min(trialSolStack[,cols-2]) 
				minY <- min(trialSolStack[,cols-1])
				maxX <- max(trialSolStack[,cols-2]) 
				maxY <- max(trialSolStack[,cols-1])
			}
			plot(trialSolStack[,cols-2], trialSolStack[,cols-1], pch=trialSolStack[,cols], col=trialSolStack[,cols], xlim= c(minX,maxX), ylim=c(minY,maxY),main='Pareto Front')
			readline()
		}
        }
	#sols <- getFitenssStack (trialStacks, fitnessFunc, inputFunc,constraints=FALSE,constrFunc , penalty)
        #min_f <- min( sols )
	#arg_min_f = unique ( trialStacks[ which ( sols == min_f ),] )
	#output
	#list( min = min_f , arg_min = arg_min_f)
	solStack <- cbind(apply(trialStacks, 1, inputFunc1), apply(trialStacks, 1, inputFunc2))
	trialSolStack <- cbind(trialStacks, solStack)
	trialSolStack <- getMultObjTable (trialSolStack)
	trialSolStack
}

#Answer <- GeneOpt( inputFunc=inputFunction, fitnessFunc =fitness, iterations=1000, crossProb = 0.8, mutate_prob=0.001, count=100, varCount = 2, range = matrix(c(-1,-1,1,1),2,2), typeNS = "tournament", constraints=TRUE, constrFunc=constraints, penalty=1 , typeCross = "wholeArithmetic", typeMutate = "realMutation", beta=1 )

#print(Answer)

getMultObjTable <- function(trialSolStack) {
	multObjTable <- NULL
	rank=1
	repeat {
		solStack <- trialSolStack[ , c( ncol(trialSolStack)-1,ncol(trialSolStack) ) ]
		rows <- nrow(solStack)
		toBeRemoved <- NULL
		for (i in 1:rows) {
			check=0
			for (j in 1:rows) {
				if (i!=j) { 
					if ( (length( unique(solStack[i,] > solStack[j,]) ) == 1) )  {
						if (unique(solStack[i,] > solStack[j,]) == TRUE)  { 
							check=1
							break
						}
					}
				}
			}
			if (check==0) {
				multObjTable <- rbind( multObjTable, c(trialSolStack[i,] , rank) )
				toBeRemoved <- c(toBeRemoved,i)
			}
		}
		solStack <- solStack[-toBeRemoved,]
		trialSolStack <- trialSolStack[-toBeRemoved,]
		rank = rank + 1
		if ( is.vector(solStack) || length(toBeRemoved) == rows) {
			break
		}
	}
	if (is.vector(solStack)) {
		solStack <- t(as.matrix(solStack))
		trialSolStack <- t(as.matrix(trialSolStack))
	}
	multObjTable <- rbind( multObjTable, cbind(trialSolStack,rep(rank,nrow(trialSolStack) )) )
	multObjTable
}

getCrowdDistance <- function(solStack) {
	if (is.vector(solStack))
		return (Inf)
	rows <- nrow(solStack)
	cols <- ncol(solStack)
	crowdDist <- NULL
	minX <- min(solStack[,1])
	maxX <- max(solStack[,1])
	for (i in 1:rows) {
		if ( (solStack[i,1] != minX) & (solStack[i,1] != maxX) ) {
			dist <- solStack[i,1] - solStack[,1]
			posInd <- which( dist == min(dist[dist>0]) )[1]
			negInd <- which( dist == max(dist[dist<0]) )[1]
			num <- 0
			for (j in 2:cols) 
				num <- num + abs(solStack[posInd,j] - solStack[negInd,j])
			num <- num + abs(solStack[posInd,1] - solStack[negInd,1])
			crowdDist <- c(crowdDist, num)
		}
		else 
			crowdDist <- c(crowdDist, Inf)
	}
	return (crowdDist)
}

answer <- GeneOpt( inputFunc1=inputFunction1, inputFunc2=inputFunction2,  fitnessFunc =fitness, iterations=10, crossProb = 0.8, mutate_prob=0.001, count=100, varCount = 2, range = matrix(c(-1,-1,1,1),2,2), typeNS = "tournament", constraints=FALSE, constrFunc=constraints, penalty=1 , typeCross = "wholeArithmetic", typeMutate = "realMutation", beta=1, display=TRUE ) 
