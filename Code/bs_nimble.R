# in BUGS code, to calculate the vector of basis matrix values for a given biomass, pass that biomass in as 'u_given', pass in the vector of u values for the knots and pass in N0,N1,N2,N3 of correct length - you can do this simply by providing N0,N1,N2,N3 as part of the 'constants' argument given to the 'nimbleModel' function

#p is the degree of the basis function. p = 3 is a cubic spline

#The following website helped with the creation of this function http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-basis.html

bs_nimble <-  nimbleFunction(
    run = function(u_given = double(0), u = double(1), N0 = double(1), N1 = double(1),
     N2 = double(1), N3 = double(1)) {
        returnType(double(1))
        
        if(u_given < u[1] | u_given > u[length(u)]){
        	print("u_given outside boundary")
        }
        
        for(i in 1:(length(u)-1))
            N0[i] = 0
        for(i in 1:(length(u)))
            N1[i] = 0
        for(i in 1:(length(u)+1))
            N2[i] = 0
        for(i in 1:(length(u)+2))
            N3[i] = 0
            
        for(i in 1:(length(u)-1)){
        	if(u_given >= u[i] & u_given < u[i+1]){
        		N0[i] = 1
        	} else {
        		N0[i] = 0
        	}
        }

        N1[1] = ((u[2] - u_given ) / (u[2] - u[1]) )* N0[1]
        
        for(i in 1:(length(N1)-1)){
            p = 1
            first = ((u_given - u[i]) / (u[i+p] - u[i])) * N0[i]
            first[is.na(first)]<-0
            second = ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N0[i+1]
            second[is.na(second)]<-0
            N1[i+1] = first + second
        }
        
        N1[is.na(N1)] <- 0 
        
        N2[1] = ((u[2] - u_given) / (u[2] - u[1])) * N1[1]
        N2[2] = ((u_given - 0) / (u[2] - 0)) * N1[1] + ((u[3] - u_given) / (u[3] - u[1])) * N1[2]
        
        for(i in 1:(length(N2)-2)){
            p = 2
            u.keep1 = u[i+p]
            u.keep1[is.na(u.keep1)] <- u[length(u)]
            first = ((u_given - u[i]) / (u.keep1- u[i])) * N1[i+1]
            first[is.na(first)]<-0
            u.keep = u[i+p+1]
            u.keep[is.na(u.keep)] <- u[length(u)]
            second = ((u.keep - u_given) / (u.keep - u[i+1])) * N1[i+2]
            second[is.na(second)]<-0
            N2[i+2] = first + second
        }
        
        N3[1] = ((u[2] - u_given) / (u[2] - 0) )* N2[1]
        N3[2] = ((u_given - 0) / (u[2] - 0)) * N2[1] + ((u[3] - u_given) / (u[3] - 0) )* N2[2]
        N3[3] = ((u_given - 0) / (u[3] - 0)) * N2[2] + ((u[4] - u_given) / (u[4] - u[1]) )* N2[3]
        
        for(i in 1:(length(N3)-3)){
            p = 3
            u.keep1 = u[i+p]
            u.keep1[is.na(u.keep1)] <- u[length(u)]
            first = ((u_given - u[i]) / (u.keep1- u[i])) * N2[i+2]
            first[is.na(first)]<-0
            u.keep = u[i+p+1]
            u.keep[is.na(u.keep)] <- u[length(u)]
            second = ((u.keep - u_given) / (u.keep - u[i+1])) * N2[i+3]
            second[is.na(second)]<-0
            N3[i+3] = first + second           
        }

        return(N3)
    })

# example

### This code calculates the cubic spline for one value at a time. Use loop to calculate more than one vector.
u <- seq(0,200,length.out=4) #u is the number of knots must be greater than 5 so length.out is an integer greater than 5. Must also surround all possible values for u_given
u_given <- 46

bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
bs(x = u_given, knots = u[2:(length(u)-1)], Boundary.knots = c(0,200),intercept=TRUE,degree=3)
splineDesign(knots=u, x=u_given, ord = 4, outer.ok=TRUE) #test against this function

### plot to think about how the splines work
knots <- u
x <- seq(min(knots)-1, max(knots)+1, length.out = 501)
bb <- splineDesign(knots, x = x, outer.ok = TRUE)

plot(range(x), c(0,1), type = "n", xlab = "x", ylab = "",
     main =  "B-splines - sum to 1 inside inner knots")
mtext(expression(B[j](x) *"  and "* sum(B[j](x), j == 1, 6)), adj = 0)
abline(v = knots, lty = 3, col = "light gray")
abline(v = knots[c(4,length(knots)-3)], lty = 3, col = "gray10")
lines(x, rowSums(bb), col = "gray", lwd = 2)
matlines(x, bb, ylim = c(0,1), lty = 1)
rug(x)

#### example nimble use
code <- nimbleCode({
   u_given ~ dunif(0,200)
   Z[1:7] <- bs_nimble(u_given, u=u[1:5], N0 = N0[1:4],N1 = N1[1:5], N2 = N2[1:6], N3 = N3[1:7]) #Z size depends on number of knots length(Z) = length(u) - 4
                                        #You must write out indexing for u and Z
})

m <- nimbleModel(code, inits = list(u_given = u_given),
                 constants = list(u = u, N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2))))

                     
spec.pred <- configureMCMC(m, thin = 10, print = TRUE)
spec.pred$addMonitors(c('Z')) 
Rmcmc.pred <- buildMCMC(spec.pred)

cm <- compileNimble(m)
Cmcmc.pred <- compileNimble(Rmcmc.pred, project = m) #Error in cModel

Cmcmc.pred$run(10000)
samples.pred <- as.matrix(Cmcmc.pred$mvSamples)

head(samples.pred)

boxplot(samples.pred)
hist(samples.pred)

