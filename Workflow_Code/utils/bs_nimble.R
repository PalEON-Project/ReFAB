# in BUGS code, to calculate the vector of basis matrix values for a 
# given biomass, pass that biomass in as 'u_given', pass in the vector 
# of u values for the knots and pass in N0,N1,N2,N3 of correct length - 
# you can do this simply by providing N0,N1,N2,N3 as part of the '
# constants' argument given to the 'nimbleModel' function

# p is the degree of the basis function. p = 3 is a cubic spline

# The following website helped with the creation of this function 
# http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-basis.html

bs_nimble <-  nimbleFunction(
    run = function(u_given = double(0), u = double(1), N0 = double(1), N1 = double(1),
     N2 = double(1), N3 = double(1)) {
        returnType(double(1))
        
        #if(u_given < u[1] | u_given > u[length(u)]){
        #	print("u_given outside boundary")
        #	print(u_given)
        #}
        
        for(i in 1:(length(u)-1))
            N0[i] <- 0
        for(i in 1:(length(u)))
            N1[i] <- 0
        for(i in 1:(length(u)+1))
            N2[i] <- 0
        for(i in 1:(length(u)+2)) #change if you increase number of knots
            N3[i] <- 0
            
        for(i in 1:(length(u)-1)){
        	if(u_given >= u[i] & u_given < u[i+1]){
        		N0[i] <- 1
        	} else {
        		N0[i] <- 0
        	}
        }

        N1[1] <- ((u[2] - u_given ) / (u[2] - u[1]) )* N0[1]
        
        for(i in 1:(length(N1)-1)){
            p <- 1
            if((i+p)>length(u)){
            	u_keep1 <- u[length(u)]
            }else{
            	u_keep1 <- u[i+p]
            }
            	
            if(i > length(u)){
            	u_keep2 <- 0
            }else{
            	u_keep2 <- u[i]
            }	
            
            if((u_given - u_keep2)==0 | (u_keep1- u_keep2)==0 | N0[i]==0){
            	first <- 0
            } else {
            	first <- ((u_given - u_keep2) / (u_keep1- u_keep2)) * N0[i]
            }
            
            if((i+p+1)>length(u)){
            	u_keep <- 0
            }else{
            	u_keep <- u[i+p+1]
            }
            
            if((i+1)>length(u)){
            	u_keep3 <- 0
            }else{
            	u_keep3 <- u[i+1]
            }	
            
            if((i+1)>length(N0)){
            	Nsub <- 0
            }else{
            	Nsub <- N0[i+1]
            }
            
            if((u_keep - u_given) == 0 | (u_keep - u_keep3)==0 | Nsub==0){
            	second <- 0
            } else {
            	second <- ((u_keep - u_given) / (u_keep - u_keep3)) * Nsub
            }
   
            N1[i+1] <- first + second
        }
        
        N2[1] <- ((u[2] - u_given) / (u[2] - u[1])) * N1[1]
        N2[2] <- ((u_given - 0) / (u[2] - 0)) * N1[1] + ((u[3] - u_given) / (u[3] - u[1])) * N1[2]
        
        for(i in 1:(length(N2)-2)){
            p <- 2
            
            if((i+p)>length(u)){
            	u_keep1 <- u[length(u)]
            }else{
            	u_keep1 <- u[i+p]
            	}
            	
            if(i>length(u)){
            	u_keep2 <- u[length(u)]
            }else{
            	u_keep2 <- u[i]
            }	
            
            if((i+1)>length(N1)){
            	Nsub <- 0
            }else{
            	Nsub <- N1[i+1]
            }
            
            if((u_given - u_keep2)==0 | (u_keep1- u_keep2)==0 | Nsub==0){
            	first <- 0
            } else {
            	first <- ((u_given - u_keep2) / (u_keep1- u_keep2)) * Nsub
            }
            
            if((i+p+1)>length(u)){
            	u_keep <- u[length(u)]
            }else{
            	u_keep <- u[i+p+1]
            }
            
            if((i+1)>length(u)){
            	u_keep3 <- u[length(u)]
            }else{
            	u_keep3 <- u[i+1]
            }	
            
            if((i+2) > length(N1)){
            	Nsub1 <- 0
            }else{
             Nsub1 <- N1[i+2]
            }
            
            if((u_keep - u_given) ==0 | (u_keep - u_keep3)==0 | Nsub1==0){
            	second <- 0
            } else {
            	second <- ((u_keep - u_given) / (u_keep - u_keep3)) * Nsub1
            }
            
            
            N2[i+2] <- first + second
        }
        
        N3[1] <- ((u[2] - u_given) / (u[2] - 0) )* N2[1]
        N3[2] <- ((u_given - 0) / (u[2] - 0)) * N2[1] + ((u[3] - u_given) / (u[3] - 0) )* N2[2]
        
        
        
        N3[3] <- ((u_given - 0) / (u[3] - 0)) * N2[2] + ((u[length(u)] - u_given) / (u[length(u)] - u[1]) )* N2[3]
        
        for(i in 1:(length(N3)-3)){
            p <- 3
            
            if((i+p)>length(u)){
            	u_keep1 <- u[length(u)]
            }else{
            	u_keep1 <- u[i+p]
            }
            
            if(i>length(u)){
            	u_keep2 <- u[length(u)]
            }else{
            	u_keep2 <- u[i]
            }	
            
            if((i+2)>length(N2)){
            	Nsub <- 0
            }else{
             Nsub <- N2[i+2]
            }
            
            if((u_given - u_keep2)==0 | (u_keep1- u_keep2)==0 | Nsub==0){
            	first <- 0
            } else {
            	first <- ((u_given - u_keep2) / (u_keep1- u_keep2)) * Nsub
            }    
            
           if((i+p+1)>length(u)){
            	u_keep <- u[length(u)]
            }else{
            	u_keep <- u[i+p+1]
            }
            
            if((i+1)>length(u)){
            	u_keep3 <- u[length(u)]
            }else{
            	u_keep3 <- u[i+1]
            }	
            
            if((i+3)> length(N2)){
            	Nsub1 <- 0
            }else{
            	Nsub1 <- N2[i+3]
            }
            
            if((u_keep - u_given) ==0 | (u_keep - u_keep3)==0 | Nsub1==0){
            	second <- 0
            } else {
            	second <- ((u_keep - u_given) / (u_keep - u_keep3)) * Nsub1
            }
            
            N3[i+3] <- first + second           
        }

        return(N3)
    })

# example

if(FALSE){

### This code calculates the cubic spline for one value at a time. Use loop to calculate more than one vector.
u <- c(0, quantile(biomass,c(.25,.5,.75)),200) #u is the number of knots
u_given <- 196

bs_nimble(u_given, u=u, N0 = rep(0, 2),N1 = rep(0, 3), N2 = rep(0,4), N3 = rep(0, 5))
bs(x = u_given, knots = u[2:(length(u)-1)], Boundary.knots = c(0,200),intercept=TRUE,degree=3)
splineDesign(knots=u, x=u_given, ord = 2, outer.ok=TRUE) #test against this function

### plot to think about how the splines work
knots <- seq(0,200,length.out=10)
x <- seq(min(knots)-1, max(knots)+1, length.out = 501)
bb <- splineDesign(knots, x = x, outer.ok = TRUE,ord=4)

plot(range(x), c(0,1), type = "n", xlab = "x", ylab = "",
     main =  "B-splines - sum to 1 inside inner knots")
mtext(expression(B[j](x) *"  and "* sum(B[j](x), j == 1, 6)), adj = 0)
abline(v = knots, lty = 3, col = "light gray")
abline(v = knots[c(4,length(knots)-3)], lty = 3, col = "gray10")
lines(x, rowSums(bb), col = "gray", lwd = 2)
matlines(x, bb, ylim = c(0,1), lty = 1)
rug(x)

Z.knots.check <- matrix(NA, u[length(u)], length(u) + 2)

for(i in 1:u[length(u)]){
	u_given <- i
	Z.knots.check[i,] = bs_nimble(u_given, u=u, N0 = rep(0, (length(u)-1)),N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2)))
}

#### Plot basis functions ####
plot(Z.knots.check[,1],xlim=c(0,u[length(u)]),pch=19,ylim=c(0,1),xlab="Biomass")
for(i in 2:ncol(Z.knots.check)){
	points(Z.knots.check[,i],col=i,pch=19)
}
abline(v=u,lwd=2)
title("Basis Functions")

#### example nimble use
code <- nimbleCode({
   u_given ~ dunif(1,150)
   Z[1:5] <- bs_nimble(u_given, u=u[1:3], N0[1:2], N1[1:3], N2[1:4], N3[1:5]) #Z size depends on number of knots length(Z) = length(u) - 4
                                        #You must write out indexing for u and Z
})

m <- nimbleModel(code, inits = list(u_given = u_given),
                 constants = list(u = u, N0 = rep(0, (length(u)-1)), N1 = rep(0, (length(u))), N2 = rep(0, (length(u)+1)), N3 = rep(0, (length(u)+2))))

                   
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

}

