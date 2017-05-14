# in BUGS code, to calculate the vector of basis matrix values for a given biomass, pass that biomass in as 'u_given', pass in the vector of u values for the knots and pass in N0,N1,N2,N3 of correct length - you can do this simply by providing N0,N1,N2,N3 as part of the 'constants' argument given to the 'nimbleModel' function

#p is the degree of the basis function. p = 3 is a cubic spline

#The following website helped with the creation of this function http://www.cs.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-basis.html

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
        for(i in 1:5) #change if you increase number of knots
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


