library(splines)
load("biomass.Rdata")
Z.save = matrix(0,length(biomass),5)

for(g in 1:length(biomass)){
u_given = biomass[g]
u = c(rep(range(biomass)[1],4),quantile(biomass,c(.5)),rep(range(biomass)[2],4))

if(biomass[g] < u[5]){
N_0 = c(0,0,0,1,0,0,0,0)	
}

if(biomass[g] > u[5]){
N_0 = c(0,0,0,0,1,0,0,0)	
}

N_1 = rep(0,7)
N_2 = rep(0,6)
N_3 = rep(0,5)

for(i in 1:length(N_1)){
	p = 1
	if(N_0[i]==0 && N_0[i+1]==0){
			N_1[i] = 0
	}
	if(N_0[i]!=0 && N_0[i+1]==0){
			N_1[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N_0[i]
	}
	if(N_0[i]==0 && N_0[i+1]!=0){
			N_1[i] = ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N_0[i+1]
	}

}

for(i in 1:length(N_2)){
	p = 2
	if(N_1[i]==0 && N_1[i+1]==0){
			N_2[i] = 0
	}
	if(N_1[i]!=0 && N_1[i+1]!=0){
			N_2[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N_1[i] +
	         ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N_1[i+1]
	}
	if(N_1[i]!=0 && N_1[i+1]==0){
			N_2[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N_1[i]
	}
	if(N_1[i]==0 && N_1[i+1]!=0){
			N_2[i] = ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N_1[i+1]
	}

}

for(i in 1:length(N_3)){
	p = 3
	if(N_2[i]==0 && N_2[i+1]==0){
			N_3[i] = 0
	}
	if(N_2[i]!=0 && N_2[i+1]!=0){
			N_3[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N_2[i] +
	         ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N_2[i+1]
	}
	if(N_2[i]!=0 && N_2[i+1]==0){
			N_3[i] = ((u_given - u[i]) / (u[i+p] - u[i])) * N_2[i]
	}
	if(N_2[i]==0 && N_2[i+1]!=0){
			N_3[i] = ((u[i+p+1] - u_given) / (u[i+p+1] - u[i+1]) )* N_2[i+1]
	}

}


Z.save[g,] = N_3
}

Z.knots = bs(biomass,degree = 3,knots=quantile(biomass,c(.5)),intercept=TRUE)

round(Z.save - Z.knots) #should be all zeros

