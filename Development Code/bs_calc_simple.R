
u_given = 1.5
u = c(1,1,1,1,2,3,3,3,3)

N_0 = c(0,0,0,1,0,0,0,0)
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


Z.test1 = bs(c(1.5,2.5),knots = c(2),Boundary.knots = c(1,3),degree = 1)
Z.test2 = bs(c(1.5,2.5),knots = c(2),Boundary.knots = c(1,3),degree = 2)
Z.test3 = bs(c(1.5,2.5),knots = c(2),Boundary.knots = c(1,3),degree = 3)

if(u_given == 2.5){
 Z.test1[2,]==N_1[4:5]
 Z.test2[2,]==N_2[3:5]
 Z.test3[2,]==N_3[2:5]	
}

if(u_given == 1.5){
 Z.test1[1,]==N_1[4:5]
 Z.test2[1,]==N_2[3:5]
 Z.test3[1,]==N_3[2:5]	
}


