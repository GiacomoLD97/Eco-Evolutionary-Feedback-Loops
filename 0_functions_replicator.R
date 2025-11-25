
# given a set of parameters, fit a random model with perfect intransitivity vs. hierarchy
run_single_sim <- function(j, seqmat, C1, C2, max_count){
	
	# set the seed
	randv<-sample.int(1e9, 1)
	set.seed(randv)
	
	# initialize the output
	my_tibb <- tibble()
	
	# get the parameters for this fit
	my_seqmat <- seqmat %>% slice(j)
	
	# sample a rho and p
	rho <- my_seqmat$r
	p <- my_seqmat$p
	
	n<- my_seqmat$n #sample(3:10, 1)
	m<-replicate(n, my_seqmat$m) #sample(minn:maxn, 1, replace=T))
	
	# get the submatrices
	my_C1 <- C1[[(my_seqmat$m-1)/2]]
	my_C2 <- C2[[(my_seqmat$m-1)/2]]
	
	# get phenotypes for each
	membership<-get_membership(m)
	
	# construct the matrices
	H<-H0<-construct_H_matrix(m, rho=rho)
	Q<-Q0<-construct_Q_matrix(m, p=p, range_p = 0, range_q = 0, ignore_min = F)
	B<-make_species_sum_matrix(m)
	
	
	# initialize with the intransitive vs. hierarchical matrix
	Hin <- Htr <- H
	for(i in 1:n){
		Htr[membership == i, membership == i] <- my_C1
		Hin[membership == i, membership == i] <- my_C2
		diag(Hin) <- diag(Htr) <- 0.5
	}
	
	# get the species-level matrix
	Hsum <- get_Hsum(H, membership)
	
	# convert to payoff matrices
	P <- H - t(H)
	Psum <- Hsum - t(Hsum)
	Pin <- Htr - t(Hin)
	Ptr <- Hin - t(Htr)
	
	# run the simulations, first time intransitive second time hierarchical
	for(z in 1:2){
		
		if(z == 1){
			Huse <- Hin
		}else{
			Huse <- Htr
		}
		
		# get the starting conditions
		xst<-find_optimal_strategy(Huse,verbose=F)+0.01
		xst<-xst/sum(xst)
		
		# do the initial integration to get past transients
		resm<-integrate_dynamics(x0=xst,pars=list(H=Huse,Q=Q),maxtime=500,lengthtime=500,maxsteps=1e6,method="lsoda",thresh=THRESH_sim,func_use=zero_sum_mutator)
		
		#now integrate for a bit to get oscillations
		resm<-integrate_dynamics(x0=resm[nrow(resm),-1],pars=list(H=Huse,Q=Q),maxtime=10000,lengthtime=1000,maxsteps=1e6,method="lsoda",thresh=THRESH_sim,func_use=zero_sum_mutator)
		soln<-apply(resm,2,mean)[-1]
		
		# continually integrate and average until we find the solution
		count<-0
		skip<-FALSE
		soln_old <- soln+100
		while((any(abs(100*(point_dxdt(soln,Huse,Q) / soln)[soln > THRESH])>0.0001) & any(point_dxdt(soln,Huse,Q) > THRESH_prune)) | any(abs(soln_old - soln) > 1e-7)){
			soln_old <- soln
			if(count==max_count){
				break()
			}
			count<-count+1
			resm<-integrate_dynamics(x0=soln,pars=list(H=Huse,Q=Q),maxtime=1000,lengthtime=500,maxsteps=1e6,method="lsoda",thresh=THRESH_sim,func_use=zero_sum_mutator)
			soln<-apply(resm,2,mean)[-1]
			alive<-soln>THRESH_prune
			alive_tab<-table(membership[alive])
		}
		
		
		# prune the system and get the eigenvalue
		H1<-Huse[alive,alive]
		Q1<-Q[alive,alive]
		soln1<-soln[alive]
		soln1<-soln1/sum(soln1)
		ev1<-Re(eigen(Jacobian_H(soln1,as.matrix(Q1),as.matrix(H1)))$values[-1])
		membership1 <- membership[alive]
		n_final <- length(unique(membership1))
		
		
		# calculate the metrics
		if(n_final<3){
			my_bab_sum <- 0
		}else{
			Hsum1 <- get_Hsum(H1, membership1)
			my_bab_sum <- calc_babington(Hsum1)
		}

		my_tibb <- my_tibb %>% bind_rows(tibble(id = randv, sim = j, iter = kind + offset,
												count = count,
												model = ifelse(z == 1, "intran", "tran"),
												n_init = n,
												n_final = n_final,
												m = mean(m),
												rho = rho,
												p = p,
												init_bab_sum = calc_babington(Hsum),
												final_bab_sum = my_bab_sum))
		
	}
	return(my_tibb)
}


# the main mutator equation, which collapses to the standard eqn when Q=I
zero_sum_mutator<- function(time, x, params){

	with(as.list(params), {

		x[x<THRESH]<-0

		x<-x/sum(x)

		dxdt<-2 * Q%*%diag(as.numeric(x))%*%H%*%(x) - x

		return(list(dxdt))
	})
}

# calculate the instantaneous rate of change. used for checking the tolerance solution
point_dxdt<- function(x,H,Q){

	dxdt<-2*Q%*%diag(x)%*%H%*%x - x

	return(dxdt)
}

# # calculate the Jacobian
Jacobian_H<-function(x,Q,H){
	D<-nrow(Q)
	if(D>1){
		2*Q%*%(diag(x)%*%(H)+diag(as.numeric((H)%*%x)))-diag(D)
	}
	else{
		return(as.matrix(1))
	}
}

# integrate the dynamics given a specified function
integrate_dynamics<- function(x0=NULL,pars,func_use,maxtime=1000,maxsteps=10000,lengthtime=100,thresh=1e-7,method="ode45",rtol=NULL,atol=NULL){

	if(is.null(x0)){
		x0<-runif(nrow(pars$H))
	}

	x0<-x0/sum(x0)

	times <- seq(0, maxtime, length=lengthtime)

	pars$THRESH<-thresh

	if(!is.null(rtol) & !is.null(atol)){
		out <- as.matrix(ode(x0, times, func_use, pars,method=method, maxsteps = maxsteps,rtol=rtol,atol=atol))
	}
	else{
		out <- as.matrix(ode(x0, times, func_use, pars,method=method, maxsteps = maxsteps))
	}

	out<-out[!is.na(rowSums(as.matrix(out[,2:ncol(out)]))),]

	out[out<thresh]<-0
	return(out)
}

# create a complementary matrix
make.comp<-function(X){
	X[lower.tri(X)]<-1-t(X)[lower.tri(t(X))]
	return(X)
}

# this creates a matrix which sums all of values within species
make_species_sum_matrix<-function(m){

	nm<-sum(m)
	n<-length(m)

	B<-matrix(0,n,nm)

	B[1,1:m[1]]<-1
	end<-m[1]
	for(i in 2:n){
		start<-end+1
		end<-start+m[i]-1
		B[i,start:end]<-1
	}

	return(B)
}

# solving for optimal solution using linear programming, assuming uncorrelated phenotypes

find_optimal_strategy<- function(H,verbose=TRUE,time_limit=5000){
    n <- dim(H)[1]
    f.obj <- rep(1, n)
    In<-diag(n)
    f.con <- H
    f.rhs <- rep(1, n)
    f.dir <-rep("<=", n)
    z<-Rglpk_solve_LP(obj=f.obj,mat=f.con,dir=f.dir,rhs=f.rhs,max=TRUE,control=list(tm_limit=time_limit,presolve=TRUE,verbose=verbose))
    return(z$solution / sum(z$solution))
}

# calculate a membership vector from the m=(m1,m2,...mn) vector
get_membership<-function(m){

	## create an expanded vector of memberships
	membership<-NULL
	n<-length(m)
	for(i in 1:n){
		membership<-c(membership,rep(i,m[i]))
	}

	return(membership)
}

# calculate the Babington measure of intransitivity
calc_babington <- function(myH, weighted = FALSE){
	if(nrow(myH) < 3){
		return(0)
	}
	if(!weighted){
		myH[myH > 0.5] <- 1
		myH[myH < 0.5] <- 0
	}
	diag(myH) <- NA
	n <- nrow(myH)
	dmax <- ifelse(n%%2 == 0, (n^3 - 4*n)/24, (n^3 - n)/24)
	d <- choose(n, 3) - sum(rowSums(myH, na.rm = T) *(rowSums(myH, na.rm = T) - 1)/2)
	return(round(d/dmax, 3))
}

# get the species-level H matrix
get_Hsum <- function(H, membership){
	usp <- unique(membership)
	n <- length(usp)
	Hsum <- matrix(0, n, n)
	diag(H) <- 0.5
	for(i in 1:n){
		for(j in i:n){
			Hsum[i,j] <- mean(H[membership == usp[i], membership == usp[j]])
			Hsum[j,i] <- 1 - Hsum[i,j]
		}
	}
	diag(Hsum) <- 0.5
	return(Hsum)
}


# give it a membership vector, a rho between zero and one, and should it have 0.5 between species?
construct_H_matrix<-function(m,rho=0,return_zipped=F,rho_within=NULL){

	n<-length(m)
	nm<-sum(m)


	if(rho<0 | rho>1){
		stop("rho must be between zero and one")
	}

	if(is.null(rho_within)){
		rho_within<-rho
	}

	# H is the full phenotype matrix, Hz is the zipped-up average matrix
	H<-Hz<-NULL

	# adjust the resource vector so as to maintain the correlation but also have the target probability
	for(i in 1:n){
		rowmat<-zrow<-NULL

		# note that this calcualtes values for all entries, not just the upper triangle
		# then calls make.comp to make sure that H+t(H)=1
		for(j in 1:n){

			# select the average value for this block
			rho_use<-rho

			# first see if we are in the diagonal, if so use rho_within
			if(i==j){
				rho_use<-rho_within
				s1<-0.5
			}
			else{
				s1<-runif(1)
			}

			#get the maximum dist from the boundary
			mv<-max(1-s1,s1)

			# sample uniformly, with the bounds scaled by rho_use
			subblock<-matrix(runif(m[i]*m[j],max(0,s1-mv*(1-rho_use)),min(1,s1+mv*(1-rho_use))),m[i],m[j])

			# bind the new columns to the previous ones
			rowmat<-cbind(rowmat,subblock)
 			# get the zipped up group average
 			zrow<-cbind(zrow,matrix(mean(subblock),1,1))
		}

		# bind the new rows to the previous ones
		H<-rbind(H,rowmat)
		Hz<-rbind(Hz,zrow)
	}


	# make them complementary
	H<-as.matrix(make.comp(H))
	Hz<-as.matrix(make.comp(Hz))

	diag(H)<-diag(Hz)<-.5

	if(return_zipped){
		return(list(H=H,Hz=Hz))
	}
	else{
		return(H)
	}

}


# construct a Q matrix with variable entries if desired
construct_Q_matrix<-function(m,p,range_p=0.1,range_q=0.1,min_above=0.05,max_p=0.95,ignore_min=FALSE){

	membership<-get_membership(m)
	max_p<-max(max_p,p)

	n<-length(m)
	nm<-sum(m)

	# set the lower bound for the diagonals
	min_vals<-rep(0,nm)
	if(!ignore_min){
		min_vals<-rep(1/m,m)
	}

	lower<-rep(p,nm)-range_p/2
	low_lim<-min_vals+min_above
	lower[lower<low_lim]<-low_lim[lower<low_lim]

	# set the upper bounds
	upper<-lower+range_p
	upper[upper>max_p]<-max_p
	lower[lower>max_p]<-max_p

	# check to see if there is only one species
	upper[rep(m,m)==1]<-lower[rep(m,m)==1]<-1

	# now sample the diagonals from a uniform dist
	Q<-matrix(0,nm,nm)
	for(i in 1:nm){
		Q[i,i]<-ifelse(upper[i]-lower[i]==0,upper[i],runif(1,lower[i],upper[i]))
	}

	# warming message in case a condition was missed
	if(any(diag(Q)<min_vals) & !ignore_min){
		stop("diagonal is less than 1/m")
	}

	# now add the off diagonals
	for(i in 1:nm){

		if(diag(Q)[i]<1){

			# relatives of the diagonal species
			relatives<-membership==membership[i] & (1:nm)!=i

			# equal weighted q
			qvals<-(1-diag(Q)[i])/(rep(m,m)[i]-1)

			if(range_q>0){
				qvals<-2
				count<-0
				# sample then normalize
				while(any(qvals>diag(Q)[i]) & !ignore_min){
					if(count>20){
						qvals<-(1-diag(Q)[i])/(rep(m,m)[i]-1)
						break()
					}
					qvals<-runif(sum(relatives),max(0,qvals-range_q),qvals+range_q)
					qvals<-qvals/sum(qvals)*(1-diag(Q)[i])
					count<-count+1
				}
			}

			if(any(qvals>diag(Q)[i]) & !ignore_min){
				stop("One of the off diagonals is larger than the diagonal")
			}

			Q[relatives,i]<-qvals
		}
	}

	if(!all(round(colSums(Q),12)==1)){
		stop("Q is not stochastic, something went wrong...")
	}


	return(Q)
}
