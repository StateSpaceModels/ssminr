SIR_inputs <- list(
	input(name="N",description="population size", value=1e7, tag="pop_size"),
	input(name="S", description="initial number of susceptible", tag="remainder"),
	input(name="R", description="initial number of recovered", value=0),
	input(name="d_infectious", description="infectious period", value=10),
	input(name="I", description="initial number of infectious", prior=unif(1,1000), value=10), 
	input(name="vol", description="vol on beta", prior=unif(0,1), value=0.1), 
	input(name="R0", description="basic reproduction number", prior=unif(0,50), value= 2), 
	input(name="rho", description="reporting rate", prior=truncnorm(mean=0.5,sd=0.1,a=0, b=1)),
	input(name="phi", description="overdispersion", prior=unif(0,1)),
	# input(name="beta", description="effective contact rate",transformation="R0/d_infectious"),
	input(name="beta", description="effective contact rate",transformation="R0/d_infectious", sde=diffusion(drift="0", volatility="vol", transformation="log")),
	# input(name="beta", description="effective contact rate",transformation="R0/d_infectious", sde=diffusion(drift="-0.1", volatility="vol*beta")),
	# input(name="beta", description="effective contact rate",transformation="R0/d_infectious", sde=diffusion(drift="-exp(-(t-25))/(1+exp(-(t-25)))^2", volatility="vol", transformation="log")),
	input(name="gamma", description="recovery rate",transformation="1/d_infectious")
	)

SIR_reactions <- list(
	reaction(from="S", to="I", description="infection", rate="beta*I/N", accumulators="incidence"),
	reaction(from="I", to="R", description="recovery", rate="gamma")
	)

SIR_observations <- list(
	poisson_obs(state="incidence", reporting="rho")
	)

data(ebola_2014)
data <- liberia2 %>% gather(time_series, value, -date)

Erlang_shapes <- c(I=1)


# the model will be created in the default temporary directory. Change the path to "wherever/you/want".
# dir_model <- tempdir()
# dir_model <- path.expand("~/Desktop")
dir_model <- "/Users/Tonton/work/presentations/talks/2016_02_Princeton/"

my_ssm <- new_ssm(
	model_path=file.path(dir_model,"SIR_ssm"),
	pop="Liberia",
	data=data,
	start_date=min(liberia1$date) - 7, # start model integration 7 days before the first observation
	inputs=SIR_inputs,
	reactions=SIR_reactions,
	observations=SIR_observations,
	erlang_shapes = Erlang_shapes
	)

if(0){
	print(my_ssm)
	plot_data(my_ssm)
	plot_model(my_ssm, collapse_erl=FALSE)

# # Have fun.. 
	my_ssm %>% simul(approx="psr", n_parts=100, dt=0.1) %>% plot_X(stat="median", hat=c(0.95, 0.5))
	my_ssm %>% kalman(hat = TRUE) %>% plot_X(stat="median", hat = c(0.5, 0.95))

# my_ssm_fit_ode <- my_ssm %>% simplex(iter=1000) %>% print %>% pmcmc(iter=1000) %>% print %>% to_tracer %>% plot_X
my_ssm_fit_sde <- my_ssm %>% ksimplex(iter=1000)

 # %>% kmcmc(iter=10000, n_traj = 500) %>% plot_X(stat="median", hat=c(0.95, 0.5))
# my_ssm_fit_psr <- my_ssm_fit_sde %>% pmcmc(id=1, approx="psr", n_parts=100, iter=1000, n_thread="max")
	my_ssm_fit_psr <- my_ssm %>% pmcmc(id=0, approx="psr", n_parts=300, iter=1000, n_thread="max")

	my_ssm_fit_psr %>% print
	my_ssm_fit_psr %>% plot_X(stat="median", hat=c(0.95, 0.5))
	my_ssm_fit_psr %>% to_tracer

# # # LHS example
my_ssm_lhs <- my_ssm %>% do_lhs(n=200, do="ksimplex", trace=FALSE, iter=100, prior=TRUE) %>% get_max_lhs %>% print

# # my_ssm_mcmc <- my_ssm_lhs %>% pmcmc(iter=100000, n_traj=1000, switch=50, eps_switch=100) 
# # my_ssm_mcmc <- my_ssm_mcmc %>% to_tracer
}




