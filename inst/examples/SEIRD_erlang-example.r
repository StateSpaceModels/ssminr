SEIRD_inputs <- list(
	input(name="N",description="population size", value=1e7, tag="pop_size"),
	input(name="S", description="initial number of susceptible", tag="remainder"),
	input(name="R", description="initial number of recovered", value=0),
	input(name="D", description="initial number of dead", value=0),
	input(name="d_incubation", description="incubation period", value=10),
	input(name="d_onset2death", description="time from onset to death", value=7),
	input(name="d_onset2recovery", description="time from onset to recovery", value=15),
	input(name="cfr", description="case fatality ratio", value=0.7),
	input(name="I", description="initial number of infectious", prior=unif(1,1000)), 
	input(name="R0", description="basic reproduction number", prior=unif(0,5)), 
	input(name="rho", description="reporting rate", prior=truncnorm(mean=0.5,sd=0.1,a=0, b=1)),
	input(name="phi", description="overdispersion", prior=unif(0,1)),
	input(name="vol", description="volatility on log(beta)", prior=unif(0,1)),
	input(name="d_infectious", description="average infectious period",transformation="d_onset2death*cfr + d_onset2recovery*(1-cfr)"),
	input(name="E", description="initial number incubating",transformation="I*d_incubation/d_infectious"),
	input(name="I_D", description="initial number of infectious deamed to die", transformation="I*cfr"), 
	input(name="I_R", description="initial number of infectious deamed to recover", transformation="I*(1-cfr)"), 
	input(name="beta", description="effective contact rate",transformation="R0/d_infectious", sde=diffusion(volatility="vol", transformation="log")),
	input(name="epsilon", description="onset of infection rate",transformation="1/d_incubation"),
	input(name="gamma", description="recovery rate",transformation="1/d_onset2recovery"),
	input(name="mu", description="death rate",transformation="1/d_onset2death")
	)

SEIRD_reactions <- list(
	reaction(from="S", to="E", description="infection", rate="beta*(I_R + I_D)/N", keywords="transmission"),
	reaction(from="E", to=c("I_D"="cfr","I_R"="1-cfr"), description="onset of infectiosity", rate="epsilon", accumulators="incidence"),
	reaction(from="I_R", to="R", description="recovery", rate="gamma"),
	reaction(from="I_D", to="D", description="death", rate="mu")
	)

Erlang_shapes <- c(E=2, I_D=4, I_R=3)

SEIRD_observations <- list(
	discretized_normal_obs(state="incidence", reporting="rho", overdispersion="phi")
	)

data(ebola_2014)
data <- liberia1 %>% gather(time_series, value, -date)

# the model will be created in the default temporary directory. Change the path to "wherever/you/want".
# dir_model <- tempdir()
dir_model <- path.expand("~/Desktop")

my_ssm <- new_ssm(
	model_path=file.path(dir_model,"SEIRD_erlang"),
	pop="Liberia",
	data=data,
	start_date=min(liberia1$date) - 7, # start model integration 7 days before the first observation
	inputs=SEIRD_inputs,
	reactions=SEIRD_reactions,
	observations=SEIRD_observations,
	erlang_shapes=Erlang_shapes
	)

# plot_model(my_ssm, display="network",  collapse_erl=FALSE)
plot_model(my_ssm, display="diag", collapse_erl=FALSE)

# # Have fun.. 

# my_ssm_fit_ode <- my_ssm %>% pmcmc(iter=100000)
# my_ssm_fit_ode <- my_ssm %>% simplex(iter=1000) %>% pmcmc(iter=100000)
# my_ssm_fit_sde <- my_ssm_fit_ode %>% ksimplex(iter=1000) %>% kmcmc(iter=1000)
# my_ssm_fit_psr <- my_ssm_fit_sde %>% pmcmc(id=1, approx="psr", n_parts=100, iter=1000, n_thread="max")

# # LHS example

# my_ssm_lhs <- my_ssm %>% do_lhs(n=20, do="simplex", trace=FALSE, iter=100, prior=TRUE) %>% get_max_lhs

# print(my_ssm_lhs)

# my_ssm_mcmc <- my_ssm_lhs %>% simplex(iter=1000) %>% pmcmc(iter=1000) %>% to_tracer

# plot_data(my_ssm_lhs)

# my_ssm_mcmc %>% plot_X(stat="median", hat=c(0.95, 0.5))





