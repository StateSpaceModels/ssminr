SEIR_inputs <- list(
	input(name="N",description="population size", value=1e7, tag="pop_size"),
	input(name="S", description="initial number of susceptible", tag="remainder"),
	input(name="R", description="initial number of recovered", value=0),
	input(name="d_incubation", description="incubation period", value=10),
	input(name="d_infectious", description="infectious period", value=10),
	input(name="I", description="initial number of infectious", prior=unif(1,1000)), 
	input(name="R0", description="basic reproduction number", prior=unif(0,5)), 
	input(name="rho", description="reporting rate", prior=unif(0,1)),
	input(name="phi", description="overdispersion", prior=unif(0,1)),
	input(name="vol", description="volatility on log(beta)", prior=unif(0,1)),
	input(name="E", description="initial number incubating",transformation="I*d_incubation/d_infectious"),
	input(name="beta", description="effective contact rate",transformation="R0/d_infectious"),
	input(name="epsilon", description="onset of infection rate",transformation="1/d_incubation"),
	input(name="gamma", description="recovery rate",transformation="1/d_infectious")
	)

SEIR_reactions <- list(
	reaction(from="S", to="E", description="infection", rate="beta*I/N", accumulators="incidence"),
	reaction(from="E", to="I", description="onset of infectiosity", rate="epsilon"),
	reaction(from="I", to="R", description="recovery", rate="gamma")
	)

SEIR_observations <- list(
	poisson_obs(state="incidence", reporting="rho")
	)

data(ebola_2014)

# the model will be created in the default temporary directory. Change the path to "wherever/you/want".
dir_model <- tempdir()
# dir_model <- path.expand("~/Desktop")

build_ssm(
	model=file.path(dir_model,"SEIR_ssm"),
	pop_name="Liberia",
	data=liberia1,
	start_date=min(liberia1$date) - 7, # start model integration 7 days before the first observation
	inputs=SEIR_inputs,
	reactions=SEIR_reactions,
	observations=SEIR_observations
	)
