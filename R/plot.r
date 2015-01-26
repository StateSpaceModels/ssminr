plot_simu <- function(model) {

	if(0){
		model <- "SIR_ssm"
	}

	df_traj <- read.csv(file.path(dir_ssm,model,"simu","X_0.csv"))

	df_data <- read.csv(file.path(dir_ssm,model,"data","data.csv"))	

	df_traj_plot <- gather(df_traj, state, value, -date, -index) %>% mutate(date=as.Date(date))
	df_data_plot <- gather(df_data, state, value, -date) %>% mutate(date=as.Date(date))

	p <- ggplot(data=df_traj_plot, aes(x=date, y=value)) + facet_wrap(~state, scales="free_y")
	p <- p + geom_line(aes(group=index), col="red", alpha=0.5)
	p <- p +geom_point(data=df_data_plot)
	print(p)


}