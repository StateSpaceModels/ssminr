# define a diffusion on log(beta) with volatility parameter "vol"	
input(name="beta", description="effective contact rate",transformation="R0/d_infectious", sde=diffusion(volatility="vol", transformation="log"))
