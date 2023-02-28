functions {
  // Changepoint model
  real elec(int occ, real t, real e0_oc, vector tau_oc, real hslope_oc, 
         real cslope_oc, real e0_un, real tau_un, real hslope_un, 
         real cslope_un, real k) {
      real dly = 0;
  	  if (occ == 1) {
  	      dly = e0_oc;
  	      if (t < tau_oc[1]) {
  	          dly += hslope_oc*(t - tau_oc[1]);
  	      }
  	      else if (t > tau_oc[2]) {
  	         dly += cslope_oc*(t - tau_oc[2]);
  	      }
  	  }
  	 else {
  	     dly = e0_un;
  	      if (t < tau_un) {
  	          dly += hslope_un*(t - tau_un);
  	      }
  	      else {
  	        dly += cslope_un*(t - tau_un);
  	      }
  	 }
  	 return(dly);
  	 }
}
data {
  // This block declares all data which will be passed to the Stan model.
  int<lower=0> N_pre;        // number of days in pre-retrofit period
  vector[N_pre] t_pre;      // pre-retrofit daily avg temperature
  int occ_pre[N_pre];       // occupied = 1, unoccupied = 0
  vector[N_pre] e_pre;      // outcome energy vector (pre-retrofit)
  int<lower=0> N_post;       // number of days in post-retrofit period
  vector[N_post] t_post;     // post-retrofit daily avg temperature
  int occ_post[N_post];      // categorical variable for the week ends
  vector[N_post] e_post;     // outcome energy vector (reporting)
}

parameters {
  // This block declares the parameters of the model. There are 11 parameters plus the error scale sigma
  real e0_oc;           // baseline consumption, occupied days
  ordered[2] tau_oc;    // heating and cooling breakpoints, occupied days
  real hslope_oc;       // heating slope, occupied days
  real cslope_oc;       // cooling slope, occupied days
  real e0_un;           // baseline consumption, unoccupied days
  real tau_un;          // Single breakpoint for unoccupied days
  real hslope_un;       // heating slope, unoccupied days
  real cslope_un;       // cooling slope, unoccupied days
  real k;               // serial correlation parameter
  real<lower=0> sigma;  // error scale
}

model {
  // Assigning prior distributions on some parameters
  e0_oc ~ normal(800,150);
  hslope_oc ~ normal(-40,15);
  cslope_oc ~ normal(40,15);
  tau_oc ~ normal(8,5);
  tau_un ~ normal(10,5);
  e0_un ~ normal(400,150);
  hslope_un ~ normal(-40,15);
  cslope_un ~ normal(40,15);
  k ~ normal(0.6,1);
  // Observational model
  for (n in 2:N_pre) {
    e_pre[n] ~ normal(k*(e_pre[n-1]-elec(occ_pre[n-1], t_pre[n-1], 
                  e0_oc, tau_oc, hslope_oc, cslope_oc, e0_un, tau_un, hslope_un, cslope_un, k)) + 
                  elec(occ_pre[n], t_pre[n], 
                  e0_oc, tau_oc, hslope_oc, cslope_oc, e0_un, tau_un, hslope_un, cslope_un, k), sigma);
  }
}

generated quantities {
  vector[N_post] e_post_pred;
  real eactpost = 0;
  real savings = 0;
  real emodpost = 0;
  vector[N_post] resid;
  resid[1] = e_post[1] - elec(occ_post[1], t_post[1], 
             e0_oc, tau_oc, hslope_oc, cslope_oc, e0_un, tau_un, hslope_un, cslope_un, k);
  for (n in 2:N_post) {
  	resid[n] = k*resid[n-1] + normal_rng(0,sigma);
  }
  for (n in 1:N_post) {
  	e_post_pred[n] = elec(occ_post[n],t_post[n], 
  	   e0_oc, tau_oc, hslope_oc, cslope_oc, e0_un, tau_un, hslope_un, cslope_un, k) + resid[n];
  } 
  for (n in 1:N_post) {
	  eactpost += e_post[n];
	  emodpost += e_post_pred[n];
  }
  savings = emodpost - eactpost;
}
