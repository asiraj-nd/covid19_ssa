##########################################################################

" 
  This file is part of the Preprint https://doi.org/10.1101/2020.04.07.20053421.
  
  Written by Amir S Siraj
  February 2020
##########################################################################

  This script sets up pomp model c++ snippets in R function formats
  using string parcing approach
  Written by: Amir S Siraj
  February 2020
  This work is part of the study 
##########################################################################
- Specific tasks
            # setup a stochastic SEIR model
            # declare pomp initial condition function
            # declare pomp process model function
            # declare pomp measurement model function
            # create the pomp object function
            # setup simulation function
            
#================================================================
"
########################################################################################
################################### POMP Delarations  States ###########################
#time.in=given.time
declare.data = function (time.in){
  eval(parse(text = paste('
                          data.temp = data.frame(
                          time = ', paste ('1:', time.in, sep=''), ', ', paste('
                                                                               YNx', ' = NA,
                                                                               YEx', ' = NA,
                                                                               YIx', ' = NA,
                                                                               YRx', ' = NA,
                                                                               YCx', ' = NA',  sep = '', collapse = ' '),
                          ')', sep = '')))

  return(data.temp)
}

########################################################################################
################################### POMP Process  declaration #########################
declare.rproc = function() {

  # parse the declaration string
  rproc.string.var = paste('
                           double N', ';
                           double E', ';
                           double I', ';
                           double R', ';
                           double C', ';
                           double rates[', 5 , '];
                           double trans[', 5 , '];
                           double xvec[', 1000, '];
                           ', sep = '')

  # Initializing (R values in R to values in c++)
  # We need to do this because c++ doesn't allow for Array declaration in statenames...

  rproc.string.assign = paste( paste('
                              N = Nx', ';
                              E = Ex', ';
                              I = Ix', ';
                              R = Rx', ';
                              C = Cx', ';', sep ='', collapse = ''),
                               paste('xvec[', 0:999, '] = xvec_', 0:999, ';', sep ='', collapse = ''),
                               sep = '', collapse = '')

  rproc.string.trans = '
  // pre processor
    double S = 0;
    double Ieff = 0;
    // time point
    int ll = floor(t); int uu = ceil(t);

  //if (t<10)  printf(\"N = %.2f, R = %.2f , I = %.2f \\n\", N, R, I); // to debug uncomment this

    S = N - E - I - R;

    double a = 0.66666667;
    double b = 0.4;
    double c = 0.5;

    double u=(double)rand();
    double d=(double)RAND_MAX;
    double div= u/d;
    int pari=  (int) (div*1000);

    double thisxi = xvec [pari];
    double thisphi = thisxi * 2;

    // make sure phi >= gamma
    if (thisphi < gamma) thisphi = gamma;

    //if (t<10)  printf(\"N = %.2f, S = %.2f , E = %.2f , I = %.2f , R = %.2f   \\n\",  N,  S, E, I, R); // to debug uncomment this

    // Compute the rates
    rates[0] = beta * I/ N ;  // transmission from S
    rates[1] = thisxi;  // transition to E=>I
    rates[2] = theta * thisphi ; // recovery I through quarantine I=>R
    rates[3] = (1-theta) * gamma;  // infection from incubation period  I=>R

  // if (t<10)  printf(\" beta = %.2f, xi = %.2f , phi = %.2f, theta = %.2f, gamma = %.2f  \\n\", beta, thisxi, thisphi,theta,gamma); // to debug uncomment this

  // Compute the transitions
  trans[0] = rpois(rates[0] * S * dt);
  reulermultinom(1, E, &rates[1], dt, &trans[1]);
  reulermultinom(1, I, &rates[2], dt, &trans[2]);
  reulermultinom(1, I, &rates[3], dt, &trans[3]);

  // Bookkeeping
  /*
  We need to keep track of all the transitions in the E, I, R, and C compartments
  */
  E += trans[0] - trans[1] ;
  I += trans[1] - trans[2] - trans[3]  ;
  R += trans[2] + trans[3] ;
  C += trans[0];'
  
  
  # Finally we set the values back to the global definitions of statenames
  rproc.string.reassign = paste('
                                Nx', '= N;','
                                Ex', '= E;','
                                Ix', '= I;','
                                Rx', '= R;','
                                Cx', '= C;',
                                sep ='', collapse = '')

  # We declare the Csnippet function by concatenating the four parced strings 
  
  rproc = Csnippet(
    paste(rproc.string.var, rproc.string.assign, rproc.string.trans, rproc.string.reassign, sep = '\n'))
  return(rproc)
}



########################################################################################
################################### POMP Measurement declaration #######################

# set up rmeasurement function ==================================

declare.rmeas = function(){
  rmeas.string = paste('
                       YNx', ' = Nx', ';
                       YEx', ' = Ex', ';
                       YIx', ' = Ix', ';
                       YRx', ' = Rx', ';
                       YCx', ' = Cx', ';', sep = '', collapse = ' ')
  rmeas = Csnippet(rmeas.string)
  return(rmeas)
}

########################################################################################
################################### POMP initial conditions  declaration ###############

# set up initital values of state variables =====================
declare.init = function(N.init.in, E.init.in, I.init.in, R.init.in, C.init.in) {

  init.string = paste('
                      Nx', ' = ', N.init.in, ';
                      Ex', ' = ', E.init.in, ';
                      Ix', ' = ', I.init.in, ';
                      Rx', ' = ', R.init.in, ';
                      Cx',  ' = ', C.init.in, ';',sep = '', collapse = ' ')
  states_init = Csnippet(init.string)
  return(states_init)
  #return(init.string) to debug
}

########################################################################################
################################### declare POMP object ################################
# set up pomp object ============================================
declare.pomp = function(data.in, start.time.in, rproc.in, rmeas.in,
                        states_init.in, times.name.in = 'time', delta.t.in = 1/6){
  pomp.obj.temp = pomp(
    data = data.in,
    times = times.name.in,
    t0 = start.time.in,
    rprocess = discrete_time(rproc.in, delta.t.in),
    rmeasure = rmeas.in,
    rinit = states_init.in,
    statenames = c( 'Nx', 'Ex', 'Ix', 'Rx', 'Cx'),
    paramnames = c('beta', 'theta', 'phi', 'xi', 'tau', 'gamma',
                   paste('xvec_', 0:999, sep = ''))
  )
  return(pomp.obj.temp)
}

########################################################################################
################################### POMP simulate function #############################

simulate.pomp = function(pomp.in, beta.in, theta.in, phi.in, xi.in, tau.in, gamma.in, xi.vec.in, nsim.in){
  eval(parse(
    text = paste(
      'params.true.temp = c(beta = ', beta.in, ', theta = ', theta.in, ', phi = ', phi.in,
      ', xi = ', xi.in, ', tau = ', tau.in,
      c(paste(', xvec_', 0:999,'=', xi.vec.in, sep = '', collapse = '')),
      ', gamma = ', gamma.in, ')', sep = '')))
  sim.obj = simulate(
    pomp.in,
    params = params.true.temp,
    nsim=nsim.in)

  return(list(sim = sim.obj, params_true = params.true.temp))
}


