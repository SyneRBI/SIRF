function obj_fun = make_Poisson_loglikelihood(ad, model)
%     Selects the objective function based on the acquisition data and acquisition
%     model types.
    obj_fun = mStir.PoissonLogLikelihoodWithLinearModelForMeanAndProjData();
    obj_fun.set_acquisition_data(ad);
end