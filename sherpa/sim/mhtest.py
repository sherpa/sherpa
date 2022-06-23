import numpy as np

from sherpa.sim.mh import dmvnorm
from sherpa.utils.err import ParameterErr
from sherpa.astro.io import io_opt
from sherpa.utils import NoNewAttributesAfterInit, print_fields
from sherpa.stats import Cash, CStat, WStat
from sherpa.astro import ui

if io_opt == 'crates_backend':
    from pycrates import read_file, get_keyval, copy_colvals, get_col_names
else:
    from astropy.io import fits

__all__ = ('mhTest',)

class PCA:

    fields = ('filename', 'emethod', 'cmethod', 'errext', 'colnames',)
    
    def __init__(self, filename):
        self.filename = filename
        self.emethod = None
        self.cmethod = None
        self.errext = None
        self.bias = None
        self.colnames = None
        self.component = None
        self.fvariance = None
        self.eigenval = None
        self.eigenvec = None
        return

    def __str__(self):
        print('bias = ', self.bias)
        print('component = ', self.component)
        print('fvariance = ', self.fvariance)
        print('eigenval = ', self.eigenval)
        print('eigenvec = ', self.eigenvec.shape, self.eigenvec)
        return print_fields(self.fields, vars(self))

class PCAcrates(PCA):

    def __init__(self, filename, blockname):
        PCA.__init__(self, filename)
        pc = read_file(filename)
        self.emethod = get_keyval(pc, 'EMETHOD')
        self.cmethod = get_keyval(pc, 'CMETHOD')
        self.errext = get_keyval(pc, 'ERREXT')
        self.bias = copy_colvals(pc, 'BIAS')
        ext2 = read_file(filename + '[' + blockname + ']')
        self.colnames = get_col_names(ext2)
        self.component = copy_colvals(ext2, 'COMPONENT')
        self.fvariance = copy_colvals(ext2, 'FVARIANCE')
        self.eigenval = copy_colvals(ext2, 'EIGENVAL')
        self.eigenvec = copy_colvals(ext2, 'EIGENVEC')
        return

class PCAastropy(PCA):

    def __init__(self, filename, blockname):
        self.filename = filename
        hdl = fits.open(filename)
        self.emethod = hdl[0].header['EMETHOD']
        self.cmethod = hdl[0].header['CMETHOD']
        self.errext = hdl[0].header['ERREXT']
        self.bias = hdl[1].data['BIAS']
        self.colnames = hdl[blockname].columns.names
        self.component = hdl[blockname].data['COMPONENT']
        self.fvariance = hdl[blockname].data['FVARIANCE']
        self.eigenval = hdl[blockname].data['EIGENVAL']
        self.eigenvec = hdl[blockname].data['EIGENVEC']
        
# def dmvnorm(x, mu, sigma, log = True):
#     if np.min(np.linalg.eigvalsh(sigma)) <= 0 :
#         raise RuntimeError("Error: sigma is not positive definite")
#     if np.max(np.abs(sigma - sigma.T)) >= 1.0e-9 :
#         raise RuntimeError("Error: sigma is not symmetric")
        
#     logdens = - mu.size / 2.0 * np.log(2 * np.pi) - \
#         1 / 2.0 * np.log(np.linalg.det(sigma)) - 1 / 2.0 * \
#         np.dot(x - mu, np.dot(np.linalg.inv(sigma), x - mu))
#     if log:
#         return logdens
#     else:
#         dens = exp(logdens)
#         return dens

class mhTest(NoNewAttributesAfterInit):

    def __init__(self):
        NoNewAttributesAfterInit.__init__(self)

        
    def dmvt(self, x, mu, sigma, df, log = True, norm = False):

        if np.min(np.linalg.eigvalsh(sigma)) <= 0 :
            raise RuntimeError("Error: sigma is not positive definite")
        if np.max(np.abs(sigma - sigma.T)) >= 1.0e-9 :
            raise RuntimeError("Error: sigma is not symmetric")

        p = len(mu)
        logdens_unnorm = - 0.5 * np.log(np.linalg.det(sigma)) - \
            (df + p) / 2.0 * np.log(df + \
            np.dot(x - mu, np.dot(np.linalg.inv(sigma), x - mu)))
        if log :
            if norm:
                logdens = logdens_unnorm + \
                    gammln((df + p) / 2.) - gammln(df / 2.) - (p / 2.) * \
                    np.log(np.pi) + (df / 2.) * np.log(df)
                return logdens
            else:
                return logdens_unnorm        
        else:
            if norm:
                logdens = logdens_unnorm + gammln((df + p) / 2.) - \
                    gammln(df / 2.) - (p / 2.) * np.log(np.pi) + (df / 2.) * \
                    np.log(df)
                dens = np.exp(logdens)
                return dens
            else:
                dens_unnorm = np.exp(logdens_unnorm)
                return dens_unnorm

    # def get_parameter_info(id = None):
    #     """Returns the parameter information needed for calling mht.

    #     This routine will call covariance() if needed, but not
    #     fit().

    #     For now only works with a single dataset, and requires
    #     the Cash statistic.
    #     """

    #     if id == None:
    #         idval = ui.get_default_id()
    #     else:
    #         idval = id

    #     #
    #     # Employ a lot of safety checks
    #     #
    #     fr = ui.get_fit_results()
    #     if fr == None:
    #         raise RuntimeError("No fit results available!")
    #     if len(fr.datasets) != 1:
    #         msg = "Fit is for multiple datasets (%d) which we do not support!" % \
    #             fr.datasets
    #         raise RuntimeError(msg)
    #     if fr.datasets[0] != idval:
    #         msg = "Fit results are for dataset %s, not %s" % (fr.datasets[0], \
    #                                                           idval)
    #         raise RuntimeError(msg)
    #     if fr.statname != "cash":
    #         msg = "Fit was run using statistic = %s rather than cash!" % \
    #             fr.statname
    #         raise RuntimeError(msg)
    #     if not fr.succeeded:
    #         # Should use standard sherpa logging
    #         msg = "fit to dataset %d did not complete successfully:\n%s" % \
    #             (idval, fr.message)
    #         warnings.warn(msg)        

    #     cr = ui.get_covar_results()
    #     if cr == None or len(cr.datasets) != 1 or cr.datasets[0] != idval:
    #         # Should use standard sherpa logging
    #         print("Running covariance for dataset %d" % idval)
    #         _, f = ui._session._get_fit(idval)
    #         cr = f.est_errors()

    #     if cr.statname != "cash":
    #         msg = "Covariance was run using statistic = %s rather than cash!" % \
    #             cr.statname
    #         raise RuntimeError(msg)

    #     if len(fr.parnames) != len(cr.parnames):
    #         msg = "Number of parameters used in fit (%d)" % len(fr.parnames)
    #         msg += "and covariance (%d) analysis do not agree!\n" % \
    #             len(cr.parnames)
    #         raise RuntimeError(msg)
    #     for (p1, p2) in zip(fr.parnames, cr.parnames):
    #         if p1 != p2:
    #             msg = "Order of fit (%s)" % p1
    #             msg += "covariance parameters (%s) does not match" % p2
    #             raise RuntimeError(msg)
    #     for (pname, v1, v2) in zip(fr.parnames, fr.parvals, cr.parvals):
    #         if v1 != v2:
    #             msg = "Value of fit / covariance parameters does not match"
    #             msg += "for parameter %s: %g vs %g" % (pname, v1, v2)
    #             raise RuntimeError(msg)

    #     if isinstance(cr.extra_output, (np.ndarray,)) is False:
    #         msg = "get_covar_results has no .extra_output or it is None"
    #         raise RuntimeError(msg)

    #     #
    #     # Store the information, we explicitly copy all items to avoid
    #     # problems if fit / covariance are run again. This is done by converting
    #     # all tuples to numpy arrays, even for strings, and is actually
    #     # not needed.
    #     #
    #     out = {
    #         "dataset":  idval, 
    #         "npars":    len(fr.parnames), 
    #         "parnames": np.asarray(fr.parnames), 
    #         "parvals":  np.asarray(fr.parvals), 
    #         "parmins":  np.asarray(cr.parmins), 
    #         "parmaxes": np.asarray(cr.parmaxes), 
    #         "sigma":    cr.sigma, 
    #         "covar":    cr.extra_output.copy(), 
    #         "statval":  fr.statval
    #         }
    #     return out

    # Simulate an arf
    def sim_arf_alt(self, pca, defspecresp, rrin, rrsig, n = 100):
        ncomp = len(pca.component)
        ncomp = min(n, ncomp)
        newarf = pca.bias + defspecresp
        rrout = rrin + rrsig * np.random.standard_normal(ncomp)
        #print rrout
        for i in range(ncomp):
            tmp = rrout[i] * pca.eigenval[i] * pca.eigenvec[i, ]
            newarf = newarf + tmp
        return newarf, rrout

    def sim_arf(self, pca, defspecresp, n = 100): # Simulate an arf
        ncomp = len(pca.component)
        ncomp = min(n, ncomp)
        newarf = pca.bias + defspecresp
        rrout = np.random.standard_normal(ncomp)
        for i in range(ncomp):
            tmp = rrout[i] * pca.eigenval[i] * pca.eigenvec[i, ]
            newarf = newarf +  tmp
        return newarf, rrout

    def _set_par_vals(self, parnames, parvals):
        """Sets the paramaters to the given values"""

        for (parname, parval) in zip(parnames, parvals):
            ui.set_par(parname, parval)


    def __call__(self, myfit, parnames, mu, sigma, num_iter, df, PCAfname, \
                 fixARFname, ids, improvedBayes=False, num_within=10, \
                 num_subiter=10, p_M=.5, comp=8, p_M_arf=.5, sd_arf=.1, \
                 thin=1, scale=1, seed=123):

        """
        p_M is mixing proportion of MH draws in the mixture of MH and Metropolis parameter draws
        p_M = 0, all m draws; p_M = 1, all mh draws
        p_M_arf is mixing proportion of MH draws in the mixture of MH and Metropolis arf draws
        p_M_arf = 0, all m draws; p_M_arf = 1, all mh draws
        """

        # _, myfit = ui._session._get_fit(ids[0], ids[1:])        
        if not isinstance(myfit.stat, (Cash, CStat, WStat)):
            raise ValueError("Fit statistic must be cash, cstat or " +
                             "wstat, not %s" % fit.stat.name)
        
        if np.min(np.linalg.eigvalsh(sigma)) < 0.0:
            raise RuntimeError("Error: sigma is not positive semi-definite")

        fr = ui.get_fit_results()
        if ids == None:
            ids = fr.datasets
            nids = len(ids)
            if nids != 1:
                msg = "Fit for multiple datasets (%d) is not supported" % nids
                raise RuntimeError(msg)

        np.random.seed(seed)

    ########################
    ### Fixed ARF Method ###
    ########################

        iterations = np.zeros((num_iter + 1, len(mu)))
        statistics = np.zeros((num_iter + 1, 1))   

        current = np.copy(mu)
        self._set_par_vals(parnames, current)

        statistics[0] = - 0.5 * ui.calc_stat(*ids)
        iterations[0] = np.copy(mu)

        zero_vec = np.zeros(len(mu))

        for i in range(1, num_iter + 1, 1):

            ##########################################
            ### DRAW parameters GIVEN arf and data ###
            ##########################################  

            current = iterations[i - 1]
            q = np.random.chisquare(df, 1)[0]
            u = np.random.uniform(0, 1, 1)
            if u > p_M:

                while True:
                    try:
                        proposal = iterations[i - 1] + \
                            np.random.multivariate_normal(zero_vec, sigma) / \
                            np.sqrt(q / df)
                        self._set_par_vals(parnames, proposal)
                        break
                    except ParameterErr:
                        pass

                stat_temp = - 0.5 * ui.calc_stat(*ids)
                alpha = np.exp(stat_temp - statistics[i - 1])

            else:

                #MH jumping rule
                while True: 
                    try:
                        proposal = mu + \
                            np.random.multivariate_normal(zero_vec, sigma) / \
                            np.sqrt(q / df)
                        self._set_par_vals(parnames, proposal)
                        break
                    except ParameterErr:
                        pass

                stat_temp = - 0.5 * ui.calc_stat(*ids)
                alpha = \
                    np.exp(stat_temp + self.dmvt(current, mu, sigma, df) - \
                           statistics[i - 1] - self.dmvt(proposal, mu, sigma, \
                                                         df))

            u = np.random.uniform(0, 1, 1)
            if u <= alpha:
                ###print("accept new para")
                iterations[i] = np.copy(proposal)
                statistics[i] = np.copy(stat_temp)
            else:
                ###print("reject new para")
                iterations[i] = np.copy(iterations[i - 1])
                statistics[i] = np.copy(statistics[i - 1])

        result1 = np.hstack((statistics, iterations))

    ########################
    ### Pragmatic Method ###
    ########################

        iterations = np.zeros((num_iter + 1, len(mu)))
        statistics = np.zeros((num_iter + 1, 1))

        current = np.copy(mu)
        self._set_par_vals(parnames, current)

        stat_temp = - 0.5 * ui.calc_stat(*ids)
        statistics[0] = stat_temp

        # Setting for PCA
        defarf = ui.get_arf()
        defresp = ui.get_arf().specresp
        if io_opt == 'crates_backend':
            mypca = PCAcrates(PCAfname, 'PCACOMP')
        else:
            mypca = PCAastropy(PCAfname, 'PCACOMP')
        ncomp = len(mypca.component)
        ncomp = min(ncomp, comp)

        arfcomp = np.zeros((num_iter + 1, ncomp))

        iterations[0] = np.copy(mu) 
        zero_vec = np.zeros(len(mu))

        new_rr = np.zeros(ncomp)
        new_arf = defresp

        num_subiterations = num_subiter

        if improvedBayes == False:
            num_within = 1

        for i in range(1, num_iter + 1, 1):

            ##########################################
            ### DRAW arf GIVEN parameters and data ###
            ##########################################

            defarf.specresp = new_arf

            if i%num_within == 0:

                new_arf, new_rr = self.sim_arf(mypca, defresp, n = ncomp)
                defarf = ui.get_arf()
                defarf.specresp = new_arf
                stat_temp = - 0.5 * ui.calc_stat(*ids)

                # run silent run deep
                # _, f = ui._session._get_fit(ids[0], ids[1:])
                # fit_result =  f.fit()
                fit_result = myfit.fit()
                covar = myfit.est_errors()
                mu = np.asarray(fit_result.parvals)
                # fit()
                # covariance()
                # tmp = get_parameter_info()
                # mu = np.copy(tmp["parvals"])
                # covar = np.copy(tmp["covar"])

            arfcomp[i] = np.copy(new_rr)

            ##########################################
            ### DRAW parameters GIVEN arf and data ###
            ##########################################  


            current = iterations[i - 1]
            self._set_par_vals(parnames, current)
            stat_prev = - 0.5 * ui.calc_stat(*ids)

            for ii in range(num_subiterations):

                q = np.random.chisquare(df, 1)[0]
                u = np.random.uniform(0, 1, 1)

                if u > p_M:
                        #Metropolis jumping rule
                    while True:
                        try:
                            proposal = current + \
                                np.random.multivariate_normal(zero_vec, \
                                                              sigma) / \
                                                              np.sqrt(q / df)
                            self._set_par_vals(parnames, proposal)
                            break
                        except ParameterErr:
                            pass

                    stat_temp = - 0.5 * ui.calc_stat(*ids)
                    alpha = np.exp(stat_temp - stat_prev)

                else:
                    #MH jumping rule
                    while True: 
                        try:
                            proposal = mu + \
                                np.random.multivariate_normal(zero_vec, \
                                                              sigma) / \
                                                              np.sqrt(q / df)
                            self._set_par_vals(parnames, proposal)
                            break
                        except ParameterErr:
                            pass

                    stat_temp = - 0.5 * ui.calc_stat(*ids)
                    alpha = np.exp(stat_temp + \
                                   self.dmvt(current, mu, sigma, df) - stat_prev - \
                                   self.dmvt(proposal, mu, sigma, df))

                u = np.random.uniform(0, 1, 1)
                if u <= alpha:
                    ###print("accept new para")
                    iterations[i] = np.copy(proposal)
                    statistics[i] = np.copy(stat_temp)
                    current = proposal
                    stat_prev = stat_temp
                else:
                    ###print("reject new para")
                    iterations[i] = np.copy(current)
                    statistics[i] = np.copy(stat_prev)

        result2 = np.hstack((statistics, iterations, arfcomp))

    ########################
    ### FullyBays Method ###
    ########################

        target = np.hstack((iterations, arfcomp))
        target = target.transpose()
        bigsigma = np.cov(target)
        bigmean = target.mean(1)

        mu1 = bigmean[:len(mu)]
        sig11 = bigsigma[:len(mu), :len(mu)]
        sig12 = bigsigma[:len(mu), len(mu):]
        sig22_inv = np.linalg.inv(bigsigma[len(mu):, len(mu):])
        con_var = scale * (sig11 - np.dot(np.dot(sig12, sig22_inv), \
                                          sig12.transpose()))

        startindex = 1

        stat_prev = statistics[startindex]
        current = iterations[startindex]
        arfcomp_current = arfcomp[startindex]

        iterations = np.zeros((num_iter + 1, len(mu)))
        statistics = np.zeros((num_iter + 1, 1))
        arfcomp = np.zeros((num_iter + 1, ncomp))

        iterations[0] = current
        statistics[0] = stat_prev
        arfcomp[0] = arfcomp_current

        con_prob_prev = 100

        mu0  = np.repeat(0, ncomp)
        sig0 = np.diag(np.repeat(1, ncomp))
        prior_prev = dmvnorm(arfcomp_current, mu0, sig0)
        # oldarf = arf

        # index = 0
        for i in range(1, num_iter + 1, 1):
            for j in range(0, thin, 1):
                u = np.random.uniform(0, 1, 1)
                if u > p_M_arf:

                    if np.min(np.linalg.eigvalsh(con_var)) < 0.0:
                        raise RuntimeError("Error: con_var is not positive semi-definite")

                    ##Mh jumping rule
                    new_arf, arfcomp_proposal = self.sim_arf(mypca, defresp, n = ncomp)
                    defarf = ui.get_arf()
                    defarf.specresp = new_arf
                    while True:
                        try:
                            tmp_mean = mu1 + np.dot(np.dot(sig12, sig22_inv), \
                                                    arfcomp_proposal.transpose())
                            proposal = np.random.multivariate_normal(tmp_mean, \
                                                                     con_var)
                            self._set_par_vals(parnames, proposal)
                            break
                        except ParameterErr:
                            pass

                    con_prob_proposal = \
                        dmvnorm(proposal,
                                mu1 + np.dot(np.dot(sig12, sig22_inv), \
                                             arfcomp_proposal.transpose()), \
                                con_var)
                    prior_proposal = dmvnorm(arfcomp_proposal, mu0, sig0)
                    stat_proposal = - 0.5 * ui.calc_stat(*ids)
                    alpha = np.exp(stat_proposal - stat_prev + con_prob_prev - \
                                   con_prob_proposal)

                else:
                    #M jumping rule
                    new_arf, arfcomp_proposal = \
                        self.sim_arf_alt(mypca, defresp, arfcomp_current, sd_arf, \
                                    n = ncomp)
                    defarf = ui.get_arf()
                    defarf.specresp = new_arf
                    while True:
                        try:
                            tmp_mean = mu1 + np.dot(np.dot(sig12, sig22_inv), \
                                                    arfcomp_proposal.transpose())
                            proposal = \
                                np.random.multivariate_normal(tmp_mean, con_var)
                            self._set_par_vals(parnames, proposal)
                            break
                        except ParameterErr:
                            pass

                    con_prob_proposal = \
                        dmvnorm(proposal, \
                                mu1 + np.dot(np.dot(sig12, sig22_inv), \
                                             arfcomp_proposal.transpose()), \
                                con_var)
                    prior_proposal = dmvnorm(arfcomp_proposal, mu0, sig0)
                    stat_proposal = - 0.5 * ui.calc_stat(*ids)
                    alpha = np.exp(stat_proposal - stat_prev + con_prob_prev - \
                                   con_prob_proposal + prior_proposal - prior_prev)

                u = np.random.uniform(0, 1, 1)
                if u <= alpha:
                    #print("accept new para")
                    current = proposal
                    arfcomp_current = arfcomp_proposal
                    stat_prev = stat_proposal
                    con_prob_prev = con_prob_proposal
                    prior_prev = prior_proposal

            iterations[i] = np.copy(current)
            statistics[i] = np.copy(stat_prev)
            arfcomp[i] = np.copy(arfcomp_current)

        result3 = np.hstack((statistics, iterations, arfcomp))
        result = np.hstack((result1, result2, result3))

        return result
#
# https://github.com/astrostat/pyblocxs
#
