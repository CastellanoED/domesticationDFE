import numpy as np, scipy.integrate
import dadi
from dadi.DFE import *
from dadi import Numerics, Integration, PhiManip, Spectrum

def Domestication_demography(params, ns, pts):

	new1,new2,ne1d,ne2d,Tbot,T3,m = params

	xx = Numerics.default_grid(pts)
	
	phi = PhiManip.phi_1D(xx) # phi for the eq. ancestral population
	phi = PhiManip.phi_1D_to_2D(xx, phi)
	phi = Integration.two_pops(phi, xx, Tbot, new1, ne1d, m12=m)
	phi = Integration.two_pops(phi, xx, T3,   new2, ne2d)

	fs = Spectrum.from_phi(phi, ns, (xx,xx))
	return fs

def Domestication_complex_demography(params, ns, pts):

    new1,new2,new3,ne1d,ne2d,ne3d,T1,T2,T3,md1,md2,md3,mw1,mw2,mw3 = params

    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx) # phi for the eq. ancestral population
    
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T1, new1, ne1d, m12=md1, m21=mw1) #time and size of the domesticated pop bottleneck 
    phi = Integration.two_pops(phi, xx, T2, new2, ne2d, m12=md2, m21=mw2) #time and size of the domesticated pop recovery 
    phi = Integration.two_pops(phi, xx, T3, new3, ne3d, m12=md3, m21=mw3) #time and size of the domesticated pop recovery 

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

def Domestication_flexible_demography(params, ns, pts):
    Tpre, nuPre, Tdiv, nu1div, nu2div, T1F, T2F, nu1F, nu2F, mw2d, md2w = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    # Change prior to divergence
    phi = Integration.one_pop(phi, xx, Tpre, nuPre)

    # Divergence 
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    # Time before either population changes size again
    phi = Integration.two_pops(phi, xx, Tdiv, nu1div, nu2div, m12=mw2d, m21=md2w)

    if T1F > T2F:
        # If population 1 changes 1st
        phi = Integration.two_pops(phi, xx, T1F-T2F, nu1F, nu2div, m12=mw2d, m21=md2w)
        phi = Integration.two_pops(phi, xx, T2F, nu1F, nu2F, m12=mw2d, m21=md2w)
    else:
        # If population 2 changes 1st
        phi = Integration.two_pops(phi, xx, T2F-T1F, nu1div, nu2F, m12=mw2d, m21=md2w)
        phi = Integration.two_pops(phi, xx, T1F, nu1F, nu2F, m12=mw2d, m21=md2w)

    fs = Spectrum.from_phi(phi, ns, (xx,xx)) 
    return fs
Domestication_flexible_demography.__param_names__ = ['Tpre', 'nuPre', 'Tdiv', 'nu1div', 'nu2div', 'T1F', 'T2F', 'nu1F', 'nu2F', 'mw2d', 'md2w']

if __name__ == "__main__":
    ns, pts = (10,10), 20
    test = Domestication_flexible_demography((0.2,2,0.1,0.2,0.4,0.3,0.2,0.9,0.1,1,1), ns, pts)

def Domestication_flexible_demography_cache(params, ns, pts):
    Tpre, nuPre, Tdiv, nu1div, nu2div, T1F, T2F, nu1F, nu2F, mw2d, md2w, gamma1, gamma2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx, gamma=gamma1)

    # Change prior to divergence
    phi = Integration.one_pop(phi, xx, Tpre, nuPre, gamma=gamma1)

    # Divergence 
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    # Time before either population changes size again
    phi = Integration.two_pops(phi, xx, Tdiv, nu1div, nu2div, m12=mw2d, m21=md2w, gamma1=gamma1, gamma2=gamma2)

    if T1F > T2F:
        # If population 1 changes 1st
        phi = Integration.two_pops(phi, xx, T1F-T2F, nu1F, nu2div, m12=mw2d, m21=md2w, gamma1=gamma1, gamma2=gamma2)
        phi = Integration.two_pops(phi, xx, T2F, nu1F, nu2F, m12=mw2d, m21=md2w, gamma1=gamma1, gamma2=gamma2)
    else:
        # If population 2 changes 1st
        phi = Integration.two_pops(phi, xx, T2F-T1F, nu1div, nu2F, m12=mw2d, m21=md2w, gamma1=gamma1, gamma2=gamma2)
        phi = Integration.two_pops(phi, xx, T1F, nu1F, nu2F, m12=mw2d, m21=md2w, gamma1=gamma1, gamma2=gamma2)

    fs = Spectrum.from_phi(phi, ns, (xx,xx)) 
    return fs
Domestication_flexible_demography_cache.__param_names__ = ['Tpre', 'nuPre', 'Tdiv', 'nu1div', 'nu2div', 'T1F', 'T2F', 'nu1F', 'nu2F', 'mw2d', 'md2w', 'gamma1', 'gamma2']

def Domestication_one_pop(params, ns, pts):
    newpars = list(params[:-1]) + [params[-1], params[-1]]
    return Domestication_flexible_demography_cache(newpars, ns, pts)

def Vourlaki_mixture(params, ns, s1, s2, theta, pts):
    """
    Inference model from Vourlaki et al.

    params = alpha, beta, ppos_wild, gamma_pos, pchange, pchange_pos
    """

    alpha, beta, ppos_wild, gamma_pos, pchange, pchange_pos = params

    # We'll scale by theta at the end of the script, so set theta=1 here.
    # Case in which gamma is negative and equal in the two pops
    m5 = s1.integrate([alpha, beta], None, PDFs.gamma, 1, None)
    # Case in which gamma is negative and indepenent in the two pops
    m6 = s2.integrate([alpha, beta], None, PDFs.biv_ind_gamma, 1, None,
                      exterior_int=True)

    # Cases in which gamma is positive in both pops.
    # For simplicity, using a single gamma value, rather than exponential dist
    #  that was simulated.
    try:
        m2 = s2.spectra[s2.gammas == gamma_pos,
                        s2.gammas == gamma_pos][0]
        # Even if selection coefficient changed to a different positive value,
        #   we just model them as equal.
        m3 = m2
    except IndexError:
        raise IndexError('Failed to find requested gamma_pos={0:.4f} '
                         'in cached spectra. Was it included '
                         'in additional_gammas during '
                         'cache generation?'.format(gamma_pos))

    # Cases in which gamma is positive in one pop and negative in the other.
    weights = PDFs.gamma(-s2.neg_gammas, [alpha, beta])

    # Case in which pop1 is positive and pop2 is negative.
    pos_neg_spectra = np.squeeze(s2.spectra[s2.gammas==gamma_pos, :len(s2.neg_gammas)])
    m4 = np.trapz(weights[:,np.newaxis,np.newaxis]*pos_neg_spectra,
                       s2.neg_gammas, axis=0)

    # Case in which pop2 is positive and pop1 is negative.
    neg_pos_spectra = np.squeeze(s2.spectra[:len(s2.neg_gammas), s2.gammas==gamma_pos])
    m7 = np.trapz(weights[:,np.newaxis,np.newaxis]*neg_pos_spectra,
                       s2.neg_gammas, axis=0)

    # Contributions to m4 and m7 for gammas that aren't covered by our cache
    # Probability toward gamma=0 that is not covered by our cache
    weight_neu, err = scipy.integrate.quad(PDFs.gamma, 0, -s2.neg_gammas[-1], args=[alpha, beta])
    # Probability toward gamma=-inf that is not covered by our cache
    weight_del, err = scipy.integrate.quad(PDFs.gamma, -s2.neg_gammas[0], np.inf, args=[alpha, beta])

    # In both cases, use the most neutral or deleterious spectra simulated
    m4 += pos_neg_spectra[0]*weight_del
    m4 += pos_neg_spectra[-1]*weight_neu
    m7 += neg_pos_spectra[0]*weight_del
    m7 += neg_pos_spectra[-1]*weight_neu

    # Weights for various parts of distribution checked against Figure S0
    fs = m5*(1-ppos_wild)*(1-pchange) +\
        m6*(1-ppos_wild)*pchange*(1-pchange_pos) +\
        m7*(1-ppos_wild)*pchange*pchange_pos +\
        m2*ppos_wild*(1-pchange) +\
        m3*ppos_wild*pchange*pchange_pos +\
        m4*ppos_wild*pchange*(1-pchange_pos)

    return theta*fs

def Vourlaki_mixture_uniq_ppos(params, ns, s1, s2, theta, pts):
    """
    Inference model from Vourlaki et al.

    params = alpha, beta, ppos_wild, gamma_pos, pchange
    """

    alpha, beta, ppos_wild, gamma_pos, pchange = params

    # We'll scale by theta at the end of the script, so set theta=1 here.
    # Case in which gamma is negative and equal in the two pops
    m5 = s1.integrate([alpha, beta], None, PDFs.gamma, 1, None)
    # Case in which gamma is negative and indepenent in the two pops
    m6 = s2.integrate([alpha, beta], None, PDFs.biv_ind_gamma, 1, None,
                      exterior_int=True)

    # Cases in which gamma is positive in both pops.
    # For simplicity, using a single gamma value, rather than exponential dist
    #  that was simulated.
    try:
        m2 = s2.spectra[s2.gammas == gamma_pos,
                        s2.gammas == gamma_pos][0]
        # Even if selection coefficient changed to a different positive value,
        #   we just model them as equal.
        m3 = m2
    except IndexError:
        raise IndexError('Failed to find requested gamma_pos={0:.4f} '
                         'in cached spectra. Was it included '
                         'in additional_gammas during '
                         'cache generation?'.format(gamma_pos))

    # Cases in which gamma is positive in one pop and negative in the other.
    weights = PDFs.gamma(-s2.neg_gammas, [alpha, beta])

    # Case in which pop1 is positive and pop2 is negative.
    pos_neg_spectra = np.squeeze(s2.spectra[s2.gammas==gamma_pos, :len(s2.neg_gammas)])
    m4 = np.trapz(weights[:,np.newaxis,np.newaxis]*pos_neg_spectra,
                       s2.neg_gammas, axis=0)

    # Case in which pop2 is positive and pop1 is negative.
    neg_pos_spectra = np.squeeze(s2.spectra[:len(s2.neg_gammas), s2.gammas==gamma_pos])
    m7 = np.trapz(weights[:,np.newaxis,np.newaxis]*neg_pos_spectra,
                       s2.neg_gammas, axis=0)

    # Contributions to m4 and m7 for gammas that aren't covered by our cache
    # Probability toward gamma=0 that is not covered by our cache
    weight_neu, err = scipy.integrate.quad(PDFs.gamma, 0, -s2.neg_gammas[-1], args=[alpha, beta])
    # Probability toward gamma=-inf that is not covered by our cache
    weight_del, err = scipy.integrate.quad(PDFs.gamma, -s2.neg_gammas[0], np.inf, args=[alpha, beta])

    # In both cases, use the most neutral or deleterious spectra simulated
    m4 += pos_neg_spectra[0]*weight_del
    m4 += pos_neg_spectra[-1]*weight_neu
    m7 += neg_pos_spectra[0]*weight_del
    m7 += neg_pos_spectra[-1]*weight_neu

    # Weights for various parts of distribution checked against Figure S0
    fs = m5*(1-ppos_wild)*(1-pchange) +\
        m6*(1-ppos_wild)*pchange*(1-ppos_wild) +\
        m7*(1-ppos_wild)*pchange*ppos_wild +\
        m2*ppos_wild*(1-pchange) +\
        m3*ppos_wild*pchange*ppos_wild +\
        m4*ppos_wild*pchange*(1-ppos_wild)

    return theta*fs


def trivial_fs(params, ns, pts): 
    return dadi.Spectrum([[0, 0.5], [0.5, 0]])

def test_Domestication_nocrash():
    """
    Simple test that the code doesn't crash.
    """
    # Test no crashes with our custom demographic model
    demo_params = [2,4,1,0.06,0.01,0.01,0.1,0.1]
    ns, pts_l = [3,3], [10,15,20]

    s1 = Cache1D(demo_params, ns, Domestication_one, pts=pts_l, gamma_pts=3, gamma_bounds=(1e-2, 1), additional_gammas=[0.1])
    s2 = Cache2D(demo_params, ns, Domestication, pts=pts_l, gamma_pts=3, gamma_bounds=(1e-2, 1), additional_gammas=[0.1])

    sel_params = 0.2, 10, 0.025, 0.1, 0.25, 0.25
    fs = Vourlaki_mixture(sel_params, None, s1, s2, 1.0, None)

def test_normalization():
    """
    Simple test that our integration is normalized for different input parameters.
    """
    ns, pts_l = [2,2], [1]
    s1 = Cache1D([], ns, trivial_fs, pts=pts_l, gamma_pts=100,
                 gamma_bounds=(1e-4, 2000), additional_gammas=[10])
    s2 = Cache2D([], ns, trivial_fs, pts=pts_l, gamma_pts=100,
                 gamma_bounds=(1e-4, 2000), additional_gammas=[10])

    # No gamma changes, no positive component: ppos_wild=0.0, pchange=0.0
    fs = Vourlaki_mixture([1, 10, 0, 10, 0, 0], None, s1, s2, 1.0, None)
    assert(np.allclose(fs.sum(), 1, atol=0.01))

    # No gamma changes, with positive component: ppos_wild=0.5, pchange=0.0
    fs = Vourlaki_mixture([1, 10, 0.5, 10, 0, 0], None, s1, s2, 1.0, None)
    assert(np.allclose(fs.sum(), 1, atol=0.01))

    # Substantial gamma changes, no positive component: ppos_wild=0.0, pchange=0.5
    fs = Vourlaki_mixture([1, 10, 0.0, 10, 0.5, 0], None, s1, s2, 1.0, None)
    assert(np.allclose(fs.sum(), 1, atol=0.01))

    # Substantial gamma changes, with positive component: ppos_wild=0.5, pchange=0.5
    fs = Vourlaki_mixture([1, 10, 0.5, 10, 0.5, 0], None, s1, s2, 1.0, None)
    assert(np.allclose(fs.sum(), 1, atol=0.01))

    # Substantial gamma changes, with positive component, some changing to positive:
    #  ppos_wild=0.5, pchange=0.5, pchange_pos=0.5
    fs = Vourlaki_mixture([1, 10, 0.5, 10, 1.0, 0.5], None, s1, s2, 1.0, None)
    assert(np.allclose(fs.sum(), 1, atol=0.01))

def visual_validation():
    """
    Simulations with a simple model and extreme parameter values, to check
    functionality and build intuition.
    """
    demo_params = [2,2,0.2,0] #nu1, nu2, T, m
    ns, pts_l = [10,10], [30,35,40]

    s1 = Cache1D(demo_params, ns, DemogSelModels.split_mig_sel_single_gamma, pts=pts_l, 
                 gamma_pts=30, gamma_bounds=(1e-4, 200), mp=True, additional_gammas=[10])
    s2 = Cache2D(demo_params, ns, DemogSelModels.split_mig_sel, pts=pts_l, 
                 gamma_pts=30, gamma_bounds=(1e-4, 200), mp=True, additional_gammas=[10])

    alpha, beta, gamma_pos = 0.2, 10, 10

    import matplotlib.pyplot as plt
    fig = plt.figure(10, figsize=(10,8))
    fig.clear()

    axw = fig.add_subplot(3,2,5)
    axd = fig.add_subplot(3,2,6)

    # alpha, beta, ppos_wild, gamma_pos, pchange, pchange_pos
    fs = Vourlaki_mixture([alpha, beta, 0, gamma_pos, 0, 0], None, s1, s2, 1.0, None)
    ax = fig.add_subplot(3,2,1); dadi.Plotting.plot_single_2d_sfs(fs, ax=ax)
    ax.set_title('ppos_wild=0.0, pchange=0.0, pchange_pos=0.0')
    ax.text(0.02,0.98,'A', fontsize='x-large', va='top', transform=ax.transAxes)
    fsw, fsd = fs.filter_pops(tokeep=[1]), fs.filter_pops(tokeep=[2])
    axw.semilogy(fsw, '-o', label='A'); axd.semilogy(fsd, '-o', label='A')

    # alpha, beta, ppos_wild, gamma_pos, pchange, pchange_pos
    fs = Vourlaki_mixture([alpha, beta, 1.0, gamma_pos, 0, 0], None, s1, s2, 1.0, None)
    ax = fig.add_subplot(3,2,2); dadi.Plotting.plot_single_2d_sfs(fs, ax=ax)
    ax.set_title('ppos_wild=1.0, pchange=0.0, pchange_pos=0.0')
    ax.text(0.02,0.98,'B', fontsize='x-large', va='top', transform=ax.transAxes)
    fsw, fsd = fs.filter_pops(tokeep=[1]), fs.filter_pops(tokeep=[2])
    axw.semilogy(fsw, '-o', label='B'); axd.semilogy(fsd, '-o', label='B')

    # alpha, beta, ppos_wild, gamma_pos, pchange, pchange_pos
    fs = Vourlaki_mixture([alpha, beta, 0.0, gamma_pos, 1.0, 0], None, s1, s2, 1.0, None)
    ax = fig.add_subplot(3,2,3); dadi.Plotting.plot_single_2d_sfs(fs, ax=ax)
    ax.set_title('ppos_wild=0.0, pchange=1.0, pchange_pos=0.0')
    ax.text(0.02,0.98,'C', fontsize='x-large', va='top', transform=ax.transAxes)
    fsw, fsd = fs.filter_pops(tokeep=[1]), fs.filter_pops(tokeep=[2])
    axw.semilogy(fsw, '-o', label='C'); axd.semilogy(fsd, '-o', label='C')

    # alpha, beta, ppos_wild, gamma_pos, pchange, pchange_pos
    # This makes me think I got m4, m7 backwards
    fs = Vourlaki_mixture([alpha, beta, 0.0, gamma_pos, 1.0, 0.9], None, s1, s2, 1.0, None)
    ax = fig.add_subplot(3,2,4); dadi.Plotting.plot_single_2d_sfs(fs, ax=ax)
    ax.set_title('ppos_wild=0.0, pchange=1.0, pchange_pos=0.9')
    ax.text(0.02,0.98,'D', fontsize='x-large', va='top', transform=ax.transAxes)
    fsw, fsd = fs.filter_pops(tokeep=[1]), fs.filter_pops(tokeep=[2])
    axw.semilogy(fsw, '-o', label='D'); axd.semilogy(fsd, '-o', label='D')

    axw.set_title('wild'); axd.set_title('domesticate')
    axw.legend(); axd.legend()

    fig.tight_layout()
    fig.savefig('validation.pdf')
    plt.show()

if __name__ == "__main__":
    test_Domestication_nocrash()
    test_normalization()
    visual_validation()



