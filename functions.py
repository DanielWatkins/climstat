def robust_rank_order_test(x,y):
    """From Lazante 1996: The robust rank-order test is similar in nature 
    to the WIlcoxon-Mann-Whitney test in that it is also a test for 
    equality of medians based on two samples; by contrast the rank-order 
    test does not assume that the parent populations have equal variance. 
    It has essentially the same power as the Wilcoxon-Mann-Whittney test. 
    It is resistant [to outliers] because it is based on ranks.
    
    x and y should both be one-dimensional Numpy arrays. They do not need 
    to have the same length.
    
    If nx > 12 or ny > 12, a two-tailed test using a normal probability 
    distribution can be used, otherwise see table K of appendix 1 in 
    Siegel and Castellan (1988).
    
    Returns the p-value of the test statistic 
    z = 0.5*(nx*NX_bar - ny*NY_bar)/np.sqrt(NX_bar*NY_bar + s2_NX + s2_NY)
    """
    
    import numpy as np
    import scipy.stats as stats
    nx = len(x)
    ny = len(y)

    if nx + ny < 13:
        print('Warning: sample size too small for two-tailed test with' +
              'normal distribution')

    # Compute ranks of the combined sample
    xy = np.concatenate((x,y))
    xy_ranks = stats.rankdata(xy)
    x_ranks = xy_ranks[0:nx]
    y_ranks = xy_ranks[nx:]

    # For each observation, calculate the number of lower-ranked
    # observations in the other sample.
    NX = np.zeros_like(x)
    NY = np.zeros_like(y)

    for i in range(nx):
        NX[i] = np.sum(y_ranks < x_ranks[i])

    for i in range(ny):
        NY[i] = np.sum(x_ranks < y_ranks[i])

    NX_bar = np.mean(NX)
    NY_bar = np.mean(NY)

    s2_NX = np.sum((NX-NX_bar)**2)
    s2_NY = np.sum((NY-NY_bar)**2)

    z = 0.5*(nx*NX_bar - ny*NY_bar)/np.sqrt(NX_bar*NY_bar + s2_NX + s2_NY)

    
    p = stats.norm.cdf(-np.abs(z),loc=0,scale=1)
    
    return p

def single_change_point_test(x):
    """
    From Lazante 1996: The change-point test presented here is used to determine if, and 
    locate a point in the time series at which, the median changes. The 
    test is based on summing the ranks of the data from the beginning to 
    each point in the series; each raw sum is adjusted by the amount 
    expected on average, which is linearly proportional to the point in 
    the time series. The maximum of the adjusted for significance. 
    
    The single-point change test is from Siegel and Castellan 1988.
    
    x should be a time series.
    Returns the p value and the location of the change point.
    
    If the sample size is too low, then the returned p value is 1.
    """

    n = len(x)
    R = stats.rankdata(x)
    SR = np.cumsum(R)
    SA = np.zeros_like(SR)
    for i in range(n):
        SA[i] = np.abs((2*SR[i] - i*(n+1)))
        
    n1 = np.argmax(SA)
    W = SR[n1]
    n2 = n - n1
    
    # If too short, return with very high p value
    if (n1 <= 10) | (n2 <= 10):
        return n1,1
    
    Wcrit = n1*(n+1)/2
    sW = np.sqrt(n1*n2*(n+1)/12)
    if W < Wcrit:
        delta = 0.5
    elif W > Wcrit:
        delta = -0.5
    else:
        delta = 0
        
    z = (W - Wcrit + delta)/sW
    
    if (n1 <= 10) | (n2 <= 10):
        print('Warning: sample size too small!')
        
    p = stats.norm.cdf(-np.abs(z),loc=0,scale=1)
    
    return n1+1,p



def change_point_snr(x, cp):
    """This measure is the the ratio of the variance associated with the
    discontinuity of the resistant mean ('signal') to the resistant 
    variance that remains after the discontinuity has been removed ('noise')
    In order to have both high resistance and efficiency, biweight estimators
    are used here.
    
    x = segment
    cp = potential changepoint (e.g. segment 1 = x[0:cp], segment 2 = x[cp:])
    
    Some considerations: A trend noise variance could be computed as well,
    using resistant regression. Minimum RDN is subjective, Lazante mentions 
    that values of 0.05 to 0.1 are good but sometimes even 0.3 to 0.5 may work
    better. Furthermore, the utility of RDN can be degraded as the length of 
    segments gets larger, so using the 50 points to each side of the change
    point may be a good idea.
    """
    
    XL = biweight_mean(x[0:cp])
    XR = biweight_mean(x[cp:])
    
    nl = len(x[0:cp])
    nr = len(x[cp:])
    n = len(x)
    
    X = (nl*XL + nr*XR)/n
    
    # Variance associated with the discontinuity
    sD2 = (nl*(XL-X)**2 + nr*(XR-X)**2)/(n-1)
    
    # Compute noise variance
    segl = x[0:cp] - XL
    segr = x[cp:] - XR
    sN2 = biweight_std(np.concatenate((segl,segr)))**2
    
    return sD2/sN2    
    
    
def biweight_mean(x,c=7.5):
    """Described in Hoaglin et al (1983). The biweight estimate is a weighted 
    average such that weighting decreases away from the centre of the 
    distribution. The parameter c determines the critical distance from the 
    center; values further than which are weighted as zero.
    
    Returns the biweight mean. """

    # compute median
    M = np.nanmedian(x)
    
    # compute median absolute deviation
    MAD = np.nanmedian(np.abs(x-M))
    
    # compute weights
    u = (x - M)/(c*MAD)
    
    # perform censoring
    u[np.abs(u)>=1.0] = 0
    
    # biweight mean
    X = M + np.nansum((x-M)*(1-u**2)**2)/np.nansum((1-u**2)**2)
    
    return X

def biweight_std(x,c=7.5):
    """Described in Hoaglin et al (1983). The biweight estimate is a weighted 
    average such that weighting decreases away from the centre of the 
    distribution. The parameter c determines the critical distance from the 
    center; values further than which are weighted as zero.
    
    Returns the biweight standard deviation. """

    # compute median
    M = np.nanmedian(x)
    
    # compute median absolute deviation
    MAD = np.nanmedian(np.abs(x-M))
    
    # compute weights
    u = (x - M)/(c*MAD)
    
    # perform censoring
    u[np.abs(u)>=1.0] = 0
    
    # biweight standard deviation
    S = np.sqrt(len(x)*np.sum((x-M)**2*(1-u**2)**4))/np.abs(np.sum((1-u**2)*(1-5*u**2)))
    
    return S
   
def pseudo_std(x, method='symmetric'):
    """Convenience function to compute pseudo-standard deviation following 
    a few methods. If method=asymmetric, compute a separate std for the
    upper and lower half of the distribution."""
    from scipy.stats import iqr

    if method=='symmetric':
        return iqr(x)/1.349
    
    if method=='asymmetric':
        lower = iqr(x,rng=(25,50))*2
        upper = iqr(x,rng=(50,75))*2
        return (lower,upper)

def multiple_change_point_test(x,alpha=0.05, rdn=0.075, min_dist=10):
    """Iterative procedure to search for multiple change-points, based on 
    application of the single change point test followed by adjustment 
    of the series. Adjustment continues as long as the significance of each 
    new change-point is less than an a priori specified level.
    x = time series (as numpy array)
    alpha = Fraction such that 1-alpha is the significance level (0 to 1).
    rdn = minimum RDN to consider a change point valid.
    min_dist = minimum distance a change point must be from the endpoints 
    of a segment.
    
    
    
    """
    
    def update_lists(x,new_cp,cp_list,check_seg,min_dist):
        """Make new list with all unique cps found so far,
        then update the check_seg list, requiring that segments 
        at least be 2xmin_dist in length in order to be searched for 
        change points."""

        new_cp_list = np.unique(np.concatenate((new_cp,cp_list)))

        # Indices of cps from previous iteration
        idx_old = np.where([cp in cp_list for cp in new_cp_list])[0][:-1]

        new_check_seg = np.ones(len(new_cp_list)-1)*True
        new_check_seg[idx_old] = check_seg

        new_seg_list = make_segment_list(x,new_cp_list)
        for i in range(len(new_seg_list)):
            if len(new_seg_list[i]) <= 2*min_dist:
                new_check_seg[i] = False

        return new_check_seg, new_seg_list, new_cp_list

    def make_segment_list(x,cp_list):
        """Make a list of segments and normalize by subtracting the median.
        x is the original sequence, and cp_list is the sorted list of change
        points including the endpoints."""
        import numpy as np
        n = len(cp_list)
        seg_list = []
        for i in range(0,n-1):
            seg_list.append(x[cp_list[i]:cp_list[i+1]])
            seg_list[i] = seg_list[i] - np.median(seg_list[i])

        return seg_list

    
    n = len(x)
    
    cp_list = np.array([0,n])
    check_seg = np.array([True])
    seg_list = make_segment_list(x,cp_list)
    num_seg = len(seg_list)
    count = 0
    while any(check_seg):
        count += 1
        new_cp = []
        for i in range(len(seg_list)):
            seg = seg_list[i]
            check = check_seg[i]
            left_idx = cp_list[i]
            if check:
                # Get candidate cp
                cp,p = single_change_point_test(seg)
                # If candidate is significant and far enough from both edges
                if np.all([cp > min_dist, 
                           len(seg)-cp > min_dist, 
                           p <= alpha, 
                           change_point_snr(seg, cp) >= rdn]):
                    
                    # Check whether difference is significant
                    
                    p = robust_rank_order_test(seg[0:cp],seg[cp:])
          
                    if p < alpha:
                        new_cp.append(cp+cp_list[i])
           
        # After going through the full list of segments,
        # we've only kept change points that divide the original
        # data into sections with significantly different medians.
        # We've already set check_seg to false for change points that
        # are too close to edges, don't have high enought snr, and 
        # are not significant. Now we update the cp_list and check_seg.
        
        if len(new_cp) > 0:
            new_cp = np.array(new_cp)
            check_seg, seg_list, cp_list = update_lists(x,new_cp,cp_list,check_seg,min_dist)
        else:
            check_seg = [False]
       
    
    return cp_list  
