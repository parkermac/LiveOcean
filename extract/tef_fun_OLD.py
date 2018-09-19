def tef_integrals(fn):
    # choices
    plot_profiles = True # plot profiles when certain conditions occur
    tidal_average = False # which king of time filtering
    NS_fact = 50 # range to look at for relative maxima (50 means 2% of the salt range)
    Crit_fact = 5 # was 50, Used for dropping small extrema (5 means 20% of max transport)
    nlay_max = 4 # maximum allowable number of layers to process
    Q_crit = 100 # mask transport layers with |Q| smaller than this (m3 s-1)
    
    # load results
    tef_dict = pickle.load(open(fn, 'rb'))
    tef_q = tef_dict['tef_q']
    tef_qs = tef_dict['tef_qs']
    sbins = tef_dict['sbins']
    smax = sbins.max()
    qnet = tef_dict['qnet']
    fnet = tef_dict['fnet']
    ot = tef_dict['ot']
    td = (ot - ot[0])/86400
    NS = len(sbins)
    #print('NS = %d' % (NS))

    # low-pass
    if tidal_average:
        # tidal averaging
        tef_q_lp = zfun.filt_godin_mat(tef_q)
        tef_qs_lp = zfun.filt_godin_mat(tef_qs)
        qnet_lp = zfun.filt_godin(qnet)
        fnet_lp = zfun.filt_godin(fnet)
        pad = 36
    else:
        # nday Hanning window
        nday = 5
        nfilt = nday*24
        tef_q_lp = zfun.filt_hanning_mat(tef_q, n=nfilt)
        tef_qs_lp = zfun.filt_hanning_mat(tef_qs, n=nfilt)
        qnet_lp = zfun.filt_hanning(qnet, n=nfilt)
        fnet_lp = zfun.filt_hanning(fnet, n=nfilt)
        pad = int(np.ceil(nfilt/2))

    # subsample
    tef_q_lp = tef_q_lp[pad:-(pad+1):24, :]
    tef_qs_lp = tef_qs_lp[pad:-(pad+1):24, :]
    td = td[pad:-(pad+1):24]
    qnet_lp = qnet_lp[pad:-(pad+1):24]
    fnet_lp = fnet_lp[pad:-(pad+1):24]

    # # find integrated TEF quantities
    # # alternate method using cumulative sum of the transport
    # # to identify the salinity dividing inflow and outflow
    # # RESULT: this way is not sensitive to the number of
    # # salinity bins.
    # #
    # start by making the low-passed flux arrays sorted
    # from high to low salinity
    rq = np.fliplr(tef_q_lp)
    rqs = np.fliplr(tef_qs_lp)
    sbinsr = sbins[::-1]
    # then form the cumulative sum (the function Q(s))
    qcs = np.cumsum(rq, axis=1)
    nt = len(td)

    # new version to handle more layers
    from scipy.signal import argrelextrema
    Imax = argrelextrema(qcs, np.greater, axis=1, order=int(NS/NS_fact))
    Imin = argrelextrema(qcs, np.less, axis=1, order=int(NS/NS_fact))

    Q = np.zeros((nt, nlay_max))
    QS = np.zeros((nt, nlay_max))

    crit = np.nanmax(np.abs(qcs)) / Crit_fact

    for tt in range(nt):
        # we use these masks because there are multiple values for a given day
        maxmask = Imax[0]==tt
        minmask = Imin[0]==tt
        imax = Imax[1][maxmask]
        imin = Imin[1][minmask]
        # drop extrema indices which are too close to the ends
        if len(imax) > 0:
            mask = np.abs(qcs[tt,imax] - qcs[tt,0]) > crit
            imax = imax[mask]
        if len(imax) > 0:
            mask = np.abs(qcs[tt,imax] - qcs[tt,-1]) > crit
            imax = imax[mask]
        if len(imin) > 0:
            mask = np.abs(qcs[tt,imin] - qcs[tt,0]) > crit
            imin = imin[mask]
        if len(imin) > 0:
            mask = np.abs(qcs[tt,imin] - qcs[tt,-1]) > crit
            imin = imin[mask]
        ivec = np.sort(np.concatenate((np.array([0]), imax, imin, np.array([NS]))))
        nlay = len(ivec)-1
    
        # combine non-alternating layers
        qq = np.zeros(nlay)
        qqs = np.zeros(nlay)
        jj = 0
        for ii in range(nlay):
            qlay = rq[tt, ivec[ii]:ivec[ii+1]].sum()
            qslay = rqs[tt, ivec[ii]:ivec[ii+1]].sum()
            if ii == 0:
                qq[0] = qlay
                qqs[0] = qslay
            else:
                if np.sign(qlay)==np.sign(qq[jj]):
                    qq[jj] += qlay
                    qqs[jj] += qslay
                    nlay -= 1
                else:
                    jj += 1
                    qq[jj] = qlay
                    qqs[jj] = qslay

        if nlay == 1:
            if qq[0] >= 0:
                Q[tt,1] = qq[0]
                QS[tt,1] = qqs[0]
            elif qq[0] < 0:
                Q[tt,2] = qq[0]
                QS[tt,2] = qqs[0]
        elif nlay ==2:
            if qq[0] >= 0:
                Q[tt,1] = qq[0]
                QS[tt,1] = qqs[0]
                Q[tt,2] = qq[1]
                QS[tt,2] = qqs[1]
            elif qq[0] < 0:
                Q[tt,1] = qq[1]
                QS[tt,1] = qqs[1]
                Q[tt,2] = qq[0]
                QS[tt,2] = qqs[0]
        elif nlay ==3:
            if qq[0] >= 0:
                Q[tt,1] = qq[0]
                QS[tt,1] = qqs[0]
                Q[tt,2] = qq[1]
                QS[tt,2] = qqs[1]
                Q[tt,3] = qq[2]
                QS[tt,3] = qqs[2]
            elif qq[0] < 0:
                Q[tt,0] = qq[0]
                QS[tt,0] = qqs[0]
                Q[tt,1] = qq[1]
                QS[tt,1] = qqs[1]
                Q[tt,2] = qq[2]
                QS[tt,2] = qqs[2]
        elif nlay ==4:
            if qq[0] >= 0:
                err_str = 'Backwards 4: td=%5.1f  nlay = %d' % (td[tt],nlay)
                print(err_str)
                Q[tt,:] = np.nan
                QS[tt,:] = np.nan
                # and plot a helpful figure
                if plot_profiles:
                    profile_plot(rq, sbinsr, smax, qcs, ivec, err_str, tt, NS)

            elif qq[0] < 0:
                Q[tt,0] = qq[0]
                QS[tt,0] = qqs[0]
                Q[tt,1] = qq[1]
                QS[tt,1] = qqs[1]
                Q[tt,2] = qq[2]
                QS[tt,2] = qqs[2]
                Q[tt,3] = qq[3]
                QS[tt,3] = qqs[3]
        else:
            print('- Excess layers: td=%5.1f  nlay = %d' % (td[tt],nlay))
            Q[tt,:] = np.nan
            QS[tt,:] = np.nan
        
    # form derived quantities
    Q[np.abs(Q)<Q_crit] = np.nan
    S = QS/Q
    
    return Q, S, QS, qnet_lp, fnet_lp, td
    
def profile_plot(rq, sbinsr, smax, qcs, ivec, err_str, tt, NS):
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(121)
    ax.plot(rq[tt,:], sbinsr)
    ax.set_ylim(smax,0)
    ax.grid(True)
    ax.set_ylabel('Salinity')
    ax.set_xlabel('q(s)')
    #
    ax = fig.add_subplot(122)
    ax.plot(qcs[tt,:], sbinsr)
    # print(ivec)
    Ivec = ivec.copy()
    Ivec[Ivec==NS] = NS-1
    ax.plot(qcs[tt,Ivec], sbinsr[Ivec],'*k')
    ax.set_ylim(smax,0)
    ax.grid(True)
    ax.set_xlabel('Q(s)')
    ax.set_title(err_str)