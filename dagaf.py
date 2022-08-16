def dagaf(S, nIMF=2, chi=1.6, eType1='d', eType2='d', TH1=20.0, TH2=0.001, \
    sFig=True):
     '''
     Data-adaptive Gaussian average filtering on signal S.
     Designed by Yue-Der Lin and Kai-Chun Liu, Taiwan.

     Reference: https://doi.org/10.1016/J.BSPC.2021.103104

     Inputs:
          S   The signal to be decomposed (an array).
       nIMF   Number of IMF expected to decompose (default = 2).
        chi   The parameter that influence the window length (default = 1.6).
              The value should be between 1.1 and 3.
     eType1   The extention type for the left boundary point:
              'p' (periodical) repeats the pattern outside the boundaries,
              'c' (constant) extends outside the boundaries with the last
                  values achieved at the boundaries (default),
              'r' (reflection) extends the signal symmetrical with respect
                  to the vertical lines over the boundary points.
              'd' (double-symmetric reflection) extends the signal firstly
                  symmetrical with respect to the vertical lines over the
                  boundary point and and next symmetrical with respect to
                  horizontal line over the boundary point.
     eType2   The extention type for the right boundary point:
              'p' (periodical) repeats the pattern outside the boundaries,
              'c' (constant) extends outside the boundaries with the last
                  values achieved at the boundaries (default),
              'r' (reflection) extends the signal symmetrical with respect
                  to the vertical lines over the boundary points.
              'd' (double-symmetric reflection) extends the signal firstly
                  symmetrical with respect to the vertical lines over the
                  boundary point and and next symmetrical with respect to
                  horizontal line over the boundary point.
        TH1    Threshold value for signal to residual energy ratio,
               the 1st decomposition stop criteria (default = 20), and is
               computed as 10*log10(||S(t)||^2/||S_k(t)||^2), where
                           S(t) is the original signal,
                           S_k(t) is the residual of the kth IMF.
        TH2    Threshold value for convergence check,
               the 2nd decomposition stop criteria (default = 0.001), and is
               computed as ||{S_{i-1}(t)-S_{i}(t)||^2/||S_{i}(t)||^2 at i-th
               iteration.
       sFig    True = Show the figures, False = Do not show the figures
               (default = True).

     Output:
        imf    Matrix containg all of the IMF, an array of (nIMF, length of S).
        res    S - (summation of imf), an array of (length of S, ).

     Usage example:
        import numpy as np
        import matplotlib.pyplot as plt
        s = np.random.random(1000)
        # Plot the data.
        plt.plot(s)
        from dagaf import dagaf
        imf, res = dagaf(s, 4, 1.6, 'd', 'd', 20, 0.001, True)
     '''
     # Import modules:
     import matplotlib
     import matplotlib.pyplot as plt
     import numpy as np
     from numpy import linalg as LA
     from scipy import signal
     from scipy.signal import argrelmax
     from scipy.signal import argrelmin
     from scipy import ndimage
     from scipy.ndimage import gaussian_filter1d
     import math

     # Initialization:
     L = np.size(S)
     S1 = S
     imf = np.empty((L))
     matplotlib.rcParams['font.family'] = 'Times'

     #   energyRatio: pre-defined threshold for signal to residual energy ratio.
     energyRatio = 10   # Sifting as energyRatio < TH1:
     #   rTolerancxe: pre-defined threshold for convergence check.
     rTolerance = 1     # Sifting as rTolerance > TH2:

     # DAGAF:
     #   ind: index for IMF matrix, indicating the ind-th row.
     ind = 0
     while (ind < nIMF):
        energyRatio = 10*math.log10(LA.norm(S,2)/LA.norm(S1,2))
        Max = argrelmax(S1)
        nMax = np.size(Max)
        Min = argrelmin(S1)
        nMin = np.size(Min)
        nMaxMin = nMax + nMin
        mask = int(2*np.floor(chi*L/nMaxMin))

        if (energyRatio > TH1) and (nMaxMin <= 2):
            break

        # Sifting process initialization:
        S2 = S1

        # Sifting process:
        rsigPrev = S2

        if (nMaxMin > 2) and (mask < L/2.0):
            gwLength = 2*mask+1
            H = signal.windows.gaussian(gwLength, std=mask/4.0728, sym=True)
            W = H/np.sum(H)

            # Generate new pattern according to eType1 and eType2:
            # ### For left boundary point:
            if eType1 == 'p':
                St_tmp = np.concatenate((S2[-(mask+1):-1], S2), axis=None)
            elif eType1 == 'r':
                St_tmp = np.concatenate((S2[mask:0:-1], S2), axis=None)
            elif eType1 == 'c':
                St_tmp = np.concatenate((S2[0]*np.ones((mask)), S2), axis=None)
            else: # eType == 'd'
                xt1 = 2*S2.mean() - S2[mask:0:-1]
                St_tmp = np.concatenate((xt1, S2), axis=None)

            ### For right boundary point:
            if eType2 == 'p':
                St = np.concatenate((St_tmp, S2[1:mask+1]), axis=None)
            elif eType2 == 'r':
                St = np.concatenate((St_tmp, S2[-2:-(mask+2):-1]), axis=None)
            elif eType2 == 'c':
                St = np.concatenate((St_tmp, S2[-1]*np.ones((mask))), axis=None)
            else: # eType2 == 'd'
                xt2 = 2*S2.mean() - S2[-2:-(mask+2):-1]
                St = np.concatenate((St_tmp, xt2), axis=None)

            # Filtering:
            ave = np.empty((L))
            for i in range(L):
                ave[i] = np.dot(W, St[i:i+gwLength])

        else:
            break

        if sFig == True:
            # Plot the iterative procedure for each imf:
            l = np.arange(L)
            fig1 = plt.figure(1)
            # plt.suptitle('Signal and Average')
            plt.subplot(nIMF, 1, ind+1)
            plt.plot(l, S2, l, ave)
            plt.ylabel('Step {}'.format(ind+1))
            fig1.subplots_adjust(hspace = 0.5)

            ### For last sub-figure:
            if ind == nIMF-1:
                plt.xlabel('Samples')
            else :
                pass

            ### Save the figure:
            plt.savefig('Fig_Signal_Ave.jpg', dpi=1200, transparent=True)
        else:
            pass


        # The resulted S2 is one imf.
        S2 = S2 - ave

        # Residual tolerance:
        rTolerance = (LA.norm(rsigPrev-S2,2)/LA.norm(S1,2))**2

        # The second convergence check:
        if (rTolerance < TH2):
            break

        if ind == 0:
            imf = S2
        else:
            imf = np.vstack((imf,S2))

        S1 = S1 - S2

        ind = ind + 1

     if imf.ndim > 1:
        (m,_) = imf.shape
        res = S

        for i in range(m):
            res = res - imf[i,:]

            if sFig == True:
               # Plot each imf:
               fig2 = plt.figure(2)
               # plt.suptitle('IMFs')
               plt.subplot(m,1,i+1)
               plt.plot(imf[i,:])
               plt.ylabel('IMF {}'.format(i+1))
               fig2.subplots_adjust(hspace = 0.5)

               ### For last sub-figure:
               if i == m-1:
                  plt.xlabel('Samples')
               else:
                  pass

               ### Save the figure:
               plt.savefig('Fig_imf.jpg', dpi=1200, transparent=True)
            else:
               pass

        if sFig == True:
           # Plot signal S and the residual:
           l = np.arange(L)
           plt.figure(3)
           plt.plot(l, S, label='Signal')
           plt.plot(l, res, label='Residue')
           plt.xlabel('Samples')
           plt.legend(loc=3)
           plt.savefig('Fig_Signal_Res.jpg', dpi=1200, transparent=True)
        else:
           pass

        return imf, res

     else:
        imf = S
        res = np.zeros(L)
        if sFig == True:
           plt.figure()
           plt.plot(imf, label='Signal')
           plt.plot(res, label='Residue')
           plt.xlabel('Samples')
           plt.legend(loc=3)
           plt.savefig('Fig_Signal_Res.jpg', dpi=1200, transparent=True)
           print('Signal is a simple mode, decomposition is not necessary.')
        else:
           pass

        return imf, res
    