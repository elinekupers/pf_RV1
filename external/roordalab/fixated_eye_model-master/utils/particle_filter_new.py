# Code for a particle filter
#
# A particle filter lets us sample the distribution of hidden states
#    given the observed state in a Markov chain where the hidden
#    state is continuous valued.
# The file defines class prototypes to be implemented to run the filter

import abc
import numpy as np
import matplotlib.pyplot as plt


def gauss_pdf(X, s):
    """
    Gaussian PDF with diagonal covariance
    s = vector of size (D_H,) giving standard deviations for each dimension
    X - input of dimension (N_S, D_H)
    Note this makes use of numpy broadcasting by comparing last axes first
       so it still works when N_S = D_H
    Returns probabilities for each sample
    """
    return np.prod(np.exp(- 0.5 * X ** 2 / s ** 2)
                   / np.sqrt(2 * np.pi * s ** 2), axis=1)


class InitPropDist():
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def __init__(self, D_H, D_O, *args):
        """
        Create an object for the proposal distribution for the first
            hidden state.
        D_H - dimension of hidden state
        D_O - dimension of observed state
        *args - other arguments to define the distribution
        """
        self.D_H = D_H
        self.D_O = D_O

    @abc.abstractmethod
    def prob(self, X0, Y0):
        """
        Return proposal probability q(X0|Y0)
        X0 - Hidden Initial State - (N_samples, D_H)
        Y0 - Initial Observation - (D_O,)
        """
        if not self.D_H == X0.shape[1]:
            raise ValueError('Hidden State Dimension mismatch')
        if not self.D_O == Y0.shape[0]:
            raise ValueError('Observed state dimension mismatch')

    @abc.abstractmethod
    def sample(self, Y0, N_S):
        """
        Generate N_S samples from the proposal distribution q(X0|Y0)
        Y0 - Initial observation - (D_O,)
        N_S - number of samples to return
        Returns samples as rows
        """
        if not self.D_O == Y0.shape[0]:
            raise ValueError('Observed state dimension mismatch')
        # return np.zeros((N_S, self.dim))


class GaussIPD(InitPropDist):

    def __init__(self, D_H, D_O, sigmas):
        """
        Create a simple Gaussian proposal distribution for the
            initial state:
        q(X_0|Y_0) = N(0, sigmas)
        sigmas - standard deviation for each variable in form (D_H,)
        D_H - dimension of hidden state
        D_O - dimension of observed state
        """
        InitPropDist.__init__(self, D_H, D_O)
        self.sigmas = sigmas

    def prob(self, X0, Y0):
        """
        Return the proposal distribution probability q(X0|Y0)
        X0 - Initial hidden state
        Y0 - Initial observed state
        """
        InitPropDist.prob(self, X0, Y0)
        return gauss_pdf(X0, self.sigmas)

    def sample(self, Y0, N_S):
        """
        Sample from proposal distribution given one sample Y0
        Y0 - Initial observation - (D_O,)
        N_S - number of samples to return
        Returns samples as rows
        """
        InitPropDist.sample(self, Y0, N_S)
        return (self.sigmas * np.random.randn(N_S, self.D_H))


class TransPropDist():
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def __init__(self, D_H, D_O, *args):
        """
        Object for the proposal distribution for transitions between
            hidden states:
        q(Xc|Yc, Xp),
        Xc - current hidden state
        Yc - current observed state
        Xp - previous hidden state
        """
        self.D_H = D_H
        self.D_O = D_O

    @abc.abstractmethod
    def prob(self, Xc, Yc, Xp):
        """
        Return value of q(Xc|Yc, Xp)
        Xc - current state - (N_S, N_H)
        Xp - previous state - (N_S, N_P)
        Yc - current observed state - (N_O,)
        """
        if not Xc.shape[0] == Xp.shape[0]:
            raise ValueError('Xc and Xp must have the same number of samples')
        if not self.D_H == Xc.shape[1]:
            raise ValueError('Xc[0] must have dimension self.D_H')
        if not self.D_H == Xp.shape[1]:
            raise ValueError('Xp[0] must have dimension self.D_H')
        if not self.D_O == Yc.shape[0]:
            raise ValueError('Yp must have shape (self.D_O,)')

    @abc.abstractmethod
    def sample(self, Yc, Xp):
        """
        Return samples Xc ~ q(Xc|Yc, Xp)
        Note Xp has rows that are samples and columns that are the
            dimensions of the hidden state
        """
        if not Xp.shape[1] == self.D_H:
            raise ValueError('Xc[0] must have shape self.D_H')
        if not self.D_O == Yc.shape[0]:
            raise ValueError('Yc must have shape (N_O,)')


class GaussTPD(TransPropDist):

    def __init__(self, D_H, D_O, sigmas, A=None):
        """
        Simple gaussian transition proposal distribution:
            q(Xc|Yc, Xp) ~ N(mu = A * Xp, sigmas)
        sigmas - standard deviations for the proposal
        D_H - dimension of hidden state of the markov model
        A - transition matrix, assume X_t's are column vectors and
             Xtp1 = A xt
        """
        TransPropDist.__init__(self, D_H, D_O)
        if A is None:
            A = np.eye(self.D_H)
        if not (self.D_H == A.shape[0] and self.D_H == A.shape[1]):
            raise ValueError('Wrong Shape for Transition matrix, A')
        self.sigmas = sigmas
        self.A = A

    def prob(self, Xc, Yc, Xp):
        """
        Returns the probability N(mu = Xp, sigma)(Xc)
        """
        return gauss_pdf(Xc - np.dot(Xp, np.transpose(self.A)), self.sigmas)

    def sample(self, Yc, Xp):
        """
        Samples from the proposal distribution
        Yc - Current observed state
        Xp - previous hidden state, each particle comes in as a row
        """
        return (np.random.randn(*Xp.shape) * self.sigmas
                + np.dot(Xp, np.transpose(self.A)))


# In addition to proposal distributions, we have potentials that give us
#   the probabilities of the HMM transition and the emission probabilities

class InitialPotential:
    __metaclass__ = abc.ABCMeta
    """
    Abstract base class for the prior on the initial hidden state of the HMM,
    p(X0)
    """
    @abc.abstractmethod
    def __init__(self, D_H):
        self.D_H = D_H

    @abc.abstractmethod
    def prob(self, X0):
        """
        X0 - samples of current hidden state in form (N_S, D_H)
        """
        if not X0.shape[1] == self.D_H:
            raise ValueError('Dimension of input does not' +
                             'match hidden state dimension')


class GaussIP(InitialPotential):

    def __init__(self, D_H, sigmas):
        """
        Simple Prior for the initial hidden state
            p(X0) = N(0, sigma)
        sigma - standard deviation
        D_H - dimension of hidden state
        """
        InitialPotential.__init__(self, D_H)
        self.sigmas = sigmas

    def prob(self, X0):
        InitialPotential.prob(self, X0)
        return gauss_pdf(X0, self.sigmas)


class TransPotential():
    __meta__ = abc.ABCMeta
    """
    Class prototype for the transition probability potential:
    p(Xc|Xp), Xc - current state, Xp - previous state
    D_H - dimension of hidden input
    """
    @abc.abstractmethod
    def __init__(self, D_H):
        self.D_H = D_H

    @abc.abstractmethod
    def prob(self, Xc, Xp):
        """
        Return probability p(Xc|Xp)
        """
        if not Xc.shape[0] == Xp.shape[0]:
            raise TypeError('Number of samples mismatch in input')
        if not Xc.shape[1] == self.D_H:
            raise TypeError('Xc[0] has incorrect dimension')
        if not Xp.shape[1] == self.D_H:
            raise TypeError('Xp[0] has incorrect dimension')


class GaussTP(TransPotential):

    """
    Gaussian transition probability potential:
    p(Xc|Xp) = N(mu = Xp, sigma)
    """

    def __init__(self, D_H, sigmas, A=None):
        """
        sigma - standard deviation
        D_H - dimension of hidden state
        A = transition matrix, so that
            Xc = A Xp where Xc,Xp are column vectors
        """
        TransPotential.__init__(self, D_H)
        self.sigmas = sigmas
        if A is None:
            self.A = np.eye(self.D_H)
        else:
            self.A = A

    def prob(self, Xc, Xp):
        """
        Return p(Xc|Xp)
        """
        TransPotential.prob(self, Xc, Xp)
        return gauss_pdf(Xc - np.dot(Xp, np.transpose(self.A)),
                         self.sigmas)


class LikelihoodPotential():
    __meta__ = abc.ABCMeta
    """
    Class prototype for likelihood function p(Yc|Xc)
    Yc - current observed state
    Xc - current hidden state
    """
    @abc.abstractmethod
    def __init__(self, D_H, D_O):
        """
        D_H - dimension of the hidden state
        D_O - dimension of the observed state
        """
        self.D_H = D_H
        self.D_O = D_O

    @abc.abstractmethod
    def prob(self, Yc, Xc):
        """
        Return p(Yc|Xc)
        Xc - hidden state of dimension (N_S, D_H)
        Yc - observed state of dimension (D_O,)
        """
        if not self.D_H == Xc.shape[1]:
            raise ValueError('Dimension of Hidden State Mismatch')
        if not self.D_O == Yc.shape[0]:
            raise ValueError('Dimension of Observed State Mismatch')


class GaussLP(LikelihoodPotential):

    """
    Simple gaussian likelihood p(Yc|Xc) = N(mu = Xc, _sigma)
    """

    def __init__(self, D_H, D_O, sigmas):
        LikelihoodPotential.__init__(self, D_H, D_O)
        if not self.D_H == self.D_O:
            raise ValueError('Input and output dimension must match')
        self.sigmas = sigmas

    def prob(self, Yc, Xc):
        LikelihoodPotential.prob(self, Yc, Xc)
        return gauss_pdf(Yc - Xc, self.sigmas)


class ParticleFilter:

    """
    Actual class for the particle filter
    Stores the potentials and proposal distributions to run the HMM
    N_T - number of time steps
    N_P - number of particles
    dim - number of dimensions of hidden state
    After calling run with some output data, the class stores the
        weights and the associated particles in the arrays
    self.XS - weights in form (N_T, N_P, dim)
    self.WS - weights in the form (N_T, N_P)
    Note XS[i], WS[i] correspond to samples and weights from the
        following distribution:
        p(X_i|Y_0,Y_1,...,Y_i)
    """

    def __init__(self, ipd, tpd, ip, tp, lp, Y, N_P):
        """
        Create a Particle Filter Object
        ipd - Initial Proposal Distribution
        tpd - Transition Proposal Distribution
        ip - Initial Potential
        tp - Transition Potential
        lp - Likelihood Potential
        Y - matrix of observed values of the form (N_T, other dims)
        N_P - number of particles to use in the particle filter
        """
        self.ipd = ipd
        self.tpd = tpd
        self.ip = ip
        self.tp = tp
        self.lp = lp
        self.D_H = self.ipd.D_H
        self.D_O = self.ipd.D_O

        self.check_dims()

        self.Y = Y
        self.N_T = self.Y.shape[0]
        self.N_P = N_P

        self.XS = np.zeros((self.N_T, self.N_P, self.D_H)).astype('float32')
        self.WS = np.zeros((self.N_T, self.N_P)).astype('float32')

        self.t = 0  # Current time of the particle filter
        self.resample_times = []

        self.means = np.zeros((self.N_T, self.D_H))
        self.sdevs = np.zeros((self.N_T, self.D_H))

    def check_dims(self):
        """
        Checks that the dimensions of the proposal distributions
            and the probabilities are consistent
        """
        if not self.D_H == self.ipd.D_H:
            raise ValueError('IPD D_H mismatch')
        if not self.D_H == self.tpd.D_H:
            raise ValueError('TPD D_H mismatch')
        if not self.D_H == self.ip.D_H:
            raise ValueError('IP D_H mismatch')
        if not self.D_H == self.tp.D_H:
            raise ValueError('TP D_H mismatch')
        if not self.D_H == self.lp.D_H:
            raise ValueError('LP D_H mismatch')

        if not self.D_O == self.ipd.D_O:
            raise ValueError('IPD D_O mismatch')
        if not self.D_O == self.tpd.D_O:
            raise ValueError('TPD D_O mismatch')
        if not self.D_O == self.lp.D_O:
            raise ValueError('LP D_O mismatch')

    def reset(self):
        """
        Resets the particle filter
        """
        self.XS[:, :] = 0.
        self.WS[:, :] = 0.
        self.t = 0
        self.resample_times = []
        self.means[:, :] = 0.
        self.sdevs[:, :] = 0.

    def resample_idx(self, W):
        """
        Implements the systematic resampling method
        Given a normalized weight matrix W for one time step
            (shape (N_P,) ), return
            indices corresponding to a resampling
        Eg. if we have our samples X, then
        X[idx] gives us a resampling of those samples
        """
        if np.abs(np.sum(W) - 1) > 0.01:
            raise ValueError('Weights for resampling not normalized')
        if not W.shape == (self.N_P,):
            raise ValueError('Bad dimension for weights to be resampled')

        u = np.random.random() / self.N_P
        U = u + np.arange(self.N_P) / (1. * self.N_P)
        S = np.cumsum(W)

        i = 0
        j = 0
        idx = np.zeros(self.N_P, dtype='int')

        while(j < self.N_P and i < self.N_P):
            if (U[j] <= S[i]):
                idx[j] = i
                j = j + 1
            else:
                i = i + 1
        return idx

    def advance(self):
        """

        Advances the particle filter one timestep

        Performs a resampling if the effective sample size becomes
            less than have the number of particles
        """

        if (self.t == 0):
            self.XS[0] = self.ipd.sample(self.Y[0], self.N_P)
            self.WS[0] = (self.ip.prob(self.XS[0])
                          * self.lp.prob(self.Y[0], self.XS[0])
                          / self.ipd.prob(self.XS[0], self.Y[0])
                          )
            self.WS[0] = self.WS[0] / np.sum(self.WS[0])
            self.t = self.t + 1

        elif (self.t < self.N_T):
            i = self.t
            self.XS[i] = self.tpd.sample(self.Y[i], self.XS[i - 1])
            self.WS[i] = (self.WS[i - 1]
                          * self.lp.prob(self.Y[i], self.XS[i])
                          * self.tp.prob(self.XS[i], self.XS[i - 1])
                          / self.tpd.prob(self.XS[i], self.Y[i], self.XS[i-1])
                          )
            self.WS[i] = self.WS[i] / np.sum(self.WS[i])
            if (np.sum(self.WS[i] ** 2) ** (-1) < self.N_P / 2.):
                idx = self.resample_idx(self.WS[i])
                self.XS[i] = self.XS[i, idx]
                self.WS[i] = 1. / self.N_P
                self.resample_times.append(i)

            self.t = self.t + 1
        else:
            raise StopIteration("Particle filter has already run to" +
                                "completion + \n" +
                                "use the reset function to do another run")

    def calculate_means_sdevs(self):
        """
        Calculates the means and standard deviations
        """
        self.means = np.sum(self.XS *
                            self.WS.reshape(self.N_T, self.N_P, 1),
                            axis=1)
        self.sdevs = np.sqrt(np.sum(self.XS ** 2 *
                                    self.WS.reshape(self.N_T, self.N_P, 1),
                                    axis=1) - self.means ** 2)

    def plot(self, X, DT, show=True):
        """
        Generate a plot of the estimated hidden state and the real hidden state
            as a function of time
        X - actual hidden state in the form (N_T, num hidden states)
        DT - timestep size
        """

        self.calculate_means_sdevs()

        plt.subplot(1, 2, 1)
        plt.fill_between(DT * np.arange(self.N_T),
                         self.means[:, 0] - self.sdevs[:, 0],
                         self.means[:, 0] + self.sdevs[:, 0],
                         alpha=0.5, linewidth=1.)
        plt.plot(DT * np.arange(self.N_T), self.means[:, 0], label='estimate')
        plt.plot(DT * np.arange(self.N_T), X[:, 0], label='actual')
        plt.xlabel('Time (s)')
        plt.ylabel('Relative position (pixels)')
        plt.legend()
        plt.title('Estimated versus actual position given the image for dim0')

        plt.subplot(1, 2, 2)
        plt.fill_between(DT * np.arange(self.N_T),
                         self.means[:, 1] - self.sdevs[:, 1],
                         self.means[:, 1] + self.sdevs[:, 1],
                         alpha=0.5, linewidth=1.)
        plt.plot(DT * np.arange(self.N_T), self.means[:, 1], label='estimate')
        plt.plot(DT * np.arange(self.N_T), X[:, 1], label='actual')
        plt.xlabel('Time (s)')
        plt.ylabel('Relative position (pixels)')
        plt.legend()
        plt.title('Estimated versus actual position given the image for dim1')

        if show:
            plt.show()


def main():
    # A quick example: consider a HMM with gaussian transitions, and
    #   gaussian noise added onto the outputs
    s1 = 0.01  # Hidden state transition standard deviation
    s2 = 0.03  # Noise for observed state
    N_T = 100  # Number of time steps
    N_P = 50  # Number of particles
    D_H = 2  # Dimension of hidden state
    D_O = D_H  # Dimension of output state

    # Generate some data according to this model
    X = np.zeros((N_T, D_H))
    Y = np.zeros((N_T, D_H))
    for i in range(D_H):
        X[:, i] = np.cumsum(s1 * np.random.randn(N_T))
    Y = np.random.randn(N_T, D_H) * s2 + X

    # Create the appropriate proposal distributions and potentials
    ipd = GaussIPD(D_H, D_O, s1)
    tpd = GaussTPD(D_H, D_O, s1)
    ip = GaussIP(D_H, s1)
    tp = GaussTP(D_H, s1)
    lp = GaussLP(D_H, D_O, s2)

    pf = ParticleFilter(ipd, tpd, ip, tp, lp, Y, N_P)
    for _ in range(N_T):
        pf.advance()

#    pf.plot(X, 0.001)


if(__name__ == '__main__'):
    main()
