import numpy as np
import cPickle as pkl
import sys
import os
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt


def inner_product(p1, l1x, l1y,
                  p2, l2x, l2y, var):
    """
    Calculate the inner product between two images with the representation:
    p1 -> array of pixels
    l1x -> array of pixel x coordinate
    l2x -> array of pixel y coordinate
    each pixel is surrounded by a gaussian with variance var
    """
    N = l1x.shape[0]
    l1x = l1x.reshape(N, 1)
    l2x = l2x.reshape(1, N)

    l1y = l1y.reshape(N, 1)
    l2y = l2y.reshape(1, N)

    coupling = np.exp(-((l1x - l2x) ** 2 +
                        (l1y - l2y) ** 2) / (4 * var))

    return np.einsum('i,j,ij->', p1, p2, coupling)


def SNR(p1, l1x, l1y, p2, l2x, l2y, var):
    """
    Using the inner product defined above, calculates the SNR
        between two images given in sum of gaussian form
    See inner product for definitinos of variables
    Note the first set of pixels and pixel locations is
        considered to be the ground truth
    """
    ip12 = inner_product(p1, l1x, l1y, p2, l2x, l2y, var)
    ip11 = inner_product(p1, l1x, l1y, p1, l1x, l1y, var)
    ip22 = inner_product(p2, l2x, l2y, p2, l2x, l2y, var)

    return ip11 / (ip11 + ip22 - 2 * ip12)


class DataAnalyzer:

    def __init__(self, data):
        """
        data - dictionary containing information about run
        Loads in data from file and saves parameters in class
        """

        self.data = data

        self.DT = self.data['DT']
        self.N_T = self.data['N_T']

        self.XR = self.data['XR'][0]
        self.YR = self.data['YR'][0]
        self.S_gen = self.data['S_gen']
        self.Var = self.data['Var'][0]

        self.blur_sdev = float(np.sqrt(self.Var))

        self.N_itr = self.data['N_itr']
        self.DC_gen = self.data['DC_gen']
        self.DC_infer = self.data['DC_infer']
        self.L_I = self.data['L_I']
        self.LAMBDA = self.data['LAMBDA']
        self.R = self.data['R']
        self.L_N = self.data['L_N']

        # Convert retinal positions to grid
        XS = self.data['XS']
        YS = self.data['YS']

        XS, YS = np.meshgrid(XS, YS)
        self.XS = XS.ravel()
        self.YS = YS.ravel()

        self.N_itr = self.data['N_itr']
        self.var = self.data['Var'][0]

    @classmethod
    def fromfilename(cls, filename):
        """
        filename - filename containing data file from EM run
        Initialize the class using a filename to a pkl file
            that contains the datadict
        """
        data = pkl.load(open(filename, 'rb'))
        return cls(data)

    def SNR_idx(self, q):
        """
        Calculates the SNR between the ground truth and the
            image produced after iteration q of the EM
        Note that we shift the image estimate by the average
            amount that the path estimate was off the true path
            (There is a degeneracy in the representation that this
            fixes. )
        """
        XYR_est = self.data['EM_data'][q]['path_means']
        S_est = self.data['EM_data'][q]['image_est']
        t = self.data['EM_data'][q]['time_steps']

        XR_est = XYR_est[:, 0]
        YR_est = XYR_est[:, 1]

        dx = np.mean(self.XR[0:t] - XR_est[0:t])
        dy = np.mean(self.YR[0:t] - YR_est[0:t])

        self.dx = dx
        self.dy = dy

        return SNR(self.S_gen.ravel(), self.XS, self.YS,
                   S_est.ravel(), self.XS + dx, self.YS + dy,
                   self.var)

    def SNR_list(self):
        """
        Returns a list giving the SNR after each iteration
        """
        return [self.SNR_idx(q) for q in range(self.N_itr)]

    def time_list(self):
        """
        Returns a list of the times for each EM iteration in ms
        """
        return self.N_T * (np.arange(self.N_itr) + 1) / self.N_itr * 1000 * self.DT

    def plot_path_estimate(self, q, d):
        """
        Plot the actual and estimated path generated from EM on iteration q
        d - dimension to plot (either 0 or 1)
        """
        est_mean = self.data['EM_data'][q]['path_means']
        est_sdev = self.data['EM_data'][q]['path_sdevs']

        if (d == 0):
            path = self.XR
            label = 'Hor.'
            dxy = self.dx
        elif (d == 1):
            path = self.YR
            label = 'Ver.'
            dxy = self.dy
        else:
            raise ValueError('d must be either 0 or 1')

        t = self.data['EM_data'][q]['time_steps']

        plt.fill_between(self.DT * np.arange(self.N_T),
                         est_mean[:, d] - est_sdev[:, d],
                         est_mean[:, d] + est_sdev[:, d],
                         alpha=0.5, linewidth=1.)
        plt.plot(self.DT * np.arange(self.N_T),
                 est_mean[:, d], label='estimate')
        plt.plot(self.DT * np.arange(self.N_T),
                 path, label='actual')
        plt.xlabel('Time (s)')
        plt.ylabel('Relative position (pixels)')
        # plt.legend()
        plt.title(label + ' Pos., shift = %.2f' % dxy)

    def plot_velocity_estimate(self, q, d):
        """
        Plot the estimate of the velocity by the EM algorithm
        q - EM iteration number
        d - dimension to plot (0 or 1)
        """
        if not self.data['motion_prior'] == 'VelocityDiffusion':
            raise RuntimeError('No velocity for this motion prior')

        est_mean = self.data['EM_data'][q]['path_means']
        est_sdev = self.data['EM_data'][q]['path_sdevs']

        if (d == 0):
            label = 'Hor.'
        elif (d == 1):
            label = 'Ver.'
        else:
            raise ValueError('d must be either 0 or 1')

        t = self.data['EM_data'][q]['time_steps']

        d = d + 2  # Correct index for est_mean

        plt.fill_between(self.DT * np.arange(self.N_T),
                         est_mean[:, d] - est_sdev[:, d],
                         est_mean[:, d] + est_sdev[:, d],
                         alpha=0.5, linewidth=1.)
        plt.plot(self.DT * np.arange(self.N_T),
                 est_mean[:, d], label='estimate')
        plt.xlabel('Time (s)')
        plt.ylabel('Velocity (pixels/sec)')
        # plt.legend()

    def plot_dynamic_vars(self, q):
        """
        Plots all of the dynamic variables (x, y, vx, vy)
        as estimated after iteration q
        """
        if not self.data['motion_prior'] == 'VelocityDiffusion':
            raise RuntimeError('Run has no velocity estimate')

        plt.subplot(2, 2, 1)
        self.plot_path_estimate(q, 0)

        plt.subplot(2, 2, 2)
        self.plot_path_estimate(q, 1)

        plt.subplot(2, 2, 3)
        self.plot_velocity_estimate(q, 0)

        plt.subplot(2, 2, 4)
        self.plot_velocity_estimate(q, 1)

    def plot_EM_estimate(self, q):
        """
        Creates a full plot of the estimated image, actual image,
            estimate and actual paths

        q - index of EM iteration to plot
        """

#        if (q >= self.N_itr):
#            raise ValueError('Iteration index, q, is too large')

        m1 = np.min(self.S_gen)
        m2 = np.max(self.S_gen)

        vmin = -0.1 * (m2 - m1) + m1
        vmax = m2 + 0.1 * (m2 - m1)

        n_time_steps = self.data['EM_data'][q]['time_steps']
        t_ms = self.DT * n_time_steps * 1000.

        fig = plt.figure(figsize=(12, 8))
        fig.suptitle('EM Reconstruction after t = %0.2f ms for DC_gen = %0.2f, LAMBDA = %0.2f' % (
            t_ms, self.DC_gen, self.LAMBDA))

        # Note actual image is convolved with a gaussian during the simulation
        #   even though the image saved has not has this happen yet

        plt.subplot(2, 3, 2)
        self.plot_spikes(n_time_steps - 1, moving_average=True, mode='ON')

        plt.subplot(2, 3, 3)
        self.plot_spikes(n_time_steps - 1, moving_average=True, mode='OFF')

        plt.subplot(2, 3, 4)
        plt.title('Estimated Image, S = DA:\n SNR = %.2f'
                  % self.SNR_idx(q))
        S_est = self.data['EM_data'][q]['image_est']

        plt.imshow(gaussian_filter(S_est, self.blur_sdev),
                   cmap=plt.cm.gray, interpolation='nearest',
                   vmin=vmin, vmax=vmax)
        plt.colorbar()

        plt.subplot(2, 3, 1)
        self.plot_base_image()

        plt.subplot(2, 3, 5)
        self.plot_path_estimate(q, 0)

        plt.subplot(2, 3, 6)
        self.plot_path_estimate(q, 1)

        plt.tight_layout()
        plt.subplots_adjust(top=0.9)

    def plot_base_image(self):
        """
        Plot the original image that generates the data

        Parameters
        ----------
        sdev: float
            blurring of image
        """

        m1 = np.min(self.S_gen)
        m2 = np.max(self.S_gen)

        vmin = -0.1 * (m2 - m1) + m1
        vmax = m2 + 0.1 * (m2 - m1)

        plt.title('Stationary Object in the World')
        plt.imshow(gaussian_filter(self.S_gen, self.blur_sdev),
                   cmap=plt.cm.gray, interpolation='nearest',
                   vmin=vmin, vmax=vmax)
        plt.colorbar()

    def save_images(self):
        for q in range(self.N_itr):
            plt.clf()
            self.plot_EM_estimate(q)
            plt.savefig('img%d.png' % (100 + q))

    def compute_spike_moving_average(self, tau=0.005):
        """
        Computes the exponential moving average of the spikes
        tau - time constant of moving average, should be a multiple of self.DT
        Saves an array self.Rav -> EMA of firing rate for each neuron
        """
        rho = 1 - self.DT / tau
        Rav = np.zeros_like(self.R)

        Rav[:, 0] = self.R[:, 0]
        for i in range(1, self.N_T):
            Rav[:, i] = rho * Rav[:, i - 1] + (1 - rho) * self.R[:, i]

        self.Rav = Rav / self.DT

    def plot_spikes(self, t, moving_average=False, mode='ON'):
        """
        Plots the spiking profile at timestep number t
        t - timestep number to plot
        """
        if t > self.N_T:
            raise ValueError('time does not go past a certain time')
        if not moving_average:
            s = self.R[:, t]
        else:
            try:
                self.Rav
            except AttributeError:
                self.compute_spike_moving_average()
            s = self.Rav[:, t]

        nn = self.L_N ** 2 * 2
        if mode == 'OFF':
            spikes = s[0: nn / 2]
        elif mode == 'ON':
            spikes = s[nn / 2: nn]
        else:
            raise ValueError('mode must be ON or OFF')

        vmin = 0.
        vmax = 200.

        if moving_average:
            st = 'Spike ExpMA'
        else:
            st = 'Spikes'

        plt.imshow(spikes.reshape(self.L_N, self.L_N),
                   interpolation='nearest', cmap=plt.cm.gray_r,
                   vmin=vmin, vmax=vmax)
        plt.title('{} for {} Cells at time {}'.format(st, mode, t))
#        plt.colorbar()

    def plot_firing_rates(self, t, mode='ON'):
        """
        Plot the firing rates for each neuron.

        Note: The visualization makes the most sense when the RFs of the
            neurons form a rectangular grid
        """
        frs = self.data['FP'][0] / self.DT
        nn = self.L_N ** 2 * 2
        if mode == 'OFF':
            fr = frs[0: nn / 2, t]
        elif mode == 'ON':
            fr = frs[nn / 2: nn, t]
        else:
            raise ValueError('mode must be ON or OFF')

        plt.imshow(fr.reshape(self.L_N, self.L_N),
                   interpolation='nearest',
                   cmap=plt.cm.gray,
                   vmin=0, vmax=100.)
        # t_str = ('lambda(t) (Hz) for {} Cells'.format(mode))
        # plt.title(t_str)

    def plot_fr_and_spikes(self, t):
        """
        Plot the base image, firing rate, and exponential moving
            average of the spikes at time t
        t : int
            Timestep to plot
        """
        plt.figure(figsize=(10, 8))

        plt.subplot(2, 2, 1)
        self.plot_base_image()

        plt.subplot(2, 2, 2)
        self.plot_firing_rates(t, mode='ON')
        plt.title('Retinal Image')

        # Spikes
        plt.subplot(2, 2, 3)
        self.plot_spikes(t, mode='ON', moving_average=True)

        plt.subplot(2, 2, 4)
        self.plot_spikes(t, mode='OFF', moving_average=True)

    def plot_RFs(self):
        """
        Creates a plot of the receptive fields of the neurons
        """
        self.XE = self.data['XE']
        self.YE = self.data['YE']
#        self.IE = self.data['IE']
        self.Var = self.data['Var']
        std = np.sqrt(np.mean(self.Var))
        fig = plt.gcf()
        ax = plt.gca()
        ax.set_xlim((np.min(self.XE), np.max(self.XE)))
        ax.set_ylim((np.min(self.YE), np.max(self.YE)))
        for xe, ye in zip(self.XE, self.YE):
            circ = plt.Circle((xe, ye), std, color='b', alpha=0.4)
            fig.gca().add_artist(circ)

    def plot_tuning_curves(self, baseline_rate=10.):
        """
        Create a plot showing the tuning curves of the neurons
        """
        x = np.arange(0, 1, 0.01)
        l0 = self.data['L0']
        l1 = self.data['L1']
        y_on = np.exp(np.log(l0) + x * np.log(l1/l0))
        y_off = np.exp(np.log(l0) + (1 - x) * np.log(l1/l0))
        plt.plot(x, y_on, label='ON')
        plt.plot(x, y_off, label='OFF')
        plt.plot(x, baseline_rate + 0 * x, '--')
        plt.xlabel('Stimulus intensity')
        plt.ylabel('Firing Rate (Hz)')
        plt.title('Firing rate as a function of Stimulus Intensity')
        plt.legend()


if __name__ == '__main__':
    fn = sys.argv[1]
    da = DataAnalyzer.fromfilename(fn)

# plt.plot(SNRs)
# plt.show()
