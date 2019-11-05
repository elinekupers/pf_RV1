from EM_burak_class import *
from scipy.signal import convolve2d


emb = EMBurak(_DT=0.005, _DC=50.)
emb.gen_data()
emb.calculate_inner_products()


R = emb.R

L_N = emb.L_N
N_N = emb.N_N
N_T = emb.N_T
DT = emb.DT

a = 10.
F = np.exp(- np.arange(0, 3 * a + 1) / (3 * a))
F = F / np.sum(F)

m = F.shape[0]


R1 = np.zeros_like(R)

for i in range(R1.shape[0]):
    for t in range(N_T):
        if (t + m < N_T):
            R1[i, t:t + m] += R[i, t] * F
        else:
            R1[i, t:N_T] += R[i, t] * F[0: N_T - t]


for t in range(50):
    print t
    plt.figure().suptitle('t = ' + str(1000 * DT * t) + ' ms')  # + '\n' +
    #                        'Lambda0,Lambda1 = (' + str(emb.L0) + ',' + str(emb.L1) + ')'
    #                        )

    plt.subplot(2, 2, 1)
    plt.imshow(emb.Ips[0, :, t].reshape(L_N, L_N),
               interpolation='nearest', cmap=plt.cm.gray_r)
    plt.title('Retinal Projection of the Image')

    plt.subplot(2, 2, 2)
    plt.imshow(R[:, t].reshape(L_N, L_N),
               interpolation='nearest', cmap=plt.cm.gray_r)
    plt.title('Spikes')

    plt.subplot(2, 2, 3)
    plt.imshow(R1[:, t].reshape(L_N, L_N),
               interpolation='nearest', cmap=plt.cm.gray_r)
    plt.title('Exponential Moving Average of Spikes')
    plt.savefig('img/retinal_projection/img' + str(1000 + t) + '.png', dpi=75)


#mean = np.mean(R, axis = 1)
# plt.imshow(mean.reshape(L_N, L_N),
#           interpolation = 'nearest', cmap = plt.cm.gray_r)
# plt.show()
