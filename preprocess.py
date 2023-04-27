import pandas as pd
import numpy as np
import matplotlib as mpl
import seaborn as sns
from scipy.stats import spearmanr
from collections import Counter
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.fftpack import fft, ifft


# normalize
def norm(arr):
    ax = np.max(arr, axis=0)
    m = np.min(arr, axis=0)
    r = ax - m
    return (arr - m) / r


# read_folder
def distance(a, b):
    dist = [(a[i] - b[i]) ** 2 for i in range(len(a))]
    return sum(dist) ** 0.5


def angle(s1, a, s2):
    va = s1 - a
    vb = s2 - a
    ab = float(np.sum(va * vb))
    abl = (float(np.sum(va ** 2)) ** 0.5) * (float(np.sum(vb ** 2)) ** 0.5)
    cos_theta = ab / abl
    theta = (np.arccos(cos_theta) / np.pi) * 180
    return theta


def pca(X, k):  # k is the components you want
    # mean of each feature
    n_samples, n_features = X.shape
    mean = np.array([np.mean(X[:, i]) for i in range(n_features)])
    # normalization
    norm_X = X - mean
    # scatter matrix
    scatter_matrix = np.dot(np.transpose(norm_X), norm_X)
    # Calculate the eigenvectors and eigenvalues
    eig_val, eig_vec = np.linalg.eig(scatter_matrix)
    eig_pairs = [(np.abs(eig_val[i]), eig_vec[:, i]) for i in range(n_features)]
    # sort eig_vec based on eig_val from highest to lowest
    eig_pairs.sort(reverse=True)
    # select the top k eig_vec
    feature = np.array([ele[1] for ele in eig_pairs[:k]])
    # get new data
    data = np.dot(norm_X, np.transpose(feature))
    return data


def dist_hm(arr):
    matrix = []
    for i in arr:
        line = []
        for j in arr:
            line.append(distance(i, j))
        matrix.append(line)
    return matrix


def get_angle(arr):
    ang = []
    for i in range(len(arr) - 2):
        a = arr[i]
        b = arr[i + 1]
        c = arr[i + 2]
        ang.append(angle(a, b, c))
    return ang


def rm_zero(x, y):
    dt = list(zip(x, y))
    dt = [n for n in dt if n[1] != 0]
    x, y = zip(*dt)
    return np.asarray(x), np.asarray(y)


def de_trend(mat):
    def fun(x, a, b, c):
        return c * np.exp(-a * x) + b

    time_series = np.linspace(0, mat.shape[1] - 1, num=mat.shape[1])
    new = []
    for data in mat:
        # remove zero point
        x, y = rm_zero(list(time_series), list(data))
        f, err = curve_fit(fun, x, y, p0=[0.003, 0.5, 0.5], maxfev=5000)
        curve = fun(time_series, f[0], f[1], f[2])
        new.append(curve)
    new = np.asarray(new)
    mat = mat / new
    return mat, new


def cor(mat, n):
    score = []
    pos = []
    neg = []
    null = []
    name = list(mat.columns)
    for i in name:
        for j in name:
            if i != j and i.isdigit()==False and j.isdigit()==False:
                s = spearmanr(mat[i], mat[j])[0]
                score.append(s)
                if s < -n:
                    neg.append((i, j))
                elif s > n:
                    pos.append((i, j))
                else:
                    null.append((i, j))
    return score, pos, neg, null


def de_noise(mat):
    mat = mat
    a = np.array([0.1] * 10)
    new = []
    for i in mat:
        new.append(np.convolve(i, a, mode='valid'))
    return np.asarray(new)


def elect(df, n):  # mat = norm(time x neuron)
    ind = list(df.var().sort_values()[-n:].index)
    mat = pd.DataFrame()
    for i in ind:
        mat = pd.concat([mat, df[i]], axis=1)
    return mat


def draw():
    ax1 = plt.subplot(5, 1, 1)
    plt.plot(fret_hm[0].T)
    ax2 = plt.subplot(5, 1, 2)
    plt.plot(fret_hm_dn[0].T)
    ax3 = plt.subplot(5, 1, 3)
    plt.plot(curve[0].T)
    ax3 = plt.subplot(5, 1, 4)
    plt.plot(fret_hm_dt[0].T)
    ax4 = plt.subplot(5, 1, 5)
    plt.plot(fret_hm_nm[0])
    plt.show()


def draw_pca(pca_d):
    dt = pca_d.T
    x = dt[0]
    y = dt[1]
    z = dt[2]
    plt.axes(projection='3d')
    plt.plot(x, y, z)
    plt.show()


def fft_filter(mat):
    new = []
    for i in mat:
        line = []
        n = np.abs(fft(i))
        for j in n:
            if 10 <= j <= 30:
                line.append(j)
            else:
                line.append(0)
        new.append(ifft(line))
    return np.asarray(new)


def fft_filter(mat):
    new = []
    for i in mat:
        line = []
        n = np.abs(fft(i))
        for j in n:
            if 10 <= j <= 100:
                line.append(j)
            else:
                line.append(0)
        line = ifft(line)
        line = line[50:len(line) - 50]
        line = line[:int(len(line) / 2)]
        new.append(line)
    return abs(np.asarray(new))

def correlation(cor_mat):
    lst = []
    for i in cor_mat:
        lst.extend(i)
    return dict(Counter(lst))


def histplot(data,label):
    dt = data.flatten()
    dt = set(dt)
    dt = [n for n in dt if n != 0]
    sns.histplot(dt,stat='probability',label=label)


def main(ch1_p,ch2_p):
    ch1 = pd.read_csv(ch1_p, index_col=0)
    ch2 = pd.read_csv(ch2_p, index_col=0)
    name = ch1.columns
    # or can manually remove error values
    # ch1 = ch1[120:]
    # ch2 = ch2[120:]
    ch1_d = ch1.values
    ch2_d = ch2.values
    fret_hm = ch2_d / ch1_d
    fret_hm = fret_hm.T
    fret_hm_dn = de_noise(fret_hm)
    fret_hm_dt, curve = de_trend(fret_hm_dn)
    fret_hm_nm = norm(fret_hm_dt.T)
    fret = pd.DataFrame(fret_hm_nm, columns=name)
    return fret
