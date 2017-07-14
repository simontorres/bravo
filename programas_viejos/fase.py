import numpy as np

def calculo_fase(mag, date, er, per, T0):
    fase2, tt2, err, t = [], [], [], []
    for i in range(len(mag)):
        fa2 = ((float(date[i]) * 1.0 - T0) / per) - int((float(date[i]) - T0) / per)
        if fa2 > 0:
            fase2.append(fa2)
            tt2.append(fa2 + 1.0)
            t.append(date[i])
            err.append(er[i])
        else:
            fase2.append(fa2 + 1)
            tt2.append(fa2 + 2)
            t.append(date[i])
            err.append(er[i])
    fase2 = np.array(fase2)
    tt2 = np.array(tt2)
    err = np.array(err)
    re2 = np.concatenate((fase2, tt2), axis=0)
    mag3 = np.concatenate((mag, mag), axis=0)
    t2 = np.concatenate((t, t), axis=0)
    err2 = np.concatenate((err, err), axis=0)

    return re2, mag3, t2, err2