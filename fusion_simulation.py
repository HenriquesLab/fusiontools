import numpy as np
import matplotlib.pyplot as plt



khf = 3
kfuse = 0.1
Nvirions = 1000
Nefc = 100
nfuse =3
#set the output directory
outfilename = ""



#sample randomly from an exponentially decaying function
def sample_k(k,N):

    out = np.random.exponential(1/k, size=N)
    return out

def get_times(khf,kfuse,Nvirions,Nefc, nfuse):


    themis_actual = np.zeros(Nvirions)
    tfuses_actual = np.zeros(Nvirions)
    ttotals_actual = np.zeros(Nvirions)

    for i in range(Nvirions):

        #for each virus calculate the hemifusion and full fusion time as the time for nfuse EFCs to activate
        themis = sample_k(khf, Nefc)
        themi_actual = min(themis)
        tfuses = sample_k(kfuse, Nefc)
        tfuse_actual = min(tfuses)
        total_sorted = np.sort(themis+tfuses)
        ttotal_actual = total_sorted[nfuse-1]
        ttotal_actual = min(themis+tfuses)

        themis_actual[i] = themi_actual
        tfuses_actual[i] = tfuse_actual
        ttotals_actual[i] = ttotal_actual

    return (themis_actual,ttotals_actual)

def getcumfreqvals(data, xmax, xstep):
    vals = []
    for x in np.arange(0,xmax,xstep):
        v = np.sum(data < x)
        vals.append(v)
    vals = np.array(vals)
    return vals

def plot_cumfreq(datas, xmax):
    colors = "bgrcmykw"
    color_index = 0
    for data in datas:
        values, base = np.histogram(data, bins=40)
        if(max(base)>xmax):
            xlim = np.argmax(base>xmax)
        else:
            xlim = len(base)-1
        cumulative = np.cumsum(values)
        plt.plot(base[:xlim], cumulative[:xlim], c=colors[color_index])
        color_index+=1
    plt.show()

if __name__ == "__main__":
    times = get_times(khf,kfuse,Nvirions,Nefc,nfuse)
    np.savetxt(outfilename, times[1])
