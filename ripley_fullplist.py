#this script calculates Ripley's L function from SMLM localisations. It can be modified to only use pre-defined ROIs
# to speed calculation or restrict it to particular areas
#particle list should be in csv format as generated from ThunderSTORM

import numpy as np

#set the locations of the particle list and for the output

plist_path = ""
out_path = ""
dstep = 5
dmax = 30
area = 80000

def get_ripl_l(pts_roi, dstep, dmax, area):
    n = len(pts_roi)
    ds = np.arange(dstep, dmax + dstep, dstep)
    l2 = []
    #A = (max(pts_roi[:, 0]) - min(pts_roi[:, 0])) * (max(pts_roi[:, 1]) - min(pts_roi[:, 1]))
    A = area
    fac = A / (np.pi * n * (n - 1))
    difs = np.zeros((n, n))
    for p in range(n):
            difs[p][p+1:] = np.linalg.norm(pts_roi[p+1:] - pts_roi[p], axis=1)
    difs = difs.reshape(n*n)
    difsnz = difs[difs > 0]
    for d in ds:
        l2d = 2*sum(difsnz < d)
        l2d *= fac
        l2.append(l2d)
    l = np.sqrt(np.array(l2))
    return ds, l


def get_pts_roi(pts, roi):
    pts_roi = []
    for pt in pts:
        if roi[1] > pt[0] > roi[0] and roi[3] > pt[1] > roi[2]:
            pts_roi.append(pt)
    return np.array(pts_roi)


def extract_xy(rois4):
    rois = np.zeros((len(rois4), 4))
    for r in range(len(rois4)):
        rois[r][0] = min(rois4[r][1][:, 1])
        rois[r][1] = max(rois4[r][1][:, 1])
        rois[r][2] = min(rois4[r][1][:, 0])
        rois[r][3] = max(rois4[r][1][:, 0])
    return rois


def roi_l_image(roi, pts_roi, pxsize, d):

    l_im = np.zeros((int(roi[1]-roi[0])/pxsize + 1, int(roi[3]-roi[2])/pxsize + 1))
    x = 0
    for px in np.arange(roi[0],roi[1],pxsize):
        y=0
        for py in np.arange(roi[2],roi[3],pxsize):
            l_im[x,y] = sum(np.linalg.norm(pts_roi - np.array([px,py]), axis=1) < d)
            y += 1
        x += 1
    return l_im

#use this to compare distribution to a simulated random uniform distribution with npts localisations
def sim_uniform(npts, x0,x1,y0,y1):
    pts = np.column_stack((np.random.uniform(x0,x1,npts), np.random.uniform(y0,y1,npts)))
    return pts

if __name__ == '__main__':

    #import all the localisation data
    with open(plist_path, 'r') as f:
        data = np.genfromtxt(f, delimiter=",", names=True)
    #extract the x,y values and form an array
    pts = np.zeros((len(data), 2))
    for d, p in zip(data, pts):
        p[0] = d['x_nm']
        p[1] = d['y_nm']
    #calculate the L function
    ds, l = get_ripl_l(pts, dstep, dmax)

    l = np.array(l)
    print l
    np.savetxt(out_path, l, delimiter=",")




