import os

import netCDF4
import numpy as np

import matplotlib.pyplot as plt

import pygmt

import csv

import os

data_folder = os.path.dirname(__file__)

def showDrifters():
    folder = f"{data_folder}/world3d_txt"

    fig = pygmt.Figure()

    fig.coast(region="d",
              projection="Q12c",
              land="lightgray",
              water="white",
              borders="1/0.5p",
              shorelines="1/0.5p",
              frame="ag",
              )

    count = 0
    for file_path in os.listdir(folder):
        with open(f"{folder}/{file_path}", "r") as f:
            data = np.array(list(csv.reader(f, delimiter=" "))).astype(float)
            
            lat, lon = worldToLL(data[:,0], data[:,1], data[:,2])

            fig.plot(region="d",
                     projection="Q12c",
                     x=lon,
                     y=lat,
                     pen="darkblue")
            print(count)
            count += 1
            # if count > 100:
            #     break
    fig.show()

def llToWorld(lat,lon):
    lat, lon = np.deg2rad(lat), np.deg2rad(lon)
    R = 6371000 # radius of the earth
    x = R * np.cos(lat) * np.cos(lon)
    y = R * np.cos(lat) * np.sin(lon)
    z = R *np.sin(lat)
    return np.array([x,y,z])

def worldToLL(x, y, z):
    R = 6371000 # radius of the earth
    lat = np.arcsin(z / R) / np.pi * 180
    lon = np.arctan2(y, x) / np.pi * 180
    return np.array([lat, lon])

withTransform = False

def _no_outlier(lat1, long1, lat2, long2, lat3, long3, factor):
    global withTransform
    if withTransform:
        dist_outlier = np.linalg.norm(llToWorld(float(lat2),float(long2)) - llToWorld(float(lat1),float(long1)))
        dist_neighbo = np.linalg.norm(llToWorld(float(lat3),float(long3)) - llToWorld(float(lat1),float(long1)))
        return (dist_outlier < factor * dist_neighbo)
    else:
        dist_outlier = np.linalg.norm(np.array([float(lat2),float(long2)]) - np.array([float(lat1),float(long1)]))
        dist_neighbo = np.linalg.norm(np.array([float(lat3),float(long3)]) - np.array([float(lat1),float(long1)]))
        return (dist_outlier < factor * dist_neighbo)

def no_outlier(lat, long, i, factor):
    return _no_outlier(lat[i-1],long[i-1],lat[i],long[i],lat[i+1],long[i+1],factor)


def has_gap(lat, long, i, factor, absolute_value_max, absolute_value_min):
    global withTransform
    if(withTransform):
        dist_to_new = np.linalg.norm(
            llToWorld(float(lat[i]),float(long[i])) - llToWorld(float(lat[i-1]),float(long[i-1])))
        dist_to_next = np.linalg.norm(
            llToWorld(float(lat[i+1]),float(long[i+1])) - llToWorld(float(lat[i]),float(long[i])))
        return (((dist_to_new > absolute_value_max)))
    else:
        dist_to_new = np.linalg.norm(np.array([float(lat[i]),float(long[i])])-np.array([float(lat[i-1]),float(long[i-1])]))
        dist_to_next = np.linalg.norm(np.array([float(lat[i+1]),float(long[i+1])])-np.array([float(lat[i]),float(long[i])]))
        return (((dist_to_new>absolute_value_max)))

def convertToWorld3d():
    global withTransform
    dist = 10000 if withTransform else 5
    counter = 1
    folder = f"{data_folder}/world3d_txt"
    hasWrittenYet = False

    # remove existing files
    existingFiles = [ f for f in os.listdir(folder) if f.endswith(".txt") ]
    for f in existingFiles:
        counter += 1
        os.remove(os.path.join(folder, f))
    print(f"Removed {counter} existing files from {folder}")

    counter = 1
    with open(f"{data_folder}/drifter_6hour_qc_8951_b49f_cfce.csv") as csvfile:
        spamreader = csv.reader(csvfile,delimiter=",")
        data = []
        curID = -1
        cur = []
        i = 0
        for row in spamreader:
            i+=1
            if i <= 2:
                print(row)
                continue
            id = row[0]
            loc = [float(row[2]),float(row[3])]

            if curID == -1:
                curID = id

            if id == curID:
                cur.append(loc)

            else:
                data.append(np.array(cur))
                curID = id
                cur = [loc]
        for d in data:
            lats = d[:,0]
            longs = d[:,1]
            newLat = np.array([lats[0]])
            newLong = np.array([longs[0]])
            for i in range(1, lats.shape[0] - 1):
                #filter into newLat and newLong
                if lats[i] != lats[i - 1] and lats[i] != 'nan' and lats[i] != '0.0' and lats[i] != '90.0' and lats[i] != '-90.0' and \
                        lats[i] != '180.0' and lats[i] != '-180.0':
                    if _no_outlier(newLat[-1],newLong[-1],lats[i],longs[i],lats[i+1],longs[i+1],5):
                        newLat = np.append(newLat,lats[i])
                        newLong = np.append(newLong,longs[i])
            lat = newLat[1:]
            long = newLong[1:]
            if(lat.shape[0]<2):
                continue
            if hasWrittenYet:
                counter += 1
                hasWrittenYet = False
            for i in range(1, lat.shape[0] - 1):
                if lat[i] != lat[i - 1] and lat[i] != 'nan' and lat[i] != '0.0' and lat[i] != '90.0' and lat[i] != '-90.0' and \
                        lat[i] != '180.0' and lat[i] != '-180.0':
                    if no_outlier(lat, long, i, 5) and no_outlier(lat, long, i-1, 5):
                        if has_gap(lat, long, i, 10, dist, 0.25):
                            if hasWrittenYet:
                                counter += 1
                                hasWrittenYet = False
                        with open(f"{folder}/{counter}_drifter.txt", "a") as f:
                            xyz = llToWorld(float(lat[i]),float(long[i]))
                            f.write(str(xyz[0])  + ' ' + str(xyz[1])+ ' ' + str(xyz[2]) + '\n')
                            #f.write(str(lat[i]) + ' ' + str(long[i]) + '\n')
                            hasWrittenYet = True
            i = lat.shape[0] - 1
            if lat[i] != lat[i - 1] and lat[i] != 'nan' and lat[i] != '0.0' and lat[i] != '90.0' and lat[i] != '-90.0' and lat[
                i] != '180.0' and lat[i] != '-180.0':
                with open(f"{folder}/{counter}_drifter.txt", "a") as f:
                    xyz = llToWorld(float(lat[i]),float(long[i]))
                    f.write(str(xyz[0])  + ' ' + str(xyz[1])+ ' ' + str(xyz[2]) + '\n')
                    #f.write(str(lat[i]) + ' ' + str(long[i]) + '\n')
                    hasWrittenYet = True
            print(f"Count: {counter}")

convertToWorld3d()
# showDrifters()