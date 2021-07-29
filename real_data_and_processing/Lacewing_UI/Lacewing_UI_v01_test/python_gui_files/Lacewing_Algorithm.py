from Lacewing_Func import binary_file_read
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve
from skimage.measure import label as label_connected_graph
import pandas as pd
from pathlib import Path
import matplotlib as mpl
# mpl.use("Qt5Agg")
from sklearn.covariance import EmpiricalCovariance, MinCovDet
from tqdm.auto import tqdm
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

def binary_to_numpy(data, f_start, f_end):
    """Returns a 3-D object containing experiment data, as well as temperature and chemical values for each time step

    Parameters
    ----------
    data : list
        Flattened data from experiment
    f_start : int
        First frame to be considered
    f_end : int, optional
        Last frame to be considered

    Returns
    -------
    time_vect : list
        Sampling times in miliseconds
    temp_vect : list
        A list of average temperatures of the solution
    frame_3d : 3-dimensional array
        Representing the data as [x,y,time]
    """
    #Each frame contains 4372 datapoints, f_end is the number of frames
    f_end = len(data)//4372
    #Number of timesteps being considered
    n_time = f_end - f_start + 1
    #Initialize 3-D object to hold the value of each pixel for the duration being considered.
    #x=78, number of rows; y=56, number of columns; z=number of time steps.
    frame_3d = np.zeros((78, 56, n_time))
    #List to hold chemical pixel values
    time_vect = np.zeros(n_time)
    #List to hold temperature values
    temp_vect = np.zeros(n_time)

    #Loop through each time instance (n) to be considered, and keep track of the loop count (i)
    for i, n in enumerate(range(f_start - 1, f_end)):
        #Extract the nth frame from the raw data
        pack = data[n * 4372:(n + 1) * 4372]
        #Convert the first 4368 data points into numpy array
        frame_1d = np.array(pack[:4368])
        #Reshape the data into 78x56 to follow the shape of the array, use fortran indexing
        frame_2d = frame_1d.reshape(78, 56, order='F')
        #The 3rd axis (z) of the 3d object represents time, set one time slice equal to the reshaped data
        frame_3d[:, :, i] = frame_2d
        #Keep track of the 43678h Value (time elapsed*10)
        time_vect[i] = pack[4368] / 10
        #Keep track of the 4370th value (the average temperature)
        temp_vect[i] = pack[4370]
    #return average chemical output, 3d object, average temperature
    return time_vect, frame_3d, temp_vect


def split_into_n_wells(frame_3d, n=2):
    assert n == 2, "Only implemented 2 wells so far!"
    n_time = frame_3d.shape[2]
    top_well = frame_3d[:39, :, :].reshape(-1, n_time, order='C').T
    bot_well = frame_3d[-39:, :, :].reshape(-1, n_time, order='C').T
    return top_well, bot_well


def filter_by_vref(X, v_thresh=70):
    return (np.diff(X[:10, :], axis=0) > v_thresh).any(axis=0)


def filter_by_vrange(X, v_range=(100, 900)):
    return (X < v_range[1]).all(axis=0) & (X > v_range[0]).all(axis=0)


def time_to_index(times, time_vect):
    indices = []
    for time in times:
        indices.append( np.argmin(np.abs(time_vect - time)) )
    return indices


def find_loading_time(time_vect, X, bounds=(600, 900), viz=False):

    search_start, search_end = time_to_index(bounds, time_vect)

    X_mean = np.mean(X, axis=1)
    X_mean_diff = np.diff(X_mean)

    loading_index = np.argmax(X_mean_diff[search_start:search_end]) + search_start + 1
    loading_time = time_vect[loading_index]

    if viz:
        fig, ax = plt.subplots(3, 1)
        fig.suptitle('Finding Loading Time...')

        ax[0].set(title='Active Chemical Pixels, ACP')
        ax[0].plot(time_vect, X)

        ax[1].set(title='Mean(ACP)')
        ax[1].plot(time_vect, X_mean)
        ax[1].axvline(time_vect[search_start], color='C1')
        ax[1].axvline(time_vect[search_end], color='C1')
        ax[1].axvline(loading_time, color='C2')

        ax[2].set(title='Diff(Mean(ACP))')
        ax[2].plot(time_vect[1:], X_mean_diff)
        ax[2].axvline(time_vect[search_start], color='C1')
        ax[2].axvline(time_vect[search_end], color='C1')
        ax[2].axvline(loading_time, color='C2')

        plt.tight_layout()
        plt.show()

    return loading_index, loading_time


def find_settled_time(time_vect, X, viz=False):
    settled_idx = 10
    return settled_idx, time_vect[settled_idx]


def split_chem_and_temp(arr):
    """Returns the chemical and temperature pixel readings given a 3-D array of [x,y,time]

    Parameters
    ----------
    arr : 3-dimensional array
        [x,y,time] for this particular experiment

    Returns
    -------
    arr_temp : 3-dimensional array
        Value of temperature pixels at every time step
    arr_chem : 3-dimensional array
        Value of chemical pixels at every time step. Temperature pixels have been replaced by the average of the
        chemical pixels around them.
    """
    # Every 3 pixels, starting from index 1, is a temperature pixel (in both dimensions)
    arr_temp = arr[1::3, 1::3, :]
    #Create a convolution matrix, 0.125 per cell
    mask = np.ones((3, 3, 1)) / 8
    #Set the middle of the mask to 0. The output of this mask is the average of the 8 pixels around the middle.
    mask[1, 1, 0] = 0
    #Convolve the mask with the array keeping the output size tme same as the input object.
    av_3d = convolve(arr, mask, mode='same')
    #Copy the original object
    arr_chem = arr.copy()
    #Replace temperature pixels with the average of those around them
    arr_chem[1::3, 1::3, :] = av_3d[1::3, 1::3, :]
    #return the 2 new arrays
    return arr_temp, arr_chem


def cleanup_pixels(idx_active):

    idx_2d = idx_active.reshape(-1, 56)
    labels_2d = label_connected_graph(idx_2d, background=0)

    values, counts = np.unique(labels_2d.reshape(-1), return_counts=True)
    ind_back = values == 0
    values, counts = values[~ind_back], counts[~ind_back]

    ind = np.argmax(counts)
    max_idx = values[ind]

    return (labels_2d == max_idx).reshape(-1)


def cleanup_by_mahalanobis(X):
    robust_cov = EmpiricalCovariance().fit(X)
    mahal_robust_cov = robust_cov.mahalanobis(X)
    idx_good = mahal_robust_cov < np.percentile(mahal_robust_cov, 60)
    return idx_good



def largest_cluster(model):
    values, counts = np.unique(model.labels_, return_counts=True)
    return values[np.argmax(counts)]

def count_pixels_in_polygon(x_coord, y_coord, polygon):
    count = 0
    for x, y in zip(x_coord, y_coord):
        point = Point(x, y)
        count += polygon.contains(point)
    return count


def center_cluster(model, idx_active, viz=False):
    y1, y2, = 20, 56
    x1, x2 = 9, 29
    polygon = Polygon([(x1, y1), (x1, y2), (x2, y2), (x2, y1)])

    best_count = -1
    for label in np.unique(model.labels_):
        temp = idx_active.copy()
        temp[temp] = model.labels_ == label
        x_coord, y_coord = np.where(temp.reshape(-1, 56))
        if viz:
            plt.plot(*polygon.exterior.xy)
            plt.plot(x_coord, y_coord, "x")
            plt.xlim((0, 39))
            plt.ylim((0, 56))
            plt.show()

        count = count_pixels_in_polygon(x_coord, y_coord, polygon)
        if count > best_count:
            best_count = count
            best_cluster = label

    return best_cluster



# idx_active [0, 0, 0, 1, 1, 1, 1, 0]

# 1st cluster [0, 0, 0, 1, 1, 0, 0, 0]

# 2nd cluster [0, 0, 0, 0, 0, 1, 1, 0]


def cleanup_by_kmeans(well, idx_active, n_max=6, method="center"):
    X = well[:, idx_active].T
    X = X / np.linalg.norm(X, axis=1).reshape(-1, 1)

    best_score = -1
    for n_clusters in tqdm(range(2, n_max+1)):
        model = KMeans(n_clusters=n_clusters)
        model.fit(X)
        score = silhouette_score(X, model.labels_)
        if score > best_score:
            best_score = score
            best_n_clusters = n_clusters

    print(f"The number of clusters found: {best_n_clusters}")
    model = KMeans(n_clusters=best_n_clusters)
    model.fit(X)
    if method == "largest":
        idx_cluster = largest_cluster(model)
    elif method == "center":
        idx_cluster = center_cluster(model, idx_active, viz=False)
    else:
        raise NotImplemented(f"The method {method} has not been implemented!")

    idx_active[idx_active] = idx_cluster == model.labels_

    return idx_active


if __name__ == "__main__":
    #Experiments to analyse
    #experiments = ["6_27", "1_28", "2_14","Testt"]
    experiments = ["6_27"]
    for experiment in experiments:
        #Get .bin files from the relevant experiment folder (glob returns a list)
        filename = [i for i in Path("../data_files", experiment).glob("*.bin")][0]
        #Read binary file, fucntion returns a 1-D array of pixel values
        #Every 4372 values is one time step (3 seconds), data len is divisible by 4372
        data = binary_file_read(filename)
        #Initialize variables holding start and stop frames
        f_start, f_end = 1, 976  # This is set by Lei somewhere
        #Convert the list-representation of data to a useful 3D object, and lists of average temperatures and times
        time_vect, frame_3d, temp_vect = binary_to_numpy(data, f_start, f_end)
        #Extract shape of 3d object
        n_rows, n_cols, n_time = frame_3d.shape
        #Seperate checmical and temperature pixels
        arr_temp, arr_chem = split_chem_and_temp(frame_3d)

        top_well, bot_well = split_into_n_wells(arr_chem, n=2)

        top_well_temp, bot_well_temp = split_into_n_wells(arr_temp, n=2)

        find_active_pixels = lambda x: filter_by_vref(x) & filter_by_vrange(x)

        idx_top_active, idx_bot_active = [find_active_pixels(i) for i in [top_well, bot_well]]

        idx_top_active, idx_bot_active = cleanup_pixels(idx_top_active), cleanup_pixels(idx_bot_active)

        idx_top_active = cleanup_by_kmeans(top_well, idx_top_active)
        idx_bot_active = cleanup_by_kmeans(bot_well, idx_bot_active)

        loading_idx_top, loading_t_top = find_loading_time(time_vect, top_well[:, idx_top_active],
                                                           bounds=(600, 900), viz=False)

        loading_idx_bot, loading_t_bot = find_loading_time(time_vect, bot_well[:, idx_bot_active],
                                                           bounds=(600, 900), viz=False)

        settled_idx_top, settled_t_top = find_settled_time(time_vect[loading_idx_top:],
                                                           temp_vect[loading_idx_top:], viz=False)

        settled_idx_bot, settled_t_bot = find_settled_time(time_vect[loading_idx_bot:],
                                                           temp_vect[loading_idx_bot:], viz=False)

        processing_idx_top = loading_idx_top + settled_idx_top
        processing_idx_bot = loading_idx_bot + settled_idx_bot

        x_top, Y_top = time_vect[processing_idx_top:], top_well[processing_idx_top:, idx_top_active]
        x_bot, Y_bot = time_vect[processing_idx_bot:], bot_well[processing_idx_bot:, idx_bot_active]

        # temp_top, temp_bot = top_well_temp[processing_idx_top:], bot_well_temp[processing_idx_bot:]
        temp_top, temp_bot = top_well_temp, bot_well_temp

        N_axes = 1
        fig, ax = plt.subplots(N_axes, 2, figsize=(10, 3*N_axes))
        #fig.suptitle(f"Experiment Name: {experiment}")
        ax_top, ax_bot = ax[0], ax[1]#ax[:, 0], ax[:, 1]
        i = 0
        control = []
        pos = []

        for ax, x, Y, idx_active, label, temp in [(ax_top, x_top, Y_top, idx_top_active, 'Top', temp_top),
                                            (ax_bot, x_bot, Y_bot, idx_bot_active, 'Bot', temp_bot)]:
            # ax[0].set_title(f'{label} Well')
            # ax[0].imshow(idx_active.reshape(-1, 56))
            #
            # ax[1].set(title='Temp Pixels, TP')
            # ax[1].plot(time_vect, temp)
            #
            # ax[2].set(title='Active Chemical Pixels, ACP')
            # ax[2].plot(x, Y)
            #
            Y_bs = Y - np.mean(Y[:3, :], axis=0)
            # ax[3].set(title='ACP - ACP[0]')
            # ax[3].plot(x, Y_bs)

            Y_bs_smooth = pd.DataFrame(Y_bs)#.rolling(30).mean().values
            ax.set(title='Output of Pixels vs Time')
            ax.plot(x, Y_bs_smooth)
            ax.plot(x, np.mean(Y_bs_smooth, axis=1), lw=2, color='k', label='Mean')
            ax.legend()

            # Y_bs_smooth_diff = np.diff(Y_bs_smooth, axis=0)
            # ax[5].set(title='')
            # ax[5].plot(x[1:], Y_bs_smooth_diff)
            # ax[5].plot(x[1:], np.mean(Y_bs_smooth_diff, axis=1), lw=2, color="k", label='Mean')
            # ax[5].legend()
            #
            # ax[1].get_shared_x_axes().join(*ax[1:])
            if i == 0:
                control = Y_bs_smooth
                control_x = x
            else:
                sample = Y_bs_smooth
                sample_x = x
            i+=1

        plt.tight_layout()
        plt.show()

        # plt.plot(x, control)
        # plt.show()
        #
        import pickle
        with open("control_signals.txt", "wb") as fp:   #Pickling
            pickle.dump(control, fp)
        with open("control_time.txt", "wb") as fp1:   #Pickling
            pickle.dump(control_x, fp1)
        with open("sample_signals.txt", "wb") as fp2:
            pickle.dump(sample, fp2)
        with open("sample_time.txt", "wb") as fp3:
            pickle.dump(sample_x, fp3)
