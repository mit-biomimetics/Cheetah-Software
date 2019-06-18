import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import os

# Plot configuration
PLOT_VERTICALLY = 0
PLOT_HORIZONTALLY = 1

# number of figures in this plot
num_figures = 5

def create_figures(subfigure_width=480, subfigure_height=600, starting_figure_no=1, starting_col_index = 0, starting_row_index=0, plot_configuration=PLOT_HORIZONTALLY):
    figure_number = starting_figure_no
    col_index = starting_col_index
    row_index = starting_row_index

    file_path = os.getcwd() + "/../sim_data/"

    ## read files
    front_contact = np.genfromtxt(file_path+'front_contact.txt', delimiter=None, dtype=(float))
    hind_contact = np.genfromtxt(file_path+'hind_contact.txt', delimiter=None, dtype=(float))

    Fr_des = np.genfromtxt(file_path+'Fr_des.txt', delimiter=None, dtype=(float))
    Fr_result = np.genfromtxt(file_path+'Fr_result.txt', delimiter=None, dtype=(float))

    fr_Fr_des = Fr_des[:, 0:3];
    fl_Fr_des = Fr_des[:, 3:6];
    hr_Fr_des = Fr_des[:, 6:9];
    hl_Fr_des = Fr_des[:, 9:12];

    fr_Fr_result = Fr_result[:, 0:3];
    fl_Fr_result = Fr_result[:, 3:6];
    hr_Fr_result = Fr_result[:, 6:9];
    hl_Fr_result = Fr_result[:, 9:12];

    data_x = np.genfromtxt(file_path+'time.txt', delimiter='\n', dtype=(float))

    scale = 0.1;
    st_idx = 0
    end_idx = len(data_x) - 10;
    data_x = data_x[st_idx:end_idx]

    axes = plt.gca()
    ################################################
    ##  POSITION
    ## plot front right foot position
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('front right foot pos')
    for i in range(3):
        ax1 = plt.subplot(3, 1, i+1)
        plt.plot(data_x, fr_Fr_des[st_idx:end_idx,i], "r-")
        plt.plot(data_x, fr_Fr_result[st_idx:end_idx,i], "b-")
        plt.plot(data_x, scale * front_contact[st_idx:end_idx], "c-")
        plt.plot(data_x, scale * hind_contact[st_idx:end_idx], "k-")
        plt.grid(True)
    plt.xlabel('time (sec)')
    ## increment figure number and index
    figure_number += 1
    if plot_configuration == PLOT_HORIZONTALLY:
        col_index += 1
    elif plot_configuration == PLOT_VERTICALLY:
        row_index +=1

    ## plot front left foot position
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('front left foot pos')
    for i in range(3):
        ax1 = plt.subplot(3, 1, i+1)
        plt.plot(data_x, fl_Fr_des[st_idx:end_idx,i], "r-")
        plt.plot(data_x, fl_Fr_result[st_idx:end_idx,i], "b-")
        plt.plot(data_x, scale * front_contact[st_idx:end_idx], "c-")
        plt.plot(data_x, scale * hind_contact[st_idx:end_idx], "k-")
        plt.grid(True)
    plt.xlabel('time (sec)')
    ## increment figure number and index
    figure_number += 1
    if plot_configuration == PLOT_HORIZONTALLY:
        col_index += 1
    elif plot_configuration == PLOT_VERTICALLY:
        row_index +=1

    ## plot hind right foot position
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('hind right foot pos')
    for i in range(3):
        ax1 = plt.subplot(3, 1, i+1)
        plt.plot(data_x, hr_Fr_des[st_idx:end_idx,i], "r-")
        plt.plot(data_x, hr_Fr_result[st_idx:end_idx,i], "b-")
        plt.plot(data_x, scale * front_contact[st_idx:end_idx], "c-")
        plt.plot(data_x, scale * hind_contact[st_idx:end_idx], "k-")
        plt.grid(True)
    plt.xlabel('time (sec)')
    ## increment figure number and index
    figure_number += 1
    if plot_configuration == PLOT_HORIZONTALLY:
        col_index += 1
    elif plot_configuration == PLOT_VERTICALLY:
        row_index +=1

    ## plot hind left foot position
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('hind left foot pos')
    for i in range(3):
        ax1 = plt.subplot(3, 1, i+1)
        plt.plot(data_x, hl_Fr_des[st_idx:end_idx,i], "r-")
        plt.plot(data_x, hl_Fr_result[st_idx:end_idx,i], "b-")
        plt.plot(data_x, scale * front_contact[st_idx:end_idx], "c-")
        plt.plot(data_x, scale * hind_contact[st_idx:end_idx], "k-")
        plt.grid(True)
    plt.xlabel('time (sec)')
    ## increment figure number and index
    figure_number += 1
    if plot_configuration == PLOT_HORIZONTALLY:
        col_index += 1
    elif plot_configuration == PLOT_VERTICALLY:
        row_index +=1


if __name__ == "__main__":
    create_figures()
    plt.show()
