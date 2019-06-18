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
    jpos_cmd = np.genfromtxt(file_path+'jpos_cmd.txt', delimiter=None, dtype=(float))
    jvel_cmd = np.genfromtxt(file_path+'jvel_cmd.txt', delimiter=None, dtype=(float))
    jacc_cmd = np.genfromtxt(file_path+'jacc_cmd.txt', delimiter=None, dtype=(float))

    config = np.genfromtxt(file_path+'config.txt', delimiter=None, dtype=(float))
    qdot = np.genfromtxt(file_path+'qdot.txt', delimiter=None, dtype=(float))

    x = np.genfromtxt(file_path+'time.txt', delimiter='\n', dtype=(float))

    st_idx = 0
    end_idx = len(x)- 10
    x = x[st_idx:end_idx]

    axes = plt.gca()

    ## plot position
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('front leg pos')
    for i in range(6):
        ax1 = plt.subplot(6, 1, i+1)
        plt.plot(x, jpos_cmd[st_idx:end_idx,i], "r-")
        plt.plot(x, config[st_idx:end_idx,i+6], "b-")
        plt.grid(True)
    plt.xlabel('time (sec)')
    ## increment figure number and index
    figure_number += 1
    if plot_configuration == PLOT_HORIZONTALLY:
        col_index += 1
    elif plot_configuration == PLOT_VERTICALLY:
        row_index +=1

    ## plot position
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('hind leg pos')
    for i in range(6):
        ax1 = plt.subplot(6, 1, i+1)
        plt.plot(x, jpos_cmd[st_idx:end_idx,i+6], "r-")
        plt.plot(x, config[st_idx:end_idx,i+12], "b-")
        plt.grid(True)
    plt.xlabel('time (sec)')
    ## increment figure number and index
    figure_number += 1
    if plot_configuration == PLOT_HORIZONTALLY:
        col_index += 1
    elif plot_configuration == PLOT_VERTICALLY:
        row_index +=1

    ## plot velocity
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('front leg vel')
    for i in range(6):
        ax1 = plt.subplot(6, 1, i+1)
        plt.plot(x, jvel_cmd[st_idx:end_idx,i], "r-")
        plt.plot(x, qdot[st_idx:end_idx,i+6], "b-")
        plt.grid(True)
    plt.xlabel('time (sec)')
    ## increment figure number and index
    figure_number += 1
    if plot_configuration == PLOT_HORIZONTALLY:
        col_index += 1
    elif plot_configuration == PLOT_VERTICALLY:
        row_index +=1

    ## plot velocity
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('hind leg vel')
    for i in range(6):
        ax1 = plt.subplot(6, 1, i+1)
        plt.plot(x, jvel_cmd[st_idx:end_idx,i+6], "r-")
        plt.plot(x, qdot[st_idx:end_idx,i+12], "b-")
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
