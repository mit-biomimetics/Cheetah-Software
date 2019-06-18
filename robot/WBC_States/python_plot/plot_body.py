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
    data_body = \
            np.genfromtxt(file_path+'body_pos.txt', delimiter=None, dtype=(float))
    cmd_body_ori_rpy = \
            np.genfromtxt(file_path+'cmd_body_ori_rpy.txt', delimiter=None, dtype=(float))
    body_ori_rpy = \
            np.genfromtxt(file_path+'body_ori_rpy.txt', delimiter=None, dtype=(float))
    body_vel = \
            np.genfromtxt(file_path+'body_vel.txt', delimiter=None, dtype=(float))
    body_acc = \
            np.genfromtxt(file_path+'body_acc.txt', delimiter=None, dtype=(float))
    data_config = \
            np.genfromtxt(file_path+'config.txt', delimiter=None, dtype=(float))
    data_qdot = \
            np.genfromtxt(file_path+'qdot.txt', delimiter=None, dtype=(float))

    data_x = np.genfromtxt(file_path+'time.txt', delimiter='\n', dtype=(float))

    st_idx = 0
    end_idx = len(data_x)
    data_x = data_x[st_idx:end_idx]
    # PHASE MARKER #
    data_phse = np.genfromtxt(file_path+'phase.txt', delimiter=None, dtype=(float))
    # get phase.txt data #
    phseChange = []
    for i in range(0,len(data_x)-1):
            if data_phse[i] != data_phse[i+1]:
                phseChange.append(i - st_idx)
            else:
                pass
    axes = plt.gca()

    ## plot position
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('body pos')
    for i in range(1,4,1):
        ax1 = plt.subplot(3, 1, i)
        plt.plot(data_x, data_body[st_idx:end_idx,i-1], "b-")
        plt.plot(data_x, data_config[st_idx:end_idx,i+2], "r-")
        plt.grid(True)
        for j in phseChange:
            # phase line
            plt.axvline(x=data_x[j],color='indigo',linestyle='-')
            # phase number
            plt.text(data_x[j],ax1.get_ylim()[1],'%d'%(data_phse[j]),color='indigo')
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
    fig.canvas.set_window_title('body vel')
    for i in range(1,4,1):
        ax1 = plt.subplot(3, 1, i)
        plt.plot(data_x, body_vel[st_idx:end_idx,i-1], "b-")
        plt.plot(data_x, data_qdot[st_idx:end_idx,i+2], "r-")
        plt.grid(True)
        for j in phseChange:
            # phase line
            plt.axvline(x=data_x[j],color='indigo',linestyle='-')
            # phase number
            plt.text(data_x[j],ax1.get_ylim()[1],'%d'%(data_phse[j]),color='indigo')
    ## increment figure number and index
    figure_number += 1
    if plot_configuration == PLOT_HORIZONTALLY:
        col_index += 1
    elif plot_configuration == PLOT_VERTICALLY:
        row_index +=1
 
    ## plot acceleration
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))    
    fig.canvas.set_window_title('body acc')
    for i in range(1,4,1):
        ax1 = plt.subplot(3, 1, i)
        plt.plot(data_x, body_acc[st_idx:end_idx,i-1], "b-")
        plt.grid(True)
        for j in phseChange:
            # phase line
            plt.axvline(x=data_x[j],color='indigo',linestyle='-')
            # phase number
            plt.text(data_x[j],ax1.get_ylim()[1],'%d'%(data_phse[j]),color='indigo')
    ## increment figure number and index
    figure_number += 1
    if plot_configuration == PLOT_HORIZONTALLY:
        col_index += 1
    elif plot_configuration == PLOT_VERTICALLY:
        row_index +=1

    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))    
    fig.canvas.set_window_title('body ori')
    for i in range(3):
        ax1 = plt.subplot(3, 1, i+1)
        plt.plot(data_x, cmd_body_ori_rpy[st_idx:end_idx,i], "r-")
        plt.plot(data_x, body_ori_rpy[st_idx:end_idx,i], "b-")
        plt.grid(True)
        for j in phseChange:
            # phase line
            plt.axvline(x=data_x[j],color='indigo',linestyle='-')
            # phase number
            plt.text(data_x[j],ax1.get_ylim()[1],'%d'%(data_phse[j]),color='indigo')
    ## increment figure number and index
    figure_number += 1
    if plot_configuration == PLOT_HORIZONTALLY:
        col_index += 1
    elif plot_configuration == PLOT_VERTICALLY:
        row_index +=1
 
    
if __name__ == "__main__":
    create_figures()
    plt.show()
