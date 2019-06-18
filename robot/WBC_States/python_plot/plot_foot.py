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
    fr_pos = np.genfromtxt(file_path+'fr_pos.txt', delimiter=None, dtype=(float))
    fl_pos = np.genfromtxt(file_path+'fl_pos.txt', delimiter=None, dtype=(float))
    hr_pos = np.genfromtxt(file_path+'hr_pos.txt', delimiter=None, dtype=(float))
    hl_pos = np.genfromtxt(file_path+'hl_pos.txt', delimiter=None, dtype=(float))

    cmd_fr_pos = np.genfromtxt(file_path+'cmd_fr_pos.txt', delimiter=None, dtype=(float))
    cmd_fl_pos = np.genfromtxt(file_path+'cmd_fl_pos.txt', delimiter=None, dtype=(float))
    cmd_hr_pos = np.genfromtxt(file_path+'cmd_hr_pos.txt', delimiter=None, dtype=(float))
    cmd_hl_pos = np.genfromtxt(file_path+'cmd_hl_pos.txt', delimiter=None, dtype=(float))

    fr_vel = np.genfromtxt(file_path+'fr_vel.txt', delimiter=None, dtype=(float))
    fl_vel = np.genfromtxt(file_path+'fl_vel.txt', delimiter=None, dtype=(float))
    hr_vel = np.genfromtxt(file_path+'hr_vel.txt', delimiter=None, dtype=(float))
    hl_vel = np.genfromtxt(file_path+'hl_vel.txt', delimiter=None, dtype=(float))

    cmd_fr_vel = np.genfromtxt(file_path+'cmd_fr_vel.txt', delimiter=None, dtype=(float))
    cmd_fl_vel = np.genfromtxt(file_path+'cmd_fl_vel.txt', delimiter=None, dtype=(float))
    cmd_hr_vel = np.genfromtxt(file_path+'cmd_hr_vel.txt', delimiter=None, dtype=(float))
    cmd_hl_vel = np.genfromtxt(file_path+'cmd_hl_vel.txt', delimiter=None, dtype=(float))


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

    ################################################
    ##  POSITION
    ## plot front right foot position
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('front right foot pos')
    for i in range(3):
        ax1 = plt.subplot(3, 1, i+1)
        plt.plot(data_x, cmd_fr_pos[st_idx:end_idx,i], "r-")
        plt.plot(data_x, fr_pos[st_idx:end_idx,i], "b-")
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

    ## plot front left foot position
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('front left foot pos')
    for i in range(3):
        ax1 = plt.subplot(3, 1, i+1)
        plt.plot(data_x, cmd_fl_pos[st_idx:end_idx,i], "r-")
        plt.plot(data_x, fl_pos[st_idx:end_idx,i], "b-")
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

    ## plot hind right foot position
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('hind right foot pos')
    for i in range(3):
        ax1 = plt.subplot(3, 1, i+1)
        plt.plot(data_x, cmd_hr_pos[st_idx:end_idx,i], "r-")
        plt.plot(data_x, hr_pos[st_idx:end_idx,i], "b-")
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

    ## plot hind left foot position
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('hind left foot pos')
    for i in range(3):
        ax1 = plt.subplot(3, 1, i+1)
        plt.plot(data_x, cmd_hl_pos[st_idx:end_idx,i], "r-")
        plt.plot(data_x, hl_pos[st_idx:end_idx,i], "b-")
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

    ################################################
    ##  VELOCITY
    ## plot front right foot velocity
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('front right foot vel')
    for i in range(3):
        ax1 = plt.subplot(3, 1, i+1)
        plt.plot(data_x, cmd_fr_vel[st_idx:end_idx,i], "r-")
        plt.plot(data_x, fr_vel[st_idx:end_idx,i], "b-")
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

    ## plot front left foot velocity
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('front left foot vel')
    for i in range(3):
        ax1 = plt.subplot(3, 1, i+1)
        plt.plot(data_x, cmd_fl_vel[st_idx:end_idx,i], "r-")
        plt.plot(data_x, fl_vel[st_idx:end_idx,i], "b-")
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

    ## plot hind right foot velocity
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('hind right foot vel')
    for i in range(3):
        ax1 = plt.subplot(3, 1, i+1)
        plt.plot(data_x, cmd_hr_vel[st_idx:end_idx,i], "r-")
        plt.plot(data_x, hr_vel[st_idx:end_idx,i], "b-")
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

    ## plot hind left foot velocity
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('hind left foot vel')
    for i in range(3):
        ax1 = plt.subplot(3, 1, i+1)
        plt.plot(data_x, cmd_hl_vel[st_idx:end_idx,i], "r-")
        plt.plot(data_x, hl_vel[st_idx:end_idx,i], "b-")
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


if __name__ == "__main__":
    create_figures()
    plt.show()
