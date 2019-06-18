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

    file_path = os.getcwd() + "/../test_data/"

    ## read files
    bezier_pos = np.genfromtxt(file_path+'bezier_pos.txt', delimiter=None, dtype=(float))
    bezier_vel = np.genfromtxt(file_path+'bezier_vel.txt', delimiter=None, dtype=(float))
    t = np.genfromtxt(file_path+'bz_time.txt', delimiter='\n', dtype=(float))

    bspline_pos = np.genfromtxt(file_path+'bspline_pos.txt', delimiter=None, dtype=(float))
    bspline_vel = np.genfromtxt(file_path+'bspline_vel.txt', delimiter=None, dtype=(float))
    bspline_acc = np.genfromtxt(file_path+'bspline_acc.txt', delimiter=None, dtype=(float))
    bs_t = np.genfromtxt(file_path+'bs_time.txt', delimiter='\n', dtype=(float))

    numeric_bezier_vel = np.copy(bezier_pos)
    bz_length = len(t)
    for i in range(bz_length-1):
        numeric_bezier_vel[i,:] = (bezier_pos[i+1,:] - bezier_pos[i,:])/(t[i+1]-t[i])
    
    numeric_bspline_vel = np.copy(bspline_pos)
    numeric_bspline_acc = np.copy(bspline_pos)
    bs_length = len(bs_t)
    for i in range(bs_length-1):
        numeric_bspline_vel[i,:] = (bspline_pos[i+1,:] - bspline_pos[i,:])/(bs_t[i+1]-bs_t[i])
        numeric_bspline_acc[i,:] = (bspline_vel[i+1,:] - bspline_vel[i,:])/(bs_t[i+1]-bs_t[i])
 
    ## plot bezier pos
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('bezier_pos')
    for i in range(3):
        ax1 = plt.subplot(3, 1, i+1)
        plt.plot(t, bezier_pos[:,i], "b-")
        plt.grid(True)
    plt.xlabel('time (sec)')
    ## increment figure number and index
    figure_number += 1
    if plot_configuration == PLOT_HORIZONTALLY:
        col_index += 1
    elif plot_configuration == PLOT_VERTICALLY:
        row_index +=1

    ## plot bezier_vel
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('bezier_vel')
    for i in range(3):
        ax1 = plt.subplot(3, 1, i+1)
        plt.plot(t, bezier_vel[:,i], "b-")
        plt.plot(t[:-1], numeric_bezier_vel[:-1,i], "r-")
        plt.grid(True)
    plt.xlabel('time (sec)')
    ## increment figure number and index
    figure_number += 1
    if plot_configuration == PLOT_HORIZONTALLY:
        col_index += 1
    elif plot_configuration == PLOT_VERTICALLY:
        row_index +=1

    ## plot bspline pos
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('bspline_pos')
    for i in range(3):
        ax1 = plt.subplot(3, 1, i+1)
        plt.plot(bs_t, bspline_pos[:,i], "b-")
        plt.grid(True)
    plt.xlabel('time (sec)')
    ## increment figure number and index
    figure_number += 1
    if plot_configuration == PLOT_HORIZONTALLY:
        col_index += 1
    elif plot_configuration == PLOT_VERTICALLY:
        row_index +=1

    ## plot bspline_vel
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('bspline_vel')
    for i in range(3):
        ax1 = plt.subplot(3, 1, i+1)
        plt.plot(bs_t, bspline_vel[:,i], "b-")
        plt.plot(bs_t[:-1], numeric_bspline_vel[:-1,i], "r-")
        plt.grid(True)
    plt.xlabel('time (sec)')
    ## increment figure number and index
    figure_number += 1
    if plot_configuration == PLOT_HORIZONTALLY:
        col_index += 1
    elif plot_configuration == PLOT_VERTICALLY:
        row_index +=1

    ## plot bspline_acc
    fig = plt.figure(figure_number)
    plt.get_current_fig_manager().window.wm_geometry(str(subfigure_width) + "x" + str(subfigure_height) +  "+" + str(subfigure_width*col_index) + "+" + str(subfigure_height*row_index))
    fig.canvas.set_window_title('bspline_acc')
    for i in range(3):
        ax1 = plt.subplot(3, 1, i+1)
        plt.plot(bs_t, bspline_acc[:,i], "b-")
        plt.plot(bs_t[:-1], numeric_bspline_acc[:-1,i], "r-")
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
