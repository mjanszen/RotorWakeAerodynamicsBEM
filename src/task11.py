# Task 11 - k :
# Operational point of the airfoils


# Idea:
# Plot polars 
# Plot aoa over the span for tsr 6 and 10 
# Discuss where out of range of optimum clcd 
# Discuss when we have most power and see whether there is correlation to the optimum working point 


from BEM import BEM
import numpy as np
import pandas as pd
from helper_functions import Helper
import matplotlib.pyplot as plt
helper = Helper()

####################### Polars

# read in airfoil data

def plot_polar(df):
    """
    Plot the polars of the given dataframe 
    """
    fig,axs = plt.subplots(3,1,figsize=[6,6])
    
    # CL over alpha
    axs[0].plot(df.alpha, df.cl)
    axs[1].plot(df.alpha, df.cd)
    axs[2].plot(df.alpha, df.cl / df.cd)
    
    axs[0].grid()
    axs[0].set_xlabel(r"Angle of attack $\alpha$ [°]")
    axs[0].set_ylabel(r"Lift coefficient $C_l$ [ ]")
    
    axs[1].grid()
    axs[1].set_xlabel(r"Angle of attack $\alpha$ [°]")
    axs[1].set_ylabel(r"Drag coefficient $C_l$ [ ]")

    axs[2].grid()
    axs[2].set_xlabel(r"Angle of attack $\alpha$ [°]")
    axs[2].set_ylabel("Ratio of lift and drag coefficients \n" + r"$C_l / C_d$ [ ]")

    maxid = df.idxmax()
    
    axs[0].scatter(df.alpha[maxid.clcd_ratio], df.cl[maxid.clcd_ratio], color="red")
    axs[1].scatter(df.alpha[maxid.clcd_ratio], df.cd[maxid.clcd_ratio], color="red")
    axs[2].scatter(df.alpha[maxid.clcd_ratio], df.clcd_ratio[maxid.clcd_ratio], color="red")

    plt.show()
    fig.savefig("../results/t11_polars.png",bbox_inches="tight")

def task11():
    polar_file = "../data/polar.xlsx"
    df_polar = pd.read_excel(polar_file,skiprows=3)
    #bem = BEM(data_root="../data",
    #          file_airfoil="polar.xlsx")

    df_polar["clcd_ratio"] = df_polar.cl/df_polar.cd 
    
    plot_polar(df_polar) 
    # get id of the maximum 
     


    breakpoint()


    # AOA over span with horizontal line for optimum cdcl ratio and maximum lift and 

    print("Task 11 - Done")



if __name__ =="__main__":
    task11()
