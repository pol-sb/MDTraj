import os
import sys
import time

try:
    import matplotlib.pyplot as plt
    import numpy as np
except ModuleNotFoundError:
    print(
        "\n[!] ERROR: One or more python module(s) have not been found. Please,"
        " install the following modules and then run again: "
    )
    with open("./src/requirements.txt", "r") as f:
        for line in f:
            print(f"\t{line.strip()}")

    print(
        "Check the README.md file for instructions in how to install the"
        " required python modules.\n"
    )
    sys.exit(1)


class MDSimulation:
    """
    Object that contains the results of a Molecular Dynamics simulation.
    """

    RDF_NAME = "rdf"
    TMP_NAME = "temp"
    PRE_NAME = "pressure"
    ENE_NAME = "energy"
    F_EXTENS = ".dat"
    OUT_PATH = "./output/"

    def __init__(self):

        # TODO: Remove this when I can confirm this is not needed anymore
        # Making the density available in the entire Class
        # self.dens = density

        # Matplotlib stylesheet which includes colorblind-friendly colors.
        plt.style.use("seaborn-colorblind")

        # Making the plot line thinner.
        plt.rc("lines", linewidth=1)

        plt.rc("font", family="serif")
        plt.rcParams.update({"font.size": 14})
        plt.tick_params(labelsize=14)
        plt.rcParams["axes.linewidth"] = 0.5
        #plt.style.use("grayscale")
        plt.grid()

        self._gather_results()

    def generate_plots(self):

        """
        Generate plots using the MD simulation results.
        The plots that can be generated consist of:

        - Kinetic Energy
        - Potential Energy
        - Total Energy
        - All of the energies in the same plot
        - Radial distribution function
        - Mean square displacement

        The plots will be stored in a new folder named results_%d-%m_%H-%M,
        which will be created in the current working directory.

        """

        # Creating the result folder to store the plots
        self._result_folder()

        # Plotting
        self._plot_temperature()
        self._plot_pressure()
        self._plot_energies()
        # self._plot_rdf()
        print(f"Done - Result plots saved in ./output/{self.fold_name}")
        # self._plot_msd()

    def _gather_results(self):

        """
        Gathers the results from MD simulations in the CWD and loads their
        values into numpy arrays. These arrays are made available in the
        MDSimulation class with the following objects:

        - self.temp_results
        - self.energ_results
        - self.press_results
        - self.rdf_results
        - self.msd_results

        """

        # TODO: Remove this later.
        # Filename treatment to be able to insert the density into the
        # filename, so results can be worked with iteratively.
        # den_str = str(self.dens).replace(".", "")
        # rdf = self.RDF_NAME + f"-{den_str}" + self.F_EXTENS
        # ene = self.ENE_NAME + f"-{den_str}" + self.F_EXTENS
        # msd = self.MSD_NAME + f"-{den_str}" + self.F_EXTENS

        temp = self.OUT_PATH + self.TMP_NAME + self.F_EXTENS
        ene = self.OUT_PATH + self.ENE_NAME + self.F_EXTENS
        press = self.OUT_PATH + self.PRE_NAME + self.F_EXTENS
        rdf = self.OUT_PATH + self.RDF_NAME + self.F_EXTENS

        # Loading the different results obtained in the MD simulation
        # into several arrays.
        self.temp_results = np.loadtxt(temp)
        self.energ_results = np.loadtxt(ene)
        self.press_results = np.loadtxt(press)
        #self.rdf_results = np.loadtxt(rdf)
        # self.msd_results = np.loadtxt(msd)

    def _result_folder(self):

        """
        Creates a timestamped folder in the output folder to store the plots.
        """

        # Getting the current time
        current_time = time.strftime("%d-%m_%H-%M")

        # Timestamped folder name
        self.fold_name = f"plots_{current_time}"

        # Creating the folder if the folder does not already exist
        if self.fold_name not in os.listdir(self.OUT_PATH):
            os.mkdir(self.OUT_PATH + self.fold_name)

    def _plot_rdf(self):

        """
        Plots the radial distribution function using the results gathered
        from the MD simulation.
        The plots are stored in the results_%d-%m_%H-%M directory.
        """

        # Plotting the RDF
        plt.plot(self.rdf_results[:, 0], self.rdf_results[:, 1])
        plt.xlim(self.rdf_results[0, 0], self.rdf_results[-2, 0])
        plt.title("Radial Distribution Function")
        plt.xlabel("r [Å]")
        plt.ylabel("g(r)")
        plt.grid()

        # Preparing a filename for the RSD plot image.
        filename = (
            self.OUT_PATH + self.fold_name + "/" + self.RDF_NAME + "plot.png"
        )

        # Saving the image and clearing the current plot.
        plt.savefig(filename, bbox_inches='tight')
        plt.clf()

    def _plot_temperature(self):

        """
        Plots the temperature as a function of time using the results gathered
        from the MD simulation.
        The plots are stored in the ./output/results_%d-%m_%H-%M directory.
        """

        # Plotting the Temperature
        plt.plot(self.temp_results[:, 0], self.temp_results[:, 1])
        plt.xlim(self.temp_results[0, 0], self.temp_results[-1, 0])
        #plt.title("Temperature")
        plt.xlabel("Iteration")
        plt.ylabel(r"$k_BT-1$")

        # Preparing a path for the T plot image.
        filename = (
            self.OUT_PATH + self.fold_name + "/" + self.TMP_NAME + "plot.png"
        )

        # Saving the image and clearing the current plot.
        plt.tight_layout()
        plt.savefig(filename, dpi=400)
        plt.clf()

    def _plot_pressure(self):

        """
        Plots the pressure as a function of time using the results gathered
        from the MD simulation.
        The plots are stored in the ./output/results_%d-%m_%H-%M directory.
        """

        # Plotting the pressure
        plt.plot(self.press_results[:, 0], self.press_results[:, 1])
        plt.xlim(self.press_results[0, 0], self.press_results[-1, 0])
        plt.title("Pressure")
        plt.xlabel("Time [r.u.]")
        plt.ylabel("Pressure [r.u.]")
        plt.grid()

        # Preparing a path for the P plot image.
        filename = (
            self.OUT_PATH + self.fold_name + "/" + self.PRE_NAME + "plot.png"
        )

        # Saving the image and clearing the current plot.
        plt.savefig(filename, bbox_inches='tight')
        plt.clf()

    def _plot_msd(self):

        """
        Plots the mean square displacement using the results gathered
        from the MD simulation.
        The plots are stored in the results_%d-%m_%H-%M directory.
        """

        # Plotting the MSD
        plt.plot(self.msd_results[:, 0], self.msd_results[:, 1])
        plt.xlim(self.msd_results[0, 0], self.msd_results[-1, 0])
        plt.title("Mean Square Displacement")
        plt.xlabel("r [Å]")
        plt.ylabel("g(r)")
        plt.grid()

        # Preparing a filename for the MSD plot image.
        filename = (
            self.OUT_PATH + self.fold_name + "/" + self.MSD_NAME + "plot.png"
        )

        # Saving the image and clearing the current plot.
        plt.savefig(filename, bbox_inches='tight')
        plt.clf()

    def _plot_energies(self):

        """
        Prepares several plots for the energies using the results gathered
        from the MD simulation.
        The plots are stored in the results_%d-%m_%H-%M directory.
        The following plots are generated:

        - Kinetic Energy
        - Potential Energy
        - Total Energy
        - All of the energies in the same plot
        """

        # Preparing filenames for the plot images
        n_pot = self.OUT_PATH + self.fold_name + "/ene-pot" + "plot.png"
        n_kin = self.OUT_PATH + self.fold_name + "/ene-kin" + "plot.png"
        n_tot = self.OUT_PATH + self.fold_name + "/ene-tot" + "plot.png"
        n_all = self.OUT_PATH + self.fold_name + "/ene-all" + "plot.png"

        # Names for the x and y axis
        x_lab = "Time [r.u.]"
        y_lab = r"$\frac{Energy}{Num. Part.}\ \ [r.u.]$"

        # Preparing the plot for the potential energy
        plt.tight_layout()
        plt.plot(self.energ_results[:, 0], self.energ_results[:, 1])
        plt.xlim(self.energ_results[0, 0], self.energ_results[-1, 0])
        plt.title(f"Potential energy")
        plt.xlabel(x_lab)
        plt.ylabel(y_lab)
        plt.grid()

        # Saving the image and clearing the current plot.
        plt.savefig(n_pot, bbox_inches='tight')
        plt.clf()

        # Preparing the plot for the kinetic energy
        plt.plot(self.energ_results[:, 0], self.energ_results[:, 2])
        plt.xlim(self.energ_results[0, 0], self.energ_results[-1, 0])
        plt.title(f"Kinetic energy")
        plt.xlabel(x_lab)
        plt.ylabel(y_lab)
        plt.grid()

        # Saving the image and clearing the current plot.
        plt.savefig(n_kin, bbox_inches='tight')
        plt.clf()

        # Preparing the plot for the total energy
        plt.plot(self.energ_results[:, 0], self.energ_results[:, 3])
        plt.xlim(self.energ_results[0, 0], self.energ_results[-1, 0])
        plt.title(f"Total energy")
        plt.xlabel(x_lab)
        plt.ylabel(y_lab)
        plt.grid()

        # Saving the image and clearing the current plot.
        plt.savefig(n_tot, bbox_inches='tight')
        plt.clf()

        # Plotting all of the energies at once
        plt.plot(
            self.energ_results[:, 0],
            self.energ_results[:, 1],
            label="Potential Energy",
        )
        plt.plot(
            self.energ_results[:, 0],
            self.energ_results[:, 2],
            label="Kinetic Energy",
        )
        plt.plot(
            self.energ_results[:, 0],
            self.energ_results[:, 3],
            label="Total Energy",
        )

        plt.title(f"All energies")
        plt.legend()
        plt.xlabel(x_lab)
        plt.ylabel(y_lab)
        plt.xlim(self.energ_results[0, 0], self.energ_results[-1, 0])
        plt.grid()

        # Saving the image and clearing the current plot.
        plt.savefig(n_all, bbox_inches='tight')
        plt.clf()


# def get_densities() -> set:

#     """
#     Gathers values for densities used in MD calculations from the filenames
#     of '.dat' MD result files.

#     Returns
#     -------
#     set
#         A set containing all of the unique densities found in the CWD.
#     """

#     # Gathering all of the possible result files
#     dat_files = [fil for fil in os.listdir() if fil.endswith(".dat")]
#     dens_list = []

#     for fil in dat_files:

#         # This try block is to filter non-compatible files, as incompatible
#         # files will raise a ValueError and they can just be ignored.
#         try:
#             # Getting the index of the hyphen in the filename
#             ini = fil.index("-")

#             # Getting the index of the dot in the filename
#             dot = fil.index(".")

#             # Slicing the filename to get the density as a string
#             dens = fil[ini + 1 : dot]

#             # Getting the middle part of the density string to be able to
#             # add a dot
#             mid = round(len(dens) / 2)
#             dens = float(dens[:mid] + "." + dens[mid:])

#             # Saving the formatted density in a list
#             dens_list.append(dens)

#         # The index method raises an error if the string is not found. We want
#         # to ignore files that don't have the desired string.
#         except ValueError:
#             pass

#     # Converting the list into a set to remove repeated densities.
#     dens_list = set(dens_list)

#     return dens_list


if __name__ == "__main__":
    sim = MDSimulation()
    sim.generate_plots()

    # TODO: Previous iteration of the code. Delete this later.
    # dens_list = get_densities()
    # for dens in dens_list:
    #    sim = MDSimulation(dens)
    #    sim.generate_plots()
