import os
import time

import matplotlib.pyplot as plt
import numpy as np


def get_densities() -> set:

    """
    Gathers values for densities used in MD calculations from the filenames
    of '.dat' MD result files.

    Returns
    -------
    set
        A set containing all of the unique densities found in the CWD.
    """

    # Gathering all of the possible result files
    dat_files = [fil for fil in os.listdir() if fil.endswith(".dat")]
    dens_list = []

    for fil in dat_files:

        # This try block is to filter non-compatible files, as incompatible
        # files will raise a ValueError and they can just be ignored.
        try:
            # Getting the index of the hyphen in the filename
            ini = fil.index("-")

            # Getting the index of the dot in the filename
            dot = fil.index(".")

            # Slicing the filename to get the density as a string
            dens = fil[ini + 1 : dot]

            # Getting the middle part of the density string to be able to
            # add a dot
            mid = round(len(dens) / 2)
            dens = float(dens[:mid] + "." + dens[mid:])

            # Saving the formatted density in a list
            dens_list.append(dens)

        # The index method raises an error if the string is not found. We want
        # to ignore files that don't have the desired string.
        except ValueError:
            pass

    # Converting the list into a set to remove repeated densities.
    dens_list = set(dens_list)

    return dens_list


class MDSimulation:
    """
    Object that contains the results of a Molecular Dynamics simulation.
    """

    RDF_NAME = "radial_dist"
    MSD_NAME = "mean_square_displ"
    ENE_NAME = "energies_per_atom"
    F_EXTENS = ".dat"

    def __init__(self, density: float):

        # Making the density available in the entire Class
        self.dens = density

        # Matplotlib stylesheet which includes colorblind-friendly colors.
        plt.style.use("seaborn-colorblind")

        # Making the plot line thinner.
        plt.rc("lines", linewidth=0.1)

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
        self._plot_energies()
        self._plot_rdf()
        # self._plot_msd()

    def _gather_results(self):

        """
        Gathers the results from MD simulations in the CWD and loads their
        values into numpy arrays. These arrays are made available in the
        MDSimulation class with the following objects:

        - self.rdf_results
        - self.energ_results
        - self.msd_results

        """

        # Filename treatment to be able to insert the density into the
        # filename, so results can be worked with iteratively.
        den_str = str(self.dens).replace(".", "")
        rdf = self.RDF_NAME + f"-{den_str}" + self.F_EXTENS
        ene = self.ENE_NAME + f"-{den_str}" + self.F_EXTENS
        msd = self.MSD_NAME + f"-{den_str}" + self.F_EXTENS

        # Loading the different results into an array
        self.rdf_results = np.loadtxt(rdf)
        self.energ_results = np.loadtxt(ene)
        # self.msd_results = np.loadtxt(msd)

    def _result_folder(self):

        """
        Creates a timestamped folder in the CWD to store plots.
        """

        # Getting the current time
        current_time = time.strftime("%d-%m_%H-%M")

        # Timestamped folder name
        self.fold_name = f"results_{current_time}"

        # Creating the folder if the folder does not already exist
        if self.fold_name not in os.listdir():
            os.mkdir(self.fold_name)

    def _plot_rdf(self):

        """
        Plots the radial distribution function using the results gathered
        from the MD simulation.
        The plots are stored in the results_%d-%m_%H-%M directory.
        """

        # Plotting the RDF
        plt.plot(self.rdf_results[:, 0], self.rdf_results[:, 1])
        plt.title("Radial Distribution Function")
        plt.xlabel("r [Å]")
        plt.ylabel("g(r)")

        # Preparing a filename for the RSD plot image.
        den_str = str(self.dens).replace(".", "")
        filename = (
            self.fold_name + "/" + self.RDF_NAME + f"_{den_str}_" + "plot.png"
        )

        # Saving the image and clearing the current plot.
        plt.savefig(filename)
        plt.clf()

    def _plot_msd(self):

        """
        Plots the mean square displacement using the results gathered
        from the MD simulation.
        The plots are stored in the results_%d-%m_%H-%M directory.
        """

        # Plotting the MSD
        plt.plot(self.msd_results[:, 0], self.msd_results[:, 1])
        plt.title("Mean Square Displacement")
        plt.xlabel("r [Å]")
        plt.ylabel("g(r)")

        # Preparing a filename for the MSD plot image.
        den_str = str(self.dens).replace(".", "")
        filename = (
            self.fold_name + "/" + self.MSD_NAME + f"_{den_str}_" + "plot.png"
        )

        # Saving the image and clearing the current plot.
        plt.savefig(filename)
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
        den_str = str(self.dens).replace(".", "")
        n_pot = self.fold_name + "/" + "ene-pot" + f"_{den_str}_" + "plot.png"
        n_kin = self.fold_name + "/" + "ene-kin" + f"_{den_str}_" + "plot.png"
        n_tot = self.fold_name + "/" + "ene-tot" + f"_{den_str}_" + "plot.png"
        n_all = self.fold_name + "/" + "ene-all" + f"_{den_str}_" + "plot.png"

        # Names for the x and y axis
        x_lab = "Time"
        y_lab = "Energy"

        # Preparing the plot for the potential energy
        plt.plot(self.energ_results[:, 0], self.energ_results[:, 1])
        plt.title(f"Potential energy for ρ = {self.dens}")
        plt.xlabel(x_lab)
        plt.ylabel(y_lab)

        # Saving the image and clearing the current plot.
        plt.savefig(n_pot)
        plt.clf()

        # Preparing the plot for the kinetic energy
        plt.plot(self.energ_results[:, 0], self.energ_results[:, 2])
        plt.title(f"Kinetic energy for ρ = {self.dens}")
        plt.xlabel(x_lab)
        plt.ylabel(y_lab)

        # Saving the image and clearing the current plot.
        plt.savefig(n_kin)
        plt.clf()

        # Preparing the plot for the total energy
        plt.plot(self.energ_results[:, 0], self.energ_results[:, 3])
        plt.title(f"Total energy for ρ = {self.dens}")
        plt.xlabel(x_lab)
        plt.ylabel(y_lab)

        # Saving the image and clearing the current plot.
        plt.savefig(n_tot)
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

        plt.title(f"All energies for ρ = {self.dens}")
        plt.legend()
        plt.xlabel(x_lab)
        plt.ylabel(y_lab)

        # Saving the image and clearing the current plot.
        plt.savefig(n_all)
        plt.clf()


if __name__ == "__main__":
    dens_list = get_densities()
    for dens in dens_list:
        sim = MDSimulation(dens)
        sim.generate_plots()
