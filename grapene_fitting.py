
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

class RamanData:

    def __init__(self, path):

        self.path = path

        data_set = self.__get_data(path)
        self.raman_data = self.__data_to_lists(data_set)
    
    def __get_data(self, path):

        with open(path) as f:
            file_data = f.readlines()

        return file_data

    def __data_to_lists(self, data_set):

        data_lists = []

        for data_set in data_set:
            data_list = data_set.replace("\n", "").split("\t")[2:]
            data_list = [float(i) for i in data_list]
            data_lists.append(data_list)

        return data_lists
    
    def plot_single_raman(self, single_no):

        plt.plot(self.raman_data[0], self.raman_data[single_no])
        plt.xlabel("Raman shift / cm$^{-1}$")
        plt.ylabel("Intensity")
        plt.show()

    def plot_all_raman(self):

        for raman in self.raman_data[1:]:
            plt.plot(self.raman_data[0], raman)
            plt.xlabel("Raman shift / cm$^{-1}$")
            plt.ylabel("Intensity")
            plt.show()


class RamanAnalysis(RamanData):

    def __init__(self, path):
        super().__init__(path)

        self.cleaned_data = self.__clean_data(self.raman_data)

        self.peaks_all = []
        self.peaks_D = []
        self.peaks_G = []
        self.peaks_2D = []

        self.ratio_D_G = []
        self.ratio_2D_G = []

    def __clean_data(self, data_set):

        cleaned_set = []
        for raman in range(1, len(data_set)):
            raman_set = data_set[raman]
            for point in range(len(raman_set) - 1):
                if raman_set[point]/raman_set[point + 1] > 5:
                    pass

    def __lorentzian(self, x, a, b, mu, sigma):
        return a + (b*sigma)/(np.pi*((x-mu)**2 + sigma**2))
    
    def single_peak_fit(self, iteration, low_freq, high_freq, guess):

        intensity_array = np.array(self.raman_data[0])

        diff_array_low = np.absolute(intensity_array - low_freq)
        low_bound = diff_array_low.argmin()
        
        diff_array_high = np.absolute(intensity_array - high_freq)
        high_bound = diff_array_high.argmin()

        popt, pcov = curve_fit(
            self.__lorentzian,
            self.raman_data[0][low_bound:high_bound],
            self.raman_data[iteration][low_bound:high_bound],
            p0 = guess
        )

        return popt, pcov
    
    def fit_three_peaks(self, instance, 
                        low_freq_D=1250, high_freq_D=1450, guess_D=[30, 1500, 1300, 100],
                        low_freq_G=1550, high_freq_G=1650, guess_G=[30, 1500, 1600, 10],
                        low_freq_2D=2600, high_freq_2D=2800, guess_2D=[30, 1500, 2700, 10]):
        
        peak_D = self.single_peak_fit(instance, low_freq_D, high_freq_D, guess_D)
        peak_G = self.single_peak_fit(instance, low_freq_G, high_freq_G, guess_G)
        peak_2D = self.single_peak_fit(instance, low_freq_2D, high_freq_2D, guess_2D)

        three_peak_results = [peak_D, peak_G, peak_2D]

        return three_peak_results
    
    def fit_all_peaks(self,
                      low_freq_D=1250, high_freq_D=1450, guess_D=[30, 1500, 1300, 100],
                      low_freq_G=1550, high_freq_G=1650, guess_G=[30, 1500, 1600, 10],
                      low_freq_2D=2600, high_freq_2D=2800, guess_2D=[30, 1500, 2700, 10]):
        
        peak_results = []
        for raman in range(1, len(self.raman_data)):
            three_peaks = self.fit_three_peaks(
                                            raman,
                                            low_freq_D, high_freq_D, guess_D,
                                            low_freq_G, high_freq_G, guess_G,
                                            low_freq_2D, high_freq_2D, guess_2D
                                        )
            peak_results.append(three_peaks)

        self.peaks_all = peak_results
        self.peaks_D = [peak[0] for peak in peak_results]
        self.peaks_G = [peak[1] for peak in peak_results]
        self.peaks_2D = [peak[2] for peak in peak_results]

        return peak_results
    
    def peak_ratio(self, peak_1, peak_2,  plot = False):

        if not self.peaks_all:
            self.fit_all_peaks()

        ratio = []
        for raman in range(len(self.peaks_all)):

            if peak_1 == "D":
                parameters_1 = self.peaks_D[raman][0]
            elif peak_1 == "G":
                parameters_1 = self.peaks_G[raman][0]
            elif peak_1 == "2D":
                parameters_1 = self.peaks_2D[raman][0]

            if peak_2 == "D":
                parameters_2 = self.peaks_D[raman][0]
            elif peak_2 == "G":
                parameters_2 = self.peaks_G[raman][0]
            elif peak_2 == "2D":
                parameters_2 = self.peaks_2D[raman][0]

            intensity_1 = self.__lorentzian(
                parameters_1[2],
                parameters_1[0],
                parameters_1[1],
                parameters_1[2],
                parameters_1[3],
            )

            intensity_2 = self.__lorentzian(
                parameters_2[2],
                parameters_2[0],
                parameters_2[1],
                parameters_2[2],
                parameters_2[3],
            )

            ratio.append(intensity_1/intensity_2)

        if peak_1 == "D" and peak_2 == "G":
            self.ratio_D_G = ratio
        elif peak_1 == "2D" and peak_2 == "G":
            self.ratio_2D_G = ratio

        if plot:
            plt.hist(ratio, bins=30)
            plt.xlabel("Ratio")
            plt.ylabel("Count")
            plt.title(f"Ratio of {peak_1} to {peak_2}")
            plt.show()

        return ratio

    def ratio_D_to_G(self, plot = False):

        ratio = self.peak_ratio("D", "G", plot)

        return ratio
    
    def ratio_2D_to_G(self, plot = False):

        ratio = self.peak_ratio("2D", "G", plot)

        return ratio
    
    def plot_single_fit(self, peak, instance, width = 200):

        if not self.peaks_all:
            self.fit_all_peaks()

        if peak == "D":
            parameters = self.peaks_D[instance][0]
        elif peak == "G":
            parameters = self.peaks_G[instance][0]
        elif peak == "2D":
            parameters = self.peaks_2D[instance][0]

        low_freq = parameters[2] - width/2
        high_freq = parameters[2] + width/2

        intensity_array = np.array(self.raman_data[0])

        diff_array_low = np.absolute(intensity_array - low_freq)
        low_bound = diff_array_low.argmin()
        
        diff_array_high = np.absolute(intensity_array - high_freq)
        high_bound = diff_array_high.argmin()

        fit_data = self.__lorentzian(
            self.raman_data[0][low_bound:high_bound],
            parameters[0],
            parameters[1],
            parameters[2],
            parameters[3]
        )

        plt.plot(self.raman_data[0][low_bound:high_bound], self.raman_data[instance + 1][low_bound:high_bound])
        plt.plot(self.raman_data[0][low_bound:high_bound], fit_data)
        plt.xlabel("Raman shirt / cm$^{-1}$")
        plt.ylabel("Intensity")
        plt.title(f"{peak} peak Raman fit number {instance + 1}")
        plt.show()

    def plot_all_fits(self):

        if not self.peaks_all:
            self.fit_all_peaks()

        for raman in range(len(self.peaks_all)):
            self.plot_single_fit("D", raman)
            self.plot_single_fit("G", raman)
            self.plot_single_fit("2D", raman)

    
    def plot_D_mu(self):

        if not self.peaks_D:
            self.fit_all_peaks()

        mu_D = [peak[0][2] for peak in self.peaks_D]

        plt.hist(mu_D)
        plt.title("Count of D peak mean")
        plt.xlabel("Peak mean / cm$^{-1}$")
        plt.ylabel("Count")
        plt.show()

    def plot_G_mu(self):

        if not self.peaks_G:
            self.fit_all_peaks()

        mu_G = [peak[0][2] for peak in self.peaks_G]

        plt.hist(mu_G)
        plt.title("Count of G peak mean")
        plt.xlabel("Peak mean / cm$^{-1}$")
        plt.ylabel("Count")
        plt.show()

    def plot_2D_mu(self):

        if not self.peaks_2D:
            self.fit_all_peaks()

        mu_2D = [peak[0][2] for peak in self.peaks_2D]

        plt.hist(mu_2D)
        plt.title("Count of 2D peak mean")
        plt.xlabel("Peak mean / cm$^{-1}$")
        plt.ylabel("Count")
        plt.show()



test_1 = RamanData(r'C:\Users\mr5014\Documents\PhD\VNA Code\Raman\Matt\graphene_03.txt')
#test_1.plot_single_raman(23)

test = RamanAnalysis(r'C:\Users\mr5014\Documents\PhD\VNA Code\Raman\Matt\graphene_03.txt')
#result = test.fit_all_peaks()
#print(result[0])
#test.plot_D_mu()
#test.plot_G_mu()
#test.plot_2D_mu()
#test.plot_all_fits()

result_D = test.ratio_D_to_G(True)
print(result_D)
result_2D = test.ratio_2D_to_G(plot=True)
print(result_2D)

plt.hist([1.0740417533676154, 1.0142977907534785, 1.2021336944179075, 1.0616699076822835, 1.1220236567384883, 1.0960294081154351, 1.0203233454543794, 0.9527631925722418, 1.0852143062066526, 1.0739953675031368, 1.0127894144485532, 1.1221366990358161, 0.9979706818918804, 0.8861951483264627, 0.9112539608576421, 1.0221152591709335, 1.0100600018064867, 1.0382101049681882, 1.0441043322353558, 0.9891394841904894, 1.004276891920068, 1.1401236615761867, 0.9283584917162361, 1.1403288865876395], bins=30)
plt.xlabel("Ratio")
plt.ylabel("Count")
plt.title(f"Ratio of D to G")
plt.show()

#print(test.raman_data)
#test.plot_all_raman()
