import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# %% Half Life Indium Decay
filename = "C:/Users/lando/PycharmProjects/PHYS-270L/RadioactivityLab-PHYS270/indiumData.txt"
time, counts = np.genfromtxt(filename, delimiter='\t', usecols=(0, 1), skip_header=7, unpack=True)


def expfit(x, a, b):
    return a * np.exp(-b * x)


indiumfit, indiumfitunc = curve_fit(expfit, time, counts, p0=(300, 0.693 / 60))
t12 = np.log(2) / indiumfit[1]

lambdaerr = np.sqrt(np.diag(indiumfitunc)[1])
t12err = t12 * lambdaerr / indiumfit[1]
t12lit = 54  # find source for lab report
t12percenterror = abs(t12 - t12lit) / t12lit * 100

print(f'The measured half life of Indium-116m is {t12:.3f} +/- {t12err:.3f} minutes.')
print(f'The percent error of this measurement is {t12percenterror:.03}%')

# *indiumfit = indiumfit[0], indiumfit[1]
plt.figure(figsize=(7, 5), dpi=300)
plt.plot(time, counts, 'o', color='blue', label='Measured Data')
plt.plot(time, expfit(time, *indiumfit), '-', color='red', label='Fitted Data')
plt.xlabel('Time (Min)')
plt.ylabel('Activity (Counts)')
plt.title('Indium Decay')
plt.show()

# %% Shielding
shieldingfilename = "C:/Users/lando/PycharmProjects/PHYS-270L/RadioactivityLab-PHYS270/cesiumData.txt"
thickness, shieldingcounts = np.genfromtxt(shieldingfilename, usecols=(0, 1), skip_header=1, unpack=True)

# sort in order of ascending x-y pairs
thickness, shieldingcounts = zip(*sorted(zip(thickness, shieldingcounts)))
thickness = np.asarray(thickness)
shieldingcounts = np.asarray(shieldingcounts)

shieldingfit, shieldingfitunc = curve_fit(expfit, thickness, shieldingcounts, p0=(600, .693 / 4000))

thickness12 = np.log(2) / indiumfit[1]

shieldinglambdaerr = np.sqrt(np.diag(indiumfitunc)[1])
thickness12err = thickness12 * shieldinglambdaerr / shieldingfit[1]

print(f'The measured half density thickness of Cesium-137 is {thickness12:.3f} +/- {thickness12err:.3f} mg/cm.')

# *indiumfit = indiumfit[0], indiumfit[1]
plt.figure(figsize=(7, 5), dpi=300)
plt.plot(thickness, shieldingcounts, 'o', color='blue', label='Measured Data')
plt.plot(thickness, expfit(thickness, *shieldingfit), '-', color='red', label='Fitted Data')
plt.xlabel('Density Thickness (mg/cm^2)')
plt.ylabel('Activity (Counts)')
plt.legend()
plt.title('Ceisum-137 Shielding')
plt.show()

# %% Distance
shieldingfilename = "C:/Users/lando/PycharmProjects/PHYS-270L/RadioactivityLab-PHYS270/cesiumData.txt"
distance, distancecounts = np.genfromtxt(shieldingfilename, usecols=(2, 3), skip_header=1, max_rows=(9), unpack=True)

expdistancefit, expdistancefitunc = curve_fit(expfit, distance, distancecounts)

def r2fit(x, a, b):
    return a/x**2 + b

r2distancefit, r2distancefitunc = curve_fit(r2fit, distance, distancecounts) # error present in this line

def rfit(x, a, b):
    return a/x + b

rdistancefit, rdistancefitunc = curve_fit(rfit, distance, distancecounts)

plt.figure(figsize=(7, 5), dpi=300)
plt.plot(distance, distancecounts, 'o', color='blue', label='Measured Data')
plt.plot(distance, expfit(distance, *expdistancefit), linestyle='solid', color='red', label='Exponential Fit')
plt.plot(distance, r2fit(distance, *r2distancefit), linestyle='dotted', color='green', label='1/r^2 Fit')
plt.plot(distance, rfit(distance, *rdistancefit), linestyle='dashed', color='yellow', label='1/r Fit')
plt.xlabel('Distance from Counter (cm)')
plt.ylabel('Activity (Counts)')
plt.legend()
plt.title('Distance and Counts of Cesium-137')
plt.show()
