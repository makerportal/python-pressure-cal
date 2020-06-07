####################################
# Calibrating Pressure Transducer
# with Manual Manometer
####################################
#
#
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

plt.style.use('ggplot')

######################################################
# Data and Conversion to Pressure/Voltage
######################################################
#
rho_2 = 997.0 # density of water [kg/m^3]
rho_1 = 1.225 # density of air [kg/m^3]
g   = 9.81 # gravity [m/s^2]
h_vals   = np.divide([13.7,23.2,43.25,79.1,97.1,145.0,111.8,126.0,-17.4,-53.8,-95.2,-37.3,
                     -73.6,-115.2,-146.2,-132.6],1000.0) # height of manometer under different pressures [m]
MPS_ADC_raw_vals = [9830572.0,9887678.0,10054668.0,10292554.0,10405446.0,
                                 10726234.0,10497706.0,10608056.0,9601672.0,9351684.0,
                                 9023918.0,9466900.0,9153698.0,8907830.0,8688792.0,8784272.0] # raw ADC vals
MPS_vals = 1000.0*(5.0*np.divide(MPS_ADC_raw_vals,np.power(2.0,24.0)))/128.0 # response of 24-bit MPS20N0040D in [mV]
P_vals   = (rho_2-rho_1)*g*h_vals/1000.0 # pressure approx in [kPa]
#
######################################################
# Calibration Analysis w/statistics
######################################################
#
slope, intercept, r_value, p_value, std_err = stats.linregress(MPS_vals,P_vals)
P_predict = (slope*MPS_vals) + intercept
rmse = np.sqrt(np.mean(np.power(np.subtract(P_predict,P_vals),2.0))) # root-mean square error
mape = np.mean(np.abs(100.0*np.subtract(P_predict,P_vals)/P_vals)) # mean absolute percent error
R_sq = r_value**2 # coefficient of determination
mae  = np.mean(np.abs(np.subtract(P_predict,P_vals))) # mean absolute error
bias = np.mean(np.subtract(P_predict,P_vals))

MPS_test = np.linspace(0.0,50.0,1000) # for drawing the pressure line as a func of [mV]
ADC_test_vals = 128.0*(MPS_test/1000.0) # for finding the max pressure
P_line = slope*MPS_test+intercept

######################################################
# Calibration Results
######################################################
#
fig,ax = plt.subplots(figsize=(14,8))
l1, = plt.plot(MPS_test,P_line,color='k',linewidth=3)
l2, = ax.plot(MPS_vals,P_vals,color=plt.cm.Set1(1),linestyle='',marker='o',markersize=10,zorder=99)

ax.set_xlabel('MPS20N0040D Voltage, $V_R$ [mV]',fontsize=16)
ax.set_ylabel('Manometer Pressure, $P$ [kPa]',fontsize=16)
ax.text(41.25, -9.0, 'Fit Statistics:\n$R^2$      = '+'{0:2.2f}\nRMSE = {1:2.2f} kPa\n'.format(R_sq,rmse)+\
                    'MAPE = {0:2.2f}%\nMAE   = {1:2.2f} kPa\nBias   = {2:2.2f} kPa'.format(mape,mae,bias),
            size=16,ha="left", va="center",bbox=dict(boxstyle="round",
                       ec=(0.9, 0.9, 0.9),fc=(1.0, 1.0, 1.0),pad = 0.75))
#
######################################################
# Theory Comparison
######################################################
#
# directly from datasheet
A_F = 128.0 # amplification from HX710B
S   = 0.05/40.0 # sensitivity [V/kPa]
V_R = np.linspace(0.0,0.05,1000) # test voltages based on +25mV bias and 50mV full-scale claim
b   = 0.025 # DC offset claim from datasheet
P_theory = (V_R/S - b/S) # implementation of sensitivity and offset

# adjust datasheet params
S_adj   = (1.0/slope)/1000.0 # sensitivity [V/kPa]
b_adj   = (intercept/slope)/1000.0 # DC offset approx. from data
P_theory_adj = (V_R/S_adj + b_adj/S_adj) # implementation of sensitivity and offset

t1, = ax.plot(V_R*1000.0,P_theory,linewidth=4,color=plt.cm.Set1(2)) # raw theory
t2, = ax.plot(V_R*1000.0,P_theory_adj,linewidth=4,color=plt.cm.Set1(3),linestyle='--') # adjusted theory

# vertical lines for marking the ADC cutoff values (identifying the effective pressure range of the sensor)
b_l = ax.vlines(x=(5.0/128.0)*1000.0, ymin=np.min(P_theory), ymax=np.max(P_theory),
                linestyle='dotted',color=plt.cm.Set1(4),linewidth=4)
b_r = ax.vlines(x=0.0, ymin=np.min(P_theory), ymax=np.max(P_theory),linestyle='dotted',
                color=plt.cm.Set1(4),linewidth=4)

# legend marking each line on the plot
ax.legend([l2,l1,b_l,t1,t2],['Data','Fit: $P$ = {0:2.2f}'.format(slope)+'$V_R$ '+'{0:2.2f}'.format(intercept),
                       '0V-5V ADC Bounds','Theory Prediction','Theory Adjusted'],fontsize=16,bbox_to_anchor=(0,0,0.45,0.925))
plt.show()
