import CalibrationClasses2 as CalibClasses2
import CalibrationSettings2 as CalibSets2

print('(Question 6)')

# create a calibration object
calibration = CalibClasses2.Calibration()

# sample the posterior of the mortality probability
calibration.sample_posterior()

# Estimate of mortality probability and the posterior interval
print('Estimate of mortality probability ({:.{prec}%} credible interval):'.format(1-CalibSets2.ALPHA, prec=0),
      calibration.get_mortality_estimate_credible_interval(CalibSets2.ALPHA, 4))

import CalibrationClasses2 as CalibSupport2
import CalibrationSettings2 as CalibSets2

# initialize a calibrated model
calibrated_model = CalibSupport2.CalibratedModel('CalibrationResults2.csv')
# simulate the calibrated model
calibrated_model.simulate(CalibSets2.SIM_POP_SIZE, CalibSets2.TIME_STEPS)

# report mean and projection interval
print('Mean survival time and {:.{prec}%} projection interval:'.format(1 - CalibSets2.ALPHA, prec=0),
      calibrated_model.get_mean_survival_time_proj_interval(CalibSets2.ALPHA, deci=4))
print('The credible interval of the estimated annual mortality probability in Question 6 is narrower than that of Question 4.')
print('The projection interval of the mean survival time in Question 6 is also narrower than that of Question 5.')



