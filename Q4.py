import CalibrationClasses as CalibClasses
import CalibrationSettings as CalibSets

print('Question 4')
# create a calibration object
calibration = CalibClasses.Calibration()

# sample the posterior of the mortality probability
calibration.sample_posterior()

# Estimate of mortality probability and the posterior interval
print('Estimate of mortality probability ({:.{prec}%} credible interval):'.format(1-CalibSets.ALPHA, prec=0),
      calibration.get_mortality_estimate_credible_interval(CalibSets.ALPHA, 4))





