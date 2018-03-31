import SurvivalModelClasses as Cls
from scipy.stats import binom

Q1 = Cls.Cohort(id=1, pop_size=573, mortality_prob=0.1)
cohortOutcome = Q1.simulate(100)
percent_5_years = Q1.get_perct_patient_survived_beyond_five_years()


print(Q1.get_survival_times())
print('(Question 1) The percentage of patients survived beyond 5 years:', percent_5_years)
print('(Question 2) It would follow binomial distribution. Parameters include k, N, q.')
print('(Question 3)', binom.pmf(400, 573, 0.5, loc=0))


