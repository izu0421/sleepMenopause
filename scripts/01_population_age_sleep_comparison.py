# ============================================================================
# 01 - Population age-stratified sleep comparison
# ----------------------------------------------------------------------------
# Project : sleepMenopause - androgen decline and perimenopausal sleep
# Purpose : Compare sleep-related outcomes between women <40 and >=40 years to
#           establish the population baseline of age-related sleep degradation
#           (manuscript Figure 1).
# Methods : Mann-Whitney U and Welch's t-tests across sleep domains; ordinal
#           models (statsmodels OrderedModel) for stratified outcomes.
# Input   : all_dt_target_variables_curated_added_tests.csv
#           (NOT included - see "Data availability" in README.md)
# Output  : Group-comparison statistics for sleep variables in women.
# Depends : pandas, numpy, scipy, statsmodels
# Note    : Raw participant data are not distributed with this repository.
# ============================================================================

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu, ttest_ind
from statsmodels.miscmodels.ordinal_model import OrderedModel

# 1. Load the dataset
# Replace with your actual file path
file_name = 'all_dt_target_variables_curated_added_tests.csv'
df = pd.read_csv(file_name)

# 2. Subset for women (Gender == 0)
# Gender=0 identified as female based on lower neck size and snoring prevalence
df_women = df[df['Gender'] == 0].copy()

# 3. Define sleep-related variables for analysis
sleep_vars = {
    'In.the.last.2.weeks..please.select.how.many.times.it.has.been.difficult.for.you.to.fall.asleep': 'Sleep Latency',
    'In.the.last.2.weeks.please.select.how.difficult.it.has.been.for.you.to.stay.asleep': 'Sleep Maintenance',
    'In.the.last.2.weeks.please.select.how.many.times.you.have.woken.too.early': 'Early Waking',
    'How.SATISFIED.DISSATISFIED.are.you.with.your.CURRENT.sleep.pattern.': 'Sleep Dissatisfaction',
    'To.what.extent.do.your.sleep.problems.INTERFERE.with.your.daily.functioning.e.g..daytime.fatigue..mood..ability.to.function.at.work.daily.chores..concentration..memory..mood..etc.': 'Daytime Interference',
    'How.do.you.feel.during.the.day..Thinking.over.the.past.week..how.tired.have.you.felt.during.the.day..1...extremely.alert..full.of.energy..9...very.sleepy..great.effort.to.stay.awake..fighting.sleep.': 'Daytime Sleepiness',
    'How.NOTICEABLE.to.others.do.you.think.your.sleep.problem.is.in.terms.of.impairing.the.quality.of.your.life.': 'Problem Noticeability',
    'How.WORRIED.DISTRESSED.are.you.about.your.current.sleep.problem.': 'Worry/Distress',
    'ISI': 'ISI Score',
    'MAPI': 'MAPI Score',
    'STOPBANG': 'STOP-BANG Score'
}

# 4. Categorize by Age (Above 40 vs Below 40)
df_women['Age_Group'] = df_women['Age'].apply(lambda x: 'Above 40' if x >= 40 else 'Below 40')

# 5. Perform Group Comparisons
comparison_results = []

for original_col, clean_name in sleep_vars.items():
    if original_col in df_women.columns:
        group_above = df_women[df_women['Age_Group'] == 'Above 40'][original_col].dropna()
        group_below = df_women[df_women['Age_Group'] == 'Below 40'][original_col].dropna()
        
        if len(group_above) > 5 and len(group_below) > 5:
            # Mann-Whitney U Test (Non-parametric)
            u_stat, p_mw = mannwhitneyu(group_above, group_below, alternative='two-sided')
            
            # Welch's T-test (Parametric)
            t_stat, p_ttest = ttest_ind(group_above, group_below, equal_var=False)
            
            comparison_results.append({
                'Variable': clean_name,
                'Mean (>= 40)': round(group_above.mean(), 2),
                'Mean (< 40)': round(group_below.mean(), 2),
                'Difference': round(group_above.mean() - group_below.mean(), 2),
                'P (Mann-Whitney)': p_mw,
                'P (Welch T-test)': p_ttest
            })

# Summarize and Sort by Significance
comparison_df = pd.DataFrame(comparison_results).sort_values(by='P (Mann-Whitney)')
print("--- Comparison of Sleep Features Between Age Groups (Women) ---")
print(comparison_df.to_string(index=False))

# 6. Ordinal Logistic Regression for Sleep Latency with BMI as Covariate
# This model assesses if age impacts sleep onset after accounting for body mass
latency_col = 'In.the.last.2.weeks..please.select.how.many.times.it.has.been.difficult.for.you.to.fall.asleep'
regression_data = df_women[[latency_col, 'Age', 'BMI']].dropna()

# Ensure target is formatted as ordered integers (0, 1, 2, 3, 4)
regression_data['Latency_Ordinal'] = regression_data[latency_col].astype(int)
regression_data['Age_Above_40'] = (regression_data['Age'] >= 40).astype(int)

# Fit the Ordered Logit Model
mod = OrderedModel(regression_data['Latency_Ordinal'], 
                   regression_data[['Age_Above_40', 'BMI']], 
                   distr='logit')
res = mod.fit(method='bfgs', disp=False)

print("\n--- Ordinal Regression Results (Latency ~ Age >= 40 + BMI) ---")
print(res.summary())

# 7. Save results for documentation
comparison_df.to_csv('age_comparison_women_results.csv', index=False)
with open('ordinal_regression_latency_summary.txt', 'w') as f:
    f.write(res.summary().as_text())