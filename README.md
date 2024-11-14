# estimators-exploration-project

Stat 624 Project

One of the key aspects within the field of statistics is the estimation of desired quantities. In order to estimate
these quantities, a multitude of estimation methods exist. The goal of this project is to explore different
estimation methods and compare their effectiveness using various metrics. The exploration of different
estimation methods is important because it helps to verify their validity, while the comparison of different
estimation methods showcases the advantages and disadvantages of each relative to the others.

This project will explore and compare three estimators: the maximum likelihood estimate, the method of
moments estimate, and the Bayes estimate under squared error loss. The exploration and comparison of
these three estimators will be two-fold. First, a simulation study will be conducted for each estimator under
various settings to evaluate its effectiveness in capturing the truth. Second, the estimation methods will be
applied to a real dataset in an effort to estimate the parameters of the selected sampling distribution.

The specific dataset analyzed in this project consists of Lionel Messi’s goals scored by season in the domestic
league in which he played. Only completed seasons were considered, resulting in 18 seasons worth of data
scraped from FBref.com. Since the data is count data, the negative binomial distribution was selected for
this project as it is used to model count data and matched the domain of the data. The negative binomial
distribution was used in both the simulation study and the application to the real dataset.

The hope is that the exploration of estimators will be useful for various potential areas. Specifically, in the
case of considering a negative binomial distribution for a dataset, a more informed decision can be made about
which estimator would be most effective. Additionally, when estimating parameters of a distribution, one key
decision is whether to use a Bayesian or a frequentist framework. This project explores both methodologies
so that the reader can be informed about how parameter estimation works in each framework. Finally, the
overall consideration of estimators is important in any type of analysis. The practice of considering different
estimators and choosing the most effective one is critical, and the general methodology presented in this
project could be applied to any analysis that includes the estimation of parameters.
