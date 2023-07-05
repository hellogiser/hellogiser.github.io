---
title: Chaper2 Exploratory data analysis - Geostatistics for Natural Resources Evaluation
date: 2023-04-10 15:12:56 +0800
categories: [Geostatistics]
tags: [Geostatistics]     # TAG names should always be lowercase
---

# Exploratory data analysis

- 2.1 The univarite (one attribute at a time) distributes of categorical and continuous variables are described.
- 2.2 looks at the joint relations between pairs of colocated metal concentrations.
- 2.3 the patterns of variation of metal concentrations are described and related to those of potential sources, such as rock types and land uses.
- 2.4 spatial relations between concentrations of different metals are analyzed.
- 2.5 The main features of the Jura data set are summarized.


## Univariate description

land use and rock type, two categorical variables.

metal concentrations, continuous variables.

### Categorical variables

land uses

- Forest
- Pasture
- Meadow
- Tillage

rock types

- Argovian
- Kimmeridgian
- Sequanian
- Portlandian
- Quaternary

### Continuous variables

#### Frequency distribution

The distribution of continuous values is typically depicted by a histogram with the range of data values discretized into a specific number of classes of equal width and the ralative proportion of data within each class expressed by the height of bars.

#### Cumulative frequency distribution

Critical threshold, the tolerable maximum for healthy soils.

#### Summary statistics

mean, median, minimum, maximum, Std. deviation, Coef. of var., skewness, tolerable max.

A simpler measure of skewness would then be the difference between the mean and median of the distribution, φ' = m - M.

#### Extreme values and data transformation

Extreme values can be handled as follows:

1. Declare the extreme values erroneous and remove them.
2. Classify the extreme values into a separate statistical population.
3. Use robust statistics, which are less sensitive to extreme values.
4. Transform the data to reduce the influence of extreme values.

Such decisions should be made carefully and call for much more than a quick look at the shape of the sample histogram and the disire to make that make histogram symmetric.

The Jura data were validated earlier.

#### The impact of land use and rock type

Split the data into several subsets according to rock type and land use, to better understan the realtion between metal concentrations and environmental factors.

##### Conditional frequencies

The distribution of the z-values, given that a particular state *s<sub>k</sub>* is observed, is said to be conditional to *s<sub>k</sub>*.

##### Conditional cumulative frequencies

Corresponding proportion of data that are above or below the critical threshold.

##### Subdivision of the data set

It would make sense to consider the concentrations for each land use or rock type as a separate population.

## Bivariate description

### The scattergram

This can be displayed in a scattergram in which the components of each data pair are plotted against one another. The figure shows the scattergrams of nickel and zinc values versus the concentrations of Cd, Cu, and Pb.

### Measures of bivariate relation

linear correlation coefficient and rank correlation coefficient

## Univariate spatial description

### Location maps

In general, observations that are close to each other on the ground are also alike in metal concentration.

### The h-scattergram



### Measures of spatial continuity and variability

#### Covariance function

#### Correlogram

#### Semivariogram

#### Example

#### Remarks

### Application to indicator transforms

#### Indicator transform

#### Indicator correlogram

#### Indicator semivariogram

#### Graphical interpretation

#### Example

#### Remarks


### Spatial continuity of metal concentrations

#### Spatial anisotropy

#### Sensitivity to extreme values

#### Interpreting patterns of spatial variation

##### Semivariograms of residuals

#### Indicator semivariograms for metal concentrations

## Bivariate spatial description

### The cross h-scattergram

### Measures of spatial cross continuity/variability

#### Cross covariance function

#### Cross correlogram

#### Pseudo cross semivariogram

#### The lag effect

### The scattergram of h-increments

### Measures of joint variability

#### Cross semivariogram

#### Codispersion function

#### Example

### Application to indicator transforms

#### Indicator cross covariance function

#### Indicator cross correlogram

#### Indicator cross semivariogram

#### Example

#### Remark

### Spatial relations between metal concentrations

#### Spatial anisotropy

#### Indicator cross correlograms

## Main features of the Jura data

1. Many individual samples are contaminated with cadmium or lead, whereas a smaller proportion of samples exceeds of the tolerable maximum for copper.
2. The distribution of Cd, Cu, and Pb concentrations are positively skewed.
3. The smallest metal concentrations are measured in forest soil or on Argovian rocks, whereas soil under pasture has the largest concentrations for all metals.
4. The metals with widespread contamination are positively are related to the better sampled zinc. There is a positive relation between nickel and cadmium concentrations.
5. A small nuggest effect, a short scale (range ≈ 200 m), and a regional scale (range ≈ 1 km) of spatial variability are observed on the semivariograms of metal concentrations. The short-range structure is the major component for the three metals with widespread contamination (Cd, Cu, and Pb) and Cr. The long-range structure dominates the semivariograms of Ni and Co concentrations. The Zn semivariogram combines the two structures in approximately equal proportions.
6. The short-range structure relates to the spatial distruction of rock types and land uses in the study area. The long-range structure reflects the influence of Argovian and Kimmeridgian rock types on metal concentrations.
7. Nickel concentrations vary more continuously in the SW-NE direction, which corresponds to the preferential orientation of the underlying geologic formations. The patterns of variation of other metals are fairly similar in all directions (isotropy).
8. Small concentrations in Cd, Ni, and Zn are better connected in space than larger concentrations. This suggests the existence of homogeneous areas of small concentrations and larger zones where high and median concentrations are intermingled.
9. The metals with widespread contamination (Cd, Cu, and Pb) show a short-range cross dependence with the better sampled Zn. There is a long-range cross dependence between Cd and Ni concentrations.
