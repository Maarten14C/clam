# clam 2.3.7
* removed links to chrono website, as it has become unstable and could disappear any moment now.
* added an option 'rule' to allow more choice when extrapolating values (see ?approx). This is an internal option and shouldn't normally be important for users.

# clam 2.3.6
* now installs and uses the IntCal package to read the radiocarbon calibration curves (this to have less replication of files between packages that use the IntCal curves)
* repaired bug with Est not being the same size as smp (e.g., whenever dates were truncated or removed, either automatically or manually)
* repaired mixed.effect (however, does not use student.t any more)

# clam 2.3.5
* updated to the latest IntCal20 calibration curves

# clam 2.3.4

* Updated the code to deal with changes in how base-R deals with c() in loops, as suggested by Martin Maechler's e-mail 29 February 2020

# clam 2.3.3

* deptime.age and deptime.depth now invisibly return the deposition time values, for future manipulation such as summary(deptime.age(900))
* Updated help functions of clam and calibrate

# clam 2.3.2

* new function calBP.14C to find the mu (IntCal 14C age) belonging to a single calendar year (suggested by Andres Christen)
* removed suppression of warnings (suggested by Maris Nartiss)
* added sep argument to mix.calibrationcurves (suggested by Thomas Dye)
* calibrate now invisibly returns dat (suggested by Andres Christen) so now can be used as, e.g., cal <- calibrate(130, 30); cal
* related functions were put into separate .R files
* unnecessary references to Bacon were removed
* dummy -5 cal BP lines were removed from the IntCal13 and SHCal13 calibration curves

# clam 2.3.1

* Now a CRAN R package
* Repaired many sundry bugs
* updated the help functions

# clam 2.2

* updated calibration curves (IntCal13 for both hemispheres and the ocean, plus postbomb curves)
* corrected behaviour of est
* corrected use of alternative separator in .csv files (e.g. ";" instead of ","). Also added option to specify alternative decimal points (e.g. "," instead of ".")
* repaired bug that caused incorrect depths for plot.proxies
* repaired bug with slump, causing some dates to be assigned wrong depths
* corrected bug with drawing sample thickness when extradates applied
* new option for est, midpoint of entire calibrated distributions (above threshold, est=7). This differs slightly from the midpoint of the calibrated hpd ranges (est=3)
* minor changes to the manual, incl. more explicit citation suggestions
* added function student.t to show effect alternative calibration

# clam 2.1

* corrected wrong behaviour when using hiatus and slump, e.g., clam(hiatus=470, slump=c(120,140))
* polynomial age-modelling now correctly deals with weights (wght)
* uncertainties for age offsets can be included in calibrate()
* marine offset can now be specified when producing a mixed terrestrial/marine curve (mix.curves())
* new function glue.curves that can be used to extend SHCal04 to 50 kcal BP using IntCal09 and a specified SH offset (obsolete since release of SHCal13)
* corrected references in the manual

#clam 2.0

* Depth segments of abrupt sedimentation can now be excised with the option slump
* different colours for C14 and calendar dates
* accumulation rates are calculated and added as a fifth column to the _ages.txt file
* future ages can be avoided in the dates and the age-depth model
* option to rotate the axes, or reverse the order for the age or depth axes
* optional calibt calibration (Christen and Perez 2009)
* option to provide .txt files with ages and probs for depths
* single instead of mirrored calibrated histograms can be drawn in age-depth graphs

* New function mix.curves
* New functions pMC.age, age.pMC to convert between 14C ages and percent modern carbon
* New function proxy plot
* histograms and confidence ranges for accumulation rates can be calculated for specific depths or ages
* probability distributions are now normalised, so that precise ages have higher heights than imprecise ages. This can be turned off
* more consistent naming of dmax, dmin, etc.
* Warnings are provided in case of clearly wrong settings (e.g. type > 5)
* Internal calculations for C14 dates are now done in F14C (as in OxCal)
* core name added to age-plots
* can now leave thickness column out of .csv file (so, need 6 columns only)
* instead of standard deviation, the statistically more correct probability is used (e.g., what was formerly sdev=2 is now prob=0.95)
* use of file _depths.txt now needs to be activated explicitly
* corrected a bug with calibrating single dates on old or young extremes of the calibration curve
* corrected wrong behaviour when outliers and mixed.effect are used together
* corrected a reservoir effect bug in mixed.effect
* corrected a bug in hpds around 0 BC/AD when BCAD=TRUE
* corrected a wrong reaction to clam(BCAD=TRUE)
* corrected wrong yr.axis titles in png and pdf
* corrected a bug in title option calibrate()
* better error message when type=5 & hiatus
* enhanced treatment of plotting parameters (par)
* better rounding for ages.txt
* squashed many more small to invisible bugs
