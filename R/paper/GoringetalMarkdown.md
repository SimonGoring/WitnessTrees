---
author: 'Simon Goring *et al*.'
biblio-files: 'goringetal_references'
bibliography: 'goringetal_references.bib'
csl: 'ecology.csl'
date: '24 June, 2015'
output:
  md_document:
    variant: markdown
  pdf_document:
    pandoc_args: '-V geometry:vmargin=1in -V geometry:hmargin=1in'
title: Witness Tree paper
...

Changes in Forest Composition, Stem Density, and Biomass from the Settlement Era (1800s) to Present in the Upper Midwestern United States
=========================================================================================================================================

<!--


-->
Simon J. Goring^1^

David J. Mladenoff^2^

Charles V. Cogbill^3^

Sydne Record^3,4^

Christopher J. Paciorek^5^

Stephen T. Jackson^6^

Michael C. Dietze^7^

Andria Dawson ^5^

Jaclyn Hatala Matthes^8^

Jason S. McLachlan^9^

John W. Williams^1,10^

^1^Department of Geography, University of Wisconsin, Madison, 550 N Park
St, Madison WI 53706

^2^Department of Forest and Wildlife Ecology, University of
Wisconsin-Madison, 1630 Linden Dr, Madison WI 53706

^3^Harvard Forest, Harvard University, 324 N Main St, Petersham MA 01366

^4^Department of Biology, Bryn Mawr College, 101 North Merion Ave., Bryn
Mawr PA 19010

^5^Department of Statistics, University of California, Berkeley, 367
Evans Hall, Berkeley CA 94720

^6^Department of the Interior Southwest Climate Science Center, U.S.
Geological Survey, 1955 E. Sixth St. Tucson, AZ 85719; School of Natural
Resources and the Environment and Department of Geosciences, University
of Arizona, Tucson AZ 85721

^7^Department of Earth and Environment, Boston University, 685
Commonwealth Ave, Boston, MA 02215

^8^Department of Geography, Dartmouth College, 6017 Fairchild, Hanover,
NH 03755

^9^Department of Biological Sciences, University of Notre Dame, 100
Galvin Life Sciences Center, Notre Dame, IN 46556

^10^Center for Climatic Research, University of Wisconsin, Madison, 1225
W Dayton St, Madison WI 53706

------------------------------------------------------------------------

Abstract
--------

*EuroAmerican land use and its legacies have transformed forest
structure and composition across the United States (US). More accurate
reconstructions of historical states are critical to understanding the
processes governing past, current, and future forest dynamics. Gridded
(8x8km) estimates of pre-settlement (1800s) forests from the upper
Midwestern US (Minnesota, Wisconsin, and most of Michigan) using 19th
Century Public Land Survey (PLS) records provide relative composition,
biomass, stem density, and basal area for 26 tree genera. This mapping
is more robust than past efforts, using spatially varying correction
factors to accommodate sampling design, azimuthal censoring, and biases
in tree selection.*

*We compare pre-settlement to modern forests using Forest Inventory and
Analysis (FIA) data, with respect to structural changes and the
prevalence of lost forests, pre-settlement forests with no current
analogue, and novel forests, modern forests with no past analogs. Stem
density, basal area and biomass are higher in contemporary forests than
in settlement-era forests, but this pattern is spatially structured.
Modern biomass is higher than pre-settlement biomass in the northwest
(Minnesota and northern Wisconsin), and lower in the east, due to shifts
in species composition and, presumably, average stand age. Modern
forests are more homogeneous, and ecotonal gradients are more diffuse
today than in the past. Novel forest represent 29% of all FIA cells,
while 25% of pre-settlement forests no longer exist in a modern
context.*

*Lost forests are centered around the forests of the Tension Zone,
particularly in hemlock dominated forests of north-central Wisconsin,
and in oak-elm-basswood forests along the forest-prairie boundary in
south central Minnesota and eastern Wisconsin. Novel FIA forest
assemblages are distributed evenly across the region, and the zones
themselves are not spatially structured, representing a broad-scale
homogenization of forest composition and structure that is strongly
influenced by current and historic land use. All datasets and
open-source analytical scripts are publicly available, providing a
benchmark for future analyses, and work is underway to extend these
settlement-era reconstructions across the northeastern US forests.*

**Key Words**: euroamerican settlement, land use change, public land
survey, historical ecology, novel ecosystems, biomass, forest inventory
and analysis, ecotone, forest ecology

Introduction:
-------------

The composition, demography, and structure of forests in eastern North
America have changed continuously over the last millennium, driven by
human land use
[@ramankutty1999estimating; @thompson2013four; @munoz2014defining; @ellis2008putting; @foster1998land]
and climate variability
[@pederson2014legacy; @booth2012multi; @hotchkiss2007response; @umbanhowar2006asymmetric].
While human effects have been a component of these systems for millenia,
the EuroAmerican settlement and industrialization period have increased
anthropogenic effects by orders of magnitude
[@brugam1978pollen; @fuller1998impact; @mcandrews1988human]. Legacies of
post-settlement land use in the upper Midwest [@grossmann2008farms] and
elsewhere have been shown to persist at local and regional scales
[@foster1998land; @dupouey2002irreversible; @etienne2013searching], and
nearly all North American forests have been affected by the
intensification of land use in the past three centuries. Hence,
contemporary ecological processes in North American forests integrate
the anthropogenic impacts of the post-EuroAmerican period and natural
influences at decadal to centennial scales. For example, the natural
processes of secession, senescense and replacement of tree species in
forests may be masked, or heavily modified by historically recent land
use change.

At a regional scale many forests in the upper Midwest (*i.e.*,
Minnesota, Wisconsin and Michigan) now have decreased species richness
and functional diversity relative to forests of the pre-EuroAmerican
settlement (hereafter pre-settlement) period
[@schulte2007homogenization; @hanberry2012comparison; @li2014drivers]
due to near complete logging. For example, forests in Wisconsin are in a
state of regrowth, with an unfilled carbon sequestration potential of 69
TgC [@rhemtulla2009historical] as a consequence of these extensive land
cover conversions and subsequent partial recovery. The upper Midwestern
United States represents a unique ecological setting, with multiple
major ecotones, including the prairie-forest boundary, historic savanna,
and the Tension Zone between southern deciduous forests and northern
evergreen forests. The extent to which these ecotones have shifted, and
their extent both prior to and following EuroAmerican settlement is of
critical importance to biogeochemical and biogeophysical
vegetation-atmosphere feedbacks [@matthes_inprep], carbon sequestration
[@rhemtulla2009historical], and regional management and conservation
policy
[@radeloff2000historical; @fritschle2008reconstructing; @knoot2010state; @gimmi2013assessing].

Land use change at the local and state-level has affected both the
structure and composition of forests in the Midwestern United States
[*e.g.* @schulte2007homogenization; @hanberry2012comparison].
Homogenization and shifts in overall forest composition are evident, but
the spatial extent and structure of this effect is less well understood.
Studies in Wisconsin have shown differential patterns of change in the
mixedwood and evergreen dominated north versus the southern driftless
and hardwood south. Does this pattern of differential change extend to
Minnesota and Michigan? To what extent are land-use effects common
across the region, and where are responses ecozone-specific? Has
homogenization [*e.g.*, @schulte2007homogenization] resulted in novel
forest assemblages relative to pre-settlement baselines across the
region, and the loss of pre-settlement forest types? Are the spatial
distributions of these novel and lost forest types overlapping, or do
they have non-overlapping extents? If broad-scale reorganization is the
norm following EuroAmerican settlement, then the ecosystems that we have
been studying for the past century may indeed be novel relative to the
reference conditions of the pre-settlement era.

Modern forest structure and composition data [*e.g.*, from the United
States Department of Agriculture Forest Service's Forest Inventory and
Analysis National Program, FIA; @gray2012forest] play a ubiquitous role
in forest management, conservation, carbon accounting, and basic
research on forest ecosystems and community dynamics. These recent
surveys (the earliest FIA surveys began in the 1930s) can be extended
with longer-term historical data to understand how forest composition
has changed since EuroAmerican settlement. The Public Land Survey was
carried out ahead of mass EuroAmerican settlement west and south of Ohio
to provide for delineation and sale of the public domain beyond the
original East Coast states [@stewart1935public; @white1983history].
Because surveyors used trees to locate survey points, recording the
identity, distance, and directory of two to four trees next to each
survey marker, we can make broad-scale inferences about forest
composition and structure in the United States prior to large-scale
EuroAmerican settlement
[@almendinger1996minnesota; @liu2011broadscale; @williams2011testing; @tomscha2014historic].
In general, FIA datasets are systematically organized and widely
available to the forest ecology and modeling community, whereas most PLS
data compilations are of local or, at most, state-level extent. This
absence of widely available data on settlement-era forest composition
and structure limits our ability to understand and model the current and
future processes governing forest dynamics at broader, regional scales.
For example, distributional models of tree species often rely upon FIA
or other contemporary observational data to build species-climate
relationships that can be used to predict potential range shifts
[@iverson1998predicting; @iverson2013tree].

Here we use survey data from the original Public Lands Surveys (PLS) in
the upper Midwest to derive estimates of pre-settlement (*ca*. mid-late
1800s) forest composition, basal area, stem density, and biomass. This
work builds upon prior digitization and classification of PLS data for
Wisconsin [@manies2000testing; @schulte2002quantitative] and for parts
of Minnesota [@hanberry2012comparison; @friedman2005regional] and
Michigan Michigan (USFS-NCRS [<http://www.ncrs.fs.fed.us/gla/>]()). Most
prior PLS-based reconstructions are for individual states or smaller
extents (among others: @duren2012vegetation; @hanberry2012comparison;
@rhemtulla2009historical; @friedman2005regional] often aggregated at the
scale of regional forest zones
[@schulte2007homogenization; @hanberry2012comparison], although
aggregation may also occur at the section [@rhemtulla2009historical] or
township scale [@kronenfeld2010influence]. Our work develops new
approaches to address major challenges to PLS data, including lack of
standardization in tree species names, azimuthal censoring by surveyors,
variations in sampling design over time, and differential biases in tree
selection among different kinds of survey points within the survey
design at any point in time. The correction factors developed here are
spatially varying, allowing us to accommodate temporal and spatial
variations in surveyor methods.

We aggregate point based estimates of stem density, basal area and
biomass to an 8 x 8km grid, and classify forest types in the upper
Midwest to facilitate comparisons between FIA and PLS data. We compare
the PLS data to late-20th-century estimates of forest composition, tree
stem density, basal area and biomass. We explore how forest
homogenization has changed the structure of ecotones along two major
ecotones from southern deciduous to northern evergreen forests and to
the forest-prairie boundary. Using analog analyses, we identify lost
forests that have no close compositional counterpart today and novel
forests with no close historical analogs. This work provides insight
into the compositional and structural changes between historic and
contemporary forests, while setting the methodological foundation for a
new generation of maps and analyses of settlement-era forests in the
Eastern US.

Methods:
--------

### Public Lands Survey Data: Assembly, and Standardization

The PLS was designed to facilitate the division and sale of federal
lands from Ohio westward and south. The survey created a 1 mile^2^ (2.56
km^2^) grid (sections) on the landscape. At each section corner, a stake
was placed as the official location marker. To mark these survey points,
PLS surveyors recorded tree stem diameters, measured distances and
azimuths of the two to four trees 'closest'to the survey point and
identified tree taxa using common (and often regionally idiosyncratic)
names. PLS data thus represent measurements by hundreds of surveyors
from 1832 until 1907, with changing sets of instructions over time
(Stewart, 1979).

The PLS was undertaken to survey land prior to assigning ownership
(Stewart 1935, White 1983), replacing earlier town proprietor surveys
(TPS) used for the northeastern states
[@cogbill2002forests; @thompson2013four]. The TPS provided estimates of
relative forest composition at the township level, but no structural
attributes. The PLS produced spatially explicit point level data, with
information about tree spacing and diameter, that can be used to
estimate absolute tree density and biomass. PLS notes include tree
identification at the plot level, disturbance [@schulte2005severe] and
other features of the pre-settlement landscape. However, uncertainties
exist within the PLS and township level dataset [@bourdo1956review].

Ecological uncertainty in the PLS arises from the dispersed spatial
sampling design (fixed sampling every 1 mile), precision and accuracy in
converting surveyor's use of common names for tree species to scientific
nomenclature [@mladenoff2002narrowing], digitization of the original
survey notes, and surveyor bias during sampling
[@schulte2001original; @bourdo1956review; @manies2001assessing; @liu2011broadscale].
Estimates vary regarding the ecological significance of surveyor bias.
Terrail *et al*. [-@terrail2014early] show strong fidelity between taxon
abundance in early land surveys versus old growth plot surveys. Liu *et
al* [-@liu2011broadscale] estimate the ecological significance of some
of the underlying sources of bias in the PLS and show ecologically
significant (\>10% difference between classes) bias in species and size
selection for corner trees. However Liu *et al*. [-@liu2011broadscale]
also indicate that the true sampling error cannot be determined,
particularly because most of these historic ecosystems are largely lost
to us.

Kronenfeld and Wang [-@kronenfeld2007accounting], working with
historical land cover datasets in western New York indicate that direct
estimates of density using plotless estimators may be off by nearly 37%
due to azimuthal censoring (*i.e.*, the tendency of surveyors to avoid
trees close to cardinal directions), while species composition estimates
may be adjusted by between -4 to +6%, varying by taxon, although
Kronenfeld [-@kronenfeld2014validating] shows adjustments of less than
1%. These biases can be minimized by appropriate analytical decisions;
many efforts over the years have assessed and corrected for the biases
and idiosyncrasies in the original surveyor data
[@liu2011broadscale; @kronenfeld2007accounting; @williams2011testing; @manies2001assessing; @cogbill2015corrections; @hanberry2012adjusting; @hanberry2012comparison; @bouldin2008some; @hanberry2011spatial].
And, even given these caveats, PLS records remain the best source of
data about both forest composition and structure in the United States
prior to EuroAmerican settlement.

This analysis builds upon and merges previous state-level efforts to
digitize and database the point-level PLS data for Wisconsin, Minnesota
and the Upper Peninsula and upper third of the Lower Peninsula of
Michigan. These datasets were combined using spatial tools in R
[@team2014r; package *rgdal*: @bivand2014] to form a common dataset for
the upper Midwest (Figure 1) using the Albers Great Lakes and St
Lawrence projection (see code in Supplement 1, file:
*step\_one\_clean\_bind.R*; proj4: *+init:EPSG:3175*).

We took several steps to standardize the dataset and minimize the
potential effects of surveyor bias upon estimates of forest composition,
density, and biomass. All steps are preserved in the supplementary R
code (Supplement 1: *step\_one\_clean\_bind.R*). First, we excluded line
and meander trees (i.e. trees encountered along survey lines, versus
trees located at section or quarter corners) because surveyor selection
biases appear to have been more strongly expressed for line trees,
meander trees have non-random habitat preferences [@liu2011broadscale],
and the inherent differences in sampling design between line, meander
and corner points. We used only the closest two trees at each corner
point because the third and fourth furthest trees have stronger biases
with respect to species composition and diameter [@liu2011broadscale].
Corner points were used only if 1) there were at least two trees at a
survey point, 2) the two trees were from different quadrants (defined by
the cardinal directions), and 3) there were valid azimuths to the trees
(a defined quadrant with an angle between 0 and 90) and valid diameters
(numeric, non-zero).

Many species-level identifications used by PLS surveyors are ambiguous.
Statistical models can predict the identity of ambiguous species
[@mladenoff2002narrowing], but these models introduce a second layer of
uncertainty into the compositional data, both from the initial
surveyors' identification, and from the statistical disambiguation.
Given the regional scale of the analysis, and the inherent uncertainty
in the survey data itself, we chose to avoid this layer of taxonomic
uncertainty, and retained only genus-level identification (Supplement 2,
*Standardized Taxonomy*). In areas of open prairie or other treeless
areas, *e.g.* southwestern Minnesota, surveyors recorded distances and
bearings to 'Non Tree' objects. When points were to be located in water
bodies the point data indicates 'Water'. Points recorded "No Tree" are
considered to have been from extremely open vegetation, with an
estimated point-level stem density of 0 stems ha^-1^. We based our
estimates on terrestrial coverage, so water cells are excluded
completely. Hence, absence of trees at âNo Treeâ locations does
reduce the gridded estimates of terrestrial stem density, but absence of
trees at 'Water' locations does not.

Digitization of the original surveyor notebooks introduces the
possibility of transcription errors. The Wisconsin dataset was compiled
by the Mladenoff lab group, and has undergone several revisions over the
last two decades in an effort to provide accurate data
[@manies2000testing; @radeloff2000historical; @schulte2002quantitative; @mladenoff2002narrowing; @liu2011broadscale].
The Minnesota transcription error rate is likely between 1 and 5%, and
the treatment of azimuths to trees varies across the dataset
[@almendinger1996minnesota]. Michigan surveyor observations were
transcribed to mylar sheets overlaid on State Quadrangle maps, so that
the points were displayed geographically, and then digititized to a
point based shapefile (Ed Schools, pers. comm.; Great Lakes Ecological
Assessment. USDA Forest Service Northern Research Station. Rhinelander,
WI. <http://www.ncrs.fs.fed.us/gla/>), carrying two potential sources of
transciption error. Preliminary assessment of Southern Michigan data
indicates a transcription error rate of 3 - 6%. To reduce errors
associated with transcription across all datasets, we exclude sites for
which multiple large trees have a distance of 1 link (20.12 cm) to plot
center, trees with very large diameters (diameter at breast height -
dbh \> 100 in; 254 cm), plots where the azimuth to the tree is unclear,
and plots where the tree is at plot center but has a recorded azimuth.
All removed plots are documented in the code used for analysis
(Supplement 1: *step\_one\_clean\_bind.R*) and are commented for review.

### Data Aggregation

We binned the point data using an 64km^2^ grid (Albers Gt. Lakes St
Lawrence projection; Supplement 1: *base\_calculations.R*) to create a
dataset that has sufficient numerical power for spatial statistical
modeling and sufficient resolution for regional scale analysis
[@thurman2015composition]. This resolution is finer than the 100km^2^
gridded scale used in Freidman and Reich [-@friedman2005regional], but
coarser than township grids used in other studies
[@kronenfeld2014validating; @rhemtulla2009historical] to provide a scale
comparable to aggregated FIA data at a broader scale. Forest
compositional data is based on the number of individuals of each genus
or plant functional type (PFT) present at all points within a cell. Stem
density, basal area and biomass are averaged across all trees at all
points within the cell.

### Stem Density

Estimating stem density from PLS data is based on a plotless density
estimator using the measured distances from each survey point to the
nearest trees [@morisita1957estimation; @perrson1971robustness]. This
Morisita density estimator is then modified to minimize error due to
different sampling geometries and several known surveyor biases
[@liu2011broadscale; @kronenfeld2007accounting; @williams2011testing; @manies2001assessing; @cogbill2015corrections; @hanberry2012adjusting; @hanberry2012comparison; @bouldin2008some; @hanberry2011spatial].
Survey sampling instructions changed throughout the implementation of
the PLS in this region and differed between section and quarter section
points and between internal and external points within a township
[@liu2011broadscale; @white1983history]. Our approach allows for spatial
variation in surveyor methods by applying various spatially different
correction factors based not only on the empirical sample geometry, but
also on known surveyor biases deviating from this design
[@cogbill2015corrections].

We estimate stem density (stems m^-2^) based on a on a modified form of
the Morisita two-tree density estimator, which uses the distance-to-tree
measurements for the two closest trees at each point
[@morisita1954estimation]. Our modified form uses explicit and spatially
varying correction factors, modeled after the Cottam correction factor
[@cottam1956use], that account for variations in sampling designs over
time and among surveyors. All code to perform the analysis is included
in Supplement 1.

We estimate the basic stem density (stems m^-2^) using the point-to-tree
distances for the closest trees to each point within a defined number of
sectors around the point ([@morisita1957estimation] eqn 31.):

$\lambda \hat{m_2} = \frac{k - 1}{\pi \times n} \times \sum_{i=1}^N\frac{k}{\sum_{j=1}^{k}\left( {r_{ij}} \right)^2}$
(1)

where $\lambda$ is density ; $k$ is the number of sectors within which
trees are sampled, $N$ is the number of points over which estimates are
aggregated, $r$ is the distance of point-to-tree (as m). This estimate
can be modified by a refinement of the Cottam quadrant factors
[@morisita1954estimation; @cottam1956use] which recognizes that
different sampling designs, and the order of the distances in different
quadrants (or sectors) carry specific weights. This correction, herein
called $\kappa$, accounts for different sampling designs. When either
four quadrants or trees are sampled (point quarter design), or when two
trees in opposite semicircles (point halves design) are sampled, the
equation is accurate and $\kappa$ = 1; when the two trees are in the
nearest of two quadrants (two nearest quadrants design), $\kappa$ =
0.857; and when two trees are in quadrants on the same side of the
direction of travel (one-sided or interior half design), $\kappa$ = 2.
This parameter, in Cottam's notation [@cottam1956use], is a divisor of
the denominator above, or here, the mathematically equivalent multiplier
in the numerator of the reciprocal of the squared distances.

We further simplify the density estimate in equation (1) so that it is
calculated at each point (N=1) and for two sample trees only (k=2):

$\lambda_{M} = \frac{2}{\pi \times \sum_{j=1}^{2}{r_{j}}^2}$

Then the point values for any sampling design can be Cottam corrected
($\kappa \times \lambda_{M}$). For example, the basic Morisita equation
for two sectors assumes trees are located in opposite halves, so if the
actual design is the nearest tree in the two nearest quadrants, the
density from equation 2 will be overestimated and must be
correspondingly corrected by multiplying by $\kappa$ = 0.857.

Further corrections account for the restriction of trees to less than
the full sector ($\theta$), censoring of trees near the cardinal
azimuths ($\zeta$), and undersampling of trees smaller than a certain
diameter limit ($\phi$). These parameters are derived from analyses of
measurements of bearing angles and diameters actually observed in
surveys of witness trees within a subset of townships across the upper
Midwest.

Sector bias ($\theta$). Although the density model for two tree points
assumes that the trees are on opposite sides of a sample line (point
halves), often the actual sample is more restricted (\< 180^o^) within
the sector or is a less restricted (\> 180^o^) angle beyond the sector.
This deviation from the equation's assumption of equal distribution of
angles across the 180^o^ sector is quantified using the empirical angle
between the bearings of the two trees (pair angle). In the pair angle
frequency plot (Figure 2), the observed proportion of trees (p) within
any restricted sector divided by the proportion of that angle within the
circle ($\alpha$ is an estimate of the bias imposed by the actual
sampling (inspired by Kronenfeld & Wang [-@kronenfeld2007accounting]).
This factor ($\theta$ = p/$\alpha$) indicates bias associated with
differences in geometry of two tree samples. This parameter ($\theta$)
varies from 0.71 to 1.27, indicating sampling from effectively 253^o^ to
141^o^ sectors.

Azimuthal censoring ($\zeta$). In addition to sector bias, surveyors did
not always sample trees near the cardinal directions
[@kronenfeld2007accounting; @bouldin2008some; @hanberry2012adjusting].
This azimuthal censoring is commonly found along the line of travel on
section lines and sometimes on the perpendicular quarter-section lines.
Trees near the cardinal directions were passed over, and a replacement
was found within a more restricted angular region. The correction for
this bias is calculated following Kronenfeld and Wang
[-@kronenfeld2007accounting] in a manner similar to the sector bias. The
factor $\zeta$ is the ratio of the proportion of trees in the restricted
area (p) divided by the proportion of the complete circle ($\alpha$)
that is used. The azimuthal censoring parameter ($\zeta$) ranges from
1.03 to 1.25 indicating an equivalent to complete elimination of trees
from 10^o^ to 72^o^ azimuths adjacent to the cardinal directions.

Diameter limit ($\phi$). Examination of the diameter distributions from
settlement era surveys across the upper Midwest clearly demonstrate
witness trees less than 8 inches in diameter were undersampled
[@bouldin2008some; @liu2011broadscale; @cogbill2015corrections]. We have
confirmed this bias in our own inspection of plots of diameter frequency
in the PLS data, which show a strong mode at 8". This bias can be
accommodated by setting a diameter limit, and only calculating the
density for trees with diameters above this limit. Total density
calculated from all trees is reduced to this reference limit by simply
multiplying the total by the percentage of trees above this limit. This
effectively eliminates the smaller trees from the total and normalizes
the value of trees above this standard. The parameter ($\phi$)
represents diameter size bias is simply the percentage of trees $\geq$
8" and, in practice, ranges from 0.6 - 0.9.

Because all surveyor bias corrections are simple multipliers of the
model density and should be independent, the bias-minimized estimate of
the point density of trees $\geq$ 8" is:

$\lambda_{M corrected} = \kappa \times \theta \times \zeta \times \phi \times \lambda_{M}$
(3)

Estimates for each point *i* can be averaged for all *N* points in any
region. Correction factors are calculated separately for different
regions, years, internal versus external lines, section versus
quarter-section points, and surveyor sampling designs (Table 1). All
code to perform the analyses is included in Supplement 1 and the full
rationale for and calculation of these measures is described further in
Cogbill et al. [-@cogbill2015corrections].

### Basal Area and Biomass Estimates

Forest basal area is calculated by multiplying the point-based stem
density estimate by the average stem basal area from the reported
diameters at breast height for the closest two trees at the point (n=2).
Aboveground dry biomass (Mg ha^-1^) is calculated using the USFS FIA
tree volume and dry aboveground biomass equations for the United States
[@jenkins2004comprehensive].

Biomass equations share the basic form:

$m = Exp(\beta_{0} + \beta_{1} * \ln dbh)$

where $m$ represents stem biomass for an individual tree in kg.
$\beta_{0}$ and $\beta_{1}$ are the parameters described in Table 2 and
$dbh$ is the stem diameter at breast height (converted to cm) recorded
in the survey notes. The biomass estimates are summed across both trees
at a survey point and multiplied by the stem density calculated at that
point to produce an estimate of aboveground biomass reported in Mg
ha^-1^ [@jenkins2004comprehensive].

Matching PLSS tree genera to the species groups defined by Jenkins *et
al*. [-@jenkins2004comprehensive] is straightforward, placing the 22
genera used in this study into 9 allometric groups (Table 2). However,
all maples are assigned to the generic "Hardwood" group since separate
allometric relationships exist for soft and hard maple (Table 2).
Biomass estimates for "Non tree" survey points are assigned 0 Mg ha^-1^.

We use the stem density thresholds of Anderson and Anderson
[-@anderson1975presettlement] to discriminate prairie, savanna, and
forest.

### FIA Stem Density, Basal Area and Biomass

The United States Forest Service has monitored the nation's forests
through the FIA Program since 1929, with an annualized state inventory
system implemented in 1998 [@woudenberg2010forest]. On average there is
one permanent FIA plot per 2,428 ha of land in the United States
classified as forested. Each FIA plot consists of four 7.2m fixed-radius
subplots in which measurements are made of all trees \>12.7cm dbh
[@woudenberg2010forest]. We used data from the most recent full plot
inventory (2007-2011). The FIA plot inventory provides a median of 3 FIA
plots per cell using the 64km^2^ grid.

We calculated mean basal area (m^2^ ha^-1^), stem density (stems
ha^-1^), mean diameter at breast height (cm), and mean biomass (Mg
ha^-1^) for all live trees with dbh greater than 20.32cm (8in). Biomass
calculations used the same set of allometric regression equations as for
the PLS data [@jenkins2004comprehensive].

### Gridding and Analysing PLS and FIA Data

Spatial maps of stem density, basal area and biomass were generated by
averaging all PLS point or FIA plot estimates within a 64km^2^ raster
cell. Differences in sampling design between PLS and FIA data combined
with spatially structured forest heterogeneity will affect the
partitioning of within-cell versus between-cell variance, but not the
expected estimates. Most 64km^2^ cells have one or a few intensively
sampled FIA plots. Therefore at this scale of aggregation, the low
density of FIA plots in heterogeneous forests could result in high
within-cell variance and high between-cell variability. For the PLS
plotless (point based) estimates, stem density estimates are sensitive
to trees close to the plot center. Point-level estimates with very high
stem densities can skew the rasterized values, and it is difficult to
distinguish artifacts from locations truly characterized by high
densities. To accommodate points with exceptionally high densities we
carry all values through the analysis, but exclude the top 2.5
percentile when reporting means and standard deviations in our analysis.
PLS-based estimates are highly variable from point to point due to the
small sample size, but have low variance among 64 km^2^ raster cells due
to the uniform sampling pattern of the data. Thus within-cell variance
is expected to be high for the PLS point data, but spatial patterns are
expected to be robust at the cell level. The base raster and all
rasterized data are available as Supplement 3.

Standard statistical analysis of the gridded data, including
correlations and regression, was carried out in R [@team2014r], and is
documented in supplementary material that includes a subset of the raw
data to allow reproducibility. Analysis and presentation uses elements
from the following R packages: `cluster` [@maechler2014cluster],
`ggplot2` [@wickham2009ggplot2; @ggplot2009], `gridExtra` [@grid2012],
`igraph` [@csardi2006], `mgcv` [@wood2011], `plyr` [@wickham2011plyr],
`raster` [@hijmans2014], `reshape2` [@reshape2007], `rgdal`
[@bivand2014], `rgeos` [@bivand2014rgeos], `sp`
[@pebesma2005sp; @bivand2013], and `spdep` [@bivand2014spdep].

We identify analogs and examine differences in composition between and
within PLS and FIA datasets using Bray-Curtis dissimilarity [*vegdist*
in `vegan`; @oksanen2014vegan] for proportional composition within
raster cells using basal area measurements. For the analog analysis we
are interested only in the minimum compositional distance between a
focal cell and its nearest compositional (not spatial) neighbor. The
distribution of compositional dissimilarities within datasets indicates
forest heterogeneity within each time period, while the search for
closest analogs between datasets indicates whether contemporary forests
lack analogs in pre-settlement forests ('novel forests'), or vice versa
('lost forests'). For the analog analyses, we compute Bray-Curtis
distance between each 64km^2^ cell in either the FIA or the PLS periods
to all other cells within the other dataset (FIA to FIA, PLS to PLS),
and between datasets (PLS to FIA and FIA to PLS), retaining only the
minimum. For within era analyses (FIA - FIA and PLS - PLS), cells were
not allowed to match to themselves. We define vegetation classes for
lost and novel forests using k-medoid clustering [*pam* in `cluster`;
@maechler2014cluster]).

The differences in sampling design and scale between the PLS and FIA
datasets, described above, potentially affect between-era assessments of
compositional similarity [*e.g.*, @kronenfeld2010influence]. The effects
of differences in scale should be strongest in regions where there are
few FIA plots per 64 km^2^ cell, or where within-cell heterogeneity is
high. For the analog analyses, this effect should increase the
compositional differences between the FIA and PLS datasets. We test for
the importance of this effect on our analog analyses via a sensitivity
analysis in which we test whether dissimilarities between FIA and PLS
grid cells are affected by the number of PLS plots per cell. We find a
small effect, suggesting that our analyses are mainly sensitive to the
compositional and structural processes operating on large spatial
scales.

Results:
--------

### Data Standardization

The original PLS dataset contains 490,818 corner points (excluding line
and meander points), with 166,607 points from Wisconsin, 231,083 points
from Minnesota and 93,095 points from Michigan. Standardizing data and
accounting for potential outliers, described above, removed
approximately 1.5% points from the dataset, yielding a final total of
366,993 points with estimates used in our analysis.

Rasterizing the PLS dataset to the Albers 64km^2^ grid produces 7,939
raster cells with data. Each cell contains between 1 and 94 corner
points, with a mean of 61.8 ($\sigma$ = 15) and a median of 67 corners
(Supplement 3). Cells with a low number of points were mainly near water
bodies or along political boundaries such as the Canadian/Minnesota
border, or southern Minnesota and Wisconsin borders. Only 2.44% of cells
have fewer than 10 points per cell.

Species assignments to genera were rarely problematic. Only 18 PLS trees
were assigned to the Unknown Tree category, representing less than 0.01%
of all points. These unknown trees largely consisted of corner trees for
which taxon could not be interpreted, but for which diameter and azimuth
data was recorded. A further 0.011% of trees were assigned to the "Other
hardwood" taxon (*e.g.*, hawthorn, "may cherry", and "white thorn").

### Spatial Patterns of Settlement-Era Forest Composition: Taxa and PFTs

### Stem Density, Basal Area and Biomass

The mean stem density for the region (Figure 3a) is 153 stems ha^-1^.
Stem density exclusive of prairie is 172 stems ha^-1^ and is 216 stems
ha^-1^ when both prairie and savanna are excluded. The 95th percentile
range is 0 - 423 stems ha^-1^, and within-cell standard deviations
between 0 and 441 stems ha^-1^. Basal area in the domain (Figure 3c) has
a 95th percentile range between 0 and 63.5 m^2^ ha^-1^, a mean of 22.2
m^2^ ha^-1^, within-cell standard deviations range from 0 to 76.7 m^2^
ha^-1^. Biomass ranges from 0 to 209 Mg ha^-1^ (Figure 3d), with cell
level standard deviations between 0 and 569 Mg ha^-1^. High within-cell
standard deviations relative to mean values within cells for density,
basal area and biomass indicate high levels of heterogeneity within
cells, as expected for the PLS data, given its dispersed sampling
design.

In the PLS data, stem density is lowest in the western and southwestern
portions of the region, regions defined as prairie and savanna (Figure
3b, Table 3). When the Anderson and Anderson
[-@anderson1975presettlement] stem density thresholds (\<47 stems ha^-1^
for Savanna, Table 3) are used, the extent of area classified as savanna
is roughly equivalent to prior reconstructions
[@curtis1959vegetation; @rhemtulla2009legacies; @bolliger2004assessing]
(Figure 3b). The highest stem densities occur in north-central Minnesota
and in north-eastern Wisconsin (Figure 3a), indicating younger forests
and/or regions of lower forest productivity. The joint controls of
broad-scale climatic structuring and local hydrology controls on forest
composition and density can be seen, particularly along the Minnesota
River in south-western Minnesota, where a corridor of savanna is
sustained in a region mostly occupied by prairie (Figure 3b).

Forest structure during the settlement era can be understood in part by
examining the ratio of stem density to biomass, a measure that
incorporates both tree size and stocking. Regions in northern Minnesota
and northwestern Wisconsin have low biomass and high stem densities
(Figure 4, brown). This indicates the presence of young, small-diameter,
even-aged stands, possibly due to frequent stand-replacing fire
disturbance in the pre-EuroAmerican period or to poor edaphic
conditions. Fire-originated vegetation is supported by co-location with
fire-prone landscapes in Wisconsin [@schulte2005spatial]. High-density,
low-biomass regions also have shallower soils, colder climate, and
resulting lower productivity. High-biomass values relative to stem
density (Figure 4, green) are found in Michigan and southern Wisconsin.
These regions have higher proportions of deciduous species, with higher
tree diameters than in northern Minnesota. This is likely due to higher
soil productivity, and more temperate climates, with lower incidences of
stand-replacing disturbance in the southern and eastern portions of the
region.

Taxon composition within settlement-era forests is spatially structured
along dominant gradients from south to north (deciduous dominated to
conifer dominated forests) and from east to west (mixed wood forests to
open prairie) (Figure 5). Oak is dominant in the south of the region,
with an average composition of 21%, however, that proportion drops to 8%
when only forested cells are considered, due to its prevalence as a
monotypic dominant in the savanna and prairie. Pine shows the opposite
trend, with average composition of 14% and 17% in unforested and
forested cells respectively. Pine distributions represent three dominant
taxa, *Pinus strobus*, *Pinus resinosa* and *Pinus banksiana*. These
three species have overlapping but ecologically dissimilar
distributions, occuring in close proximity in some regions, such as
central Wisconsin, and are typically associated with sandy soils with
low water availability. Other taxa with high average composition in
forested cells include maple (10%), birch (10%), tamarack (9%) and
hemlock (8%).

For a number of taxa, proportions are linked to the total basal area
within the cell. For 4 taxa - hemlock, birch, maple and cedar - taxon
proportions are positively related to total basal area. For 17 taxa
including oak, ironwood, poplar, tamarack and elm, high proportions are
strongly associated with lower basal areas (Figures 3 and 5). This
suggests that hemlock, birch, maple and cedar occured in well-stocked
forests, with higher average dbh. These taxa are most common in Michigan
and in upper Wisconsin. Taxa with negative relationships to total basal
area (*e.g.*, spruce and tamarack) are more common in the northwestern
part of the domain.

Spruce in the PLS represents two species (*Picea glauca*, *Picea
mariana*) with overlapping distributions, but complex site preferences
(dry upland to wet-mesic for white spruce, and hydric for *P. mariana*,
although in northern Minnesota *P. mariana* also occurs on uplands).
Both cedar (*Thuja occidentalis*) and fir (*Abies balsamea*) are
mono-specific genera in this region.

Northern hardwoods, such as yellow birch and sugar maple, and beech, are
much less common in the lower peninsula of Michigan, and southern
Wisconsin, except along Lake Michigan. Birch has extensive cover in the
north, likely reflecting high pre-settlement proportions of yellow birch
(*Betula alleghaniensis*) on mesic soils, and paper birch on sandy
fire-prone soils and in northern Minnesota (birch proportions reach
upwards of 34.1% in northeastern Minnesota). Hardwoods in the southwest,
such as oak, elm, ironwood and basswood, are most typically
mono-specific groupings, with the exception of oak, which comprises 7
species (see Supplement 2). Hardwoods in the southwest are located
primarily along the savanna and southern forest margins, or in the
southern temperate deciduous forests. Finally, maple and poplar (aspen)
have a broad regional distribution, occupying nearly the entire wooded
domain. Poplar comprises four species in the region, while maple
comprises five species (Supplement 2). Both hardwood classes, those
limited to the southern portions of the region, and those with
distributions across the domain, correspond to well-defined vegetation
patterns for the region [@curtis1959vegetation]. Thus overlap among PFT
distributions (Figure 6) emerges from the changing composition within
the plant functional type from deciduous broadleaved species associated
with the southern, deciduous dominated region, to broadleafed deciduous
species associated with more northern regions in the upper Midwest.

### Structural Changes Between PLS and FIA Forests

Modern forests (FIA) have higher stem densities (146 stems ha^-1^,
$t_{1,5177}$ = 51.8, $p$ \< 0.01) and basal areas (-4.5 m^2^ ha^-1^,
$t_{1,5177}$ = -16.4, $p$ \< 0.01) than PLS forests, but overall, lower
biomass (-8.72 Mg ha^-1^, $t_{1,5177}$ = -6.55, $p$ \< 0.01) than
historical forests (Figure 7). We use only point pairs where both FIA
and PLS data occur since non-forested regions are excluded from the FIA
and as such . The similarity in biomass despite lower stem density and
total basal area in the PLS data is surprising. Two likely factors are
shifts in allometric scaling associated with changes in species
composition, or a higher mean diameter of PLS trees (Figure 7d).

The PLS has a lower overall mean diameter than the FIA ($\delta_{diam}$
= -2.9 cm, 95%CI from -17.3 to 8.18cm). FIA diameters are higher than
PLS diameters in the northwestern parts of the domain (on average 6.47
cm higher), overlapping almost exactly with regions where we have shown
low biomass-high density stands (Figure 4). At the same time, regions
with high biomass and low density stands, in northeastern Wisconsin, and
the Upper and Lower Peninsulas of Michigan, had higher average diameters
during the PLS than in the FIA, on average 3.65 cm higher. Thus we are
seeing an overal increase in tree size in the sub-boreal region and a
decrease in temperate mixedwood forests, where we find tree species with
much higher dbh:biomass ratios [@jenkins2004comprehensive]. This is
coupled with declining variance in dbh across the domain (from within
cell variance of 37.9 cm the PLS to 30.7 cm in the FIA). Thus, the
mechanism by which low density and basal area produce roughly equivalent
biomass estimates between the FIA and PLS is likely due to the
differential impacts of land clearence and subesequent forest management
in the south east vs the northwest. The loss of high biomass southern
hardwood forests is balanced by higher biomass in the northeast due to
fire supression and regeneration of hardwoods in the northwest.
Declining diameters from the PLS to FIA are most strongly associated
with higher abundances of poplar, ironwood and oak, while declining
diameters are associated with maple and hemlock, further supporting the
assertion that much of the loss in diameter, and, subsequently biomass,
is occuring in southeastern mixedwood/hardwood forest, while diameter
and biomass increases are occuring in the northwest.

Differences between FIA and PLS data in sampling design are unlikely to
be a factor; these differences are expected to affect how these datasets
sample local- to landscape-scale heterogeneity, but should not affect
the overall trends between datasets. Differences in variability
introduce noise into the relationship, but given the large number of
samples used here, the trends should be robust.

### Compositional Changes Between PLS and FIA Forests: Novel and Lost Forests

Both the PLS- and FIA-era compositional data show similar patterns of
within-dataset dissimilarity, with the highest dissimilarities found in
central Minnesota and northwestern Wisconsin. High within-PLS
dissimilarities are associated with high proportions of maple, birch and
fir while high within-FIA dissimilarities are associated with high
proportions of hemlock, cedar and fir. Dissimilarity values in the FIA
dataset are less spatially structured than in the PLSS. Moran's I for
dissimilarities within the FIA ($I_{FIA}$ = 0.198, p \< 0.001) are lower
than the dissimilarities within the PLSS ($I_{PLSS}$ = 0.496, p \<
0.001), suggesting lower spatial autocorrelation in the FIA dataset.
Cells with identical pairs represent 5.6% of the PLS cells and 7.44% of
FIA cells. Identical cells in the PLS are largely located along the
southern margin and most (69.5%) are composed entirely of oak. Cells in
the FIA with identical neighbors are composed of either pure oak
(19.4%), pure poplar (26%) or pure ash (14%).

There is a small but significant positive relationship ($F_{1,5964}$=
920, p \< 0.001) between the number of FIA plots and within-FIA
dissimilarity. The relationship accounts for 13% of total variance and
estimates an increase of $\delta_{d}$ = 0.0134 for every FIA plot within
a cell. This increase represents only 3.08% of the total range of
dissimilarity values for the FIA data. There is a gradient of species
richness that is co-linear with the number of FIA plots within a cell,
where plot number increases from open forest in the south-west to closed
canopy, mixed forest in the Upper Peninsula of Michigan. Hence,
differences in within- and between-cell variability between the PLS and
FIA datasets seem to be having only a minor effect on these
regional-scale dissimilarity analyses.

We define no-analog communities as those whose nearest neighbour is
beyond the 95%ile for dissimilarities within a particular dataset. In
the PLS dataset, forests that have no modern analogs are defined as
"lost forests", while forest types in the FIA with no past analogs are
defined as "novel forests". More than 25% of PLS sites have no analog in
the FIA dataset ('lost forests'; PLS-FIA dissmimilarity, Figure 8c),
while 29% of FIA sites have no analog in the PLS data ('novel forests';
FIA-PLS dissimilarity, Figure 8d). Lost forests show strong spatial
coherence, centered on the "Tension Zone" [@curtis1959vegetation], the
ecotone between deciduous forests and hemlock-dominated mixed forest
(Figure 5).

Lost forests are drawn from across the domain, and show strong
ecological and spatial coherence (Figure 8c). Forest classes generally
fall into five classes: Tamarack-Pine-Birch-Spruce-Poplar accounts for
28.8% of all lost forests and 7.97% of the total region. This forest
type is largely found in north eastern Minnesota, extending southward to
central Minnesota, into Wisconsin and along the Upper Peninsula of
Michigan, as well as in scattered locations on the Lower Peninsula of
Michigan (Figure 8c). This forest likely represents a mesic to hydric
forest assemblage, particularly further eastward. Modern forests
spatially overlapping this lost type are largely composed of poplar
($\bar{x}_{FIA}$ = 12%) and oak ($\bar{x}_{FIA}$ = 12%). Tamarack in
these forests has declined significantly, from 23% to only 5% in the
FIA, while Poplar has increased from 10% to 22%, resulting in forests
that look less mesic and more like early seral forests.

Cedar/juniper-Hemlock-Pine accounts for 19.8% of all lost forests and
5.49% of the total region. This forest type is found largely in
northeastern Wisconsin, and the Upper and Lower Peninsulas of Michigan.
This lost forest type has been predominantly replaced by maple, poplar,
and pine, retaining relatively high levels of cedar ($\bar{x}_{PLS}$ =
19%; $\bar{x}_{FIA}$ = 14%). The loss of hemlock is widespread across
the region, but particularly within this forest type, declining to only
3% from a pre-settlement average of 18%.

Elm-Oak-Basswood-Ironwood accounts for 19.6% of all lost forests and
5.42% of the total region. The region is centered largely within savanna
and prairie-forest margins, both in south-central Minnesota and in
eastern Wisconsin, but, is largely absent from savanna in the Driftless
area of southwestern Wisconsin. These forests were historically elm
dominated ($\bar{x}_{PLS}$ = 25%), not oak dominated savanna, as
elsewhere (particularly in the Driftless). Modern forests replacing
these stands are dominated by oak and ash, with strong components of
maple, and basswood. Elm has declined strongly in modern forests
($\bar{x}_{FIA}$ = 1%), possibly in part due to Dutch Elm Disease and
land use. The increase in ash in these forests is substantial, from
$\bar{x}_{PLS}$ = 5% to $\bar{x}_{FIA}$ = 15%.

Hemlock-Birch-Maple-Pine accounts for 19.2% of all lost forests and
5.33% of the total region. This forest type, dominant in north central
Wisconsin, was dominated by hemlock ($\bar_{x}_{PLS}$ = 26%) and what
was likely late seral yellow birch ($\bar{x}_{PLS}$ = 24%), replaced
largely by maple (from $\bar{x}_{PLS}$ = 12\$ to $\bar{x}_{FIA}$ = 27%).
Poplar increases from 1% to 13% in the FIA, again indicating a shift to
earlier seral forests in the FIA. Hemlock is almost entirely lost from
the forests, declining from 26% to 4% in the FIA.

Lastly, Beech-Maple-Hemlock accounts for 12.6% of all lost forests and
3.49% of the total region. This forest type is found exclusively on the
central, western shore of Lake Michigan and in the Lower Peninsula, in
part due to the limited geographic range of Beech in the PLS dataset
(Figure 5). Beech is almost entirely excluded from the modern forests in
this region, declining from $\bar{x}_{PLS}$ = 37% to $\bar{x}_{FIA}$ =
4%. Pine in the region increases from 9% to 16%, while maple, the
dominant taxa in the modern forests, increases from 16 - 25%.

On average lost forests contain higher proportions of ironwood (r =
0.203), beech (r = 0.2), birch (r = 0.189) and hemlock (r = 0.188) than
the average PLS forest, and lower proportions of oak (r = -0.28), poplar
(r = -0.145), and pine (r = -0.107).

The distribution of novel ecosystems (Figure 8d) is spatially diffuse
relative to the lost forest of the PLS and the forest types tend to have
fewer co-dominant taxa. FIA novel forest types also have a more uneven
distribution in proportion than the PLS lost forests. Overall, novel
forests are associated with higher proportions of maple (r = 0.02), ash
(r = 0.03) and basswood (r = -0.04), although basswood is dominant in
only one forest type (Poplar-Cedar/juniper-Maple). Novel forests are
associated with lower proportions of oak (r = -0.28), and pine (r
= -0.11). This analysis suggests that the loss of particular forest
types associated with post-settlement land use was concentrated in mesic
deciduous forests and the ecotonal transition between southern and
northern hardwood forests, while the gains in novelty were more
dispersed, resulting from an overall decline in seral age.

By far the largest novel forest type is Maple, which accounts for 37.2%
of all novel forests and 2.68% of the total region. As with all novel
forest types, this forest type is broadly distributed across the region.
This forest type is associated with co-dominant maple ($\bar{x}_{FIA}$ =
23%) and ash ($\bar{x}_{FIA}$ = 22%). Hemlock has declined significantly
across this forest type, from $\bar{x}_{PLS}$ = 24% to $\bar{x}_{FIA}$ =
4%.

Poplar-Cedar/juniper-Maple, accounts for 28.8% of all novel forests and
2.08% of the total region. The broad distributiof these novel forests
makes assigning a past forest type more difficult than for the PLS lost
forests, the distribution replaces two classes of past forest, one
dominated by oak, in southern Wisconsin and Minnesota, the other by
mixed hemlock, beech, birch and cedar forests.

Pine-Cedar/juniper-Poplar-Maple forest accounts for 16.3% of all novel
forests and 1.17% of the total region. This forest type is again broadly
distributed, and is widely distributed across the region, representing a
homogenous, early seral forest type, likely associated with more mesic
sites. Oak forest accounts for 13.3% of all novel forests and 0.96% of
the total region. This grouping again shows a pattern of broad
distribution across the region, associated with cedar/juniper
percentages near 40%, with smaller components of poplar (14%) and maple
(13%).

### Compositional Changes Between PLS and FIA Forests: Ecotone Structure

To understand how the ecotonal structure has been transformed by
post-settlement land use, we constructed two transects of the FIA and
PLS data (Figure 10a), and fitted GAM models to genus abundances along
these transects. Transect One (T1) runs from northern prairie (in
northern Minnesota) to southern deciduous savanna in southeastern
Wisconsin (left to right in Figures 10c-f), while Transect Two (T2) runs
from southern prairie in southwestern Minnesota to northern mixedwood
forest in the Upper Peninsula of Michigan (left to right in Figures
10g-j). In general, these transect analyses show: 1) significant
differences in ecotonal structure between the present and
pre-settlement, and 2) steeper ecotones in the past and more diffuse
ecotones today.

For T1, GAM models show significant differences (using AIC) between time
periods in curves for all broadleafed taxa (Figure 10b & c) and for al
needleleafed taxa (Figures 10d and e). The PLS curves show a rapid
transition in the northwest from oak to poplar dominated open forest
that then transitions to a needleleafed forest composed of pine, spruce
and tamarack, with high proportions of tamarack grading to pine further
to the south east. Tamarack and poplar proportions decline gradually
from the east, being replaced first by pine, then briefly by maple and
birch, and then. ultimately by oak as the transect grades into oak
savanna. In the FIA dataset oak and poplar in the northwest appears to
decline simultaneously, grading into needleleafed forests that are
absent from the FIA dataset in the first 100km along the transect. While
the PLS transect shows distinct vegetation types in the central partof
the transect, the FIA shows relatively constant proportions of oak,
pine, spruce, poplar and maple before pine, oak and elm increase in the
southeastern portions of the transect.

The second transect shows a similar pattern, with well defined ecotones
in the pre-settlement period(Figure 10f and h), that are largely absent
from the FIA data (Figure 10g and i). Oak forest, with a component of
elm and poplar in the southwest grades slowly to a rapid transition zone
where pine, elm, maple (first), then rapidly birch, hemlock and
tamarack, and later, spruce, increase. This region, the Tension Zone,
extends from 3 x 10^5^ to 4.5x10^5^ meters East, eventually becoming a
forest that shows co-dominance between birch, pine, maple, spruce and
tamarack, likely reflecting some local variability as a result of
topographic and hydrological factors. Missing data at the beginning of
the FIA transect reflects a lack of FIA plots in unforested regions in
the west

Contemporary forests show broader homogenization and increased
heterogeneity (evidenced by the lower within-FIA Moran's I estimates for
near-neighbor distances) at a local scale in the region. Homogenization
is evident across T1, where Bray-Curtis dissimilarity between adjacent
cells declines from the PLSS to the FIA ($\delta_{beta}$ = -0.22,
$t_{113}$ = -7.93, p\<0.001), mirroring declines in the pine barrens
between the 1950s and the present [@li2014drivers]. The PLS shows strong
differentiation in the central region of T2 where maple-pine-oak shifts
to pine-poplar-birch forest (Figure 9d). This sharp ecotone is not
apparent in the FIA data, which shows gradual and blurred changes in
species composition across the ecotone (Figure 10i). $\beta$-diversity
along T2 is lower in the FIA than in the PLSS ($\delta_{beta}$ = -0.19,
$t_{65}$=-7.34, p \< 0.01), indicating higher heterogeneity in the PLS
data at the 64 km^2^ meso-scale.

Across the entire domain, $\beta$ diversity is lower in the FIA than in
the PLS ($\delta_{\beta}$ = -0.172, $t_{1.3e7}$ = 2480, p \<0.001),
lending support to the hypothesis of overall homogenization. Differences
in sampling design between PLS and FIA data cannot explain this
homogenzation, since its effect would have been expected to increase
$\beta$-diversity along linear transects and at larger spatial scales.

Discussion
----------

Many forests of the PLS, are no longer a part of the modern landscape.
Forest types have been lost at the 64 km^2^ mesoscale, the grain of our
analysis, reflected in the loss of some forest types and the gain of
novel forest types. Ecotones in forest composition are weaker now than
in the past (Fig. 10), with clear signs of increased homogenization at
local and regional scales and decreased spatial structure in vegetation
assemblages. Decreased $\beta$ diversity along regional transects
indicates homogenization at meso-scales of 100s of km^2^, while the
overall reduction in Moran's I for dissimilarity in the FIA indicates a
regional reduction in heterogeneity on the scale of 1000s of km^2^. The
selective loss or weakening of major vegetation ecotones, particularly
in central Wisconsin, and the development of novel species assemblages
across the region. These changes are the result of land use change, both
agricultural and logging, but affect forests in contrasting ways across
the domain. Maple has become one of the most dominant taxa across the
region, while in northern Minnesota, forest biomass has increased and
species shifts have reflected increases in poplar and pine, while in
southern Wisconsin, biomass has declined, and hemlock has been lost
almost completely.

Anthropogenic shifts in forest composition over decades and centuries
seen here and elsewhere [@thompson2013four; @cogbill2002forests] are
embedded within a set of interacting systems that operate on multiple
scales of space and time [macrosystems, *sensu*
@heffernan2014macrosystems]. Combining regional historical baselines,
long term ecological studies and high frequency analyses can reveal
complex responses to climate change at local and regional scales
[@groffman2012long]. Estimates of pre-settlement forest composition and
structure are critical to understanding the processes that govern forest
dynamics because they represent a snapshot of the landscape prior to
major EuroAmerican land-use conversion. Pre-settlement vegetation
provides an opportunity to test forest-climate relationships prior to
land-use conversion and to test dynamic vegetation models in a data
assimilation framework [*e.g.*, @hartig2012connecting]. For these
reason, the widespread loss of regional forest associations common in
the PLS (Figure 8d), and the rapid rise of novel forest assemblages
(Figure 8e) have important implications for our ability to understand
ecological responses to changing climate. The loss of these forests
implies that the modern understanding of forest cover, climate
relationships, realized and potential niches and species associations
may be strongly biased toward a single state, even though 38% of the
total regional cover is novel relative to forests only two centuries
ago.

Beyond shifts in composition at a meso-scale, the broader shifts in
ecotones can strongly impact models of species responses and
co-occurence on the landscape. For example, the heterogeneity,
distribution, and control of savanna-forest boundaries
[@staver2011global] is of particular interest to ecologists and modelers
given the ecological implications of current woody encroachment on
savanna ecosystems [@ratajczak2012woody]. Declines in landscape
heterogeneity may also strongly affect ecosystem models, and predictions
of future change. Recent work using the FLUXNET tower network has shown
that energy budgets are strongly related to landscape measures of
heterogeneity [@stoy2013data]. Our data show higher levels of
heterogeneity at mesoscales during the pre-settlement era, and greater
fine scaled turnover along transects. Lower $\beta$ diversity shown here
and elsewhere [@li2014drivers] indicate increasing heterogeneity at a
very large spatial scale, and the loss of resolution along major
historical ecotones. Thus analysis of the processes governing vegetation
heterogeneity and ecotones may inadvertently and substantially
incorporate anthropogenic processes.

Methodological advances of the current work include 1) the systematic
standardization of PLS data to enable mapping at broad spatial extent
and high spatial resolution, 2) the use of spatially varying correction
factors to accommodate variations among surveyors in sampling design,
and 3) parallel analysis of FIA datasets to enable comparisons of forest
composition and structure between contemporary and historical time
periods. This approach is currently being extended to TPS and PLS
datasets across the north-central and northeastern US, with the goal of
providing consistent reconstructions of forest composition and structure
for northeastern US forests at the time of EuroAmerican forests.

All datasets and analytic codes presented here are publicly available
and open source ([<http://github.som/SimonGoring/WitnessTrees>]()), with
the goal of enabling further analyses of ecological patterns across the
region and the effects of post-settlement land use on forest composition
and structure. Our results support the consensus that robust estimates
of pre-settlement forest composition and structure can be obtained from
PLS data [*e.g.*, Wisconsin: @schulte2002quantitative; California:
@williams2011testing; Iowa: @rayburn2009integrating; Oregon:
@duren2012vegetation]. Patterns of density, basal area and biomass are
roughly equivalent to previous estimates
[@schulte2007homogenization; @rhemtulla2009historical]. Our results for
stem density are lower than those estimated by Hanberrry *et al*.
[@hanberry2012comparison] for eastern Minnesota, but density and basal
area are similar to those in the northern Lower Peninsula of Michigan
[@leahy2003comparison] and biomass estimates are in line with estimates
of aboveground carbon for Wisconsin [@rhemtulla2009historical].

These maps of settlement-era forest composition and structure can also
provide a useful calibration dataset for pollen-based vegetation
reconstructions for time periods prior to the historic record. Many
papers have used calibration datasets comprised of modern pollen samples
to build transfer functions for inferring past climates and vegetation
from fossil pollen records
[@birks2010strengths; @paciorek2009mapping; @jacques2008pre; @goring2009new].
However, modern pollen datasets are potentially confounded by recent
land use, which can alter paleoclimatic reconstructions
[@jacques2008pre]. By linking pollen and vegetation at modern and
historical periods we develop capacity to provide compositional datasets
at broader spatio-temporal scales, providing more data for model
validation and improvement. Ultimately, it should be possible to
assimilate these empirical reconstructions of past vegetation with
dynamic vegetation models in order to infer forest composition and
biomass during past climate changes. Data assimilation, however,
requires assessment of observational and model uncertainty in the data
sources used for data assimilation. Spatiotemporal models of uncertainty
are being developed for the compositional data [@thurman2015composition]
and biomass data (Feng *et al*. in prep.).

Ultimately the pre-settlement vegetation data present an opportunity to
develop and refine statistical and mechanistic models of terrestrial
vegetation that can take multiple structural and compositional forest
attributes into account. The future development of uncertainty estimates
for the data remains an opportunity that can help integrate
pre-settlement estimates of composition and structure into a data
assimilation framework to build more complete and more accurate
reconstructions of past vegetation dynamics, and to help improve
predictions of future vegetation under global change scenarios.

------------------------------------------------------------------------

**Table 1**. *Correction values based on plot level survey design using
state, year, and location within township as a basis for assignment.
Years reported represent the upper bound for each set of survey years.
Internal points are points within the township, external points are on
the township boundary; no sampling occured outside of a township
boundary so plots were limited to half of the space for internal points.
Townships are divided into Section and Quarter Sections, at most section
points andsome quarter section points, suvey instructions indicated four
trees were to be sampled, these were '2nQ' plots, wheras others surveyed
only two points in adjacent plot halves ('P' plots).*

  State   Survey Year   Internal   Section   Trees   kappa    theta   zeta   phi
  ------- ------------- ---------- --------- ------- -------- ------- ------ ------
  Wisc    1845          ext        Sec       P       2        0.82    1.14   0.89
  Wisc    1845          ext        QSec      P       1        1.29    1.11   0.89
  Wisc    1845          int        Sec       P       1        1.14    1.17   0.89
  Wisc    1845          int        QSec      P       1        1.08    1.06   0.85
  Wisc    1845          ext        Sec       2nQ     0.86     1       1.21   0.86
  Wisc    1845          ext        QSec      2nQ     0.8563   1       1.11   0.91
  Wisc    1845          int        Sec       2nQ     0.86     1       1.24   0.92
  Wisc    1845          int        QSec      2nQ     0.86     1       0.75   0
  Wisc    1907          ext        Sec       P       2        0.89    1.16   0.9
  Wisc    1907          ext        QSec      P       2        0.9     1.14   0.84
  Wisc    1907          int        Sec       P       1        1.07    1.12   0.9
  Wisc    1907          int        QSec      P       1        1.04    1.04   0.8
  Wisc    1907          ext        Sec       2nQ     0.86     1       1.13   0.99
  Wisc    1907          ext        QSec      2nQ     0.86     1       1.12   0
  Wisc    1907          int        Sec       2nQ     0.8563   1       1.24   0.83
  Wisc    1907          int        QSec      2nQ     0.8563   1       1      0
  Mich    all           ext        Sec       P       2        0.87    1.25   0.85
  Mich    all           ext        QSec      P       1        0.94    1.21   0.76
  Mich    all           int        Sec       P       1        1.27    1.24   0.85
  Mich    all           int        QSec      P       1        1.26    1.15   0.77
  Mich    all           ext        Sec       2nQ     0.86     1       1.24   0.84
  Mich    all           ext        QSec      2nQ     0.86     1       1.35   0.85
  Mich    all           int        Sec       2nQ     0.8563   1       1.26   0.84
  Mich    all           int        QSec      2nQ     0.8563   1       1.28   0.68
  Minn    1855          ext        Sec       P       2        0.71    1.19   0.67
  Minn    1855          ext        QSec      P       1        1.05    1.11   0.68
  Minn    1855          int        Sec       P       1        0.71    1.05   0.76
  Minn    1855          int        QSec      P       1        1.09    1.03   0.6
  Minn    1855          ext        Sec       2nQ     0.86     1       1.17   0.66
  Minn    1855          ext        QSec      2nQ     0.86     1       1      0.68
  Minn    1855          int        Sec       2nQ     0.8563   1       1.5    0.59
  Minn    1855          int        QSec      2nQ     0.8563   1       1      0.25
  Minn    1907          ext        Sec       P       2        0.71    1.19   0.67
  Minn    1907          ext        QSec      P       1        1.05    1.11   0.68
  Minn    1907          int        Sec       P       1        0.71    1.05   0.76
  Minn    1907          int        QSec      P       1        1.09    1.03   0.6
  Minn    1907          ext        Sec       2nQ     0.86     1       1.17   0.66
  Minn    1907          ext        QSec      2nQ     0.86     1       1      0.68
  Minn    1907          int        Sec       2nQ     0.8563   1       1.5    0.59
  Minn    1907          int        QSec      2nQ     0.8563   1       1      0.25

------------------------------------------------------------------------

**Table 2.** *Biomass parameters used for the calculation of biomass in
the pre-settlement dataset(rounded for clarity).*

  Jenkins Species Group             $\beta_{0}$   $\beta_{1}$   PalEON Taxa Included (Supp. 2)
  --------------------------------- ------------- ------------- ------------------------------------------------------------------------------------
  Aspen, Alder, Poplar, Willow      -2.20         2.38          Poplar, Willow, Alder
  Soft Maple, Birch                 -1.91         2.36          Birch
  Mixed Hardwood                    -2.48         2.48          Ash, Elm, Maple, Basswood, Ironwood, Walnut, Hackberry, Cherries, Dogwood, Buckeye
  Hard Maple, Oak, Hickory, Beech   -2.01         2.43          Oak, Hickory, Beech, Other Hardwood
  Cedar and Larch                   -2.03         2.26          Tamarack, Cedar
  Fir and Hemlock                   -2.54         2.43          Fir, Hemlock
  Pine                              -2.54         2.43          Pine
  Spruce                            -2.08         2.33          Spruce

------------------------------------------------------------------------

**Table 3**. *Forest classification scheme used in this paper for
comparison between pre-settlement forests and the Haxeltine and Prentice
[-@haxeltine1996biome3] potential vegetation classes represented in
Ramankutty and Foley [@ramankutty1999estimating]. Plant functional types
(PFTs) for grasslands (CG, grassland; Non-Tree samples in the PLS),
broad leafed deciduous taxa (BDT) and needleleaded evergreen taxa (NET)
are used, but leaf area index used in Haxeltine and Prentice
[-@haxeltine1996biome3] is replaced by stem density classes from
Anderson and Anderson [@anderson1975presettlement].*

  Forest Class          Haxeltine & Prentice Rules                         Current Study
  --------------------- -------------------------------------------------- ------------------------------------------------------------------
  Prairie               Dominant PFT CG, LAI \> 0.4                        Stem dens. \< 0.5 stem/ha
  Savanna               Dominant PFT CG, LAI \> 0.6                        1 \< Stem dens. \< 47 stems ha^-1^
  Temperate Deciduous   Dominant PFT BDT, LAI \> 2.5                       Stem dens. \> 48 stems ha^-1^, BDT \> 70% composition
  Temperate Conifer     Dominant PFT (NET + NDT), LAI \> 2.5               Stem dens. \> 47 stems ha^-1^, NET + NDT \> 70% composition
  Mixedwood             Both BDT (LAI \> 1.5) & NET (LAI \> 2.5) present   Stem dens. \> 47 stems ha^-1^, BDT & NET both \< 70% composition

------------------------------------------------------------------------

**Table 4**. *Classification proportions, patch size and edge cell
proportion for various classification schemes used with the Public Land
Survey data. All patch size estimates are in 1000s of km^2^.*

  Classification Metric                               Prairie   Savanna   Temperate Deciduous   Temperate Evergreen   Mixed Wood   Mean Patch Size (km^2^)   Edge Cells (%)
  --------------------------------------------------- --------- --------- --------------------- --------------------- ------------ ------------------------- ----------------
  Ramankutty and Foley [-@ramankutty1999estimating]   86        63        9                     50                    241          181                       27.5
  PLS Data                                            51        90        83                    149                   100          37                        72

------------------------------------------------------------------------

**Figure 1**. *The domain of the Public Land Survey investigated in this
study. The broad domain includes Minnesota, Wisconsin and the upper two
thirds of Michigan state. A 8x8km grid is superimposed over the region
to aggregate data, resulting in a total of 7940 cells containing data.*

**Figure 2**. *Correction factors for $\zeta$ in the PLS data, and the
associated distribution of azimuths for each $\zeta$ value, by panel.
High peaks represent midpoints for quadrants where azimuth is defined
as* e.g.*, NE or SW. Greater differences between cardinal directions and
other azimuths result in higher $\zeta$ values, excluding the peaked
values.*

**Figure 3.** *Total stem density (a) in the Upper Midwest, along with
forest type classification (b) based on PLS data and the stem density
thresholds defined by Anderson and Anderson
[-@anderson1975presettlement]; Table 3). Fine lines represent major
rivers. To a first order, basal area (c) and biomass (d) show similar
patterns to stem density (but see* Figure 4*).*

**Figure 4.** *Maps of the ratio between biomass and stem density as an
indicator of forest structure. Regions with high stem density to biomass
ratios (blue) indicate dense stands of smaller trees, while regions with
low stem density to biomass ratios (red) indicate larger trees with
wider spacings.*

**Figure 5**. *Forest composition (%) for the 15 most abundant tree
taxa. The scale is drawn using a square-root transform to emphasize low
abundances. Shading of the bar above individual taxon maps indicates
plant functional type assignments (dark gray: needleleafed deciduous;
light gray: needleleafed evergreen; white: broadleafed deciduous).*

**Figure 6**. *Proportional distribution of Plant Functional Types
(PFTs) in the upper Midwest from PLS data, for broadleaved deciduous
trees (BDT), needleleaved deciduous trees (NDT), and needleleaved
evergreen trees (NET). Distributions are shown as proportions relative
to total basal area, total biomass, and composition (*Figure 3*). The
grassland PFT is mapped onto non-tree cells with the assumption that if
trees were available surveyors would have sampled them.*

**Figure 7**. *The relationship between average stem density, total
basal area and biomass values in the PLS and FIA datasets. Stem density
and total basal area are higher in the FIA than in the PLS, however
average cell biomass is higher in the PLS.*

**Figure 8**. *Minimum dissimilarity maps. Distributions of minimum
(within dataset) dissimilarities during the PLS (a) and FIA (b) show
spatially structured patterns of dissimilarity, with stronger spatial
coherence for the PLS. Lost forests (c) show strong compositional and
spatial coherence, and have more taxa with percent composition \> 10%
than within Novel forests during the FIA era (d).*

**Figure 9**. *Transects (a) across the region show clear changes in the
ecotonal strength. Transect One shows shifts in broad-leafed taxon
distributions from the PLS to FIA (b and c) and in needle-leafed
distributions (d and e). Transect Two broadleaf (f and g) and needleleaf
(h and i) taxa show shifts that again appear to represent regional scale
homogenization. Ecotones in the pre-settlement era were stronger in the
past than they are in the present. Fitted curves represent smoothed
estimates across the transects using Generalized Additive Models using a
beta family.*

------------------------------------------------------------------------

![](figure/fig1_output-1.png)

------------------------------------------------------------------------

![](figure/fig2_output-1.png)

------------------------------------------------------------------------

![](figure/figX_output-1.png)

**Currently Unscripted Plot** Angle pairs for theta.

------------------------------------------------------------------------

![](figure/fig3_output-1.png)

------------------------------------------------------------------------

![](figure/fig4_output-1.png)

------------------------------------------------------------------------

![](figure/fig5_output-1.png)

------------------------------------------------------------------------

![](figure/fig6_output-1.png)

------------------------------------------------------------------------

------------------------------------------------------------------------

![](figure/fig9_output-1.png)

------------------------------------------------------------------------

![](figure/fig10_output-1.png)

------------------------------------------------------------------------
