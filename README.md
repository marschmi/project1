# Significant Differences between Particle-Associated and Free-Living Bacterial Community Composition in Freshwater Lakes
Marian L. Schmidt  

###Introduction
Though aquatic ecosystems may appear homogeneous upon first glance, there is quite a lot of heterogeneity in habitat structure, especially in dynamic environments like freshwater lakes.  Freshwater lakes have a very seasonal behavior that is driven by temperature-induced density changes in water. Since water is most dense at 3.9 degrees celsius and as it heats up or cools, it becomes less dense.  Typically, temperate freshwater lakes stratify in the winter and summer and undergo full water column mixis in spring and fall.  At the onset of stratification, surface waters and bottom waters are restricted from interacting.  Therefore, as stratification prolongs (especially in summer), there are major shifts in the chemistry and biology surface waters and bottom waters.  

An important habitat for bacteria in many aquatic ecosystems are particles.  Particles in aquatic ecosystems have been deemed “hotspots of microbial activity” (Azam, 1998; Grossart, 2010), “hotbeds for genome reshuffling” (Ganesh et al., 2014), and active sites for organic matter mineralization (Grossart, 2010).  Therefore, I was interested in analyzing particle-associated bacteria in comparison to free-living bacteria in stratified lakes.

While I have done most of the analysis of bacterial community composition,  I was interested in testing for significant differences between particle associated and free-living bacteria in the surface waters versus the bottom waters in stratified lakes.

###Research Questions and Hypotheses
1.  What bacterial phyla are the dominant members of the overall freshwater lake community that I sampled?
2.  Do particle-associated bacteria differ from free-living bacteria in **surface** waters?  
      + **Hypothesis:** Particle-associated and free-living bacteria **do not** significantly differ in surface waters.
3.  Do particle-associated bacteria differ from free-living bacteria in **bottom** waters?
      + **Hypothesis:** Particle-associated and free-living bacteria **do** significantly differ in surface waters.


###Methods
I have characterized the bacterial community composition of 10 lakes in southwestern Michigan by Michigan State University’s Kellogg Biological station that fall across a productivity gradient, which include three eutrophic (high productivity: Baker, Baseline and Wintergreen), three mesotrophic (medium productivity: Bassett, Bristol, and Payne), and four oligotrophic (low productivity: Gull, Lee, Little Long and Sixteen) lakes.  At the deepest location within each lake, I performed a lake profile using a multi-probe Hydrolab sonde for pH, temperature, conductivity, and dissolved oxygen concentration.  I collected duplicate samples of the particle-associated (20–3 μm fraction) and free-living (3-0.22 μm fraction) bacterial communities from the surface and bottom layer of the lake.  Next, I inferred bacterial community composition using high throughput tag sequencing of the V4 hypervariable region of the 16S rRNA gene (Caporaso et al., 2012). I subsequently processed the sequences using the MiSeq SOP (Kozich et al., 2013, accessed September 2014) with the Mothur program (Schloss et al., 2009) to create a shared and a taxonomy file.  I then inported my shared and taxonomy file into RStudio to combine them into taxonomic-level tables.  These tables include categorical information (e.g. particle-association, lake layer, etc), sequence abundance, and means and standard deviations for categorical information.  Here, I will be working with my phylum-level and my class level data.  

####Programming methodology:
1. Combine my phylum-level and class-level tables to include the sub-phyla of the proteobacteria.  
2. Plot the relative abundance of each bacterial phyla in all of my samples.  
3. Subset the data into the surface and bottom waters.  
4. Calculate and plot a Log2-fold odds ratio as a visualization of changes in relative abundance.  
5. Run a significance test via the non-parametric paired Wilcoxon Test.  
      + I have decided to do a paired Wilcoxon test because my sample collection methodology included a successive filtration of the same water sample through a 3 μm filter (particle-associated) and 0.22 μm (free-living) filter. 


###Results
First, I wondered what the mean relative abundance of each bacterial phyla across all my samples.  Who are the dominant members of the overall freshwater community in the lakes that I sampled?


<img src="README_files/figure-html/Read in Functions & Abundance-1.png" title="" alt="" style="display: block; margin: auto;" />

_**Figure 1:**  The relative abundance of each bacterial phyla above 0.5% in all samples.  Error bars represent standard error._  

_Bacteriodetes_ is the most common member of my freshwater lake samples as it has a mean of 23.55%.  _Cyanobacteria_ rank second with a mean abundance of 14.21% and then _Verrucomicrobia_ at 11.88%.  This goes in line with other freshwater bacterioplankton studies (Newton et al., 2011; Bizic-Ionescu et al., 2014).












#### Comparing the Surface and Bottom Water Particle-Associated versus Free-Living Bacterial Community Composition 
<img src="README_files/figure-html/sidebyside log2-fold abundance ratio-1.png" title="" alt="" style="display: block; margin: auto;" />

_**Figure 2:**  The Log2-Fold Abundance Ratio of bacterial phyla in surface and bottom waters.  Positive values indicate that a phylum is more abundant in particles while negative values indicates a phylum is more abundant in free-water.  Signifant _  

Figure 1a corresponds to the significance values presented in table 1 (below).  13 of the 25 phyla present in the surface waters of the lakes that I sampled are significantly different in particle-associated and free-living fractions.  

On the other hand, Figure 1b corresponds to table 2 (below).  20 of 33 phyla present in the bottom waters of the lakes that I sampled are significantly different in particle-associated and free-living fractions.  

####Overlapping Significance Between Surface and Bottom Waters


_Betaproteobacteria, TA18, Actinobacteria, Armatimonadetes, Candidate Division OD1, Candidate Division OP3, Cyanobacteria, NPL-UPA2,_ and _Planctomycetes_ are significantally different between particle-association and free-living in **both** surface and bottom waters.  




####Surface Water Significant results: Figure 1a 
**Table 1:** Ran 25 Paired Wilcoxon Tests with a Bonferonni Corrected P-value Cut-Off of 0.002  

Phylum | Original P-Value | Bonferroni P-Value
--- | --- | ---
Alphaproteobacteria | 0.0006103516 | 2.441406e-05
Betaproteobacteria | 0.0001220703 | 4.882813e-06
Gammaproteobacteria | 0.0302764099 | 1.211056e-03
TA18 | 0.0222541248 | 8.901650e-04
Actinobacteria | 0.0001220703 | 4.882813e-06
Armatimonadetes | 0.0041430084 | 1.657203e-04
Candidate_division_OD1 | 0.0038398001 | 1.535920e-04
Candidate_division_OP3 | 0.0141474039 | 5.658962e-04
Chlamydiae | 0.0108269213 | 4.330769e-04
Cyanobacteria | 0.0001220703 | 4.882813e-06
NPL-UPA2 | 0.0050339674 | 2.013587e-04
Planctomycetes | 0.0001220703 | 4.882813e-06
Spirochaetae | 0.0016508126 | 6.603250e-05

####Bottom Water Significant results:  Figure 1b
**Table 2: **Ran 33 Paired Wilcoxon Tests with a Bonferonni Corrected P-value Cut-Off of 0.0015152  

Phylum | Original P-Value | Bonferroni P-Value
--- | --- | ---
Betaproteobacteria | 3.204271e-04 | 9.709911e-06
Deltaproteobacteria | 2.578735e-03 | 7.814350e-05
Epsilonproteobacteria | 5.921537e-03 | 1.794405e-04
TA18 | 3.461056e-02 | 1.048805e-03
Actinobacteria | 1.525879e-05 | 4.623875e-07
Armatimonadetes | 2.164494e-02 | 6.559071e-04
BD1-5 | 4.454755e-02 | 1.349926e-03
Candidate_division_BRC1 | 5.857099e-03 | 1.774879e-04
Candidate_division_OD1 | 5.846590e-04 | 1.771694e-05
Candidate_division_OP3 | 1.090835e-03 | 3.305561e-05
Candidate_division_SR1 | 3.748593e-02 | 1.135937e-03
Candidate_division_WS3 | 5.857099e-03 | 1.774879e-04
Chloroflexi | 2.578735e-03 | 7.814350e-05
Cyanobacteria | 1.525879e-05 | 4.623875e-07
Fibrobacteres | 5.889270e-03 | 1.784627e-04
Lentisphaerae | 1.069882e-02 | 3.242066e-04
NPL-UPA2 | 7.265139e-04 | 2.201557e-05
Planctomycetes | 2.090454e-03 | 6.334709e-05
Tenericutes | 5.793045e-03 | 1.755468e-04
TM6 | 2.343773e-03 | 7.102344e-05




###Conclusions
1.  What bacterial phyla are the dominant members of the overall freshwater lake community that I sampled?
      + **Conclusion:**  5 most abundant Phyla:  _Bacteriodetes, Cyanobacteria, Verrucomicrobia, Betaproteobacteria_ and _Actinobacteria._
2.  Do particle-associated bacteria differ from free-living bacteria in **surface** waters?  
      + **Hypothesis:** Particle-associated and free-living bacteria **do not** significantly differ in surface waters.  
        + **Conclusion:** 13 of the 25 phyla in the surface waters differ significantly between particle-associated and free-living fractions.     
3.  Do particle-associated bacteria differ from free-living bacteria in **bottom** waters?
      + **Hypothesis:** Particle-associated and free-living bacteria **do** significantly differ in surface waters.  
        + **Conclusion:** 20 of the 33 phyla in the surface waters differ significantly between particle-associated and free-living fractions.    

**Overlap between the surface and bottom waters:** 9 phyla are significant in both the surface and bottom waters.  Including:_Betaproteobacteria, TA18, Actinobacteria, Armatimonadetes, Candidate Division OD1, Candidate Division OP3, Cyanobacteria, NPL-UPA2,_ and _Planctomycetes_ are significantally different between particle-association and free-living in **both** surface and bottom waters.  


###References
+ Azam, F. (1998)  Microbial control of oceanic carbon flux: the plot thickens, _Science 280(5364):_ 694-696.  
+ Bižić-Ionescu, M., Zeder, M., Ionescu, D., Orlić, S., Fuchs, B.M., Grossart, H.P., and R. Amann (2014)  Compariston of bacterial communities on limnic versus coastal marine particles reveals profound differences in colonization, _Environmental Microbiology:_ 1-15. 
+ Ganesh, S., Parris, D.J., DeLong, E.F., and F. J. Stewart. (2014) Metagenomic analysis of size-fractionated picoplankton in a marine oxygen minimum zone, _The International Society for Microbial Ecology Journal 8:_ 187 – 211.  
+ Grossart, H.P. (2010) Ecological consequences of bacterioplankton lifestyles: changes in concepts are needed, _Environmental Microbiology Reports 2:_ 706-714.  
+ Newton, R.J., Jones, S.E., Eiler, A., McMahon, K.D., and S. Bertilsson. (2011) A guide to the natural history of freshwater lake bacteria, _Microbiology and Molecular Biology Reviews 75(1):_ 14-49.  



